### ============================================================
### Distance to Nearest Waterway (line features)
### Kilimanjaro Thesis
### ============================================================
# Computes DistanceToWaterway_km: Euclidean distance from each
# cell centroid to the nearest WATERWAY LINE feature (rivers,
# streams -- OSM "waterway" tag), as distinct from
# DistanceToWater_km, which (per the existing indicator table)
# measures distance to WATER BODIES / polygons.
#
# Rationale: rivers/streams are frequently linear features not
# well captured by polygon water bodies, and are theoretically
# relevant under the Coherence concept (Ode et al. 2008) alongside
# DistanceToWater_km -- consistent with the merging of waterway
# lines into the water class already done for the categorical
# visible_water indicator.
#
# Data source: OpenStreetMap "waterway" lines (river, stream,
# canal, drain), fetched via osmdata. If you already have a local
# waterway shapefile (e.g. from HydroSHEDS or a prior OSM extract
# used for the visible_water merge step), set use_local_file <- TRUE
# and point local_waterway_path at it instead -- avoids re-querying
# the Overpass API.
### ============================================================

library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. ADJUST PATHS / OPTIONS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"

use_local_file      <- FALSE   # TRUE = read existing waterway shapefile instead of querying OSM
local_waterway_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Hydrology/waterways_kilimanjaro.gpkg"

# Bounding box for the OSM query, in WGS84 (lon/lat) -- use the
# SAME extended bounding box already used for your GEE exports,
# to guarantee consistent spatial coverage
bbox_wgs84 <- c(xmin = 36.9, ymin = -3.5, xmax = 37.9, ymax = -2.7)  # ADJUST to your actual bbox

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_dist_waterway.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/dist_waterway_results.csv"

# ---------------------------------------------------------------
# 1. LOAD GRID
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_path, quiet = TRUE)
pud_grid$cell_id <- seq_len(nrow(pud_grid))
centroids <- st_centroid(pud_grid)

# ---------------------------------------------------------------
# 2. LOAD WATERWAY LINES
# ---------------------------------------------------------------
if (use_local_file) {

  stopifnot(file.exists(local_waterway_path))
  waterways <- st_read(local_waterway_path, quiet = TRUE)

} else {

  if (!requireNamespace("osmdata", quietly = TRUE)) {
    install.packages("osmdata")
  }
  library(osmdata)

  cat("Querying OSM Overpass API for waterway lines in the study bbox...\n")
  cat("(river, stream, canal, drain -- adjust the 'waterway' value filter if needed)\n")

  osm_query <- opq(bbox = bbox_wgs84, timeout = 120) %>%
    add_osm_feature(key = "waterway", value = c("river", "stream", "canal", "drain")) %>%
    osmdata_sf()

  waterways <- osm_query$osm_lines
  stopifnot(!is.null(waterways), nrow(waterways) > 0)
  cat(sprintf("Retrieved %d waterway line features from OSM.\n", nrow(waterways)))
}

# Reproject to the same CRS as the grid (should be UTM 37S, metric,
# for meaningful distance calculations)
if (!same.crs(vect(waterways), vect(pud_grid))) {
  waterways <- st_transform(waterways, st_crs(pud_grid))
}
if (!same.crs(vect(centroids), vect(pud_grid))) {
  centroids <- st_transform(centroids, st_crs(pud_grid))
}

# Union into a single geometry for fast nearest-distance queries
waterways_union <- st_union(waterways)

# ---------------------------------------------------------------
# 3. DISTANCE CALCULATION (Euclidean, cell centroid to nearest
#    waterway line), identical logic to DistanceToWater_km
# ---------------------------------------------------------------
cat("Computing distance to nearest waterway per cell...\n")

dist_m <- st_distance(centroids, waterways_union)
dist_m <- as.numeric(dist_m)   # matrix -> vector (one column, since waterways_union is a single geometry)

pud_grid$DistanceToWaterway_km <- dist_m / 1000

# ---------------------------------------------------------------
# 4. SAVE
# ---------------------------------------------------------------
st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))
cat("\n--- Summary ---\n")
print(summary(pud_grid$DistanceToWaterway_km))

cat("\nNOTE: If the Overpass query times out or returns too few features,\n")
cat("consider using a local HydroSHEDS/OSM extract instead\n")
cat("(set use_local_file <- TRUE above).\n")
