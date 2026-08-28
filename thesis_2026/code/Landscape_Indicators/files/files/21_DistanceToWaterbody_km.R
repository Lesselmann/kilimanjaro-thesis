### ============================================================
### Distance to Nearest OSM Water Body (natural=water polygons)
### Kilimanjaro Thesis
### ============================================================
# Computes DistanceToWaterbody_OSM_km: Euclidean distance from each
# cell centroid to the nearest OSM natural=water polygon boundary
# (lakes, reservoirs, ponds -- human-digitized, not spectrally
# derived, and therefore not subject to the mountain-shadow
# misclassification documented for ESA WorldCover in steep terrain).
#
# NAMING: explicitly tagged "_OSM_" to distinguish it from the
# EXISTING WorldCover-based DistanceToWater_km. Both are kept as
# separate candidates for now -- run the scatterplot self-check on
# both, compare which shows a clearer/stronger relationship with
# avg_annual_PUD, then drop whichever is weaker before the final
# screening run (or let Spearman screening decide automatically,
# since they are likely to be highly correlated with each other).
#
# DistanceToWaterway_km (OSM waterway LINES, rivers/streams) stays
# a third, separate candidate -- conceptually distinct (lines vs.
# polygons, streams vs. standing water), not part of this comparison.
#
# IMPORTANT: only natural=='water' polygons are used here --
# NOT natural=='wetland' (a much larger, separate OSM category,
# e.g. the 145 km2 marsh polygon found in Waterbodies.gpkg, which
# corresponds conceptually to WorldCover class 90, not class 80).
### ============================================================

library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. ADJUST PATHS
# ---------------------------------------------------------------
pud_grid_path      <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
waterbodies_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Shape/Waterbodies.gpkg"   # ADJUST if moved

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_dist_waterbody_osm.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/dist_waterbody_osm_results.csv"

# ---------------------------------------------------------------
# 1. LOAD GRID
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_path, quiet = TRUE)
pud_grid$cell_id <- seq_len(nrow(pud_grid))
centroids <- st_centroid(pud_grid)

# ---------------------------------------------------------------
# 2. LOAD OSM WATER BODIES -- filter to natural=='water' ONLY
#    (excludes natural=='wetland', a different, much larger
#    category not comparable to WorldCover's "Water" class)
# ---------------------------------------------------------------
waterbodies_raw <- st_read(waterbodies_path, quiet = TRUE)

cat("Feature counts in Waterbodies.gpkg by 'natural' tag:\n")
print(table(waterbodies_raw$natural, useNA = "ifany"))

waterbodies <- waterbodies_raw %>% filter(natural == "water")
cat(sprintf("\nUsing %d 'natural=water' polygons (excluded %d wetland + %d other/NA features).\n",
            nrow(waterbodies),
            sum(waterbodies_raw$natural == "wetland", na.rm = TRUE),
            nrow(waterbodies_raw) - nrow(waterbodies) - sum(waterbodies_raw$natural == "wetland", na.rm = TRUE)))

# Reproject to match the grid CRS (should already both be UTM 37S)
if (!same.crs(vect(waterbodies), vect(pud_grid))) {
  waterbodies <- st_transform(waterbodies, st_crs(pud_grid))
}
if (!same.crs(vect(centroids), vect(pud_grid))) {
  centroids <- st_transform(centroids, st_crs(pud_grid))
}

waterbodies_union <- st_union(st_make_valid(waterbodies))

# ---------------------------------------------------------------
# 3. DISTANCE CALCULATION (Euclidean, cell centroid to nearest
#    water body polygon boundary), identical logic to
#    DistanceToWater_km / DistanceToWaterway_km
# ---------------------------------------------------------------
cat("\nComputing distance to nearest OSM water body per cell...\n")

dist_m <- st_distance(centroids, waterbodies_union)
dist_m <- as.numeric(dist_m)

pud_grid$DistanceToWaterbody_OSM_km <- dist_m / 1000

# ---------------------------------------------------------------
# 4. SAVE
# ---------------------------------------------------------------
st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))
cat("\n--- Summary ---\n")
print(summary(pud_grid$DistanceToWaterbody_OSM_km))

cat("\nNext step: merge dist_waterbody_osm_results.csv into the master table\n")
cat("alongside the EXISTING DistanceToWater_km (WorldCover-based) -- both\n")
cat("stay in as separate candidates for now. Run the scatterplot self-check\n")
cat("script to visually compare DistanceToWater_km vs. DistanceToWaterbody_OSM_km\n")
cat("against avg_annual_PUD (rho/edf/p in the caption), then decide which one\n")
cat("to drop before your next screening run -- or let Spearman screening\n")
cat("catch the redundancy automatically if you keep both.\n")
