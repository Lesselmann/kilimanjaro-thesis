### ============================================================
### Distance to Summit PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Goal: Euclidean distance from each grid cell centroid to the
# Kilimanjaro summit (Uhuru Peak) -- Imageability indicator
# (landmark distance following Lynch 1960).
### ============================================================

library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_summit_distance.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/summit_distance_results.csv"

# Uhuru Peak, Kilimanjaro summit (WGS84)
# Coordinates: 3°04'33"S, 37°21'12"E (most commonly cited value,
# consistent across several independent sources)
summit_lon <- 37.3533
summit_lat <- -3.0758

# ---------------------------------------------------------------
# 1. CREATE SUMMIT POINT AND PROJECT TO TARGET CRS
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_path, quiet = TRUE)
pud_grid$cell_id <- seq_len(nrow(pud_grid))

summit_pt_wgs84 <- st_sfc(st_point(c(summit_lon, summit_lat)), crs = 4326)
summit_pt <- st_transform(summit_pt_wgs84, crs = st_crs(pud_grid))

cat("Summit coordinate (projected):\n")
print(st_coordinates(summit_pt))

# ---------------------------------------------------------------
# 2. COMPUTE DISTANCE PER CENTROID
# ---------------------------------------------------------------
centroids <- st_centroid(pud_grid)
distances_m <- as.numeric(st_distance(centroids, summit_pt))

pud_grid$DistanceToSummit_km <- distances_m / 1000

cat(sprintf("\nDistance range: %.2f to %.2f km\n",
            min(pud_grid$DistanceToSummit_km), max(pud_grid$DistanceToSummit_km)))

# ---------------------------------------------------------------
# 3. SAVE
# ---------------------------------------------------------------
st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nSaved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(pud_grid$DistanceToSummit_km))

