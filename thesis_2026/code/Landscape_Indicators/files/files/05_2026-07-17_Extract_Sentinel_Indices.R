### ============================================================
### Sentinel-2 Indicators PER GRID CELL (NDVI, NDSI, NDVI_var)
### Kilimanjaro Thesis
### ============================================================
# Goal: Extract the three bands from your GEE composite
# (kilimanjaro_ndvi_ndsi_ephemera.tif) via zonal statistics
# (mean) onto your 228 PUD grid cells.
#
# Much simpler than the landscape metrics scripts: no patch/shape
# computation needed, just a direct mean per cell and band.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS (please adjust)
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"

# Path to your GEE export file downloaded from Google Drive
sentinel_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Indices/kilimanjaro_ndvi_ndsi_ephemera.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_sentinel.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/sentinel_results.csv"

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(sentinel_path))
sentinel <- rast(sentinel_path)

cat("Bands in the Sentinel-2 composite:\n")
print(names(sentinel))
cat(sprintf("Resolution: %.1f x %.1f m | CRS: %s\n\n", res(sentinel)[1], res(sentinel)[2], crs(sentinel, describe = TRUE)$name))

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
pud_grid$cell_id <- seq_len(nrow(pud_grid))

# Align CRS (your composite may already be in UTM 37S from GEE, since
# you requested crs='EPSG:32737' at export time -- checked/aligned
# here regardless, for safety)
if (!same.crs(vect(pud_grid), sentinel)) {
  cat("CRS differs -- reprojecting pud_grid to the Sentinel CRS...\n")
  pud_grid <- st_transform(pud_grid, crs(sentinel))
}

# ---------------------------------------------------------------
# 2. ZONAL STATISTICS PER CELL (mean per band)
# ---------------------------------------------------------------
cat("Extracting zonal means per grid cell...\n")

extracted <- terra::extract(sentinel, vect(pud_grid), fun = mean, na.rm = TRUE, ID = FALSE)

# Align column names with your existing naming scheme
names(extracted) <- gsub("NDVI_mean", "NDVI_mean", names(extracted))
names(extracted) <- gsub("NDSI", "NDSI", names(extracted))
names(extracted) <- gsub("NDVI_var", "NDVI_var", names(extracted))

extracted$cell_id <- pud_grid$cell_id

cat("\n--- Preview of extracted values ---\n")
print(head(extracted, 10))

# ---------------------------------------------------------------
# 3. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, extracted, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone: Sentinel-2 indicators extracted for %d grid cells.\n", nrow(pud_grid)))
cat(sprintf("Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("NDVI_mean", "NDSI", "NDVI_var")]))

# ---------------------------------------------------------------
# 4. PLAUSIBILITY CHECK
# ---------------------------------------------------------------
cat("\n--- Plausibility check ---\n")
cat(sprintf("NDVI_mean range: %.3f to %.3f (expected: -0.2 to 0.9)\n",
            min(pud_grid$NDVI_mean, na.rm = TRUE), max(pud_grid$NDVI_mean, na.rm = TRUE)))
cat(sprintf("NDSI range: %.3f to %.3f (expected: -0.5 to 0.8, higher values only near the summit/snow)\n",
            min(pud_grid$NDSI, na.rm = TRUE), max(pud_grid$NDSI, na.rm = TRUE)))
cat(sprintf("NDVI_var range: %.3f to %.3f (expected: 0 to ~0.3)\n",
            min(pud_grid$NDVI_var, na.rm = TRUE), max(pud_grid$NDVI_var, na.rm = TRUE)))
cat(sprintf("NAs: NDVI_mean=%d, NDSI=%d, NDVI_var=%d\n",
            sum(is.na(pud_grid$NDVI_mean)), sum(is.na(pud_grid$NDSI)), sum(is.na(pud_grid$NDVI_var))))
