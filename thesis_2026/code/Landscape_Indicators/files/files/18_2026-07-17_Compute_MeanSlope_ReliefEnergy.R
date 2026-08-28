### ============================================================
### Mean Slope + Relief Energy PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Two classic, simple terrain indicators from Table 3
# (Schirpke et al. 2013; Bishop & Hulse 1994), not yet
# operationalised previously (only Relief Diversity had been built):
#
#  - Mean Slope: average slope angle (degrees)
#  - Relief Energy: elevation difference (max-min) within the cell
#
# Both from FABDEM (bare earth, consistent with Relief Diversity --
# no mixing of terrain complexity with vegetation height).
# Technically trivial: terrain() function + simple zonal extraction,
# no per-cell loop needed.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
fab_path_utm  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/fabdem_kilimanjaro_v3_utm37s.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_slope_relief.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/slope_relief_results.csv"

buffer_m <- 100  # same buffer as the other per-cell terrain indicators

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(fab_path_utm))
fab <- rast(fab_path_utm)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), fab)) pud_grid <- st_transform(pud_grid, crs(fab))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

# Buffer polygons for the zonal extraction (consistent with the
# other terrain indicators)
pud_grid_buf <- st_buffer(pud_grid, buffer_m)

# ---------------------------------------------------------------
# 2. MEAN SLOPE (computed once over the whole area, then
#    zonally extracted -- faster than a per-cell loop)
# ---------------------------------------------------------------
cat("Computing slope raster...\n")
slope_raster <- terrain(fab, "slope", unit = "degrees")

cat("Extracting mean slope per cell...\n")
mean_slope <- terra::extract(slope_raster, vect(pud_grid_buf), fun = mean, na.rm = TRUE, ID = FALSE)
pud_grid$MeanSlope_deg <- mean_slope[, 1]

# ---------------------------------------------------------------
# 3. RELIEF ENERGY (max - min elevation per cell)
# ---------------------------------------------------------------
cat("Extracting min/max elevation per cell...\n")
min_elev <- terra::extract(fab, vect(pud_grid_buf), fun = min, na.rm = TRUE, ID = FALSE)
max_elev <- terra::extract(fab, vect(pud_grid_buf), fun = max, na.rm = TRUE, ID = FALSE)

pud_grid$ReliefEnergy_m <- max_elev[, 1] - min_elev[, 1]

# ---------------------------------------------------------------
# 4. SAVE
# ---------------------------------------------------------------
st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(pud_grid[, c("MeanSlope_deg", "ReliefEnergy_m")] %>% st_drop_geometry()))
