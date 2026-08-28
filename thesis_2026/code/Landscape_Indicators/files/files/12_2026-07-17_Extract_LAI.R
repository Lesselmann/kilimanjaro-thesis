### ============================================================
### Leaf Area Index + LAI Variability PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Extracts LAI_mean AND LAI_var from the MODIS composite
# (export_lai_modis.js) via zonal statistics (mean) onto your
# grid cells.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"

# LAI: from export_lai_modis.js (Google Drive -> GEE_exports/kilimanjaro_lai.tif)
lai_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Indices/kilimanjaro_lai.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_lai.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/lai_results.csv"

# ---------------------------------------------------------------
# 1. LOAD PUD GRID
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_path, quiet = TRUE)
pud_grid$cell_id <- seq_len(nrow(pud_grid))

# ---------------------------------------------------------------
# 2. EXTRACT LAI (both bands)
# ---------------------------------------------------------------
stopifnot(file.exists(lai_path))
cat("Loading Leaf Area Index (mean + variability)...\n")
lai <- rast(lai_path)
cat("Bands found:", paste(names(lai), collapse = ", "), "\n")

grid_lai <- pud_grid
if (!same.crs(vect(grid_lai), lai)) grid_lai <- st_transform(grid_lai, crs(lai))

lai_extracted <- terra::extract(lai, vect(grid_lai), fun = mean, na.rm = TRUE, ID = FALSE)
lai_extracted$cell_id <- pud_grid$cell_id

cat(sprintf("-> LAI_mean: range %.3f to %.3f, NAs: %d\n",
            min(lai_extracted$LAI_mean, na.rm = TRUE),
            max(lai_extracted$LAI_mean, na.rm = TRUE),
            sum(is.na(lai_extracted$LAI_mean))))
cat(sprintf("-> LAI_var: range %.3f to %.3f, NAs: %d\n",
            min(lai_extracted$LAI_var, na.rm = TRUE),
            max(lai_extracted$LAI_var, na.rm = TRUE),
            sum(is.na(lai_extracted$LAI_var))))

# ---------------------------------------------------------------
# 3. JOIN TO pud_grid, SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, lai_extracted, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(pud_grid[, c("LAI_mean", "LAI_var")] %>% st_drop_geometry()))
