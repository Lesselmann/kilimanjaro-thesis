### ============================================================
### Contagion Index PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Standard FRAGSTATS metric (lsm_l_contag) -- measures how strongly
# land cover patches are spatially "clumped" (high values) vs.
# scattered/interspersed (low values). Fills the previously entirely
# vacant Coherence concept WITHOUT requiring a self-defined class
# contrast matrix (in contrast to the Edge Contrast Index).
#
# Direct function call (lsm_l_contag), not via calculate_lsm() with a
# metric vector -- consistent with the fix applied in the other
# landscape metrics scripts (calculate_lsm() silently dropped some
# metrics when given combined metric vectors).
### ============================================================

library(terra)
library(sf)
library(dplyr)
library(landscapemetrics)
library(future)
library(future.apply)

terraOptions(threads = TRUE)
options(future.globals.maxSize = 1000 * 1024^2)

n_cores <- min(4, max(1, parallel::detectCores() - 1))
plan(multisession, workers = n_cores)
cat(sprintf("Using %d of %d available CPU cores.\n", n_cores, parallel::detectCores()))

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
lc_utm_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Landcover/kilimanjaro_worldcover.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_contagion.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/contagion_results.csv"

buffer_m <- 100

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(lc_utm_path))
lc <- rast(lc_utm_path)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), lc)) pud_grid <- st_transform(pud_grid, crs(lc))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

grid_bbox <- st_bbox(pud_grid)
crop_ext <- ext(grid_bbox["xmin"] - buffer_m - 1000, grid_bbox["xmax"] + buffer_m + 1000,
                 grid_bbox["ymin"] - buffer_m - 1000, grid_bbox["ymax"] + buffer_m + 1000)
lc <- crop(lc, crop_ext)

# ---------------------------------------------------------------
# 2. STEP A: CROP SEQUENTIALLY PER CELL
# ---------------------------------------------------------------
cat("Cropping small extracts per cell (sequentially)...\n")
cell_crops <- vector("list", nrow(pud_grid))
for (i in seq_len(nrow(pud_grid))) {
  cell_poly <- vect(pud_grid[i, ])
  cell_buf  <- buffer(cell_poly, buffer_m)
  cell_crops[[i]] <- list(id = pud_grid$cell_id[i], lc = wrap(crop(lc, cell_buf, mask = TRUE)))
  if (i %% 20 == 0) cat(sprintf("  %d/%d cells cropped\n", i, nrow(pud_grid)))
}
rm(lc)
gc()

# ---------------------------------------------------------------
# 3. COMPUTE CONTAGION INDEX IN PARALLEL
# ---------------------------------------------------------------
cat("Computing Contagion Index...\n")
contagion_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  res <- tryCatch(lsm_l_contag(lcc), error = function(e) NULL)
  data.frame(cell_id = cc$id,
             ContagionIndex = if (!is.null(res) && nrow(res) > 0) res$value[1] else NA)
}, future.seed = TRUE)

contagion_df <- bind_rows(contagion_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(contagion_df$ContagionIndex)), nrow(contagion_df)))
print(head(contagion_df, 10))

# ---------------------------------------------------------------
# 4. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, contagion_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary (0-100, high = clumped/coherent, low = fragmented) ---\n")
print(summary(pud_grid$ContagionIndex))
