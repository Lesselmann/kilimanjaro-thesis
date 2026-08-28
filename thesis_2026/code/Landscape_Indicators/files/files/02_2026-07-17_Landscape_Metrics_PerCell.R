### ============================================================
### Landscape Metrics Analysis PER GRID CELL (individual metrics)
### Kilimanjaro Thesis
### ============================================================
# Like the fast per-cell version, but EACH metric is computed
# INDIVIDUALLY (not all 9 bundled per cell), with dataframe output
# for checking after each individual metric -- so that problems
# with a metric (like the earlier ShapeIndex bug) are noticed
# immediately, instead of seeing everything only at the end.
#
# Both earlier fixes are already built in:
#  - ShapeIndex: area-weighted mean from class level
#    (since 'shape_mn' is missing at landscape level in some versions)
#  - WaterTotalPatchArea/WaterEdgeLength/ForestEcotoneDensity:
#    NA (class missing) is converted to 0 (a genuinely measured absence)
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
cat(sprintf("Using %d of %d available CPU cores (deliberately limited).\n", n_cores, parallel::detectCores()))

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
lc_utm_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Landcover/kilimanjaro_worldcover.tif"
fab_path_utm  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/fabdem_kilimanjaro_v3_utm37s.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_landscapemetrics.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/landscapemetrics_results.csv"

buffer_m <- 100

hemeroby_lookup <- data.frame(
  class = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
  hemeroby = c(1, 2, 2, 6, 7, 3, 1, 1, 1, 1, 1)
)

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(lc_utm_path), file.exists(fab_path_utm))
lc  <- rast(lc_utm_path)
dtm <- rast(fab_path_utm)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), lc)) pud_grid <- st_transform(pud_grid, crs(lc))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

grid_bbox <- st_bbox(pud_grid)
crop_ext <- ext(grid_bbox["xmin"] - buffer_m - 1000, grid_bbox["xmax"] + buffer_m + 1000,
                 grid_bbox["ymin"] - buffer_m - 1000, grid_bbox["ymax"] + buffer_m + 1000)
lc  <- crop(lc, crop_ext)
dtm <- crop(dtm, crop_ext)
cat(sprintf("Raster cropped (WorldCover: %.0f MB, FABDEM: %.0f MB in memory)\n\n",
            as.numeric(object.size(wrap(lc))) / 1e6, as.numeric(object.size(wrap(dtm))) / 1e6))

# ---------------------------------------------------------------
# 2. STEP A: CROP SEQUENTIALLY PER CELL
# ---------------------------------------------------------------
cat("Step A: Cropping small extracts per cell (sequentially)...\n")
cell_crops <- vector("list", nrow(pud_grid))
failed_cells <- character(0)
pb_step <- max(1, round(nrow(pud_grid) / 10))
for (i in seq_len(nrow(pud_grid))) {
  cell_poly <- vect(pud_grid[i, ])
  cell_buf  <- buffer(cell_poly, buffer_m)

  crop_result <- tryCatch({
    lc_crop  <- crop(lc, cell_buf, mask = TRUE)
    dtm_crop <- crop(dtm, cell_buf, mask = TRUE)
    list(id = pud_grid$cell_id[i], lc = wrap(lc_crop), dtm = wrap(dtm_crop))
  }, error = function(e) {
    message(sprintf("Cell %d (cell_id=%s): crop failed (%s) -- likely outside raster extent, skipping",
                     i, pud_grid$cell_id[i], conditionMessage(e)))
    failed_cells <<- c(failed_cells, as.character(pud_grid$cell_id[i]))
    NULL
  })

  cell_crops[[i]] <- crop_result

  if (i %% pb_step == 0) cat(sprintf("  %d/%d cells cropped\n", i, nrow(pud_grid)))
}

if (length(failed_cells) > 0) {
  cat(sprintf("\nWARNING: %d cell(s) could not be cropped (likely outside DSM/FABDEM extent):\n", length(failed_cells)))
  print(failed_cells)
  cat("These cells will have NA for all metrics in this script. Consider re-downloading\n")
  cat("a larger DSM/FABDEM extent covering the full 373-cell grid if this is unexpected.\n")
}

# Remove failed cells from further processing (keeps indices aligned via NULL entries)
cell_crops <- cell_crops[!sapply(cell_crops, is.null)]

cat("Step A done.\n\n")

rm(lc, dtm)
gc()

# =================================================================
# METRIC 1: SHDI
# =================================================================
cat("--- Metric 1: SHDI ---\n")
shdi_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  res <- tryCatch(lsm_l_shdi(lcc), error = function(e) NULL)
  data.frame(cell_id = cc$id, SHDI = if (!is.null(res) && nrow(res) > 0) res$value[1] else NA)
}, future.seed = TRUE)
shdi_df <- bind_rows(shdi_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(shdi_df$SHDI)), nrow(shdi_df)))
print(head(shdi_df, 5))
cat("\n")

# =================================================================
# METRIC 2: Number of Patches
# =================================================================
cat("--- Metric 2: NumLandElements (np) ---\n")
np_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  res <- tryCatch(lsm_l_np(lcc), error = function(e) NULL)
  data.frame(cell_id = cc$id, NumLandElements = if (!is.null(res) && nrow(res) > 0) res$value[1] else NA)
}, future.seed = TRUE)
np_df <- bind_rows(np_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(np_df$NumLandElements)), nrow(np_df)))
print(head(np_df, 5))
cat("\n")

# =================================================================
# METRIC 3: Shape Index (area-weighted from class level)
# =================================================================
cat("--- Metric 3: ShapeIndex (area-weighted) ---\n")
shape_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  # Direct function calls instead of calculate_lsm() with a combined
  # metric vector -- calculate_lsm(metric = c("shape_mn","pland"))
  # silently dropped shape_mn (a bug/quirk of the version), even
  # though both metrics work reliably individually.
  shape_res <- tryCatch(lsm_c_shape_mn(lcc), error = function(e) NULL)
  pland_res <- tryCatch(lsm_c_pland(lcc), error = function(e) NULL)
  val <- NA
  if (!is.null(shape_res) && !is.null(pland_res)) {
    merged <- merge(shape_res[, c("class", "value")], pland_res[, c("class", "value")],
                     by = "class", suffixes = c("_shape", "_pland"))
    if (nrow(merged) > 0) val <- weighted.mean(merged$value_shape, w = merged$value_pland, na.rm = TRUE)
  }
  data.frame(cell_id = cc$id, ShapeIndex = val)
}, future.seed = TRUE)
shape_df <- bind_rows(shape_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(shape_df$ShapeIndex)), nrow(shape_df)))
print(head(shape_df, 5))
cat("\n")

# =================================================================
# METRIC 4+5: Water Total Patch Area + Edge Length
# =================================================================
cat("--- Metric 4+5: WaterTotalPatchArea + WaterEdgeLength ---\n")
water_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  ca_res <- tryCatch(lsm_c_ca(lcc), error = function(e) NULL)
  te_res <- tryCatch(lsm_c_te(lcc), error = function(e) NULL)
  ca_val <- if (!is.null(ca_res)) { v <- ca_res$value[ca_res$class == 80]; if (length(v) == 0) 0 else v[1] } else 0
  te_val <- if (!is.null(te_res)) { v <- te_res$value[te_res$class == 80]; if (length(v) == 0) 0 else v[1] } else 0
  data.frame(cell_id = cc$id, WaterTotalPatchArea = ca_val, WaterEdgeLength = te_val)
}, future.seed = TRUE)
water_df <- bind_rows(water_list)
cat(sprintf("-> %d/%d cells with water present\n", sum(water_df$WaterTotalPatchArea > 0), nrow(water_df)))
print(head(water_df, 5))
cat("\n")

# =================================================================
# METRIC 6: Open Space Proportion (1 - Built-up%)
# =================================================================
cat("--- Metric 6: OpenSpaceProportion ---\n")
open_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  res <- tryCatch(lsm_c_pland(lcc), error = function(e) NULL)
  builtup <- if (!is.null(res)) { v <- res$value[res$class == 50]; if (length(v) == 0) 0 else v[1] } else 0
  data.frame(cell_id = cc$id, OpenSpaceProportion = 1 - builtup / 100)
}, future.seed = TRUE)
open_df <- bind_rows(open_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(open_df$OpenSpaceProportion)), nrow(open_df)))
print(head(open_df, 5))
cat("\n")

# =================================================================
# METRIC 7: Forest Edge Density
# =================================================================
cat("--- Metric 7: ForestEcotoneDensity ---\n")
forest_list <- future_lapply(cell_crops, function(cc) {
  library(terra); library(landscapemetrics)
  lcc <- unwrap(cc$lc)
  res <- tryCatch(lsm_c_ed(lcc), error = function(e) NULL)
  val <- if (!is.null(res)) { v <- res$value[res$class == 10]; if (length(v) == 0) 0 else v[1] } else 0
  data.frame(cell_id = cc$id, ForestEcotoneDensity = val)
}, future.seed = TRUE)
forest_df <- bind_rows(forest_list)
cat(sprintf("-> %d/%d cells with forest present\n", sum(forest_df$ForestEcotoneDensity > 0), nrow(forest_df)))
print(head(forest_df, 5))
cat("\n")

# =================================================================
# METRIC 8: Hemeroby Index
# =================================================================
cat("--- Metric 8: HemerobyIndex ---\n")
hemeroby_list <- future_lapply(cell_crops, function(cc) {
  library(terra)
  lcc <- unwrap(cc$lc)
  val <- tryCatch({
    hr <- classify(lcc, as.matrix(hemeroby_lookup))
    as.numeric(global(hr, "mean", na.rm = TRUE)[1, 1])
  }, error = function(e) NA)
  data.frame(cell_id = cc$id, HemerobyIndex = val)
}, future.seed = TRUE)
hemeroby_df <- bind_rows(hemeroby_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(hemeroby_df$HemerobyIndex)), nrow(hemeroby_df)))
print(head(hemeroby_df, 5))
cat("\n")

# =================================================================
# METRIC 9: Relief Diversity
# =================================================================
cat("--- Metric 9: ReliefDiversity ---\n")
relief_list <- future_lapply(cell_crops, function(cc) {
  library(terra)
  dtmc <- unwrap(cc$dtm)
  val <- tryCatch({
    surf3d <- surfArea(dtmc, filename = "")
    planar_area_m2 <- prod(res(dtmc))
    ratio <- surf3d / planar_area_m2
    as.numeric(global(ratio, "mean", na.rm = TRUE)[1, 1])
  }, error = function(e) NA)
  data.frame(cell_id = cc$id, ReliefDiversity = val)
}, future.seed = TRUE)
relief_df <- bind_rows(relief_list)
cat(sprintf("-> %d/%d successful\n", sum(!is.na(relief_df$ReliefDiversity)), nrow(relief_df)))
print(head(relief_df, 5))
cat("\n")

# ---------------------------------------------------------------
# 3. COMBINE ALL METRICS, JOIN TO pud_grid, SAVE
# ---------------------------------------------------------------
results_df <- shdi_df %>%
  left_join(np_df, by = "cell_id") %>%
  left_join(shape_df, by = "cell_id") %>%
  left_join(water_df, by = "cell_id") %>%
  left_join(open_df, by = "cell_id") %>%
  left_join(forest_df, by = "cell_id") %>%
  left_join(hemeroby_df, by = "cell_id") %>%
  left_join(relief_df, by = "cell_id")

pud_grid <- left_join(pud_grid, results_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone: landscape metrics computed for %d grid cells.\n", nrow(pud_grid)))
cat(sprintf("Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Overall summary ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("SHDI", "ShapeIndex", "NumLandElements", "HemerobyIndex",
                                              "WaterTotalPatchArea", "WaterEdgeLength",
                                              "ReliefDiversity", "OpenSpaceProportion", "ForestEcotoneDensity")]))
