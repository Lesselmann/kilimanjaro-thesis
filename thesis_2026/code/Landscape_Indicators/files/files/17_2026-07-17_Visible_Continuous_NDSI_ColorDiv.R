### ============================================================
### Visible NDSI + Visible Color Diversity PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Complements Visible Composition (categorical, land cover) with two
# CONTINUOUS visibility-filtered indicators:
#
#  - Visible NDSI: is the Kilimanjaro summit/glacier actually VISIBLE
#    from this point? Direct Imageability relevance.
#  - Visible Color Diversity: colour diversity ONLY within the
#    visible area.
#
# MULTI-POINT SAMPLING: As with viewshed_analysis.R and
# visible_composition_analysis.R, this version samples 5 points per
# cell (quincunx: centroid + 4 points offset 250m towards N/S/E/W)
# and averages the resulting metrics, for consistency and to better
# approximate the areal characteristics of the full 1km2 cell.
#
# NDVI is deliberately NOT recomputed with visibility filtering,
# since visible_naturalness already covers this conceptually.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
dsm_path_utm  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/kilimanjaro_dsm_cop30_utm37s.tif"
fab_path_utm  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/fabdem_kilimanjaro_v3_utm37s.tif"
sentinel_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Indices/kilimanjaro_ndvi_ndsi_ephemera.tif"
rgb_path      <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Satellite/kilimanjaro_rgb.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_visible_continuous.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/visible_continuous_results.csv"

max_distance <- 10000
observer_h   <- 1.7
target_h     <- 0
n_clusters   <- 5
subpoint_offset_m <- 250  # consistent with the other multi-point scripts

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(dsm_path_utm), file.exists(fab_path_utm),
          file.exists(sentinel_path), file.exists(rgb_path))

dsm <- rast(dsm_path_utm)
fab <- rast(fab_path_utm)
sentinel <- rast(sentinel_path)
rgb <- rast(rgb_path)

cat("Sentinel bands:", paste(names(sentinel), collapse = ", "), "\n")
cat("RGB bands:", paste(names(rgb), collapse = ", "), "\n")

ndsi_band <- sentinel[["NDSI"]]

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), dsm)) pud_grid <- st_transform(pud_grid, crs(dsm))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

centroids <- st_centroid(pud_grid)
coords    <- st_coordinates(centroids)

set.seed(42)

# ---------------------------------------------------------------
# 2. FUNCTION: VISIBLE NDSI + COLOR DIVERSITY AT ONE POINT
# ---------------------------------------------------------------
compute_visible_continuous_point <- function(pt, dsm, fab, ndsi_band, rgb, max_distance,
                                              observer_h, target_h, n_clusters) {

  na_result <- list(ndsi = NA, color_div = NA, brightness = NA)

  buf      <- st_buffer(st_sfc(st_point(pt), crs = "EPSG:32737"), max_distance)
  dsm_crop <- tryCatch(crop(dsm, vect(buf)), error = function(e) NULL)
  if (is.null(dsm_crop)) return(na_result)

  ground_h <- tryCatch(extract(fab, matrix(pt, ncol = 2))[1, 1], error = function(e) NA)
  if (!is.na(ground_h)) {
    obs_cell <- cellFromXY(dsm_crop, matrix(pt, ncol = 2))
    if (!is.na(obs_cell)) dsm_crop[obs_cell] <- ground_h
  }

  vs <- tryCatch(
    viewshed(dsm_crop, loc = pt, observer = observer_h, target = target_h),
    error = function(e) NULL
  )
  if (is.null(vs)) return(na_result)

  # ---- NDSI ----
  ndsi_crop <- crop(ndsi_band, ext(vs))
  vs_resampled_ndsi <- resample(vs, ndsi_crop, method = "near")
  ndsi_vals <- values(ndsi_crop)[values(vs_resampled_ndsi) == 1]
  ndsi_vals <- ndsi_vals[!is.na(ndsi_vals)]
  ndsi_mean <- if (length(ndsi_vals) > 0) mean(ndsi_vals) else NA

  # ---- RGB / Color Diversity ----
  rgb_crop <- crop(rgb, ext(vs))
  vs_resampled_rgb <- resample(vs, rgb_crop, method = "near")
  rgb_vals <- values(rgb_crop)
  vis_mask_rgb <- values(vs_resampled_rgb) == 1
  rgb_visible <- rgb_vals[vis_mask_rgb & complete.cases(rgb_vals), , drop = FALSE]

  color_div <- NA
  brightness <- NA
  if (nrow(rgb_visible) >= 20) {
    km <- tryCatch(kmeans(rgb_visible, centers = min(n_clusters, nrow(unique(rgb_visible))), nstart = 5),
                   error = function(e) NULL)
    if (!is.null(km)) {
      cluster_props <- table(km$cluster) / length(km$cluster)
      p <- cluster_props[cluster_props > 0]
      color_div <- -sum(p * log(p))
      brightness <- mean(rowMeans(rgb_visible))
    }
  }

  list(ndsi = ndsi_mean, color_div = color_div, brightness = brightness)
}

# ---------------------------------------------------------------
# 2b. FUNCTION: 5-POINT QUINCUNX AVERAGE PER CELL
# ---------------------------------------------------------------
compute_visible_continuous_cell_average <- function(center_x, center_y, dsm, fab, ndsi_band, rgb,
                                                      max_distance, observer_h, target_h,
                                                      n_clusters, offset_m) {
  sub_points <- list(
    c(center_x, center_y),
    c(center_x, center_y + offset_m),
    c(center_x, center_y - offset_m),
    c(center_x + offset_m, center_y),
    c(center_x - offset_m, center_y)
  )

  sub_results <- lapply(sub_points, function(p) {
    tryCatch(
      compute_visible_continuous_point(c(X = p[1], Y = p[2]), dsm, fab, ndsi_band, rgb,
                                        max_distance, observer_h, target_h, n_clusters),
      error = function(e) list(ndsi = NA, color_div = NA, brightness = NA)
    )
  })

  sub_df <- bind_rows(lapply(sub_results, as.data.frame))
  list(mean_ndsi = mean(sub_df$ndsi, na.rm = TRUE),
       mean_color_div = mean(sub_df$color_div, na.rm = TRUE),
       mean_brightness = mean(sub_df$brightness, na.rm = TRUE))
}

# ---------------------------------------------------------------
# 3. MAIN LOOP (5-point average per cell)
# ---------------------------------------------------------------
results <- vector("list", nrow(coords))

for (i in seq_len(nrow(coords))) {
  pt_x <- unname(coords[i, "X"])
  pt_y <- unname(coords[i, "Y"])
  stopifnot(!is.na(pt_x), !is.na(pt_y))

  cell_avg <- compute_visible_continuous_cell_average(pt_x, pt_y, dsm, fab, ndsi_band, rgb,
                                                        max_distance, observer_h, target_h,
                                                        n_clusters, subpoint_offset_m)

  results[[i]] <- data.frame(
    cell_id = i,
    visible_NDSI_mean = cell_avg$mean_ndsi,
    visible_ColorDiversity = cell_avg$mean_color_div,
    visible_MeanBrightness = cell_avg$mean_brightness
  )

  if (i %% 20 == 0) cat(sprintf("Cell %d/%d done\n", i, nrow(coords)))
}

results_df <- bind_rows(results)
cat(sprintf("\n-> NDSI: %d/%d successful | ColorDiversity: %d/%d successful\n",
            sum(!is.na(results_df$visible_NDSI_mean)), nrow(results_df),
            sum(!is.na(results_df$visible_ColorDiversity)), nrow(results_df)))

# ---------------------------------------------------------------
# 4. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, results_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("visible_NDSI_mean", "visible_ColorDiversity", "visible_MeanBrightness")]))
