### ============================================================
### Viewshed Analysis for PUD Grid Cells (Kilimanjaro Thesis)
### ============================================================
# Goal: Compute viewshed metrics per grid cell as predictors.
#
# MULTI-POINT SAMPLING: Rather than computing the viewshed only at
# the cell centroid, this version samples 5 points per cell (a
# "quincunx" pattern: centroid + 4 points offset 250m towards
# N/S/E/W) and averages the resulting metrics. This better
# approximates the true areal average across the full 1km2 cell,
# consistent with how the zonal-statistic-based indicators (SHDI,
# NDVI_mean, etc.) already represent the whole cell rather than a
# single point. A single centroid point remains a common, defensible
# simplification (cf. Schirpke et al. 2013's discrete viewpoint
# sampling), but this averages over that simplification's main
# weakness.
#
# FABDEM correction: at each observer point, the DSM value is
# replaced by the FABDEM ground height so that, for points located
# in forested areas, the observer does not "stand" on the tree
# canopy. The surrounding area remains unchanged DSM (forest still
# acts as a visual obstruction).
#
# Note: This script was NOT tested in an R environment (no R/terra
# available here). Please check the parameter names of
# terra::viewshed() against your installed terra version
# (?terra::viewshed), as the API may differ slightly between versions.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. ADJUST PATHS
# ---------------------------------------------------------------
pud_grid_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
dsm_path_utm    <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/kilimanjaro_dsm_cop30_utm37s.tif"
fab_path_utm    <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/fabdem_kilimanjaro_v3_utm37s.tif"
out_gpkg        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_viewshed.gpkg"
out_csv         <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/viewshed_metrics.csv"

# ---------------------------------------------------------------
# 1. LOAD DSM + FABDEM (already reprojected to UTM 37S)
# ---------------------------------------------------------------
stopifnot(file.exists(dsm_path_utm), file.exists(fab_path_utm))
dsm <- rast(dsm_path_utm)
fab <- rast(fab_path_utm)

cat("DSM CRS:", crs(dsm, describe = TRUE)$name, "| resolution:", res(dsm), "\n")
cat("FABDEM CRS:", crs(fab, describe = TRUE)$name, "| resolution:", res(fab), "\n")

# ---------------------------------------------------------------
# 2. LOAD PUD GRID AND ALIGN CRS
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), dsm)) pud_grid <- st_transform(pud_grid, crs(dsm))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

# ---------------------------------------------------------------
# 2b. COVERAGE CHECK
# ---------------------------------------------------------------
grid_bbox <- st_bbox(pud_grid)
dsm_ext   <- ext(dsm)
if (grid_bbox["xmin"] < dsm_ext[1] || grid_bbox["xmax"] > dsm_ext[2] ||
    grid_bbox["ymin"] < dsm_ext[3] || grid_bbox["ymax"] > dsm_ext[4]) {
  warning("WARNING: pud_grid extends partially beyond the DSM extent! ",
          "Edge grid cells may receive incomplete viewsheds.")
}

# ---------------------------------------------------------------
# 3. EXTRACT CENTROIDS
# ---------------------------------------------------------------
centroids <- st_centroid(pud_grid)
coords    <- st_coordinates(centroids)

# ---------------------------------------------------------------
# 4. PARAMETERS
# ---------------------------------------------------------------
max_distance <- 10000   # search radius in metres (Schirpke et al. 2013 boundary)
observer_h   <- 1.7
target_h     <- 0
near_cutoff  <- 1500     # near/middle zone boundary, metres
subpoint_offset_m <- 250 # offset of the 4 off-centre sample points from the cell centroid

# ---------------------------------------------------------------
# 5. FUNCTION: COMPUTE VIEWSHED METRICS FOR ONE POINT
# ---------------------------------------------------------------
compute_viewshed_metrics <- function(pt, dsm, fab, max_distance, observer_h, target_h, near_cutoff) {

  na_result <- list(visible_area_km2 = NA, visible_prop = NA, dist_sd_m = NA,
                     visible_area_near_km2 = NA, visible_prop_near = NA, dist_sd_near_m = NA,
                     visible_area_middle_km2 = NA, visible_prop_middle = NA, dist_sd_middle_m = NA)

  buf      <- st_buffer(st_sfc(st_point(pt), crs = "EPSG:32737"), max_distance)
  dsm_crop <- tryCatch(crop(dsm, vect(buf)), error = function(e) NULL)
  if (is.null(dsm_crop)) return(na_result)

  # FABDEM correction at this specific observer point
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

  vis_vals      <- values(vs)
  n_visible     <- sum(vis_vals == 1, na.rm = TRUE)
  n_total       <- sum(!is.na(vis_vals))
  cell_area_km2 <- prod(res(vs)) / 1e6

  vis_cells <- which(vis_vals == 1)

  result <- na_result
  result$visible_area_km2 <- n_visible * cell_area_km2
  result$visible_prop <- ifelse(n_total > 0, n_visible / n_total, NA)

  if (length(vis_cells) > 1) {
    vis_xy  <- xyFromCell(vs, vis_cells)
    dists   <- sqrt((vis_xy[, 1] - pt[1])^2 + (vis_xy[, 2] - pt[2])^2)

    result$dist_sd_m <- sd(dists, na.rm = TRUE)

    near_dists   <- dists[dists <= near_cutoff]
    middle_dists <- dists[dists > near_cutoff]

    near_zone_area_km2   <- pi * (near_cutoff / 1000)^2
    middle_zone_area_km2 <- pi * (max_distance / 1000)^2 - near_zone_area_km2

    result$visible_area_near_km2   <- length(near_dists) * cell_area_km2
    result$visible_area_middle_km2 <- length(middle_dists) * cell_area_km2
    result$visible_prop_near   <- result$visible_area_near_km2 / near_zone_area_km2
    result$visible_prop_middle <- result$visible_area_middle_km2 / middle_zone_area_km2
    result$dist_sd_near_m   <- if (length(near_dists) > 1) sd(near_dists) else NA
    result$dist_sd_middle_m <- if (length(middle_dists) > 1) sd(middle_dists) else NA
  }

  result
}

# ---------------------------------------------------------------
# 5b. FUNCTION: 5-POINT QUINCUNX WITHIN A CELL, AVERAGE THE METRICS
# ---------------------------------------------------------------
compute_viewshed_cell_average <- function(center_x, center_y, dsm, fab, max_distance,
                                           observer_h, target_h, near_cutoff, offset_m) {
  sub_points <- list(
    c(center_x, center_y),
    c(center_x, center_y + offset_m),
    c(center_x, center_y - offset_m),
    c(center_x + offset_m, center_y),
    c(center_x - offset_m, center_y)
  )

  sub_results <- lapply(sub_points, function(p) {
    tryCatch(
      compute_viewshed_metrics(c(X = p[1], Y = p[2]), dsm, fab, max_distance, observer_h, target_h, near_cutoff),
      error = function(e) list(visible_area_km2 = NA, visible_prop = NA, dist_sd_m = NA,
                                visible_area_near_km2 = NA, visible_prop_near = NA, dist_sd_near_m = NA,
                                visible_area_middle_km2 = NA, visible_prop_middle = NA, dist_sd_middle_m = NA)
    )
  })

  sub_df <- bind_rows(sub_results)
  # Average each metric across the (up to) 5 valid subpoints
  as.list(colMeans(sub_df, na.rm = TRUE))
}

# ---------------------------------------------------------------
# 6. MAIN LOOP (5-point average per cell)
# ---------------------------------------------------------------
results <- vector("list", nrow(coords))

for (i in seq_len(nrow(coords))) {
  pt_x <- unname(coords[i, "X"])
  pt_y <- unname(coords[i, "Y"])
  stopifnot(!is.na(pt_x), !is.na(pt_y))

  cell_avg <- compute_viewshed_cell_average(pt_x, pt_y, dsm, fab, max_distance,
                                             observer_h, target_h, near_cutoff, subpoint_offset_m)

  results[[i]] <- data.frame(cell_id = i,
                              visible_area_km2 = cell_avg$visible_area_km2,
                              visible_prop = cell_avg$visible_prop,
                              dist_sd_m = cell_avg$dist_sd_m,
                              visible_area_near_km2 = cell_avg$visible_area_near_km2,
                              visible_prop_near = cell_avg$visible_prop_near,
                              dist_sd_near_m = cell_avg$dist_sd_near_m,
                              visible_area_middle_km2 = cell_avg$visible_area_middle_km2,
                              visible_prop_middle = cell_avg$visible_prop_middle,
                              dist_sd_middle_m = cell_avg$dist_sd_middle_m)

  if (i %% 20 == 0) cat(sprintf("Viewshed %d/%d done\n", i, nrow(coords)))
}

viewshed_df <- bind_rows(results)

# ---------------------------------------------------------------
# 7. JOIN BACK TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, viewshed_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("Done: viewshed metrics (5-point cell average) computed for %d grid cells.\n", nrow(pud_grid)))
cat(sprintf("Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("visible_area_km2", "visible_prop", "dist_sd_m",
                                              "visible_prop_near", "dist_sd_near_m",
                                              "visible_prop_middle", "dist_sd_middle_m")]))
