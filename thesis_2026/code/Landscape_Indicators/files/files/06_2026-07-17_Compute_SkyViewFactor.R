### ============================================================
### Sky View Factor (SVF) PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Native R implementation, NO SAGA/GRASS required (solution for
# the earlier software-dependency problem).
#
# Method: Classic horizon-angle approach (cf. Dozier & Frew 1990;
# Zaksek et al. 2011, simplified variant):
#  1. From each centroid, elevation profiles along the DSM are
#     sampled in N evenly spaced compass directions (azimuths) up
#     to a search radius.
#  2. Per direction, the maximum elevation angle (horizon angle)
#     to the observer is determined -- this is the angle beyond
#     which terrain/vegetation obstructs the sky.
#  3. SVF = mean across all directions of (1 - sin(horizon_angle))
#     -- 1 = fully open sky, 0 = fully obstructed.
#
# Uses the DSM for the SURROUNDING horizon search (vegetation/buildings
# genuinely obstruct the sky), but the FABDEM (bare-earth) ground
# height AT THE OBSERVER POINT ONLY -- otherwise, in forested cells,
# the 30m DSM pixel value already includes the average canopy height
# at that location, and adding eye height on top of it would place
# the "observer" effectively on top of the tree canopy rather than
# on the ground looking up through it. This mirrors the identical
# correction already applied in viewshed_analysis.R.
#
# MULTI-POINT SAMPLING: Rather than computing SVF only at the cell
# centroid, this version samples 5 points per cell (a "quincunx"
# pattern: centroid + 4 points offset 250m towards N/S/E/W) and
# averages the resulting SVF values. A single centroid point is a
# common, defensible simplification (cf. Schirpke et al. 2013's
# discrete viewpoint sampling), but this multi-point approach better
# approximates the true areal average across the full 1km2 cell,
# consistent with how the zonal-statistic-based indicators (SHDI,
# NDVI_mean, etc.) already represent the whole cell rather than a
# single point.
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

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_svf.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/svf_results.csv"

n_directions   <- 16     # number of azimuth directions (more = more accurate, but slower)
search_radius  <- 500    # metres -- SVF is a LOCAL measure, a smaller radius than viewshed (8km) makes sense
sample_step    <- 30     # metres between sample points along each direction (~DSM resolution)
observer_h     <- 1.7
subpoint_offset_m <- 250 # distance of the 4 off-centre sample points from the cell centroid

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(dsm_path_utm), file.exists(fab_path_utm))
dsm <- rast(dsm_path_utm)
fab <- rast(fab_path_utm)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), dsm)) pud_grid <- st_transform(pud_grid, crs(dsm))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

centroids <- st_centroid(pud_grid)
coords    <- st_coordinates(centroids)

# ---------------------------------------------------------------
# 2. FUNCTION: COMPUTE SVF FOR ONE POINT
# ---------------------------------------------------------------
compute_svf <- function(pt_x, pt_y, dsm_crop, fab_crop, observer_h, n_directions, search_radius, sample_step) {

  # Ground height AT THE OBSERVER POINT from FABDEM (bare earth) --
  # not from DSM, to avoid effectively placing the observer on top
  # of the local tree canopy in forested cells.
  ground_h <- tryCatch(as.numeric(extract(fab_crop, matrix(c(pt_x, pt_y), ncol = 2))[1, 1]),
                        error = function(e) NA)
  if (is.na(ground_h)) return(NA)
  eye_z <- ground_h + observer_h

  azimuths <- seq(0, 360 - 360 / n_directions, length.out = n_directions)
  distances <- seq(sample_step, search_radius, by = sample_step)

  horizon_angles <- numeric(n_directions)

  for (a in seq_along(azimuths)) {
    az_rad <- azimuths[a] * pi / 180
    dx <- sin(az_rad)
    dy <- cos(az_rad)

    sample_x <- pt_x + dx * distances
    sample_y <- pt_y + dy * distances
    sample_pts <- cbind(sample_x, sample_y)

    # Horizon search uses the DSM (surface incl. vegetation/buildings)
    # -- trees and structures around the observer genuinely obstruct
    # the sky, so this part is intentionally NOT the bare-earth model
    elevs <- tryCatch(extract(dsm_crop, sample_pts)[, 1], error = function(e) rep(NA, length(distances)))

    valid <- !is.na(elevs)
    if (sum(valid) == 0) {
      horizon_angles[a] <- 0
      next
    }

    elev_angles <- atan2(elevs[valid] - eye_z, distances[valid])
    horizon_angles[a] <- max(c(0, elev_angles), na.rm = TRUE)  # min. 0 (no negative horizon)
  }

  svf <- mean(1 - sin(pmax(horizon_angles, 0)))
  svf
}

# ---------------------------------------------------------------
# 2b. FUNCTION: 5-POINT QUINCUNX WITHIN A CELL, AVERAGE THE SVF
# ---------------------------------------------------------------
compute_svf_cell_average <- function(center_x, center_y, dsm_crop, fab_crop,
                                      observer_h, n_directions, search_radius,
                                      sample_step, offset_m) {
  sub_points <- list(
    c(center_x, center_y),                 # centroid
    c(center_x, center_y + offset_m),      # north
    c(center_x, center_y - offset_m),      # south
    c(center_x + offset_m, center_y),      # east
    c(center_x - offset_m, center_y)       # west
  )

  sub_svfs <- sapply(sub_points, function(p) {
    tryCatch(
      compute_svf(p[1], p[2], dsm_crop, fab_crop, observer_h, n_directions, search_radius, sample_step),
      error = function(e) NA
    )
  })

  list(mean_svf = mean(sub_svfs, na.rm = TRUE),
       sd_svf = if (sum(!is.na(sub_svfs)) > 1) sd(sub_svfs, na.rm = TRUE) else NA,
       n_valid = sum(!is.na(sub_svfs)))
}

# ---------------------------------------------------------------
# 3. COMPUTE PER CELL (5-point average within each cell)
# ---------------------------------------------------------------
cat(sprintf("Computing SVF for %d cells (%d directions, %dm radius, 5-point cell average)...\n",
            nrow(coords), n_directions, search_radius))

svf_results <- vector("list", nrow(coords))
failed_cells <- character(0)

for (i in seq_len(nrow(coords))) {
  pt_x <- coords[i, "X"]
  pt_y <- coords[i, "Y"]

  # Buffer must cover both the search radius AND the +/-250m subpoint offset
  buf <- st_buffer(st_sfc(st_point(c(pt_x, pt_y)), crs = st_crs(pud_grid)),
                    search_radius + subpoint_offset_m + 50)

  crop_result <- tryCatch({
    list(dsm = crop(dsm, vect(buf)), fab = crop(fab, vect(buf)))
  }, error = function(e) {
    message(sprintf("Cell %d: crop failed (%s) -- likely outside DSM/FABDEM extent, skipping",
                     i, conditionMessage(e)))
    failed_cells <<- c(failed_cells, as.character(i))
    NULL
  })

  if (is.null(crop_result)) {
    svf_results[[i]] <- data.frame(cell_id = i, SkyViewFactor = NA,
                                    SkyViewFactor_sd_within_cell = NA,
                                    SkyViewFactor_n_valid_subpoints = 0)
    next
  }

  cell_result <- tryCatch(
    compute_svf_cell_average(pt_x, pt_y, crop_result$dsm, crop_result$fab, observer_h,
                              n_directions, search_radius, sample_step, subpoint_offset_m),
    error = function(e) list(mean_svf = NA, sd_svf = NA, n_valid = 0)
  )

  svf_results[[i]] <- data.frame(
    cell_id = i,
    SkyViewFactor = cell_result$mean_svf,
    SkyViewFactor_sd_within_cell = cell_result$sd_svf,
    SkyViewFactor_n_valid_subpoints = cell_result$n_valid
  )

  if (i %% 20 == 0) cat(sprintf("  Cell %d/%d done\n", i, nrow(coords)))
}

if (length(failed_cells) > 0) {
  cat(sprintf("\nWARNING: %d cell(s) could not be cropped (likely outside DSM/FABDEM extent): %s\n",
              length(failed_cells), paste(failed_cells, collapse = ", ")))
  cat("Consider re-downloading a larger DSM/FABDEM extent if this is unexpected.\n")
}

svf_df <- bind_rows(svf_results)
cat(sprintf("\n-> %d/%d cells successful\n", sum(!is.na(svf_df$SkyViewFactor)), nrow(svf_df)))
print(head(svf_df, 10))

# ---------------------------------------------------------------
# 4. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, svf_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone: Sky View Factor computed for %d grid cells.\n", nrow(pud_grid)))
cat(sprintf("Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary (expected: 0 = fully obstructed, 1 = fully open) ---\n")
print(summary(pud_grid$SkyViewFactor))

cat("\n--- Within-cell variability (SD across the 5 sample points) ---\n")
cat("Higher values indicate more heterogeneous cells, where a single centroid\n")
cat("point would have been a less representative choice.\n")
print(summary(pud_grid$SkyViewFactor_sd_within_cell))
