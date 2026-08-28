### ============================================================
### Visible Landscape Composition (Viewshed x WorldCover)
### Kilimanjaro Thesis
### ============================================================
# Goal: For each grid cell, not just "how much is visible"
# (visible_area_km2 etc. from viewshed_analysis.R), but "WHAT is
# visible" -- i.e. the land cover composition WITHIN the visible area.
#
# MULTI-POINT SAMPLING: As with viewshed_analysis.R, this version
# samples 5 points per cell (quincunx: centroid + 4 points offset
# 250m towards N/S/E/W) and averages the resulting composition
# metrics, for consistency and to better approximate the areal
# characteristics of the full 1km2 cell.
#
# WorldCover classes (reference):
#  10 Tree cover   20 Shrubland      30 Grassland    40 Cropland
#  50 Built-up     60 Bare/sparse    70 Snow/ice     80 Water
#  90 Herb. wetland 95 Mangroves     100 Moss/lichen
#
# Derived indicators per grid cell (all cell-averaged across 5 points):
#  - visible_naturalness : share (tree+shrub+grass+wetland+mangrove+moss) of visible area
#  - visible_water        : share of water (class 80) of visible area
#  - visible_anthropogenic: share (cropland+built-up) of visible area
#  - visible_diversity     : Shannon diversity of the visible class distribution
#  - class_XX columns      : full raw distribution (share of each class), averaged across points
#
# Note: Not tested in an R environment. Please check parameter/function
# names against your installed package versions.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. ADJUST PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"

dsm_path_utm  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/kilimanjaro_dsm_cop30_utm37s.tif"
fab_path_utm  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/fabdem_kilimanjaro_v3_utm37s.tif"
lc_utm_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Landcover/kilimanjaro_worldcover.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_visible_composition.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/visible_composition_results.csv"

max_distance <- 10000  # identical to viewshed_analysis.R (Schirpke et al. 2013 boundary) -- do not change for comparability
observer_h   <- 1.7
target_h     <- 0
near_cutoff  <- 1500    # metres -- near/middle zone boundary, identical to viewshed_analysis.R
subpoint_offset_m <- 250 # offset of the 4 off-centre sample points, consistent with viewshed_analysis.R

# WorldCover class assignment for the derived indicators
natural_classes <- c(10, 20, 30, 70, 90, 95, 100)  # tree, shrub, grass, snow/ice, wetland, mangrove, moss
anthro_classes  <- c(40, 50)                        # cropland, built-up area
water_class     <- 80

# ---------------------------------------------------------------
# 1. LOAD RASTERS (already reprojected, from earlier scripts)
# ---------------------------------------------------------------
stopifnot(file.exists(dsm_path_utm), file.exists(fab_path_utm), file.exists(lc_utm_path))

dsm <- rast(dsm_path_utm)
fab <- rast(fab_path_utm)
lc  <- rast(lc_utm_path)

# ---------------------------------------------------------------
# 2. LOAD PUD GRID
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), dsm)) {
  pud_grid <- st_transform(pud_grid, crs(dsm))
}
pud_grid$cell_id <- seq_len(nrow(pud_grid))

centroids <- st_centroid(pud_grid)
coords    <- st_coordinates(centroids)

# ---------------------------------------------------------------
# 3. HELPER FUNCTION: COMPOSITION METRICS FROM A VECTOR OF CLASSES
# ---------------------------------------------------------------
compute_composition <- function(lc_classes) {
  if (length(lc_classes) == 0) {
    return(list(naturalness = NA, water = NA, anthropogenic = NA, diversity = NA, props = NULL))
  }
  class_counts <- table(lc_classes)
  class_props  <- as.numeric(class_counts) / sum(class_counts)
  names(class_props) <- names(class_counts)

  naturalness   <- sum(class_props[names(class_props) %in% as.character(natural_classes)])
  water         <- sum(class_props[names(class_props) %in% as.character(water_class)])
  anthropogenic <- sum(class_props[names(class_props) %in% as.character(anthro_classes)])

  p <- class_props[class_props > 0]
  diversity <- -sum(p * log(p))

  list(naturalness = naturalness, water = water, anthropogenic = anthropogenic,
       diversity = diversity, props = class_props)
}

# ---------------------------------------------------------------
# 4. FUNCTION: VISIBLE COMPOSITION AT ONE POINT
# ---------------------------------------------------------------
compute_visible_composition_point <- function(pt, dsm, fab, lc, max_distance, observer_h,
                                               target_h, near_cutoff) {

  na_result <- list(naturalness = NA, water = NA, anthropogenic = NA, diversity = NA,
                     naturalness_near = NA, water_near = NA, anthropogenic_near = NA, diversity_near = NA,
                     naturalness_middle = NA, water_middle = NA, anthropogenic_middle = NA, diversity_middle = NA,
                     props = NULL)

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

  lc_crop <- crop(lc, ext(vs))
  vs_resampled <- resample(vs, lc_crop, method = "near")

  vis_vals <- values(vs_resampled)
  lc_vals  <- values(lc_crop)
  valid_mask <- vis_vals == 1 & !is.na(lc_vals)

  if (sum(valid_mask, na.rm = TRUE) == 0) return(na_result)

  vis_cell_idx <- which(valid_mask)
  vis_xy <- xyFromCell(lc_crop, vis_cell_idx)
  vis_dist <- sqrt((vis_xy[, 1] - pt[1])^2 + (vis_xy[, 2] - pt[2])^2)
  vis_lc_class <- lc_vals[vis_cell_idx]

  comp_total  <- compute_composition(vis_lc_class)
  comp_near   <- compute_composition(vis_lc_class[vis_dist <= near_cutoff])
  comp_middle <- compute_composition(vis_lc_class[vis_dist > near_cutoff])

  list(naturalness = comp_total$naturalness, water = comp_total$water,
       anthropogenic = comp_total$anthropogenic, diversity = comp_total$diversity,
       naturalness_near = comp_near$naturalness, water_near = comp_near$water,
       anthropogenic_near = comp_near$anthropogenic, diversity_near = comp_near$diversity,
       naturalness_middle = comp_middle$naturalness, water_middle = comp_middle$water,
       anthropogenic_middle = comp_middle$anthropogenic, diversity_middle = comp_middle$diversity,
       props = comp_total$props)
}

# ---------------------------------------------------------------
# 5. FUNCTION: 5-POINT QUINCUNX AVERAGE PER CELL
# ---------------------------------------------------------------
compute_visible_composition_cell_average <- function(center_x, center_y, dsm, fab, lc,
                                                       max_distance, observer_h, target_h,
                                                       near_cutoff, offset_m) {
  sub_points <- list(
    c(center_x, center_y),
    c(center_x, center_y + offset_m),
    c(center_x, center_y - offset_m),
    c(center_x + offset_m, center_y),
    c(center_x - offset_m, center_y)
  )

  sub_results <- lapply(sub_points, function(p) {
    tryCatch(
      compute_visible_composition_point(c(X = p[1], Y = p[2]), dsm, fab, lc,
                                         max_distance, observer_h, target_h, near_cutoff),
      error = function(e) NULL
    )
  })
  sub_results <- sub_results[!sapply(sub_results, is.null)]

  numeric_fields <- c("naturalness", "water", "anthropogenic", "diversity",
                       "naturalness_near", "water_near", "anthropogenic_near", "diversity_near",
                       "naturalness_middle", "water_middle", "anthropogenic_middle", "diversity_middle")

  averaged <- lapply(numeric_fields, function(f) {
    vals <- sapply(sub_results, function(r) r[[f]])
    mean(vals, na.rm = TRUE)
  })
  names(averaged) <- numeric_fields

  # Average the raw class proportions across points too (union of all classes seen)
  all_props <- lapply(sub_results, function(r) r$props)
  all_classes <- unique(unlist(lapply(all_props, names)))
  if (length(all_classes) > 0) {
    avg_props <- sapply(all_classes, function(cl) {
      vals <- sapply(all_props, function(p) if (cl %in% names(p)) p[[cl]] else 0)
      mean(vals, na.rm = TRUE)
    })
    names(avg_props) <- all_classes
    raw_row <- as.list(setNames(avg_props, paste0("class_", all_classes)))
  } else {
    raw_row <- list()
  }

  c(averaged, raw_row)
}

# ---------------------------------------------------------------
# 6. MAIN LOOP (5-point average per cell)
# ---------------------------------------------------------------
results <- vector("list", nrow(coords))

for (i in seq_len(nrow(coords))) {
  pt_x <- unname(coords[i, "X"])
  pt_y <- unname(coords[i, "Y"])
  stopifnot(!is.na(pt_x), !is.na(pt_y))

  cell_avg <- compute_visible_composition_cell_average(pt_x, pt_y, dsm, fab, lc, max_distance,
                                                         observer_h, target_h, near_cutoff, subpoint_offset_m)

  results[[i]] <- c(
    list(cell_id = i,
         visible_naturalness   = cell_avg$naturalness,
         visible_water         = cell_avg$water,
         visible_anthropogenic = cell_avg$anthropogenic,
         visible_diversity     = cell_avg$diversity,
         visible_naturalness_near   = cell_avg$naturalness_near,
         visible_water_near         = cell_avg$water_near,
         visible_anthropogenic_near = cell_avg$anthropogenic_near,
         visible_diversity_near     = cell_avg$diversity_near,
         visible_naturalness_middle   = cell_avg$naturalness_middle,
         visible_water_middle         = cell_avg$water_middle,
         visible_anthropogenic_middle = cell_avg$anthropogenic_middle,
         visible_diversity_middle     = cell_avg$diversity_middle),
    cell_avg[grepl("^class_", names(cell_avg))]
  )

  if (i %% 20 == 0) cat(sprintf("Cell %d/%d done\n", i, nrow(coords)))
}

results_df <- bind_rows(lapply(results, as.data.frame))

# ---------------------------------------------------------------
# 7. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, results_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone: visible composition (5-point cell average) computed for %d grid cells.\n", nrow(pud_grid)))
cat(sprintf("Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary (total 0-10km) ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("visible_naturalness", "visible_water",
                                              "visible_anthropogenic", "visible_diversity")]))

cat("\n--- Summary (near zone 0-1.5km) ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("visible_naturalness_near", "visible_water_near",
                                              "visible_anthropogenic_near", "visible_diversity_near")]))

cat("\n--- Summary (middle zone 1.5-10km) ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("visible_naturalness_middle", "visible_water_middle",
                                              "visible_anthropogenic_middle", "visible_diversity_middle")]))
