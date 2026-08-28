### ============================================================
### Refuge / Enclosure Proximity PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Operationalisation: share of forest (class 10) + shrubland
# (class 20) within a 250m buffer -- complements Prospect (viewshed)
# with the Refuge component of Prospect-Refuge theory (Appleton
# 1975): how much shelter/enclosure potential exists in the
# immediate vicinity of the observation point.
#
# Deliberately small radius (250m) in contrast to the 10km viewshed --
# Refuge is a local, not a landscape-wide, measure.
#
# MULTI-POINT SAMPLING: As with viewshed_analysis.R and
# compute_sky_view_factor.R, this version samples 5 points per cell
# (quincunx: centroid + 4 points offset 250m towards N/S/E/W) and
# averages the resulting RefugeIndex values, for consistency and to
# better approximate the areal average across the full 1km2 cell
# rather than relying on a single centroid estimate. Note that with
# a 250m buffer radius AND a 250m subpoint offset, the five buffers
# overlap substantially -- this is expected and acceptable here,
# since the goal is a smoother areal estimate, not five independent
# samples.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
lc_utm_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Landcover/kilimanjaro_worldcover.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_refuge.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/refuge_results.csv"

refuge_radius_m   <- 250   # local buffer -- deliberately small, see header
refuge_classes    <- c(10, 20)  # forest, shrubland
subpoint_offset_m <- 250   # offset of the 4 off-centre sample points, consistent with the other scripts

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(lc_utm_path))
lc <- rast(lc_utm_path)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), lc)) pud_grid <- st_transform(pud_grid, crs(lc))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

centroids <- st_centroid(pud_grid)
coords    <- st_coordinates(centroids)

# ---------------------------------------------------------------
# 2. FUNCTION: REFUGE INDEX AT ONE POINT
# ---------------------------------------------------------------
compute_refuge_point <- function(pt_x, pt_y, lc, refuge_radius_m, refuge_classes) {
  buf <- st_buffer(st_sfc(st_point(c(pt_x, pt_y)), crs = "EPSG:32737"), refuge_radius_m)
  lc_crop <- tryCatch(crop(lc, vect(buf), mask = TRUE), error = function(e) NULL)
  if (is.null(lc_crop)) return(NA)

  vals <- values(lc_crop)
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) return(NA)
  sum(vals %in% refuge_classes) / length(vals)
}

# ---------------------------------------------------------------
# 2b. FUNCTION: 5-POINT QUINCUNX AVERAGE PER CELL
# ---------------------------------------------------------------
compute_refuge_cell_average <- function(center_x, center_y, lc, refuge_radius_m,
                                         refuge_classes, offset_m) {
  sub_points <- list(
    c(center_x, center_y),
    c(center_x, center_y + offset_m),
    c(center_x, center_y - offset_m),
    c(center_x + offset_m, center_y),
    c(center_x - offset_m, center_y)
  )

  sub_vals <- sapply(sub_points, function(p) {
    tryCatch(compute_refuge_point(p[1], p[2], lc, refuge_radius_m, refuge_classes),
             error = function(e) NA)
  })

  list(mean_refuge = mean(sub_vals, na.rm = TRUE),
       sd_refuge = if (sum(!is.na(sub_vals)) > 1) sd(sub_vals, na.rm = TRUE) else NA)
}

# ---------------------------------------------------------------
# 3. PER CELL (5-point average)
# ---------------------------------------------------------------
cat(sprintf("Computing Refuge index for %d cells (%dm buffer, 5-point cell average)...\n",
            nrow(coords), refuge_radius_m))

refuge_results <- vector("list", nrow(coords))

for (i in seq_len(nrow(coords))) {
  pt_x <- coords[i, "X"]
  pt_y <- coords[i, "Y"]

  cell_result <- compute_refuge_cell_average(pt_x, pt_y, lc, refuge_radius_m,
                                              refuge_classes, subpoint_offset_m)

  refuge_results[[i]] <- data.frame(cell_id = i,
                                     RefugeIndex = cell_result$mean_refuge,
                                     RefugeIndex_sd_within_cell = cell_result$sd_refuge)

  if (i %% 20 == 0) cat(sprintf("  Cell %d/%d done\n", i, nrow(coords)))
}

refuge_df <- bind_rows(refuge_results)
cat(sprintf("\n-> %d/%d cells successful\n", sum(!is.na(refuge_df$RefugeIndex)), nrow(refuge_df)))
print(head(refuge_df, 10))

# ---------------------------------------------------------------
# 4. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, refuge_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone: Refuge index computed for %d grid cells.\n", nrow(pud_grid)))
cat(sprintf("Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary (0 = no forest/shrub nearby, 1 = fully enclosed) ---\n")
print(summary(pud_grid$RefugeIndex))
