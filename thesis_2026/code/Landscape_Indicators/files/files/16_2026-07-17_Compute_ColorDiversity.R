### ============================================================
### Color Diversity + Color Diversity Variability PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Methodology following Vaz et al. (2020, Conservation Letters):
# RGB pixels within each grid cell are grouped into colour clusters
# via k-means, then a Shannon-Wiener index is computed across the
# cluster distribution.
#
# NOW ALSO computes ColorDiversity_var: ColorDiversity is computed
# SEPARATELY for each of 4 seasonal RGB composites, then the
# standard deviation across the 4 seasonal values is taken --
# analogous to NDVI_var and LAI_var (an Ephemera-relevant measure
# of how much colour diversity itself fluctuates seasonally).
#
# Requires:
#   - kilimanjaro_rgb.tif (dry-season composite, for ColorDiversity/MeanBrightness)
#   - kilimanjaro_rgb_seasonal.tif (12-band, 4 seasons x R/G/B, for ColorDiversity_var)
#   both from export_rgb_bands.js, downloaded from Google Drive.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
rgb_path        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Satellite/kilimanjaro_rgb.tif"
rgb_seasonal_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Satellite/kilimanjaro_rgb_seasonal.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_colordiv.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/colordiv_results.csv"

buffer_m <- 100
n_clusters <- 5
min_pixels_for_clustering <- 20
season_names <- c("dry1", "dry2", "wet1", "wet2")

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(rgb_path), file.exists(rgb_seasonal_path))
rgb <- rast(rgb_path)
rgb_seasonal <- rast(rgb_seasonal_path)
cat("Dry-season bands:", paste(names(rgb), collapse = ", "), "\n")
cat("Seasonal bands:", paste(names(rgb_seasonal), collapse = ", "), "\n")

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), rgb)) pud_grid <- st_transform(pud_grid, crs(rgb))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

set.seed(42)

# ---------------------------------------------------------------
# 2. HELPER FUNCTION: COLOR DIVERSITY FROM A SET OF RGB VALUES
# ---------------------------------------------------------------
compute_color_diversity_from_vals <- function(vals, n_clusters, min_pixels) {
  vals <- vals[complete.cases(vals), , drop = FALSE]
  if (nrow(vals) < min_pixels) return(list(diversity = NA, brightness = NA))

  km <- tryCatch(
    kmeans(vals, centers = min(n_clusters, nrow(unique(vals))), nstart = 5),
    error = function(e) NULL
  )
  if (is.null(km)) return(list(diversity = NA, brightness = NA))

  cluster_props <- table(km$cluster) / length(km$cluster)
  p <- cluster_props[cluster_props > 0]
  diversity <- -sum(p * log(p))
  brightness <- mean(rowMeans(vals))
  list(diversity = diversity, brightness = brightness)
}

# ---------------------------------------------------------------
# 3. PER CELL: DRY-SEASON ColorDiversity/MeanBrightness (as before)
#    PLUS per-season ColorDiversity -> ColorDiversity_var
# ---------------------------------------------------------------
cat(sprintf("Computing Color Diversity (+ variability) for %d cells...\n", nrow(pud_grid)))

results <- vector("list", nrow(pud_grid))

for (i in seq_len(nrow(pud_grid))) {
  cell_poly <- vect(pud_grid[i, ])
  cell_buf  <- buffer(cell_poly, buffer_m)

  # ---- Dry-season ColorDiversity / MeanBrightness (unchanged) ----
  rgb_crop <- crop(rgb, cell_buf, mask = TRUE)
  vals <- values(rgb_crop)
  dry_result <- compute_color_diversity_from_vals(vals, n_clusters, min_pixels_for_clustering)

  # ---- Per-season ColorDiversity, for variability ----
  seasonal_crop <- crop(rgb_seasonal, cell_buf, mask = TRUE)
  seasonal_diversities <- numeric(length(season_names))

  for (s in seq_along(season_names)) {
    season_bands <- paste0(season_names[s], c("_Red", "_Green", "_Blue"))
    season_vals <- values(seasonal_crop)[, season_bands, drop = FALSE]
    season_result <- compute_color_diversity_from_vals(season_vals, n_clusters, min_pixels_for_clustering)
    seasonal_diversities[s] <- season_result$diversity
  }

  color_div_var <- if (sum(!is.na(seasonal_diversities)) >= 2) {
    sd(seasonal_diversities, na.rm = TRUE)
  } else {
    NA
  }

  results[[i]] <- data.frame(
    cell_id = i,
    ColorDiversity = dry_result$diversity,
    MeanBrightness = dry_result$brightness,
    ColorDiversity_var = color_div_var
  )

  if (i %% 20 == 0) cat(sprintf("  Cell %d/%d done\n", i, nrow(pud_grid)))
}

colordiv_df <- bind_rows(results)
cat(sprintf("\n-> ColorDiversity: %d/%d successful | ColorDiversity_var: %d/%d successful\n",
            sum(!is.na(colordiv_df$ColorDiversity)), nrow(colordiv_df),
            sum(!is.na(colordiv_df$ColorDiversity_var)), nrow(colordiv_df)))
print(head(colordiv_df, 10))

# ---------------------------------------------------------------
# 4. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, colordiv_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("ColorDiversity", "MeanBrightness", "ColorDiversity_var")]))
