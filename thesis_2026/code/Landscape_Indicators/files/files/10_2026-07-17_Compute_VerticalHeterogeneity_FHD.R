### ============================================================
### Vertical Structural Heterogeneity PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Computes the classic "Foliage Height Diversity" index
# (MacArthur & MacArthur 1961) from the ETH Global Canopy Height
# model: canopy heights are binned into ecologically meaningful
# strata, then a Shannon-Wiener index is computed ACROSS THE STRATA
# (not across land cover classes as with SHDI) -- measures how
# evenly vegetation is distributed vertically across height layers.
#
# Requires: kilimanjaro_canopy_height.tif (from export_canopy_height.js,
# downloaded from Google Drive).
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
canopy_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/DEM/kilimanjaro_canopy_height.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_fhd.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/fhd_results.csv"

buffer_m <- 100  # same buffer as the other per-cell scripts

# Height strata (metres) -- following common vegetation layering
# (ground layer, shrub layer, understory, canopy, emergent)
height_breaks <- c(-Inf, 1, 5, 15, 30, Inf)
height_labels <- c("ground_0_1m", "shrub_1_5m", "subcanopy_5_15m", "canopy_15_30m", "emergent_30m_plus")

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(canopy_path))
canopy <- rast(canopy_path)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), canopy)) pud_grid <- st_transform(pud_grid, crs(canopy))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

# ---------------------------------------------------------------
# 2. PER CELL: HEIGHT DISTRIBUTION -> FOLIAGE HEIGHT DIVERSITY
# ---------------------------------------------------------------
cat(sprintf("Computing Foliage Height Diversity for %d cells...\n", nrow(pud_grid)))

fhd_results <- vector("list", nrow(pud_grid))

for (i in seq_len(nrow(pud_grid))) {
  cell_poly <- vect(pud_grid[i, ])
  cell_buf  <- buffer(cell_poly, buffer_m)

  canopy_crop <- crop(canopy, cell_buf, mask = TRUE)
  vals <- values(canopy_crop)
  vals <- vals[!is.na(vals)]

  if (length(vals) == 0) {
    fhd_results[[i]] <- data.frame(cell_id = i, VerticalStructuralHeterogeneity = NA,
                                    MeanCanopyHeight_m = NA, MaxCanopyHeight_m = NA)
    next
  }

  strata <- cut(vals, breaks = height_breaks, labels = height_labels)
  strata_props <- table(strata) / length(vals)
  p <- strata_props[strata_props > 0]

  fhd <- -sum(p * log(p))  # Shannon-Wiener across height layers

  fhd_results[[i]] <- data.frame(
    cell_id = i,
    VerticalStructuralHeterogeneity = fhd,
    MeanCanopyHeight_m = mean(vals),
    MaxCanopyHeight_m = max(vals)
  )

  if (i %% 20 == 0) cat(sprintf("  Cell %d/%d done\n", i, nrow(pud_grid)))
}

fhd_df <- bind_rows(fhd_results)
cat(sprintf("\n-> %d/%d cells successful\n", sum(!is.na(fhd_df$VerticalStructuralHeterogeneity)), nrow(fhd_df)))
print(head(fhd_df, 10))

# ---------------------------------------------------------------
# 3. JOIN TO pud_grid AND SAVE
# ---------------------------------------------------------------
pud_grid <- left_join(pud_grid, fhd_df, by = "cell_id")

st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nDone. Saved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(st_drop_geometry(pud_grid)[, c("VerticalStructuralHeterogeneity", "MeanCanopyHeight_m", "MaxCanopyHeight_m")]))
