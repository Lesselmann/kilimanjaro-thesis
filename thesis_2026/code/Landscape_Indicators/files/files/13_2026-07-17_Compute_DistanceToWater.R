### ============================================================
### Distance to Waterbody PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Goal: Distance from each grid cell centroid to the nearest water
# pixel (WorldCover class 80) -- Coherence/Imageability indicator
# (Arriaza et al. 2004; White et al. 2010), already cited in Table 3
# but not yet operationalised.
#
# Technical approach: A SINGLE distance transform over the entire
# (cropped) area, instead of 228 individual queries -- much more
# efficient, since terra::distance() is optimised for this.
### ============================================================

library(terra)
library(sf)
library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_path <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
lc_utm_path   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/Landcover/kilimanjaro_worldcover.tif"

out_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/pud_grid_with_waterdist.gpkg"
out_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/waterdist_results.csv"

water_class <- 80
buffer_km   <- 15  # buffer around the grid bounding box, to also capture distant water bodies

# ---------------------------------------------------------------
# 1. LOAD DATA
# ---------------------------------------------------------------
stopifnot(file.exists(lc_utm_path))
lc <- rast(lc_utm_path)

pud_grid <- st_read(pud_grid_path, quiet = TRUE)
if (!same.crs(vect(pud_grid), lc)) pud_grid <- st_transform(pud_grid, crs(lc))
pud_grid$cell_id <- seq_len(nrow(pud_grid))

# ---------------------------------------------------------------
# 2. CROP TO GRID BBOX + BUFFER
# ---------------------------------------------------------------
grid_bbox <- st_bbox(pud_grid)
buffer_m <- buffer_km * 1000
crop_ext <- ext(grid_bbox["xmin"] - buffer_m, grid_bbox["xmax"] + buffer_m,
                 grid_bbox["ymin"] - buffer_m, grid_bbox["ymax"] + buffer_m)
lc_crop <- crop(lc, crop_ext)

cat(sprintf("Cropped raster: %d x %d pixels\n", ncol(lc_crop), nrow(lc_crop)))

# ---------------------------------------------------------------
# 3. WATER MASK + DISTANCE TRANSFORM (ONCE)
# ---------------------------------------------------------------
cat("Creating water mask...\n")
water_mask <- lc_crop == water_class
water_mask[water_mask == 0] <- NA  # only water cells remain as "target"

n_water_pixels <- sum(values(water_mask), na.rm = TRUE)
cat(sprintf("Number of water pixels in the cropped area: %d\n", n_water_pixels))

if (n_water_pixels == 0) {
  stop("No water bodies found in the cropped area -- increase the buffer (buffer_km).")
}

cat("Computing distance transform (may take a few minutes)...\n")
dist_raster <- distance(water_mask)
cat("Distance transform done.\n")

# ---------------------------------------------------------------
# 4. EXTRACT PER CENTROID
# ---------------------------------------------------------------
centroids <- st_centroid(pud_grid)
dist_m <- terra::extract(dist_raster, vect(centroids))[, 2]

pud_grid$DistanceToWater_km <- dist_m / 1000

cat(sprintf("\nDistance range: %.3f to %.2f km\n",
            min(pud_grid$DistanceToWater_km, na.rm = TRUE),
            max(pud_grid$DistanceToWater_km, na.rm = TRUE)))
cat(sprintf("NAs: %d\n", sum(is.na(pud_grid$DistanceToWater_km))))

# ---------------------------------------------------------------
# 5. SAVE
# ---------------------------------------------------------------
st_write(pud_grid, out_gpkg, delete_dsn = TRUE)
write.csv(st_drop_geometry(pud_grid), out_csv, row.names = FALSE)

cat(sprintf("\nSaved: %s, %s\n", out_gpkg, out_csv))

cat("\n--- Summary ---\n")
print(summary(pud_grid$DistanceToWater_km))
