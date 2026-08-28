### ============================================================
### Photo Count (n_images) PER GRID CELL
### Kilimanjaro Thesis
### ============================================================
# Computes the RAW number of images per 1km2 grid cell -- separate
# from PUD (unique user-days), for comparison purposes (e.g.
# scatterplot PUD vs. n_images, showing the "power user" reduction
# effect from Wood et al. 2013).
#
# Uses TWO input files:
#   - labels_csv: Level II classification results (file_name, pred_label)
#   - metadata_csv: Flickr metadata (Photo_ID, owner, latitude, longitude...)
# These are merged via the numeric ID embedded in file_name
# (e.g. "image_10002.jpg" -> Photo_ID 10002).
### ============================================================

library(dplyr)
library(sf)

# ---------------------------------------------------------------
# 0. PATHS -- adjust to your actual files
# ---------------------------------------------------------------
labels_csv   <- "C:/Users/Lukas/masterthesis/thesis_2026/data/Classification_Individual_CoEco_Landscape/1786350997845_final_class_individual_coeco_landscape_LINEARPROBE.csv"
metadata_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/final_flickr_dataset_with_metadata.csv"
out_csv      <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/n_images_per_grid_cell.csv"

grid_size_m <- 1000
source_crs  <- 4326
target_crs  <- 32737

aesthetic_labels <- c("landscape", "community_ecosystem")  # same as PUD "aesthetic proxy"

# ---------------------------------------------------------------
# 1. LOAD BOTH FILES
# ---------------------------------------------------------------
labels_df   <- read.csv(labels_csv)
metadata_df <- read.csv(metadata_csv)

cat(sprintf("Labels file: %d rows\n", nrow(labels_df)))
cat(sprintf("Metadata file: %d rows\n", nrow(metadata_df)))

# ---------------------------------------------------------------
# 2. EXTRACT Photo_ID FROM file_name (e.g. "image_10002.jpg" -> 10002)
# ---------------------------------------------------------------
labels_df$Photo_ID <- as.integer(gsub("image_|\\.jpg", "", labels_df$file_name))

# ---------------------------------------------------------------
# 3. MERGE LABELS WITH METADATA (coordinates, owner, date)
# ---------------------------------------------------------------
merged_df <- labels_df %>%
  inner_join(metadata_df, by = "Photo_ID")

cat(sprintf("Merged (labels + metadata): %d rows\n", nrow(merged_df)))

# ---------------------------------------------------------------
# 4. GRID ASSIGNMENT (same convention as the PUD calculation script)
# ---------------------------------------------------------------
points_sf   <- st_as_sf(merged_df, coords = c("longitude", "latitude"), crs = source_crs, remove = FALSE)
points_proj <- st_transform(points_sf, crs = target_crs)
coords      <- st_coordinates(points_proj)

merged_df$easting  <- coords[, "X"]
merged_df$northing <- coords[, "Y"]
merged_df$cell_x   <- floor(merged_df$easting / grid_size_m)
merged_df$cell_y   <- floor(merged_df$northing / grid_size_m)
merged_df$grid_id  <- paste0(merged_df$cell_x, "_", merged_df$cell_y)

# ---------------------------------------------------------------
# 5. FILTER TO THE SAME AESTHETIC PROXY CATEGORIES AS PUD
# ---------------------------------------------------------------
df_aesthetic <- merged_df %>% filter(pred_label %in% aesthetic_labels)
cat(sprintf("Aesthetic proxy images (landscape + community_ecosystem): %d\n", nrow(df_aesthetic)))

# ---------------------------------------------------------------
# 6. COUNT RAW IMAGES PER GRID CELL
# ---------------------------------------------------------------
n_images_per_cell <- df_aesthetic %>%
  group_by(grid_id, cell_x, cell_y) %>%
  summarise(n_images = n(), .groups = "drop")

n_users_per_cell <- df_aesthetic %>%
  group_by(grid_id) %>%
  summarise(n_unique_users = n_distinct(owner), .groups = "drop")

result <- n_images_per_cell %>%
  left_join(n_users_per_cell, by = "grid_id")

cat(sprintf("\nGrid cells with at least one aesthetic-proxy image: %d\n", nrow(result)))
print(summary(result[, c("n_images", "n_unique_users")]))

# ---------------------------------------------------------------
# 7. SAVE
# ---------------------------------------------------------------
write.csv(result, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_csv))
cat("Merge this with your PUD table (on 'grid_id') to compare n_images vs. PUD.\n")
