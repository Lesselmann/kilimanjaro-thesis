### ============================================================
### Merging All Indicator Tables
### Kilimanjaro Thesis
### ============================================================
# Goal: Merge the viewshed, landscape metrics, and Sentinel-2
# results (previously separate CSVs, all with the same 228 cells
# but different indicator columns) into ONE master table -- the
# basis for the subsequent collinearity screening and GAM modelling.
### ============================================================

library(dplyr)

# ---------------------------------------------------------------
# 0. PATHS (please adjust)
# ---------------------------------------------------------------
viewshed_csv        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/viewshed_metrics.csv"
landscape_csv       <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/landscapemetrics_results.csv"
sentinel_csv        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/sentinel_results.csv"
visible_csv         <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/visible_composition_results.csv"
svf_csv             <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/svf_results.csv"
summit_csv          <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/summit_distance_results.csv"
refuge_csv          <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/refuge_results.csv"
fhd_csv             <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/fhd_results.csv"
lai_csv              <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/lai_results.csv"
waterdist_csv       <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/waterdist_results.csv"
contagion_csv       <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/contagion_results.csv"
colordiv_csv        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/colordiv_results.csv"
sloperelief_csv     <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/slope_relief_results.csv"
visiblecont_csv     <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/visible_continuous_results.csv"

out_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table.csv"

# ---------------------------------------------------------------
# 1. LOAD ALL FOURTEEN TABLES
# ---------------------------------------------------------------
stopifnot(file.exists(viewshed_csv), file.exists(landscape_csv), file.exists(sentinel_csv),
          file.exists(visible_csv), file.exists(svf_csv), file.exists(summit_csv),
          file.exists(refuge_csv), file.exists(fhd_csv), file.exists(lai_csv),
          file.exists(waterdist_csv), file.exists(contagion_csv), file.exists(colordiv_csv),
          file.exists(sloperelief_csv), file.exists(visiblecont_csv))

viewshed_df     <- read.csv(viewshed_csv)
landscape_df    <- read.csv(landscape_csv)
sentinel_df     <- read.csv(sentinel_csv)
visible_df      <- read.csv(visible_csv)
svf_df          <- read.csv(svf_csv)
summit_df       <- read.csv(summit_csv)
refuge_df       <- read.csv(refuge_csv)
fhd_df          <- read.csv(fhd_csv)
lai_df          <- read.csv(lai_csv)
waterdist_df    <- read.csv(waterdist_csv)
contagion_df    <- read.csv(contagion_csv)
colordiv_df     <- read.csv(colordiv_csv)
sloperelief_df  <- read.csv(sloperelief_csv)
visiblecont_df  <- read.csv(visiblecont_csv)

cat(sprintf("Viewshed:            %d rows, %d columns\n", nrow(viewshed_df), ncol(viewshed_df)))
cat(sprintf("Landscape Metrics:   %d rows, %d columns\n", nrow(landscape_df), ncol(landscape_df)))
cat(sprintf("Sentinel-2:          %d rows, %d columns\n", nrow(sentinel_df), ncol(sentinel_df)))
cat(sprintf("Visible Composition: %d rows, %d columns\n", nrow(visible_df), ncol(visible_df)))
cat(sprintf("Sky View Factor:     %d rows, %d columns\n", nrow(svf_df), ncol(svf_df)))
cat(sprintf("Distance to Summit:  %d rows, %d columns\n", nrow(summit_df), ncol(summit_df)))
cat(sprintf("Refuge Index:        %d rows, %d columns\n", nrow(refuge_df), ncol(refuge_df)))
cat(sprintf("Foliage Height Div.: %d rows, %d columns\n", nrow(fhd_df), ncol(fhd_df)))
cat(sprintf("LAI:                 %d rows, %d columns\n", nrow(lai_df), ncol(lai_df)))
cat(sprintf("Distance to Water:   %d rows, %d columns\n", nrow(waterdist_df), ncol(waterdist_df)))
cat(sprintf("Contagion Index:     %d rows, %d columns\n", nrow(contagion_df), ncol(contagion_df)))
cat(sprintf("Color Diversity:     %d rows, %d columns\n", nrow(colordiv_df), ncol(colordiv_df)))
cat(sprintf("Slope/Relief Energy: %d rows, %d columns\n", nrow(sloperelief_df), ncol(sloperelief_df)))
cat(sprintf("Visible Continuous:  %d rows, %d columns\n\n", nrow(visiblecont_df), ncol(visiblecont_df)))

# ---------------------------------------------------------------
# 2. IDENTIFY SHARED/DUPLICATE COLUMNS
# ---------------------------------------------------------------
base_cols <- Reduce(intersect, list(names(viewshed_df), names(landscape_df), names(sentinel_df),
                                     names(visible_df), names(svf_df), names(summit_df),
                                     names(refuge_df), names(fhd_df), names(lai_df),
                                     names(waterdist_df), names(contagion_df), names(colordiv_df),
                                     names(sloperelief_df), names(visiblecont_df)))
cat("Shared base columns (kept only once):\n")
print(base_cols)
cat("\n")

landscape_only    <- landscape_df    %>% select(cell_id, setdiff(names(landscape_df), base_cols))
sentinel_only     <- sentinel_df     %>% select(cell_id, setdiff(names(sentinel_df), base_cols))
visible_only      <- visible_df      %>% select(cell_id, setdiff(names(visible_df), base_cols))
svf_only          <- svf_df          %>% select(cell_id, setdiff(names(svf_df), base_cols))
summit_only       <- summit_df       %>% select(cell_id, setdiff(names(summit_df), base_cols))
refuge_only       <- refuge_df       %>% select(cell_id, setdiff(names(refuge_df), base_cols))
fhd_only          <- fhd_df          %>% select(cell_id, setdiff(names(fhd_df), base_cols))
lai_only          <- lai_df          %>% select(cell_id, setdiff(names(lai_df), base_cols))
waterdist_only    <- waterdist_df    %>% select(cell_id, setdiff(names(waterdist_df), base_cols))
contagion_only    <- contagion_df    %>% select(cell_id, setdiff(names(contagion_df), base_cols))
colordiv_only     <- colordiv_df     %>% select(cell_id, setdiff(names(colordiv_df), base_cols))
sloperelief_only  <- sloperelief_df  %>% select(cell_id, setdiff(names(sloperelief_df), base_cols))
visiblecont_only  <- visiblecont_df  %>% select(cell_id, setdiff(names(visiblecont_df), base_cols))

# ---------------------------------------------------------------
# 3. MERGE (viewshed_df retains all base/PUD columns)
# ---------------------------------------------------------------
master_df <- viewshed_df %>%
  left_join(landscape_only, by = "cell_id") %>%
  left_join(sentinel_only, by = "cell_id") %>%
  left_join(visible_only, by = "cell_id") %>%
  left_join(svf_only, by = "cell_id") %>%
  left_join(summit_only, by = "cell_id") %>%
  left_join(refuge_only, by = "cell_id") %>%
  left_join(fhd_only, by = "cell_id") %>%
  left_join(lai_only, by = "cell_id") %>%
  left_join(waterdist_only, by = "cell_id") %>%
  left_join(contagion_only, by = "cell_id") %>%
  left_join(colordiv_only, by = "cell_id") %>%
  left_join(sloperelief_only, by = "cell_id") %>%
  left_join(visiblecont_only, by = "cell_id")

cat(sprintf("Master table: %d rows, %d columns\n\n", nrow(master_df), ncol(master_df)))
cat("Columns in the master table:\n")
print(names(master_df))

# ---------------------------------------------------------------
# 4. PLAUSIBILITY CHECK: MISSING VALUES PER COLUMN
# ---------------------------------------------------------------
cat("\n--- NA check per column ---\n")
na_counts <- sapply(master_df, function(x) sum(is.na(x)))
print(na_counts[na_counts > 0])
if (all(na_counts == 0)) cat("No NAs in the entire table.\n")

# ---------------------------------------------------------------
# 5. SAVE
# ---------------------------------------------------------------
write.csv(master_df, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_csv))

# ---------------------------------------------------------------
# 6. OVERVIEW: WHICH INDICATORS ARE NOW AVAILABLE?
# ---------------------------------------------------------------
indicator_cols <- setdiff(names(master_df), c(base_cols))
cat(sprintf("\n--- %d indicator columns now available ---\n", length(indicator_cols)))
print(indicator_cols)
