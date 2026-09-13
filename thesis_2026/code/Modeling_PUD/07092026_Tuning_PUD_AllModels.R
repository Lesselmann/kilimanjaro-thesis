### ============================================================
### Hyperparameter Tuning -- PUD only, ALL models needing tuning
### Kilimanjaro Thesis
### ============================================================
# Consolidates the tuning of every model that HAS hyperparameters
# to tune, for PUD only (spatial CV, consistent with the final
# comparison script):
#   - Random Forest: ntree, mtry, nodesize        (grid search)
#   - XGBoost:        max_depth, nrounds, learning_rate,
#                      subsample, colsample_bytree  (grid search)
#   - GAM:            k (basis dimension)           (grid search)
#
# NOT tuned here (no traditional hyperparameters):
#   - LM:  no tuning parameters
#   - GLM: no tuning parameters (Gamma family, log link fixed)
#
# PUD uses Gamma (link="log") for GLM/GAM and "reg:gamma" for
# XGBoost, since avg_annual_PUD contains no exact zeros (verified).
### ============================================================

library(dplyr)
library(mgcv)
library(randomForest)
library(xgboost)
library(sf)

set.seed(42)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_gpkg        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
master_csv           <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv"
final_predictors_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/screening_final_predictors.csv"
out_dir              <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/"

id_col       <- "cell_id"
response_var <- "avg_annual_PUD"
k_folds      <- 5

predictors_full <- read.csv(final_predictors_csv)$variable  # loaded dynamically -- currently 17 (corrected screening result)

# GAM-specific reduced predictor list (from script 04's backward
# selection on the corrected 17-predictor set) -- used for GAM k
# tuning, since the FINAL model (script 06, PUD_selected scenario)
# fits GAM on this reduced set, not the full 17. RF/XGBoost keep
# using predictors_full below, since script 06 uses the full set
# for both models in both scenarios.
pud_gam <- c("dist_sd_near_m", "visible_area_middle_km2", "ShapeIndex", "ForestEcotoneDensity",
             "ReliefDiversity", "NDSI", "NDVI_var", "visible_water_near",
             "VerticalStructuralHeterogeneity", "ColorDiversity", "ColorDiversity_var", "visible_NDSI_mean")

# ---------------------------------------------------------------
# 1. LOAD COORDINATES + DATA
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_gpkg, quiet = TRUE)
coords_df <- data.frame(cell_id = seq_len(nrow(pud_grid)), cell_x = pud_grid$cell_x, cell_y = pud_grid$cell_y)

master_df <- read.csv(master_csv) %>% left_join(coords_df, by = id_col)

stopifnot(all(master_df[[response_var]] > 0))
cat("Confirmed: avg_annual_PUD contains no zero/negative values -- Gamma family is valid.\n\n")

data_pud <- master_df[, c(id_col, "cell_x", "cell_y", response_var, predictors_full)] %>% filter(complete.cases(.))
cat(sprintf("Complete cases: %d rows\n\n", nrow(data_pud)))

# ---------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------------
compute_rmse <- function(observed, predicted) sqrt(mean((observed - predicted)^2, na.rm = TRUE))

make_spatial_folds <- function(data, k, seed) {
  set.seed(seed)
  km <- kmeans(data[, c("cell_x", "cell_y")], centers = k, nstart = 1)  # nstart=1: allows seed to vary clustering
  km$cluster
}

spatial_fold <- make_spatial_folds(data_pud, k_folds, seed = 2001)

# ---------------------------------------------------------------
# 3. TUNE RANDOM FOREST
# ---------------------------------------------------------------
cat("############## Tuning Random Forest ##############\n")
rf_grid <- expand.grid(ntree = c(300, 500, 1000), mtry = c(3, 6, 9, 12), nodesize = c(1, 5, 10))
cat(sprintf("Testing %d configurations...\n", nrow(rf_grid)))

rf_grid$CV_RMSE <- NA
for (g in seq_len(nrow(rf_grid))) {
  preds <- rep(NA, nrow(data_pud))
  for (f in sort(unique(spatial_fold))) {
    train <- data_pud[spatial_fold != f, ]; test <- data_pud[spatial_fold == f, ]
    m <- randomForest(as.formula(paste(response_var, "~", paste(predictors_full, collapse = " + "))),
                       data = train, ntree = rf_grid$ntree[g], mtry = rf_grid$mtry[g], nodesize = rf_grid$nodesize[g])
    preds[spatial_fold == f] <- predict(m, newdata = test)
  }
  rf_grid$CV_RMSE[g] <- compute_rmse(data_pud[[response_var]], preds)
  if (g %% 10 == 0) cat(sprintf("  %d/%d tested\n", g, nrow(rf_grid)))
}
rf_best <- rf_grid[order(rf_grid$CV_RMSE), ][1, ]
cat(sprintf("\nBest RF: ntree=%d, mtry=%d, nodesize=%d (CV_RMSE=%.4f)\n\n",
            rf_best$ntree, rf_best$mtry, rf_best$nodesize, rf_best$CV_RMSE))
write.csv(rf_grid, paste0(out_dir, "tuning_grid_RandomForest_PUD.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# 4. TUNE XGBOOST
# ---------------------------------------------------------------
cat("############## Tuning XGBoost ##############\n")
xgb_grid <- expand.grid(max_depth = c(2, 3, 4), nrounds = c(30, 50, 100),
                         learning_rate = c(0.03, 0.05, 0.1), subsample = c(0.7, 1.0), colsample_bytree = c(0.7, 1.0))
cat(sprintf("Testing %d configurations...\n", nrow(xgb_grid)))

X <- as.matrix(data_pud[, predictors_full]); y <- data_pud[[response_var]]
xgb_grid$CV_RMSE <- NA
for (g in seq_len(nrow(xgb_grid))) {
  preds <- rep(NA, nrow(data_pud))
  for (f in sort(unique(spatial_fold))) {
    train_idx <- which(spatial_fold != f); test_idx <- which(spatial_fold == f)
    m <- xgboost(x = X[train_idx, ], y = y[train_idx], nrounds = xgb_grid$nrounds[g], max_depth = xgb_grid$max_depth[g],
                 learning_rate = xgb_grid$learning_rate[g], subsample = xgb_grid$subsample[g],
                 colsample_bytree = xgb_grid$colsample_bytree[g], objective = "reg:gamma")
    preds[test_idx] <- predict(m, X[test_idx, ])
  }
  xgb_grid$CV_RMSE[g] <- compute_rmse(y, preds)
  if (g %% 10 == 0) cat(sprintf("  %d/%d tested\n", g, nrow(xgb_grid)))
}
xgb_best <- xgb_grid[order(xgb_grid$CV_RMSE), ][1, ]
cat(sprintf("\nBest XGBoost: max_depth=%d, nrounds=%d, learning_rate=%.2f, subsample=%.1f, colsample_bytree=%.1f (CV_RMSE=%.4f)\n\n",
            xgb_best$max_depth, xgb_best$nrounds, xgb_best$learning_rate, xgb_best$subsample, xgb_best$colsample_bytree, xgb_best$CV_RMSE))
write.csv(xgb_grid, paste0(out_dir, "tuning_grid_XGBoost_PUD.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# 5. TUNE GAM (basis dimension k) -- uses the REDUCED pud_gam list
#    (14 predictors), matching what the final model actually fits
# ---------------------------------------------------------------
cat("############## Tuning GAM (k) ##############\n")

# Apply the same sqrt-transform as script 06 for GAM's known skewed
# predictors. visible_water_near is NEW (not in the original
# diagnostic) -- check its skewness separately if needed; not
# transformed here by default (see note in script 06).
transform_cols <- c("visible_water_middle", "ForestEcotoneDensity", "LAI_var",
                     "DistanceToWaterway_km", "visible_water_near")
master_df_gam <- master_df
for (col in transform_cols) {
  if (col %in% names(master_df_gam)) master_df_gam[[col]] <- sqrt(master_df_gam[[col]])
}
data_pud_gam <- master_df_gam[, c(id_col, "cell_x", "cell_y", response_var, pud_gam)] %>% filter(complete.cases(.))
spatial_fold_gam <- make_spatial_folds(data_pud_gam, k_folds, seed = 2001)  # same seed, own data subset

k_grid <- c(3, 5, 7, 9, 12, 15)
gam_k_results <- data.frame(k = k_grid, CV_RMSE = NA)

for (i in seq_along(k_grid)) {
  # visible_water_near: linear term, not smooth (see script 06 notes -- 80% zeros)
  smooth_vars <- setdiff(pud_gam, "visible_water_near")
  smooth_part <- paste0("s(", smooth_vars, ", bs='ts', k=", k_grid[i], ")", collapse = " + ")
  linear_part <- if ("visible_water_near" %in% pud_gam) "visible_water_near" else NULL
  formula <- paste(response_var, "~", paste(c(smooth_part, linear_part), collapse = " + "))
  preds <- rep(NA, nrow(data_pud_gam))
  for (f in sort(unique(spatial_fold_gam))) {
    train <- data_pud_gam[spatial_fold_gam != f, ]; test <- data_pud_gam[spatial_fold_gam == f, ]
    m <- tryCatch(gam(as.formula(formula), data = train, family = Gamma(link = "log"), method = "REML", select = TRUE),
                  error = function(e) NULL)
    if (!is.null(m)) preds[spatial_fold_gam == f] <- predict(m, newdata = test, type = "response")
  }
  gam_k_results$CV_RMSE[i] <- compute_rmse(data_pud_gam[[response_var]], preds)
  cat(sprintf("  k=%d: CV_RMSE=%.4f\n", k_grid[i], gam_k_results$CV_RMSE[i]))
}
best_k <- gam_k_results$k[which.min(gam_k_results$CV_RMSE)]
cat(sprintf("\nBest GAM k: %d (CV_RMSE=%.4f)\n\n", best_k, min(gam_k_results$CV_RMSE, na.rm = TRUE)))
write.csv(gam_k_results, paste0(out_dir, "tuning_grid_GAM_k_PUD.csv"), row.names = FALSE)

# ---------------------------------------------------------------
# 6. SUMMARY -- ready to paste into the final model script
# ---------------------------------------------------------------
cat("\n\n############## SUMMARY: BEST HYPERPARAMETERS (PUD) ##############\n")
cat(sprintf("Random Forest: ntree=%d, mtry=%d, nodesize=%d\n", rf_best$ntree, rf_best$mtry, rf_best$nodesize))
cat(sprintf("XGBoost:       max_depth=%d, nrounds=%d, learning_rate=%.2f, subsample=%.1f, colsample_bytree=%.1f\n",
            xgb_best$max_depth, xgb_best$nrounds, xgb_best$learning_rate, xgb_best$subsample, xgb_best$colsample_bytree))
cat(sprintf("GAM:           k=%d\n", best_k))
cat("\nLM and GLM have no hyperparameters to tune.\n")
