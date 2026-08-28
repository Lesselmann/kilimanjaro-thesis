### ============================================================
### Hyperparameter Tuning (Gamma family): Random Forest, XGBoost, GAM
### Kilimanjaro Thesis
### ============================================================
# Combines the two separate tuning scripts (05_Hyperparameter_Tuning_
# RF_XGBoost.R and 05b_Hyperparameter_Tuning_GAM.R) into one, with the
# PUD model family switched from Tweedie to Gamma (log link) --
# everything else (grid definitions, CV logic, output structure) is
# UNCHANGED from the two source scripts. Concretely, this means:
#
#   1. GAM k-tuning (Part C below): family_fun for PUD changes from
#      tw (Tweedie) to a Gamma(link="log") wrapper. n_images stays on
#      Negative Binomial (nb()) -- unaffected by the PUD family choice.
#   2. XGBoost's PUD objective (Part B) changes from "reg:squarederror"
#      to "reg:gamma", consistent with the rest of the Gamma pipeline
#      (matches 06_FINAL_MultiModel_All4Scenarios.R's PUD XGBoost call).
#   3. predictors_pud_gam is swapped for the Gamma-family backward-AIC
#      result (12 predictors, from backward_selection_ALL.R / used in
#      06_FINAL) instead of the old Tweedie-derived 16-predictor set --
#      this follows directly from the family change, since the reduced
#      predictor set itself depends on which family the backward
#      elimination in script 04/backward_selection_ALL.R was run under.
#      predictors_nimg_gam is untouched (n_images was never on
#      Tweedie/Gamma to begin with).
#   4. Random Forest (Part A) needs NO change -- randomForest() has no
#      family/distribution argument, so its tuning is identical either
#      way. Included here only so all three tunings live in one script.
#
# Nothing else was touched: grid definitions (36 RF combos, 108 XGBoost
# combos, k_grid for GAM), CV fold logic, and output file structure are
# copied verbatim from the two source scripts.
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
pud_grid_gpkg <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
n_images_csv  <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/n_images_per_grid_cell.csv"
master_csv    <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv"
final_predictors_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/screening_final_predictors.csv"
out_dir       <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/"

id_col  <- "cell_id"
k_folds <- 5
predictors_full <- read.csv(final_predictors_csv)$variable

# ---------------------------------------------------------------
# 1. LOAD DATA (both response variables)
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_gpkg, quiet = TRUE)
grid_id_lookup <- data.frame(cell_id = seq_len(nrow(pud_grid)), grid_id = pud_grid$grid_id)

n_images_df <- read.csv(n_images_csv) %>%
  dplyr::select(grid_id, n_images) %>%
  left_join(grid_id_lookup, by = "grid_id")

master_df <- read.csv(master_csv) %>%
  left_join(n_images_df %>% dplyr::select(cell_id, n_images), by = id_col)
master_df$n_images[is.na(master_df$n_images)] <- 0

# Untransformed data -- used for RF/XGBoost (Parts A/B), which always
# retain all 18 predictors regardless of family.
data_pud  <- master_df[, c("avg_annual_PUD", predictors_full)] %>% filter(complete.cases(.))
data_nimg <- master_df[, c("n_images", predictors_full)]        %>% filter(complete.cases(.))

# ---------------------------------------------------------------
# 2. HELPER: run k-fold CV for one hyperparameter configuration
# ---------------------------------------------------------------
compute_rmse <- function(observed, predicted) sqrt(mean((observed - predicted)^2, na.rm = TRUE))

cv_rf <- function(data, dep, vars, ntree, mtry, nodesize, fold_id) {
  preds <- rep(NA, nrow(data))
  for (f in sort(unique(fold_id))) {
    train <- data[fold_id != f, ]; test <- data[fold_id == f, ]
    m <- randomForest(as.formula(paste(dep, "~", paste(vars, collapse = " + "))),
                       data = train, ntree = ntree, mtry = mtry, nodesize = nodesize)
    preds[fold_id == f] <- predict(m, newdata = test)
  }
  compute_rmse(data[[dep]], preds)
}

cv_xgb <- function(data, dep, vars, nrounds, max_depth, learning_rate, subsample, colsample_bytree, objective, fold_id) {
  X <- as.matrix(data[, vars]); y <- data[[dep]]
  preds <- rep(NA, nrow(data))
  for (f in sort(unique(fold_id))) {
    train_idx <- which(fold_id != f); test_idx <- which(fold_id == f)
    m <- xgboost(x = X[train_idx, ], y = y[train_idx], nrounds = nrounds, max_depth = max_depth,
                 learning_rate = learning_rate, subsample = subsample, colsample_bytree = colsample_bytree,
                 objective = objective)
    preds[test_idx] <- predict(m, X[test_idx, ])
  }
  compute_rmse(y, preds)
}

fold_id_pud  <- sample(rep(1:k_folds, length.out = nrow(data_pud)))
fold_id_nimg <- sample(rep(1:k_folds, length.out = nrow(data_nimg)))

# =================================================================
# PART A -- RANDOM FOREST (no family-dependent change)
# =================================================================
rf_grid <- expand.grid(ntree = c(300, 500, 1000), mtry = c(3, 6, 9, 12), nodesize = c(1, 5, 10))
cat(sprintf("Random Forest: testing %d configurations...\n", nrow(rf_grid)))

rf_grid$rmse_PUD <- NA
rf_grid$rmse_nimages <- NA
for (g in seq_len(nrow(rf_grid))) {
  rf_grid$rmse_PUD[g]     <- cv_rf(data_pud, "avg_annual_PUD", predictors_full, rf_grid$ntree[g], rf_grid$mtry[g], rf_grid$nodesize[g], fold_id_pud)
  rf_grid$rmse_nimages[g] <- cv_rf(data_nimg, "n_images", predictors_full, rf_grid$ntree[g], rf_grid$mtry[g], rf_grid$nodesize[g], fold_id_nimg)
  if (g %% 10 == 0) cat(sprintf("  %d/%d tested\n", g, nrow(rf_grid)))
}

rf_best_pud  <- rf_grid[order(rf_grid$rmse_PUD), ][1, ]
rf_best_nimg <- rf_grid[order(rf_grid$rmse_nimages), ][1, ]

cat(sprintf("\nBest RF for PUD:      ntree=%d, mtry=%d, nodesize=%d (RMSE=%.4f)\n",
            rf_best_pud$ntree, rf_best_pud$mtry, rf_best_pud$nodesize, rf_best_pud$rmse_PUD))
cat(sprintf("Best RF for n_images: ntree=%d, mtry=%d, nodesize=%d (RMSE=%.4f)\n",
            rf_best_nimg$ntree, rf_best_nimg$mtry, rf_best_nimg$nodesize, rf_best_nimg$rmse_nimages))

write.csv(rf_grid, paste0(out_dir, "tuning_grid_RandomForest.csv"), row.names = FALSE)

# =================================================================
# PART B -- XGBOOST (PUD objective: reg:squarederror -> reg:gamma)
# =================================================================
xgb_grid <- expand.grid(max_depth = c(2, 3, 4), nrounds = c(30, 50, 100),
                         learning_rate = c(0.03, 0.05, 0.1), subsample = c(0.7, 1.0), colsample_bytree = c(0.7, 1.0))
cat(sprintf("\nXGBoost: testing %d configurations...\n", nrow(xgb_grid)))

xgb_grid$rmse_PUD <- NA
xgb_grid$rmse_nimages <- NA
for (g in seq_len(nrow(xgb_grid))) {
  xgb_grid$rmse_PUD[g] <- cv_xgb(data_pud, "avg_annual_PUD", predictors_full, xgb_grid$nrounds[g], xgb_grid$max_depth[g],
                                  xgb_grid$learning_rate[g], xgb_grid$subsample[g], xgb_grid$colsample_bytree[g],
                                  "reg:gamma", fold_id_pud)   # CHANGED from "reg:squarederror"
  xgb_grid$rmse_nimages[g] <- cv_xgb(data_nimg, "n_images", predictors_full, xgb_grid$nrounds[g], xgb_grid$max_depth[g],
                                      xgb_grid$learning_rate[g], xgb_grid$subsample[g], xgb_grid$colsample_bytree[g],
                                      "count:poisson", fold_id_nimg)
  if (g %% 10 == 0) cat(sprintf("  %d/%d tested\n", g, nrow(xgb_grid)))
}

xgb_best_pud  <- xgb_grid[order(xgb_grid$rmse_PUD), ][1, ]
xgb_best_nimg <- xgb_grid[order(xgb_grid$rmse_nimages), ][1, ]

cat(sprintf("\nBest XGBoost for PUD:      max_depth=%d, nrounds=%d, learning_rate=%.2f, subsample=%.1f, colsample_bytree=%.1f (RMSE=%.4f)\n",
            xgb_best_pud$max_depth, xgb_best_pud$nrounds, xgb_best_pud$learning_rate, xgb_best_pud$subsample, xgb_best_pud$colsample_bytree, xgb_best_pud$rmse_PUD))
cat(sprintf("Best XGBoost for n_images: max_depth=%d, nrounds=%d, learning_rate=%.2f, subsample=%.1f, colsample_bytree=%.1f (RMSE=%.4f)\n",
            xgb_best_nimg$max_depth, xgb_best_nimg$nrounds, xgb_best_nimg$learning_rate, xgb_best_nimg$subsample, xgb_best_nimg$colsample_bytree, xgb_best_nimg$rmse_nimages))

write.csv(xgb_grid, paste0(out_dir, "tuning_grid_XGBoost.csv"), row.names = FALSE)

best_params <- data.frame(
  Response = c("PUD","PUD","n_images","n_images"),
  Model = c("RandomForest","XGBoost","RandomForest","XGBoost"),
  Params = c(
    sprintf("ntree=%d, mtry=%d, nodesize=%d", rf_best_pud$ntree, rf_best_pud$mtry, rf_best_pud$nodesize),
    sprintf("max_depth=%d, nrounds=%d, learning_rate=%.2f, subsample=%.1f, colsample_bytree=%.1f",
            xgb_best_pud$max_depth, xgb_best_pud$nrounds, xgb_best_pud$learning_rate, xgb_best_pud$subsample, xgb_best_pud$colsample_bytree),
    sprintf("ntree=%d, mtry=%d, nodesize=%d", rf_best_nimg$ntree, rf_best_nimg$mtry, rf_best_nimg$nodesize),
    sprintf("max_depth=%d, nrounds=%d, learning_rate=%.2f, subsample=%.1f, colsample_bytree=%.1f",
            xgb_best_nimg$max_depth, xgb_best_nimg$nrounds, xgb_best_nimg$learning_rate, xgb_best_nimg$subsample, xgb_best_nimg$colsample_bytree)
  )
)
write.csv(best_params, paste0(out_dir, "best_hyperparameters_summary.csv"), row.names = FALSE)
cat("\n\n############## BEST RF/XGBOOST HYPERPARAMETERS ##############\n")
print(best_params, row.names = FALSE)
cat(sprintf("\nSaved: %sbest_hyperparameters_summary.csv\n", out_dir))

# =================================================================
# PART C -- GAM BASIS DIMENSION k (PUD family: Tweedie -> Gamma)
# =================================================================
# GAM predictor sets -- n_images unchanged; PUD swapped for the
# Gamma-family backward-AIC result (12 predictors, from
# backward_selection_ALL.R / used in 06_FINAL_MultiModel_All4Scenarios.R)
# instead of the old Tweedie-derived 16-predictor set.
predictors_pud_gam <- c("visible_area_middle_km2", "dist_sd_middle_m", "ForestEcotoneDensity",
                         "ReliefDiversity", "NDSI", "NDVI_var", "visible_naturalness_middle",
                         "visible_water_middle", "VerticalStructuralHeterogeneity", "ContagionIndex",
                         "ColorDiversity", "ColorDiversity_var")

predictors_nimg_gam <- c("visible_area_middle_km2", "ForestEcotoneDensity", "ReliefDiversity",
                          "NDSI", "NDVI_var", "visible_naturalness_middle", "visible_water_middle",
                          "VerticalStructuralHeterogeneity", "LAI_var", "DistanceToWater_km",
                          "ContagionIndex", "ColorDiversity", "ColorDiversity_var")

transform_cols <- c("visible_water_middle", "ForestEcotoneDensity", "LAI_var", "DistanceToWaterway_km")

# k values to test -- unchanged from the original script
k_grid <- c(5, 7, 9, 12, 15)

master_df_gam <- master_df
for (col in transform_cols) {
  if (col %in% names(master_df_gam)) master_df_gam[[col]] <- sqrt(master_df_gam[[col]])
}

data_pud_gam  <- master_df_gam[, c("avg_annual_PUD", predictors_pud_gam)]  %>% filter(complete.cases(.))
data_nimg_gam <- master_df_gam[, c("n_images", predictors_nimg_gam)]        %>% filter(complete.cases(.))

cat(sprintf("\nComplete cases -- PUD (GAM): %d | n_images (GAM): %d\n", nrow(data_pud_gam), nrow(data_nimg_gam)))

cv_gam_k <- function(data, dep, vars, k_val, family_fun, fold_id) {
  formula <- paste(dep, "~", paste0("s(", vars, ", bs='ts', k=", k_val, ")", collapse = " + "))
  preds <- rep(NA, nrow(data))
  for (f in sort(unique(fold_id))) {
    train <- data[fold_id != f, ]; test <- data[fold_id == f, ]
    m <- tryCatch(gam(as.formula(formula), data = train, family = family_fun(), method = "REML", select = TRUE),
                  error = function(e) NULL)
    if (!is.null(m)) preds[fold_id == f] <- predict(m, newdata = test, type = "response")
  }
  compute_rmse(data[[dep]], preds)
}

# Gamma(link="log") wrapped as a zero-arg function, mirroring the tw/nb
# call style used elsewhere (family_fun() -> family object)
gamma_log <- function() Gamma(link = "log")

fold_id_pud_gam  <- sample(rep(1:k_folds, length.out = nrow(data_pud_gam)))
fold_id_nimg_gam <- sample(rep(1:k_folds, length.out = nrow(data_nimg_gam)))

cat("\nTesting k values for GAM (PUD, Gamma)...\n")
k_results_pud <- data.frame(k = k_grid, RMSE_PUD = NA)
for (i in seq_along(k_grid)) {
  k_results_pud$RMSE_PUD[i] <- cv_gam_k(data_pud_gam, "avg_annual_PUD", predictors_pud_gam, k_grid[i], gamma_log, fold_id_pud_gam)
  cat(sprintf("  k=%d: RMSE=%.4f\n", k_grid[i], k_results_pud$RMSE_PUD[i]))
}

cat("\nTesting k values for GAM (n_images, Negative Binomial)...\n")
k_results_nimg <- data.frame(k = k_grid, RMSE_nimages = NA)
for (i in seq_along(k_grid)) {
  k_results_nimg$RMSE_nimages[i] <- cv_gam_k(data_nimg_gam, "n_images", predictors_nimg_gam, k_grid[i], nb, fold_id_nimg_gam)
  cat(sprintf("  k=%d: RMSE=%.4f\n", k_grid[i], k_results_nimg$RMSE_nimages[i]))
}

best_k_pud  <- k_results_pud$k[which.min(k_results_pud$RMSE_PUD)]
best_k_nimg <- k_results_nimg$k[which.min(k_results_nimg$RMSE_nimages)]

cat(sprintf("\nBest k for PUD GAM:      k=%d (RMSE=%.4f)\n", best_k_pud, min(k_results_pud$RMSE_PUD, na.rm=TRUE)))
cat(sprintf("Best k for n_images GAM: k=%d (RMSE=%.4f)\n", best_k_nimg, min(k_results_nimg$RMSE_nimages, na.rm=TRUE)))

write.csv(k_results_pud, paste0(out_dir, "tuning_grid_GAM_k_PUD.csv"), row.names = FALSE)
write.csv(k_results_nimg, paste0(out_dir, "tuning_grid_GAM_k_nimages.csv"), row.names = FALSE)

cat(sprintf("\nSaved: tuning_grid_GAM_k_PUD.csv, tuning_grid_GAM_k_nimages.csv\n"))
cat(sprintf("\nUse s(x, bs='ts', k=%d) for PUD and s(x, bs='ts', k=%d) for n_images in script 06.\n", best_k_pud, best_k_nimg))

cat("\nNOTE: This tunes a SINGLE k applied uniformly to all smooth terms, for\n")
cat("simplicity and computational feasibility. A per-predictor k (e.g. higher k\n")
cat("only for ColorDiversity, as noted earlier in the diagnostic history) would be\n")
cat("more precise but requires manual, term-by-term adjustment rather than grid search.\n")
