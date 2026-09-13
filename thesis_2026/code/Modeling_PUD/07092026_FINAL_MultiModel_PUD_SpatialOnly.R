### ============================================================
### FINAL Multi-Model Comparison -- PUD ONLY, SPATIAL CV ONLY
### Selected predictors only (no all-17 baseline scenario)
### Kilimanjaro Thesis
### ============================================================
# Only the "selected" scenario: backward-AIC-reduced predictor
# sets per model (LM/GLM/GAM), full 17-predictor set for RF/XGBoost
# (tree-based models don't need external feature selection). Only
# spatial CV (random CV omitted). 50x repeated spatial CV.
#
# PUD family: GLM and GAM use Gamma (log link) instead of Tweedie,
# and XGBoost uses "reg:gamma" instead of "reg:squarederror" --
# avg_annual_PUD contains no zero values (verified).
#
# RF/XGBoost hyperparameters and GAM k are from the tuning results
# (script 05, spatial 5-fold CV):
#   RF ntree=500/mtry=3/nodesize=1, XGB max_depth=2/nrounds=100/
#   lr=0.05/subsample=0.7/colsample=0.7, GAM k=5
### ============================================================

library(dplyr)
library(tidyr)
library(mgcv)
library(randomForest)
library(xgboost)
library(sf)

set.seed(42)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_gpkg        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
master_csv_screening <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv"
final_predictors_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/screening_final_predictors.csv"
out_dir              <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/"

id_col       <- "cell_id"
response_var <- "avg_annual_PUD"
k_folds      <- 5
n_repeats    <- 50   # reduce to e.g. 20 for a faster preliminary run

# Tuned hyperparameters from script 05 (spatial CV, corrected predictor sets)
rf_ntree <- 1000; rf_mtry <- 3; rf_nodesize <- 1
xgb_nrounds <- 100; xgb_max_depth <- 2; xgb_learning_rate <- 0.05
xgb_subsample <- 0.7; xgb_colsample_bytree <- 1.0

# UPDATE 2: k=3 (the lowest viable value for bs='ts') was added to
# the tuning grid and confirmed as the new best (CV_RMSE=0.306,
# vs. 0.344 at k=5) -- continuing the pattern that lower flexibility
# is consistently more stable at this sample size (~298 per training
# fold). k=7+ remain unstable regardless.
gam_k <- 3

# ---------------------------------------------------------------
# 1. PREDICTOR SETS
# ---------------------------------------------------------------
predictors_full <- read.csv(final_predictors_csv)$variable  # loaded dynamically -- currently 17 (corrected screening result)

pud_lm  <- c("ReliefDiversity", "visible_water_near", "VerticalStructuralHeterogeneity",
             "ColorDiversity", "visible_NDSI_mean")
pud_glm <- c("visible_area_middle_km2", "ShapeIndex", "ReliefDiversity", "NDSI", "NDVI_var",
             "VerticalStructuralHeterogeneity", "ColorDiversity", "visible_NDSI_mean", "visible_ColorDiversity")
pud_gam <- c("dist_sd_near_m", "visible_area_middle_km2", "ShapeIndex", "ForestEcotoneDensity",
             "ReliefDiversity", "NDSI", "NDVI_var", "visible_water_near",
             "VerticalStructuralHeterogeneity", "ColorDiversity", "ColorDiversity_var", "visible_NDSI_mean")

# IMPORTANT: visible_water_near has 80% exact zeros (median=0,
# mean/median ratio = Inf). It was treated as a LINEAR term in GAM
# (not smooth), which improved but did not fully resolve instability.
# The LM run then showed catastrophic spatial-CV instability
# (Test Pseudo-R2 = -17.1, RMSE = 1.43 vs. Train RMSE = 0.28) --
# LM has no link function to bound predictions, unlike GLM/GAM, so
# an unbounded coefficient on this extreme predictor can explode
# under spatial extrapolation. visible_water_near is now added to
# transform_cols and applied to ALL models that use it (LM, GLM,
# GAM), not just GAM as before.
transform_cols <- c("visible_water_middle", "ForestEcotoneDensity", "LAI_var",
                     "DistanceToWaterway_km", "visible_water_near")

# ---------------------------------------------------------------
# 2. LOAD COORDINATES + DATA
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_gpkg, quiet = TRUE)
coords_df <- data.frame(cell_id = seq_len(nrow(pud_grid)), cell_x = pud_grid$cell_x, cell_y = pud_grid$cell_y)

master_df <- read.csv(master_csv_screening) %>% left_join(coords_df, by = id_col)

master_df_transformed <- master_df
for (col in transform_cols) {
  if (col %in% names(master_df_transformed)) master_df_transformed[[col]] <- sqrt(master_df_transformed[[col]])
}

stopifnot(all(master_df[[response_var]] > 0))
cat("Confirmed: avg_annual_PUD contains no zero/negative values -- Gamma family is valid.\n\n")

# ---------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ---------------------------------------------------------------
compute_metrics <- function(observed, predicted) {
  resid <- observed - predicted
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  mae <- mean(abs(resid), na.rm = TRUE)
  ss_res <- sum(resid^2, na.rm = TRUE)
  ss_tot <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
  c(RMSE = rmse, MAE = mae, Pseudo_R2 = 1 - ss_res / ss_tot)
}
make_spatial_folds <- function(data, k, seed) {
  set.seed(seed)
  km <- kmeans(data[, c("cell_x", "cell_y")], centers = k, nstart = 1)
  km$cluster
}
run_cv <- function(data, dep, vars, fold_id, fit_predict_fun, ...) {
  n <- nrow(data)
  test_preds <- rep(NA, n)
  train_rmse_folds <- c(); train_mae_folds <- c(); train_r2_folds <- c()

  for (f in sort(unique(fold_id))) {
    train <- data[fold_id != f, ]; test <- data[fold_id == f, ]
    result <- fit_predict_fun(train, test, dep, vars, ...)   # now returns list(train_pred=..., test_pred=...)
    test_preds[fold_id == f] <- result$test_pred

    train_metrics_fold <- compute_metrics(train[[dep]], result$train_pred)
    train_rmse_folds <- c(train_rmse_folds, train_metrics_fold["RMSE"])
    train_mae_folds  <- c(train_mae_folds,  train_metrics_fold["MAE"])
    train_r2_folds   <- c(train_r2_folds,   train_metrics_fold["Pseudo_R2"])
  }

  test_metrics <- compute_metrics(data[[dep]], test_preds)
  train_metrics <- c(RMSE = mean(train_rmse_folds), MAE = mean(train_mae_folds), Pseudo_R2 = mean(train_r2_folds))
  # train_metrics is averaged ACROSS THE 5 FOLDS' in-sample fits
  # (i.e. how well each fold's model fits the data it was trained on)

  list(test = test_metrics, train = train_metrics)
}

fit_predict_lm <- function(train, test, dep, vars) {
  m <- lm(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train)
  list(train_pred = predict(m, newdata = train), test_pred = predict(m, newdata = test))
}
fit_predict_glm_gamma <- function(train, test, dep, vars) {
  m <- glm(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train, family = Gamma(link = "log"))
  list(train_pred = predict(m, newdata = train, type = "response"), test_pred = predict(m, newdata = test, type = "response"))
}
fit_predict_gam_gamma <- function(train, test, dep, vars, k) {
  # visible_water_near has ~80% exact zeros (confirmed: Mean/Median
  # ratio = Inf, median = 0) -- too few informative (non-zero) values
  # (~73 of 373 cells) for a stable flexible smooth term. It is
  # therefore included as a LINEAR term (no s(), same as in LM/GLM),
  # while all other predictors remain flexible smooth terms.
  linear_terms <- intersect(vars, "visible_water_near")
  smooth_vars  <- setdiff(vars, "visible_water_near")

  smooth_part <- if (length(smooth_vars) > 0) paste0("s(", smooth_vars, ", bs='ts', k=", k, ")", collapse = " + ") else NULL
  linear_part <- if (length(linear_terms) > 0) paste(linear_terms, collapse = " + ") else NULL
  rhs <- paste(c(smooth_part, linear_part), collapse = " + ")

  formula <- paste(dep, "~", rhs)
  m <- gam(as.formula(formula), data = train, family = Gamma(link = "log"), method = "REML", select = TRUE)
  list(train_pred = predict(m, newdata = train, type = "response"), test_pred = predict(m, newdata = test, type = "response"))
}
fit_predict_rf <- function(train, test, dep, vars, ntree, mtry, nodesize) {
  m <- randomForest(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train,
                     ntree = ntree, mtry = mtry, nodesize = nodesize)
  list(train_pred = predict(m, newdata = train), test_pred = predict(m, newdata = test))
}
fit_predict_xgb <- function(train, test, dep, vars, nrounds, max_depth, learning_rate, subsample, colsample_bytree) {
  X_train <- as.matrix(train[, vars]); y_train <- train[[dep]]
  X_test <- as.matrix(test[, vars])
  m <- xgboost(x = X_train, y = y_train, nrounds = nrounds, max_depth = max_depth,
               learning_rate = learning_rate, subsample = subsample, colsample_bytree = colsample_bytree,
               objective = "reg:gamma")
  list(train_pred = predict(m, X_train), test_pred = predict(m, X_test))
}

# ---------------------------------------------------------------
# 4. DEFINE THE 2 PUD SCENARIOS
# ---------------------------------------------------------------
data_pud_lm   <- master_df_transformed[, c(id_col, "cell_x", "cell_y", response_var, pud_lm)]  %>% filter(complete.cases(.))
data_pud_glm  <- master_df_transformed[, c(id_col, "cell_x", "cell_y", response_var, pud_glm)] %>% filter(complete.cases(.))
data_pud_gam  <- master_df_transformed[, c(id_col, "cell_x", "cell_y", response_var, pud_gam)] %>% filter(complete.cases(.))
data_pud_full <- master_df[, c(id_col, "cell_x", "cell_y", response_var, predictors_full)] %>% filter(complete.cases(.))
# predictors_full/data_pud_full are still used below for RF/XGBoost,
# which keep the full predictor set (no backward selection needed
# for tree-based models -- see script 04/05 notes)

scenarios <- list(
  PUD_selected = list(
    models = list(
      LM  = list(fun = fit_predict_lm, args = list(), data = data_pud_lm, vars = pud_lm),
      GLM = list(fun = fit_predict_glm_gamma, args = list(), data = data_pud_glm, vars = pud_glm),
      GAM = list(fun = fit_predict_gam_gamma, args = list(k = gam_k), data = data_pud_gam, vars = pud_gam),
      RandomForest = list(fun = fit_predict_rf, args = list(ntree=rf_ntree, mtry=rf_mtry, nodesize=rf_nodesize), data = data_pud_full, vars = predictors_full),
      XGBoost = list(fun = fit_predict_xgb, args = list(nrounds=xgb_nrounds, max_depth=xgb_max_depth, learning_rate=xgb_learning_rate, subsample=xgb_subsample, colsample_bytree=xgb_colsample_bytree), data = data_pud_full, vars = predictors_full)
    )
  )
)

# ---------------------------------------------------------------
# 5. RUN BOTH SCENARIOS x 50 REPETITIONS x 5 MODELS -- SPATIAL CV ONLY
# ---------------------------------------------------------------
all_results <- data.frame()

for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  cat(sprintf("\n\n########## SCENARIO: %s ##########\n", scenario_name))

  for (rep_i in 1:n_repeats) {
    if (rep_i %% 10 == 0) cat(sprintf("  [%s] repetition %d/%d\n", scenario_name, rep_i, n_repeats))

    for (model_name in names(scenario$models)) {
      fit_info <- scenario$models[[model_name]]
      spatial_fold <- make_spatial_folds(fit_info$data, k_folds, seed = 2000 + rep_i)

      run_result <- tryCatch(do.call(run_cv, c(list(fit_info$data, response_var, fit_info$vars, spatial_fold, fit_info$fun), fit_info$args)),
                              error = function(e) list(test = c(RMSE=NA, MAE=NA, Pseudo_R2=NA), train = c(RMSE=NA, MAE=NA, Pseudo_R2=NA)))

      all_results <- rbind(all_results,
        data.frame(Scenario = scenario_name, Repetition = rep_i, Model = model_name, Set = "Test",
                   RMSE = run_result$test["RMSE"], MAE = run_result$test["MAE"], Pseudo_R2 = run_result$test["Pseudo_R2"]),
        data.frame(Scenario = scenario_name, Repetition = rep_i, Model = model_name, Set = "Train",
                   RMSE = run_result$train["RMSE"], MAE = run_result$train["MAE"], Pseudo_R2 = run_result$train["Pseudo_R2"])
      )
    }
  }
  write.csv(all_results, paste0(out_dir, "PUD_spatialCV_all_runs.csv"), row.names = FALSE)
  cat(sprintf("Scenario %s complete, intermediate results saved.\n", scenario_name))
}

# ---------------------------------------------------------------
# 6. FINAL SUMMARY
# ---------------------------------------------------------------
summary_df <- all_results %>%
  group_by(Scenario, Model, Set) %>%
  summarise(Mean_RMSE = mean(RMSE, na.rm=TRUE), SD_RMSE = sd(RMSE, na.rm=TRUE),
            Mean_MAE = mean(MAE, na.rm=TRUE), SD_MAE = sd(MAE, na.rm=TRUE),
            Mean_Pseudo_R2 = mean(Pseudo_R2, na.rm=TRUE), SD_Pseudo_R2 = sd(Pseudo_R2, na.rm=TRUE), .groups = "drop")

cat("\n\n############## FINAL SUMMARY: PUD, SPATIAL CV, TRAIN vs. TEST ##############\n")
print(as.data.frame(summary_df), row.names = FALSE)

write.csv(summary_df, paste0(out_dir, "PUD_spatialCV_summary.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %sPUD_spatialCV_summary.csv\n", out_dir))

# ---------------------------------------------------------------
# 7. OVERFITTING GAP -- Train R2 minus Test R2, per model
#    (large positive gap = model fits training data much better
#    than unseen data = sign of overfitting)
# ---------------------------------------------------------------
gap_df <- summary_df %>%
  select(Scenario, Model, Set, Mean_Pseudo_R2) %>%
  tidyr::pivot_wider(names_from = Set, values_from = Mean_Pseudo_R2) %>%
  mutate(Overfitting_Gap = Train - Test)

cat("\n\n############## OVERFITTING GAP (Train R2 - Test R2) ##############\n")
print(as.data.frame(gap_df), row.names = FALSE)
write.csv(gap_df, paste0(out_dir, "PUD_spatialCV_overfitting_gap.csv"), row.names = FALSE)
