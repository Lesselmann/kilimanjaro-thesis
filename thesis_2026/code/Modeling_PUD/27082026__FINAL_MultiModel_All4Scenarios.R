### ============================================================
### COMBINED: All 4 Scenarios in One Script
### PUD (all18 / selected) x n_images (all18 / selected)
### Kilimanjaro Thesis -- for supervisory meeting
### ============================================================
# Runs all 4 comparisons in sequence, each with 50x repeated
# random + spatial CV, for all 5 models:
#   1. PUD_all18       -- PUD, all 18 predictors, untransformed
#   2. PUD_selected     -- PUD, backward-AIC-reduced sets + GAM transform
#   3. NImages_all18    -- n_images, all 18 predictors, untransformed
#   4. NImages_selected -- n_images, backward-AIC-reduced sets
#
# PUD family: GLM and GAM use Gamma (log link) instead of Tweedie, and
# XGBoost uses "reg:gamma" instead of "reg:squarederror" -- avg_annual_PUD
# contains no zero values (verified), consistent with the rest of the
# pipeline. n_images stays on Negative Binomial / Poisson.
#
# RF/XGBoost hyperparameters and GAM k are response-specific, taken from
# the tuning results (scripts 05 and 05b, blockCV::cv_spatial 5-fold CV):
#   PUD:      RF ntree=500/mtry=3/nodesize=1, XGB max_depth=2/nrounds=100/lr=0.05/subsample=0.7/colsample=0.7, GAM k=5
#   n_images: RF ntree=300/mtry=3/nodesize=1, XGB max_depth=2/nrounds=100/lr=0.10/subsample=1.0/colsample=0.7, GAM k=15
# NOTE: the n_images GAM (k=15, Pseudo_R2 ~ -0.002 in tuning) performed
# essentially no better than predicting the mean -- clearly worse than
# RF/XGBoost for n_images (Pseudo_R2 ~ 0.12). This is kept here for a
# complete/fair 5-model comparison, but the n_images GAM result should be
# flagged and discussed as a likely non-viable model for that response,
# not silently reported alongside the others as if it were comparable.
#
# WARNING: this is a LONG run -- 4 scenarios x 50 repetitions x 5
# models x 2 CV types x 5 folds = 10,000 individual model fits.
# Expect several hours depending on your machine. Consider running
# overnight, or reduce n_repeats below for a faster preliminary check.
### ============================================================

library(dplyr)
library(mgcv)
library(randomForest)
library(xgboost)
library(MASS)
library(sf)

set.seed(42)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_gpkg        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
n_images_csv         <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/n_images_per_grid_cell.csv"
master_csv_screening <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv"
final_predictors_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/screening_final_predictors.csv"
out_dir              <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/"

id_col    <- "cell_id"
k_folds   <- 5
n_repeats <- 50   # reduce to e.g. 20 for a faster preliminary run

# Tuned hyperparameters from scripts 05 / 05b, response-specific
rf_ntree_pud <- 500; rf_mtry_pud <- 3; rf_nodesize_pud <- 1
rf_ntree_nimg <- 300; rf_mtry_nimg <- 3; rf_nodesize_nimg <- 1

xgb_nrounds_pud <- 100; xgb_max_depth_pud <- 2; xgb_learning_rate_pud <- 0.05
xgb_subsample_pud <- 0.7; xgb_colsample_bytree_pud <- 0.7

xgb_nrounds_nimg <- 100; xgb_max_depth_nimg <- 2; xgb_learning_rate_nimg <- 0.10
xgb_subsample_nimg <- 1.0; xgb_colsample_bytree_nimg <- 0.7

gam_k_pud <- 5
gam_k_nimg <- 15

# ---------------------------------------------------------------
# 1. PREDICTOR SETS
# ---------------------------------------------------------------
predictors_full <- read.csv(final_predictors_csv)$variable  # all 18

# PUD-specific reduced sets (backward AIC selection, Gamma family for GLM/GAM)
pud_lm  <- c("visible_area_middle_km2", "ReliefDiversity", "NDSI",
             "visible_naturalness_middle", "VerticalStructuralHeterogeneity",
             "ColorDiversity", "visible_ColorDiversity")
pud_glm <- c("visible_area_middle_km2", "dist_sd_middle_m", "ShapeIndex",
             "ReliefDiversity", "NDVI_var", "visible_naturalness_middle",
             "VerticalStructuralHeterogeneity", "LAI_var", "DistanceToWater_km",
             "ColorDiversity", "visible_ColorDiversity")
pud_gam <- c("visible_area_middle_km2", "dist_sd_middle_m", "ForestEcotoneDensity",
             "ReliefDiversity", "NDSI", "NDVI_var", "visible_naturalness_middle",
             "visible_water_middle", "VerticalStructuralHeterogeneity", "ContagionIndex",
             "ColorDiversity", "ColorDiversity_var")

# n_images-specific reduced sets (backward AIC selection, Negative Binomial)
nimg_lm  <- c("visible_area_middle_km2", "ReliefDiversity", "NDSI",
              "VerticalStructuralHeterogeneity", "LAI_var", "ColorDiversity")
nimg_glm <- c("dist_sd_near_m", "visible_area_middle_km2", "ReliefDiversity", "NDSI",
              "NDVI_var", "VerticalStructuralHeterogeneity", "DistanceToWater_km",
              "ContagionIndex", "ColorDiversity", "visible_ColorDiversity")
nimg_gam <- c("visible_area_middle_km2", "ForestEcotoneDensity", "ReliefDiversity",
              "NDSI", "NDVI_var", "visible_naturalness_middle", "visible_water_middle",
              "VerticalStructuralHeterogeneity", "LAI_var", "DistanceToWater_km",
              "ContagionIndex", "ColorDiversity", "ColorDiversity_var")

transform_cols <- c("visible_water_middle", "ForestEcotoneDensity", "LAI_var", "DistanceToWaterway_km")

# ---------------------------------------------------------------
# 2. LOAD COORDINATES + BASE DATA
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_gpkg, quiet = TRUE)
coords_df <- data.frame(cell_id = seq_len(nrow(pud_grid)), cell_x = pud_grid$cell_x, cell_y = pud_grid$cell_y)
grid_id_lookup <- data.frame(cell_id = seq_len(nrow(pud_grid)), grid_id = pud_grid$grid_id)

n_images_df <- read.csv(n_images_csv) %>%
  dplyr::select(grid_id, n_images) %>%
  left_join(grid_id_lookup, by = "grid_id")

master_df <- read.csv(master_csv_screening) %>%
  left_join(coords_df, by = id_col) %>%
  left_join(n_images_df %>% dplyr::select(cell_id, n_images), by = id_col)
master_df$n_images[is.na(master_df$n_images)] <- 0

master_df_transformed <- master_df
for (col in transform_cols) {
  if (col %in% names(master_df_transformed)) master_df_transformed[[col]] <- sqrt(master_df_transformed[[col]])
}

cat(sprintf("avg_annual_PUD: mean=%.3f | n_images: mean=%.2f, max=%d\n",
            mean(master_df$avg_annual_PUD, na.rm=TRUE), mean(master_df$n_images), max(master_df$n_images)))

stopifnot(all(master_df$avg_annual_PUD > 0))
cat("Confirmed: avg_annual_PUD contains no zero/negative values -- Gamma family is valid.\n")

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
  preds <- rep(NA, n)
  for (f in sort(unique(fold_id))) {
    train <- data[fold_id != f, ]; test <- data[fold_id == f, ]
    preds[fold_id == f] <- fit_predict_fun(train, test, dep, vars, ...)
  }
  compute_metrics(data[[dep]], preds)
}

fit_predict_lm <- function(train, test, dep, vars) {
  m <- lm(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train)
  predict(m, newdata = test)
}
fit_predict_glm_gamma <- function(train, test, dep, vars) {
  m <- glm(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train,
           family = Gamma(link = "log"))
  predict(m, newdata = test, type = "response")
}
fit_predict_gam_gamma <- function(train, test, dep, vars, k) {
  formula <- paste(dep, "~", paste0("s(", vars, ", bs='ts', k=", k, ")", collapse = " + "))
  m <- gam(as.formula(formula), data = train, family = Gamma(link = "log"), method = "REML", select = TRUE)
  predict(m, newdata = test, type = "response")
}
fit_predict_glm_nb <- function(train, test, dep, vars) {
  m <- tryCatch(MASS::glm.nb(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train),
                error = function(e) glm(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train, family = poisson()))
  predict(m, newdata = test, type = "response")
}
fit_predict_gam_nb <- function(train, test, dep, vars, k) {
  formula <- paste(dep, "~", paste0("s(", vars, ", bs='ts', k=", k, ")", collapse = " + "))
  m <- gam(as.formula(formula), data = train, family = nb(), method = "REML", select = TRUE)
  predict(m, newdata = test, type = "response")
}
fit_predict_rf <- function(train, test, dep, vars, ntree, mtry, nodesize) {
  m <- randomForest(as.formula(paste(dep, "~", paste(vars, collapse = " + "))), data = train,
                     ntree = ntree, mtry = mtry, nodesize = nodesize)
  predict(m, newdata = test)
}
fit_predict_xgb <- function(train, test, dep, vars, nrounds, max_depth, learning_rate, subsample, colsample_bytree, objective) {
  X_train <- as.matrix(train[, vars]); y_train <- train[[dep]]
  X_test <- as.matrix(test[, vars])
  m <- xgboost(x = X_train, y = y_train, nrounds = nrounds, max_depth = max_depth,
               learning_rate = learning_rate, subsample = subsample, colsample_bytree = colsample_bytree,
               objective = objective)
  predict(m, X_test)
}

# ---------------------------------------------------------------
# 4. DEFINE ALL 4 SCENARIOS
# ---------------------------------------------------------------
data_pud_full <- master_df[, c(id_col, "cell_x", "cell_y", "avg_annual_PUD", predictors_full)] %>% filter(complete.cases(.))
data_pud_lm   <- master_df[, c(id_col, "cell_x", "cell_y", "avg_annual_PUD", pud_lm)]  %>% filter(complete.cases(.))
data_pud_glm  <- master_df[, c(id_col, "cell_x", "cell_y", "avg_annual_PUD", pud_glm)] %>% filter(complete.cases(.))
data_pud_gam  <- master_df_transformed[, c(id_col, "cell_x", "cell_y", "avg_annual_PUD", pud_gam)] %>% filter(complete.cases(.))

data_nimg_full <- master_df[, c(id_col, "cell_x", "cell_y", "n_images", predictors_full)] %>% filter(complete.cases(.))
data_nimg_lm   <- master_df[, c(id_col, "cell_x", "cell_y", "n_images", nimg_lm)]  %>% filter(complete.cases(.))
data_nimg_glm  <- master_df[, c(id_col, "cell_x", "cell_y", "n_images", nimg_glm)] %>% filter(complete.cases(.))
data_nimg_gam  <- master_df_transformed[, c(id_col, "cell_x", "cell_y", "n_images", nimg_gam)] %>% filter(complete.cases(.))

scenarios <- list(
  PUD_all18 = list(
    dep = "avg_annual_PUD",
    models = list(
      LM  = list(fun = fit_predict_lm, args = list(), data = data_pud_full, vars = predictors_full),
      GLM = list(fun = fit_predict_glm_gamma, args = list(), data = data_pud_full, vars = predictors_full),
      GAM = list(fun = fit_predict_gam_gamma, args = list(k = gam_k_pud), data = data_pud_full, vars = predictors_full),
      RandomForest = list(fun = fit_predict_rf, args = list(ntree=rf_ntree_pud, mtry=rf_mtry_pud, nodesize=rf_nodesize_pud), data = data_pud_full, vars = predictors_full),
      XGBoost = list(fun = fit_predict_xgb, args = list(nrounds=xgb_nrounds_pud, max_depth=xgb_max_depth_pud, learning_rate=xgb_learning_rate_pud, subsample=xgb_subsample_pud, colsample_bytree=xgb_colsample_bytree_pud, objective="reg:gamma"), data = data_pud_full, vars = predictors_full)
    )
  ),
  PUD_selected = list(
    dep = "avg_annual_PUD",
    models = list(
      LM  = list(fun = fit_predict_lm, args = list(), data = data_pud_lm, vars = pud_lm),
      GLM = list(fun = fit_predict_glm_gamma, args = list(), data = data_pud_glm, vars = pud_glm),
      GAM = list(fun = fit_predict_gam_gamma, args = list(k = gam_k_pud), data = data_pud_gam, vars = pud_gam),
      RandomForest = list(fun = fit_predict_rf, args = list(ntree=rf_ntree_pud, mtry=rf_mtry_pud, nodesize=rf_nodesize_pud), data = data_pud_full, vars = predictors_full),
      XGBoost = list(fun = fit_predict_xgb, args = list(nrounds=xgb_nrounds_pud, max_depth=xgb_max_depth_pud, learning_rate=xgb_learning_rate_pud, subsample=xgb_subsample_pud, colsample_bytree=xgb_colsample_bytree_pud, objective="reg:gamma"), data = data_pud_full, vars = predictors_full)
    )
  ),
  NImages_all18 = list(
    dep = "n_images",
    models = list(
      LM  = list(fun = fit_predict_lm, args = list(), data = data_nimg_full, vars = predictors_full),
      GLM_NB = list(fun = fit_predict_glm_nb, args = list(), data = data_nimg_full, vars = predictors_full),
      GAM_NB = list(fun = fit_predict_gam_nb, args = list(k = gam_k_nimg), data = data_nimg_full, vars = predictors_full),
      RandomForest = list(fun = fit_predict_rf, args = list(ntree=rf_ntree_nimg, mtry=rf_mtry_nimg, nodesize=rf_nodesize_nimg), data = data_nimg_full, vars = predictors_full),
      XGBoost = list(fun = fit_predict_xgb, args = list(nrounds=xgb_nrounds_nimg, max_depth=xgb_max_depth_nimg, learning_rate=xgb_learning_rate_nimg, subsample=xgb_subsample_nimg, colsample_bytree=xgb_colsample_bytree_nimg, objective="count:poisson"), data = data_nimg_full, vars = predictors_full)
    )
  ),
  NImages_selected = list(
    dep = "n_images",
    models = list(
      LM  = list(fun = fit_predict_lm, args = list(), data = data_nimg_lm, vars = nimg_lm),
      GLM_NB = list(fun = fit_predict_glm_nb, args = list(), data = data_nimg_glm, vars = nimg_glm),
      GAM_NB = list(fun = fit_predict_gam_nb, args = list(k = gam_k_nimg), data = data_nimg_gam, vars = nimg_gam),
      RandomForest = list(fun = fit_predict_rf, args = list(ntree=rf_ntree_nimg, mtry=rf_mtry_nimg, nodesize=rf_nodesize_nimg), data = data_nimg_full, vars = predictors_full),
      XGBoost = list(fun = fit_predict_xgb, args = list(nrounds=xgb_nrounds_nimg, max_depth=xgb_max_depth_nimg, learning_rate=xgb_learning_rate_nimg, subsample=xgb_subsample_nimg, colsample_bytree=xgb_colsample_bytree_nimg, objective="count:poisson"), data = data_nimg_full, vars = predictors_full)
    )
  )
)

# ---------------------------------------------------------------
# 5. RUN ALL 4 SCENARIOS x 50 REPETITIONS x 5 MODELS x 2 CV TYPES
# ---------------------------------------------------------------
all_results <- data.frame()

for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  cat(sprintf("\n\n########## SCENARIO: %s ##########\n", scenario_name))

  for (rep_i in 1:n_repeats) {
    if (rep_i %% 10 == 0) cat(sprintf("  [%s] repetition %d/%d\n", scenario_name, rep_i, n_repeats))

    for (model_name in names(scenario$models)) {
      fit_info <- scenario$models[[model_name]]
      random_fold  <- {set.seed(1000 + rep_i); sample(rep(1:k_folds, length.out = nrow(fit_info$data)))}
      spatial_fold <- make_spatial_folds(fit_info$data, k_folds, seed = 2000 + rep_i)

      random_metrics  <- tryCatch(do.call(run_cv, c(list(fit_info$data, scenario$dep, fit_info$vars, random_fold, fit_info$fun), fit_info$args)),
                                   error = function(e) c(RMSE=NA, MAE=NA, Pseudo_R2=NA))
      spatial_metrics <- tryCatch(do.call(run_cv, c(list(fit_info$data, scenario$dep, fit_info$vars, spatial_fold, fit_info$fun), fit_info$args)),
                                   error = function(e) c(RMSE=NA, MAE=NA, Pseudo_R2=NA))

      all_results <- rbind(all_results,
        data.frame(Scenario = scenario_name, Repetition = rep_i, Model = model_name, CV_Type = "Random",
                   RMSE = random_metrics["RMSE"], MAE = random_metrics["MAE"], Pseudo_R2 = random_metrics["Pseudo_R2"]),
        data.frame(Scenario = scenario_name, Repetition = rep_i, Model = model_name, CV_Type = "Spatial",
                   RMSE = spatial_metrics["RMSE"], MAE = spatial_metrics["MAE"], Pseudo_R2 = spatial_metrics["Pseudo_R2"])
      )
    }
  }
  # Save intermediate results after each scenario, in case of interruption
  write.csv(all_results, paste0(out_dir, "COMBINED_all4_scenarios_all_runs.csv"), row.names = FALSE)
  cat(sprintf("Scenario %s complete, intermediate results saved.\n", scenario_name))
}

# ---------------------------------------------------------------
# 6. FINAL SUMMARY -- ALL 4 SCENARIOS SIDE BY SIDE
# ---------------------------------------------------------------
summary_df <- all_results %>%
  group_by(Scenario, Model, CV_Type) %>%
  summarise(Mean_RMSE = mean(RMSE, na.rm=TRUE), SD_RMSE = sd(RMSE, na.rm=TRUE),
            Mean_MAE = mean(MAE, na.rm=TRUE), SD_MAE = sd(MAE, na.rm=TRUE),
            Mean_Pseudo_R2 = mean(Pseudo_R2, na.rm=TRUE), SD_Pseudo_R2 = sd(Pseudo_R2, na.rm=TRUE), .groups = "drop")

cat("\n\n############## FINAL SUMMARY: ALL 4 SCENARIOS ##############\n")
print(as.data.frame(summary_df), row.names = FALSE)

write.csv(summary_df, paste0(out_dir, "COMBINED_all4_scenarios_summary.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %sCOMBINED_all4_scenarios_summary.csv\n", out_dir))

