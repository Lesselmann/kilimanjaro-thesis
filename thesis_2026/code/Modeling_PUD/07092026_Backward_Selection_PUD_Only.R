### ============================================================
### Backward Feature Selection -- PUD ONLY (LM, GLM, GAM)
### Kilimanjaro Thesis
### ============================================================
# AIC-based backward elimination (no CV -- see earlier discussion:
# AIC is an in-sample criterion by design and does not require
# train/test splitting), following the course methodology
# (Zeuss, 2026) and Burnham & Anderson (2002) for AIC-based model
# selection.
#
# Starts from the CORRECTED 17-predictor screening result
# (screening_final_predictors.csv), NOT the earlier 18-predictor
# list used in prior pipeline versions -- differences:
#   removed: visible_naturalness_middle, visible_water_middle, DistanceToWater_km
#   added:   visible_water_near, visible_NDSI_mean
#
# PUD family: GLM and GAM use Gamma (link="log") instead of Tweedie,
# since avg_annual_PUD contains no exact zero values (verified).
# GAM uses select=FALSE during elimination (to avoid interference
# with its own select=TRUE penalization, which is only used in the
# final reported model fit).
### ============================================================

library(dplyr)
library(mgcv)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
master_csv           <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv"
final_predictors_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/screening_final_predictors.csv"
out_dir              <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/"

response_var <- "avg_annual_PUD"

# ---------------------------------------------------------------
# 1. LOAD DATA (dynamically reads whichever predictor list is at
#    final_predictors_csv -- currently the corrected 17-predictor set)
# ---------------------------------------------------------------
master_df <- read.csv(master_csv)
predictor_vars <- read.csv(final_predictors_csv)$variable

cat(sprintf("Loaded %d predictors from screening result:\n", length(predictor_vars)))
print(predictor_vars)

# ---------------------------------------------------------------
# 1b. OPTIONAL: outlier removal toggle (set to TRUE to test the
#     effect of removing extreme avg_annual_PUD values). This
#     script does NOT depend on the GeoPackage/coordinates at all
#     (AIC-based backward selection uses no spatial folds), so
#     filtering rows here is fully independent of any gpkg-based
#     script (05/06) -- no cell_x/cell_y mismatch risk here.
# ---------------------------------------------------------------
REMOVE_OUTLIERS <- FALSE   # <-- toggle: TRUE = remove outliers, FALSE = keep all

if (REMOVE_OUTLIERS) {
  q1 <- quantile(master_df[[response_var]], 0.25, na.rm = TRUE)
  q3 <- quantile(master_df[[response_var]], 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr

  n_before <- nrow(master_df)
  master_df <- master_df[master_df[[response_var]] >= lower & master_df[[response_var]] <= upper, ]
  cat(sprintf("\nOutlier removal ON: %d -> %d rows (removed %d, bounds: %.3f-%.3f)\n",
              n_before, nrow(master_df), n_before - nrow(master_df), lower, upper))
} else {
  cat(sprintf("\nOutlier removal OFF: using all %d rows.\n", nrow(master_df)))
}

stopifnot(all(master_df[[response_var]] > 0))
cat("Confirmed: avg_annual_PUD contains no zero/negative values -- Gamma family is valid.\n")

model_df <- master_df[, c(response_var, predictor_vars)]
model_df <- model_df[complete.cases(model_df), ]
cat(sprintf("\nComplete cases: %d rows\n\n", nrow(model_df)))

# ---------------------------------------------------------------
# 2. BACKWARD SELECTION FUNCTIONS
# ---------------------------------------------------------------
recursive_feature_selection_lm <- function(data, dep, vars) {
  rec_fs <- lapply(seq(0, length(vars)), function(v) {
    if (v == 0) { formula <- paste(dep, " ~ ", paste(vars, collapse = " + ")); v <- "all" }
    else { formula <- paste(dep, " ~ ", paste(vars[-v], collapse = " + ")); v <- vars[v] }
    mod <- lm(as.formula(formula), data = data)
    data.frame(Variable = v, AIC = round(AIC(mod), 4))
  })
  rec_fs <- do.call("rbind", rec_fs)
  rec_fs$Diff <- rec_fs$AIC - rec_fs$AIC[1]
  print(rec_fs)
  list(table = rec_fs, min_diff = min(rec_fs$Diff[-1]))
}

recursive_feature_selection_glm_gamma <- function(data, dep, vars) {
  rec_fs <- lapply(seq(0, length(vars)), function(v) {
    if (v == 0) { formula <- paste(dep, " ~ ", paste(vars, collapse = " + ")); v <- "all" }
    else { formula <- paste(dep, " ~ ", paste(vars[-v], collapse = " + ")); v <- vars[v] }
    mod <- glm(as.formula(formula), data = data, family = Gamma(link = "log"))
    data.frame(Variable = v, AIC = round(AIC(mod), 4))
  })
  rec_fs <- do.call("rbind", rec_fs)
  rec_fs$Diff <- rec_fs$AIC - rec_fs$AIC[1]
  print(rec_fs)
  list(table = rec_fs, min_diff = min(rec_fs$Diff[-1]))
}

recursive_feature_selection_gam_gamma <- function(data, dep, vars) {
  rec_fs <- lapply(seq(0, length(vars)), function(v) {
    if (v == 0) { smooth_terms <- paste0("s(", vars, ", bs='ts')", collapse = " + "); v <- "all" }
    else { smooth_terms <- paste0("s(", vars[-v], ", bs='ts')", collapse = " + "); v <- vars[v] }
    formula <- paste(dep, " ~ ", smooth_terms)
    mod <- gam(as.formula(formula), data = data, family = Gamma(link = "log"), method = "REML", select = FALSE)
    dev_expl <- round(summary(mod)$dev.expl * 100, 2)
    mod_aic <- tryCatch(AIC(mod), error = function(e) NA)
    data.frame(Variable = v, Dev_Expl_pct = dev_expl, AIC = round(mod_aic, 4))
  })
  rec_fs <- do.call("rbind", rec_fs)
  rec_fs$Diff <- rec_fs$AIC - rec_fs$AIC[1]
  print(rec_fs)
  list(table = rec_fs, min_diff = min(rec_fs$Diff[-1]))
}

run_backward_selection <- function(fit_fun, data, dep, vars) {
  round_num <- 1
  repeat {
    cat(sprintf("\n=== Round %d (%d variables) ===\n", round_num, length(vars)))
    result <- fit_fun(data, dep, vars)
    if (any(is.na(result$table$AIC))) { cat("\nWARNING: NA AIC -- stopping.\n"); break }
    if (result$min_diff >= 0) { cat("\nNo further AIC improvement -- stopping.\n"); break }
    exclude <- result$table$Variable[which(result$table$Diff == result$min_diff)]
    cat(sprintf("Removing: %s (AIC change: %.4f)\n", exclude, result$min_diff))
    vars <- vars[vars != exclude]
    round_num <- round_num + 1
    if (length(vars) <= 1) { cat("\nOnly one variable left -- stopping.\n"); break }
  }
  vars
}

# ---------------------------------------------------------------
# 3. RUN FOR LM, GLM (Gamma), GAM (Gamma)
# ---------------------------------------------------------------
cat("\n\n#################### PUD: LM ####################\n")
pud_lm <- run_backward_selection(recursive_feature_selection_lm, model_df, response_var, predictor_vars)
cat(sprintf("\nFinal LM predictors (%d): %s\n", length(pud_lm), paste(pud_lm, collapse=", ")))

cat("\n\n#################### PUD: GLM (Gamma) ####################\n")
pud_glm <- run_backward_selection(recursive_feature_selection_glm_gamma, model_df, response_var, predictor_vars)
cat(sprintf("\nFinal GLM predictors (%d): %s\n", length(pud_glm), paste(pud_glm, collapse=", ")))

cat("\n\n#################### PUD: GAM (Gamma) ####################\n")
pud_gam <- run_backward_selection(recursive_feature_selection_gam_gamma, model_df, response_var, predictor_vars)
cat(sprintf("\nFinal GAM predictors (%d): %s\n", length(pud_gam), paste(pud_gam, collapse=", ")))

# ---------------------------------------------------------------
# 4. SAVE
# ---------------------------------------------------------------
selection_summary <- data.frame(
  Model = c("LM", "GLM_Gamma", "GAM_Gamma"),
  N_Predictors = c(length(pud_lm), length(pud_glm), length(pud_gam)),
  Predictors = c(paste(pud_lm, collapse="; "), paste(pud_glm, collapse="; "), paste(pud_gam, collapse="; "))
)
cat("\n\n############## SUMMARY ##############\n")
print(selection_summary, row.names = FALSE)

out_suffix <- ifelse(REMOVE_OUTLIERS, "_noOutliers", "_allData")
out_file <- paste0(out_dir, "backward_selection_PUD_summary", out_suffix, ".csv")
write.csv(selection_summary, out_file, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_file))
