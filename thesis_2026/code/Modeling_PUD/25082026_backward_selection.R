### ============================================================
### Recursive (Backward) AIC-Based Feature Selection -- ALL
### PUD (LM/GLM/GAM) and n_images (LM/GLM_NB/GAM_NB)
### Kilimanjaro Thesis -- following course Unit 07 methodology
### ============================================================
# Combines the three separate backward selection scripts into one:
#  - PUD:      LM, GLM (Gamma, log link), GAM (Gamma, log link, select=FALSE)
#  - n_images: LM, GLM (Negative Binomial), GAM (Neg. Binomial, select=FALSE)
#
# GAM-specific note: select=FALSE is used deliberately during backward
# elimination (for BOTH response variables), so manual elimination does
# not interact with GAM's own select=TRUE penalization -- select=TRUE
# is only used in the FINAL reported model fit (script 05).
#
# The 4 heavily right-skewed predictors (visible_water_middle,
# ForestEcotoneDensity, LAI_var, DistanceToWaterway_km) are
# sqrt-transformed ONLY for the GAM runs (both response variables),
# after being diagnosed as the cause of unstable log-link predictions
# during cross-validation.
#
# PUD family note: avg_annual_PUD contains no zero values (verified),
# so the Gamma family (log link) is used directly for GLM and GAM --
# no offset or hurdle approach needed.
### ============================================================

library(dplyr)
library(mgcv)
library(MASS)
library(sf)

set.seed(42)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
pud_grid_gpkg        <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/pud_grid_aesthetic.gpkg"
n_images_csv         <- "C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/n_images_per_grid_cell.csv"
master_csv           <- "C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv"
final_predictors_csv <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/screening_final_predictors.csv"
out_dir              <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results/"

id_col <- "cell_id"
predictor_vars <- read.csv(final_predictors_csv)$variable

# ---------------------------------------------------------------
# 1. LOAD DATA (both response variables, merged)
# ---------------------------------------------------------------
pud_grid <- st_read(pud_grid_gpkg, quiet = TRUE)
grid_id_lookup <- data.frame(cell_id = seq_len(nrow(pud_grid)), grid_id = pud_grid$grid_id)

n_images_df <- read.csv(n_images_csv) %>%
  dplyr::select(grid_id, n_images) %>%
  left_join(grid_id_lookup, by = "grid_id")

master_df <- read.csv(master_csv) %>%
  left_join(n_images_df %>% dplyr::select(cell_id, n_images), by = id_col)
master_df$n_images[is.na(master_df$n_images)] <- 0

# GAM-specific transformed version (4 skewed predictors, sqrt)
transform_cols <- c("visible_water_middle", "ForestEcotoneDensity", "LAI_var", "DistanceToWaterway_km")
master_df_gam <- master_df
for (col in transform_cols) {
  if (col %in% names(master_df_gam)) {
    cat(sprintf("Applying sqrt() transform to: %s (GAM runs only)\n", col))
    master_df_gam[[col]] <- sqrt(master_df_gam[[col]])
  }
}

model_df_pud       <- master_df[, c("avg_annual_PUD", predictor_vars)]     %>% filter(complete.cases(.))
model_df_pud_gam   <- master_df_gam[, c("avg_annual_PUD", predictor_vars)] %>% filter(complete.cases(.))
model_df_nimg      <- master_df[, c("n_images", predictor_vars)]           %>% filter(complete.cases(.))
model_df_nimg_gam  <- master_df_gam[, c("n_images", predictor_vars)]       %>% filter(complete.cases(.))

cat(sprintf("Complete cases -- PUD: %d | PUD (GAM): %d | n_images: %d | n_images (GAM): %d\n",
            nrow(model_df_pud), nrow(model_df_pud_gam), nrow(model_df_nimg), nrow(model_df_nimg_gam)))

# Sanity check: Gamma family requires strictly positive response values
stopifnot(all(model_df_pud$avg_annual_PUD > 0))
stopifnot(all(model_df_pud_gam$avg_annual_PUD > 0))
cat("Confirmed: avg_annual_PUD contains no zero/negative values -- Gamma family is valid.\n\n")

# ---------------------------------------------------------------
# 2. BACKWARD SELECTION FUNCTIONS -- ONE PER MODEL TYPE/FAMILY
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
    mod_aic <- tryCatch(AIC(mod), error = function(e) NA)
    data.frame(Variable = v, AIC = round(mod_aic, 4))
  })
  rec_fs <- do.call("rbind", rec_fs)
  rec_fs$Diff <- rec_fs$AIC - rec_fs$AIC[1]
  print(rec_fs)
  list(table = rec_fs, min_diff = min(rec_fs$Diff[-1]))
}

recursive_feature_selection_glm_nb <- function(data, dep, vars) {
  rec_fs <- lapply(seq(0, length(vars)), function(v) {
    if (v == 0) { formula <- paste(dep, " ~ ", paste(vars, collapse = " + ")); v <- "all" }
    else { formula <- paste(dep, " ~ ", paste(vars[-v], collapse = " + ")); v <- vars[v] }
    mod <- tryCatch(MASS::glm.nb(as.formula(formula), data = data),
                     error = function(e) glm(as.formula(formula), data = data, family = poisson()))
    mod_aic <- tryCatch(AIC(mod), error = function(e) NA)
    data.frame(Variable = v, AIC = round(mod_aic, 4))
  })
  rec_fs <- do.call("rbind", rec_fs)
  rec_fs$Diff <- rec_fs$AIC - rec_fs$AIC[1]
  print(rec_fs)
  list(table = rec_fs, min_diff = min(rec_fs$Diff[-1]))
}

recursive_feature_selection_gam <- function(data, dep, vars, family_fun) {
  rec_fs <- lapply(seq(0, length(vars)), function(v) {
    if (v == 0) { smooth_terms <- paste0("s(", vars, ", bs='ts')", collapse = " + "); v <- "all" }
    else { smooth_terms <- paste0("s(", vars[-v], ", bs='ts')", collapse = " + "); v <- vars[v] }
    formula <- paste(dep, " ~ ", smooth_terms)
    mod <- gam(as.formula(formula), data = data, family = family_fun(), method = "REML", select = FALSE)
    mod_aic <- tryCatch(AIC(mod), error = function(e) NA)
    data.frame(Variable = v, AIC = round(mod_aic, 4))
  })
  rec_fs <- do.call("rbind", rec_fs)
  rec_fs$Diff <- rec_fs$AIC - rec_fs$AIC[1]
  print(rec_fs)
  list(table = rec_fs, min_diff = min(rec_fs$Diff[-1]))
}

run_backward_selection <- function(fit_fun, data, dep, vars, ...) {
  round_num <- 1
  repeat {
    cat(sprintf("\n=== Round %d (%d variables) ===\n", round_num, length(vars)))
    result <- fit_fun(data, dep, vars, ...)
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
# 3. RUN ALL 6 BACKWARD SELECTIONS
# ---------------------------------------------------------------
cat("\n\n#################### PUD: LM ####################\n")
pud_lm <- run_backward_selection(recursive_feature_selection_lm, model_df_pud, "avg_annual_PUD", predictor_vars)
cat(sprintf("\nFinal PUD-LM predictors (%d): %s\n", length(pud_lm), paste(pud_lm, collapse=", ")))

cat("\n\n#################### PUD: GLM (Gamma, log link) ####################\n")
pud_glm <- run_backward_selection(recursive_feature_selection_glm_gamma, model_df_pud, "avg_annual_PUD", predictor_vars)
cat(sprintf("\nFinal PUD-GLM predictors (%d): %s\n", length(pud_glm), paste(pud_glm, collapse=", ")))

cat("\n\n#################### PUD: GAM (Gamma, log link, transformed) ####################\n")
pud_gam <- run_backward_selection(recursive_feature_selection_gam, model_df_pud_gam, "avg_annual_PUD", predictor_vars,
                                   family_fun = function() Gamma(link = "log"))
cat(sprintf("\nFinal PUD-GAM predictors (%d): %s\n", length(pud_gam), paste(pud_gam, collapse=", ")))

cat("\n\n#################### n_images: LM ####################\n")
nimg_lm <- run_backward_selection(recursive_feature_selection_lm, model_df_nimg, "n_images", predictor_vars)
cat(sprintf("\nFinal n_images-LM predictors (%d): %s\n", length(nimg_lm), paste(nimg_lm, collapse=", ")))

cat("\n\n#################### n_images: GLM (Neg. Binomial) ####################\n")
nimg_glm <- run_backward_selection(recursive_feature_selection_glm_nb, model_df_nimg, "n_images", predictor_vars)
cat(sprintf("\nFinal n_images-GLM predictors (%d): %s\n", length(nimg_glm), paste(nimg_glm, collapse=", ")))

cat("\n\n#################### n_images: GAM (Neg. Binomial, transformed) ####################\n")
nimg_gam <- run_backward_selection(recursive_feature_selection_gam, model_df_nimg_gam, "n_images", predictor_vars, family_fun = nb)
cat(sprintf("\nFinal n_images-GAM predictors (%d): %s\n", length(nimg_gam), paste(nimg_gam, collapse=", ")))

# ---------------------------------------------------------------
# 4. SAVE ALL 6 RESULTS
# ---------------------------------------------------------------
selection_summary <- data.frame(
  Response = c("PUD","PUD","PUD","n_images","n_images","n_images"),
  Model = c("LM","GLM","GAM","LM","GLM_NB","GAM_NB"),
  N_Predictors = c(length(pud_lm), length(pud_glm), length(pud_gam), length(nimg_lm), length(nimg_glm), length(nimg_gam)),
  Predictors = c(paste(pud_lm, collapse="; "), paste(pud_glm, collapse="; "), paste(pud_gam, collapse="; "),
                 paste(nimg_lm, collapse="; "), paste(nimg_glm, collapse="; "), paste(nimg_gam, collapse="; "))
)

cat("\n\n############## SUMMARY: ALL 6 BACKWARD SELECTIONS ##############\n")
print(selection_summary, row.names = FALSE)

write.csv(selection_summary, paste0(out_dir, "backward_selection_summary_ALL.csv"), row.names = FALSE)
cat(sprintf("\nSaved: %sbackward_selection_summary_ALL.csv\n", out_dir))
