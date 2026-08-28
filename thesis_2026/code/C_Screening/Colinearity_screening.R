#############################################################
# STEP 1 OF THE WORKFLOW: Pairwise Multicollinearity Screening
# (flat, NO concept groups)
#
# This script ONLY performs the pre-fit collinearity screening:
#   Step A: pairwise Spearman screening (|rho| >= 0.6 -> drop)
#   Step B: iterative VIF screening      (VIF >= 5     -> drop)
#
# It does NOT fit the GAM. The GAM fit is a separate script
# (see e.g. gam_fit.R) that reads the exported
# 'screening_final_predictors.csv' produced here as its input.
# This keeps the pre-fit filter (linear redundancy) cleanly
# separated from model fitting and post-fit diagnostics
# (concurvity - see concurvity_diagnostics.R).
#############################################################

library(car)

# ---- 0) Output directory ---------------------------------------------------
out_dir <- "C:/Users/Lukas/masterthesis/thesis_2026/data/screening/results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 1) Load data -----------------------------------------------------------
master <- read.csv("C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators/master_indicator_table_SCREENING_ONLY_4.csv")

response_col <- "avg_annual_PUD"
id_col       <- "cell_id"

df       <- master[, !(names(master) %in% c(response_col, id_col))]
response <- master[[response_col]]

spearman_threshold <- 0.6
vif_threshold       <- 5

# ---- 2) Pairwise Spearman screening (flat, no groups) ------------------------
# Returns both the retained variables AND a full log of every drop decision.

screen_spearman <- function(df, threshold = 0.6) {
  cor_mat <- cor(df, method = "spearman", use = "pairwise.complete.obs")
  diag(cor_mat) <- 0
  keep <- colnames(df)

  drop_log <- data.frame(
    step        = integer(0),
    method      = character(0),
    dropped_var = character(0),
    partner_var = character(0),
    statistic   = numeric(0),
    reason      = character(0),
    stringsAsFactors = FALSE
  )
  step_counter <- 0

  repeat {
    sub_mat <- cor_mat[keep, keep, drop = FALSE]
    max_val <- max(abs(sub_mat))
    if (max_val < threshold || length(keep) <= 1) break

    idx  <- which(abs(sub_mat) == max_val, arr.ind = TRUE)[1, ]
    var1 <- rownames(sub_mat)[idx[1]]
    var2 <- colnames(sub_mat)[idx[2]]

    mean_cor1 <- mean(abs(sub_mat[var1, ]))
    mean_cor2 <- mean(abs(sub_mat[var2, ]))
    drop_var  <- if (mean_cor1 >= mean_cor2) var1 else var2
    keep_var  <- setdiff(c(var1, var2), drop_var)

    step_counter <- step_counter + 1
    message(sprintf("Spearman drop: '%s' (rho = %.2f with '%s')",
                    drop_var, max_val, keep_var))

    drop_log <- rbind(drop_log, data.frame(
      step        = step_counter,
      method      = "spearman",
      dropped_var = drop_var,
      partner_var = keep_var,
      statistic   = max_val,
      reason      = sprintf("|rho| = %.3f >= %.2f threshold; higher mean |rho| across remaining set (%.3f vs %.3f)",
                             max_val, threshold, mean_cor1, mean_cor2),
      stringsAsFactors = FALSE
    ))

    keep <- setdiff(keep, drop_var)
  }
  list(keep = keep, drop_log = drop_log)
}

spearman_result      <- screen_spearman(df, spearman_threshold)
vars_after_spearman  <- spearman_result$keep
spearman_drop_log    <- spearman_result$drop_log

cat("\nRetained after Spearman (n =", length(vars_after_spearman), "):\n")
print(vars_after_spearman)

# ---- 3) Iterative VIF screening ----------------------------------------------
# Returns final variable set, final VIF values, AND a drop log in the same
# format as the Spearman step, so both can be combined into one master log.

screen_vif <- function(df, keep_vars, response, threshold = 5) {
  drop_log <- data.frame(
    step        = integer(0),
    method      = character(0),
    dropped_var = character(0),
    partner_var = character(0),
    statistic   = numeric(0),
    reason      = character(0),
    stringsAsFactors = FALSE
  )
  step_counter <- 0
  vif_vals <- NULL

  repeat {
    form <- as.formula(paste("response ~", paste(keep_vars, collapse = " + ")))
    mod  <- lm(form, data = data.frame(response = response,
                                       df[, keep_vars, drop = FALSE]))
    vif_vals <- car::vif(mod)
    max_vif  <- max(vif_vals)
    if (max_vif < threshold || length(keep_vars) <= 1) break

    drop_var <- names(vif_vals)[which.max(vif_vals)]
    step_counter <- step_counter + 1
    message(sprintf("VIF drop: '%s' (VIF = %.2f)", drop_var, max_vif))

    drop_log <- rbind(drop_log, data.frame(
      step        = step_counter,
      method      = "vif",
      dropped_var = drop_var,
      partner_var = NA_character_,
      statistic   = max_vif,
      reason      = sprintf("VIF = %.2f >= %.1f threshold (highest among remaining %d predictors)",
                             max_vif, threshold, length(keep_vars)),
      stringsAsFactors = FALSE
    ))

    keep_vars <- setdiff(keep_vars, drop_var)
  }
  list(final_vars = keep_vars, final_vif = vif_vals, drop_log = drop_log)
}

vif_result       <- screen_vif(df, vars_after_spearman, response, vif_threshold)
final_predictors <- vif_result$final_vars
vif_drop_log     <- vif_result$drop_log

cat("\nFinal predictor set after screening (n =", length(final_predictors), "):\n")
print(final_predictors)
cat("\nFinal VIF values:\n")
print(sort(vif_result$final_vif, decreasing = TRUE))

# ---- 4) Combine drop histories & EXPORT screening results --------------------

full_drop_log <- rbind(spearman_drop_log, vif_drop_log)
full_drop_log$step <- seq_len(nrow(full_drop_log))  # continuous numbering

final_vif_df <- data.frame(
  variable = names(vif_result$final_vif),
  vif      = as.numeric(vif_result$final_vif)
)
final_vif_df <- final_vif_df[order(-final_vif_df$vif), ]

final_predictors_df <- data.frame(
  rank     = seq_along(final_predictors),
  variable = final_predictors
)

# CSV exports (human-readable, easy to paste into thesis tables/appendix)
write.csv(full_drop_log,       file.path(out_dir, "screening_drop_history.csv"),     row.names = FALSE)
write.csv(final_predictors_df, file.path(out_dir, "screening_final_predictors.csv"), row.names = FALSE)
write.csv(final_vif_df,        file.path(out_dir, "screening_final_vif.csv"),        row.names = FALSE)

# RDS export (full R objects, for later re-analysis / GAM script input)
saveRDS(
  list(
    spearman_threshold  = spearman_threshold,
    vif_threshold        = vif_threshold,
    vars_after_spearman  = vars_after_spearman,
    final_predictors     = final_predictors,
    final_vif            = vif_result$final_vif,
    drop_log             = full_drop_log,
    n_candidates_start   = ncol(df),
    n_after_spearman     = length(vars_after_spearman),
    n_final              = length(final_predictors),
    timestamp            = Sys.time()
  ),
  file.path(out_dir, "screening_results.rds")
)

cat("\nScreening results exported to:\n -", out_dir, "\n")
cat("  - screening_drop_history.csv     (every dropped variable, step, reason)\n")
cat("  - screening_final_predictors.csv (surviving predictor set)\n")
cat("  - screening_final_vif.csv        (final VIF values)\n")
cat("  - screening_results.rds          (full R object -> input for the GAM-fit script)\n")

cat("\nDone. Next step: run the GAM-fit script, which reads\n")
cat("'screening_final_predictors.csv' (or the .rds) as its predictor set.\n")
