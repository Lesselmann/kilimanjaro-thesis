### ============================================================
### Full-Area GAM Prediction + Map Export
### Kilimanjaro Thesis -- FINAL STEP
### ============================================================
# Loads the final trained GAM model and the prepared full-area
# predictor table, generates PUD predictions for the entire study
# area, and exports both a CSV/joined GeoPackage and a raster map
# for visualisation.
#
# FIX 2026-09-28: predict.gam() returns a 1D array. terra treated the
# column as categorical, so the GeoTIFF held factor codes (0..n)
# instead of PUD values. Predictions are now coerced to numeric and
# the raster is checked + written as float.
### ============================================================

library(dplyr)
library(mgcv)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------------
# 0. PATHS
# ---------------------------------------------------------------
dir_model    <- "C:/Users/Lukas/downloads/screening_2/results/"
dir_fullarea <- "C:/Users/Lukas/masterthesis/thesis_2026/data/FullArea_Prediction/"

gam_model_rds   <- paste0(dir_model, "FINAL_model_GAM.rds")   # from script 06
predictors_csv  <- paste0(dir_fullarea, "fullarea_GAM_predictors_FINAL.csv")
grid_gpkg       <- paste0(dir_fullarea, "full_study_area_grid_1km.gpkg")

out_csv      <- paste0(dir_fullarea, "fullarea_PUD_prediction.csv")
out_gpkg     <- paste0(dir_fullarea, "fullarea_PUD_prediction.gpkg")
out_raster   <- paste0(dir_fullarea, "fullarea_PUD_prediction.tif")
out_raster_cl <- paste0(dir_fullarea, "fullarea_PUD_prediction_classes.tif")

# ---------------------------------------------------------------
# 1. LOAD THE FINAL TRAINED MODEL
# ---------------------------------------------------------------
gam_final <- readRDS(gam_model_rds)
cat("Loaded final GAM model.\n")
cat("Model formula:\n")
print(formula(gam_final))

# ---------------------------------------------------------------
# 2. LOAD FULL-AREA PREDICTOR DATA
# ---------------------------------------------------------------
full_data <- read.csv(predictors_csv)
cat(sprintf("\nLoaded %d cells (%d predictable, %d excluded due to missing predictors)\n",
            nrow(full_data), sum(full_data$can_predict), sum(!full_data$can_predict)))

# ---------------------------------------------------------------
# 3. PREDICT -- only for cells with complete predictors
# ---------------------------------------------------------------
predictable <- full_data %>% filter(can_predict)

# as.numeric(): predict.gam() returns a named 1D array -> plain numeric vector
predictable$PUD_predicted <- as.numeric(
  predict(gam_final, newdata = predictable, type = "response")
)

cat("\n--- Prediction summary ---\n")
print(summary(predictable$PUD_predicted))

# Merge predictions back into the full table (NA for excluded cells)
full_data <- full_data %>%
  left_join(predictable %>% select(cell_id, PUD_predicted), by = "cell_id")

write.csv(full_data, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_csv))

# ---------------------------------------------------------------
# 4. JOIN TO GRID GEOMETRY AND SAVE AS GEOPACKAGE
# ---------------------------------------------------------------
grid <- st_read(grid_gpkg, quiet = TRUE)
grid_pred <- grid %>%
  left_join(full_data %>% select(cell_id, PUD_predicted, can_predict), by = "cell_id")

stopifnot(is.numeric(grid_pred$PUD_predicted), !is.array(grid_pred$PUD_predicted))

old_gpkg <- paste0(out_gpkg, c("", "-wal", "-shm"))
unlink(old_gpkg[file.exists(old_gpkg)])
if (file.exists(out_gpkg)) stop("GPKG is locked -- close it in QGIS (or restart R) and run again.")
st_write(grid_pred, out_gpkg, quiet = TRUE)
cat(sprintf("Saved: %s\n", out_gpkg))

# ---------------------------------------------------------------
# 5. RASTERIZE FOR VISUALISATION
# ---------------------------------------------------------------
# remove old raster + sidecar files (an old .aux.xml would re-attach
# the broken category table when the TIF is opened)
old_files <- paste0(c(out_raster, out_raster_cl), rep(c("", ".aux.xml"), each = 2))
unlink(old_files[file.exists(old_files)])

grid_vect <- vect(grid_pred)
grid_vect$PUD_predicted <- as.numeric(grid_vect$PUD_predicted)

template   <- rast(grid_vect, resolution = 1000)
pud_raster <- rasterize(grid_vect, template, field = "PUD_predicted")
names(pud_raster) <- "PUD_predicted"

# sanity checks: must be continuous and match the prediction range
stopifnot(!is.factor(pud_raster))
rng <- global(pud_raster, c("min", "max"), na.rm = TRUE)
cat("\nRaster value range:\n"); print(rng)
stopifnot(abs(rng$min - min(predictable$PUD_predicted)) < 1e-6,
          abs(rng$max - max(predictable$PUD_predicted)) < 1e-6)

writeRaster(pud_raster, out_raster, datatype = "FLT4S", overwrite = TRUE)
cat(sprintf("Saved: %s\n", out_raster))

# ---------------------------------------------------------------
# 6. CLASSIFIED MAP (5 quantile classes)
# ---------------------------------------------------------------
# PUD predictions are strongly right-skewed (skewness ~5, log(PUD)
# ~symmetric) -> quantile classes show spatial structure across the
# whole area; hotspots fall into the top class.
v    <- values(pud_raster, na.rm = TRUE)
brks <- quantile(v, probs = seq(0, 1, 0.2))
cat("\nQuantile class breaks:\n"); print(round(brks, 4))

pud_classes <- classify(pud_raster, rcl = brks, include.lowest = TRUE)
names(pud_classes) <- "PUD_class"
writeRaster(pud_classes, out_raster_cl, datatype = "INT1U", overwrite = TRUE)
cat(sprintf("Saved: %s\n", out_raster_cl))

par(mfrow = c(1, 2))
hist(v, breaks = 50, main = "PUD prediction", xlab = "PUD")
hist(log(v), breaks = 50, main = "log(PUD prediction)", xlab = "log(PUD)")
par(mfrow = c(1, 1))

plot(pud_classes, col = hcl.colors(5, "YlOrRd", rev = TRUE),
     main = "Predicted PUD (quantile classes)")

cat("\n============================================================\n")
cat("DONE. Full-area PUD prediction complete.\n")
cat(sprintf("Predicted cells: %d | Excluded (missing predictors): %d\n",
            sum(!is.na(full_data$PUD_predicted)), sum(is.na(full_data$PUD_predicted))))
cat("============================================================\n")

# ---------------------------------------------------------------
# 7. MODEL DIAGNOSTICS
# ---------------------------------------------------------------
plot(gam_final, pages = 1, shade = TRUE)
gam.check(gam_final)

# ---------------------------------------------------------------
# 8. THESIS FIGURES: PUD vs. PREDICTORS + PARTIAL EFFECTS
# ---------------------------------------------------------------
# Both figures use the training data stored in the model object
# (gam_final$model), i.e. the cells with observed PUD, and share
# one visual style.
dir_fig <- paste0(dir_fullarea, "figures/")
dir.create(dir_fig, showWarnings = FALSE, recursive = TRUE)

# --- 8.0 Style ---------------------------------------------------
col_main  <- "#1F4E79"   # line / smoother
col_band  <- "#9DB9D5"   # confidence band
col_point <- "grey35"

# Readable axis labels -- edit freely (names not listed stay as they are)
var_labels <- c(
  dist_sd_near_m                  = "Distance SD,\nnear zone (m)",
  visible_area_middle_km2         = "Visible area,\nmiddle zone (km\u00b2)",
  ShapeIndex                      = "Shape index",
  ForestEcotoneDensity            = "Forest ecotone density",
  ReliefDiversity                 = "Relief diversity",
  NDSI                            = "NDSI",
  NDVI_var                        = "NDVI variance",
  visible_water_near              = "Visible water,\nnear zone",
  VerticalStructuralHeterogeneity = "Vertical structural\nheterogeneity",
  ColorDiversity                  = "Colour diversity",
  ColorDiversity_var              = "Colour diversity\nvariance",
  visible_NDSI_mean               = "Mean visible NDSI"
)
# Use "\n" for line breaks in long names; "\u00b2" = superscript 2
nice <- function(x) ifelse(x %in% names(var_labels), var_labels[x], x)

theme_thesis <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", hjust = 0,
                                      size = base_size, margin = margin(b = 3)),
      axis.line        = element_line(linewidth = 0.3, colour = "grey30"),
      axis.ticks       = element_line(linewidth = 0.3, colour = "grey30"),
      axis.text        = element_text(colour = "grey20", size = base_size - 1),
      axis.title       = element_text(size = base_size),
      panel.spacing.x  = unit(0.9, "lines"),
      panel.spacing.y  = unit(1.1, "lines"),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

save_fig <- function(p, name, n_panels, ncol = 3) {
  h <- 5 * ceiling(n_panels / ncol) + 1        # cm per panel row
  ggsave(paste0(dir_fig, name, ".png"), p, width = 16, height = h,
         units = "cm", dpi = 600, bg = "white")
  tryCatch(
    ggsave(paste0(dir_fig, name, ".pdf"), p, width = 16, height = h,
           units = "cm", device = cairo_pdf),
    error = function(e) ggsave(paste0(dir_fig, name, ".pdf"), p,
                               width = 16, height = h, units = "cm")
  )
  cat(sprintf("Saved: %s%s.png / .pdf\n", dir_fig, name))
}

# --- 8.1 Training data from the model ----------------------------
mf        <- gam_final$model
resp_name <- names(mf)[1]
pred_vars <- names(mf)[-1]
pred_vars <- pred_vars[!grepl("^offset\\(|^\\(", pred_vars)]   # drop offset / weights
pred_vars <- pred_vars[sapply(mf[pred_vars], is.numeric)]

train <- data.frame(PUD = as.numeric(mf[[resp_name]]), mf[pred_vars],
                    check.names = FALSE)

# --- 8.2 PUD vs. each predictor (Spearman correlation) -----------
long <- pivot_longer(train, all_of(pred_vars),
                     names_to = "variable", values_to = "value")

rho <- sapply(pred_vars, function(v)
  suppressWarnings(cor(train$PUD, train[[v]], method = "spearman",
                       use = "complete.obs")))
p_val <- sapply(pred_vars, function(v)
  suppressWarnings(cor.test(train$PUD, train[[v]], method = "spearman",
                            exact = FALSE)$p.value))
stars <- cut(p_val, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))

cor_tab <- data.frame(variable = pred_vars, rho = rho, p = p_val,
                      label = sprintf('italic(r)[s] == "%.2f%s"', rho, stars))
cor_tab <- cor_tab[order(-abs(cor_tab$rho)), ]     # strongest first
print(cor_tab[, c("variable", "rho", "p")], row.names = FALSE)
write.csv(cor_tab, paste0(dir_fig, "PUD_spearman_correlations.csv"), row.names = FALSE)

lev <- nice(cor_tab$variable)
long$variable    <- factor(nice(long$variable), levels = lev)
cor_tab$variable <- factor(nice(cor_tab$variable), levels = lev)

p_cor <- ggplot(long, aes(value, PUD)) +
  geom_point(size = 0.5, alpha = 0.35, colour = col_point, stroke = 0) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE,
              colour = col_main, fill = col_band, linewidth = 0.6) +
  geom_text(data = cor_tab, aes(x = Inf, y = Inf, label = label),
            hjust = 1.05, vjust = 1.3, size = 2.7, colour = "grey15",
            parse = TRUE, inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_x", ncol = 3) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
  coord_cartesian(ylim = c(0, max(train$PUD, na.rm = TRUE))) +
  labs(x = NULL, y = "Observed PUD (log scale)") +
  theme_thesis()

save_fig(p_cor, "Fig_PUD_vs_predictors", length(pred_vars))

# --- 8.3 Partial effect plots ------------------------------------
# plot.gam() computes the curves; we only take its numbers and draw
# them with ggplot. seWithMean = TRUE includes intercept uncertainty.
pdf(NULL)
pg <- plot(gam_final, pages = 1, seWithMean = TRUE, n = 200)
dev.off()

edf <- summary(gam_final)$s.table[, "edf"]
pv  <- summary(gam_final)$s.table[, "p-value"]

smooth_df <- do.call(rbind, lapply(seq_along(pg), function(i) {
  s <- pg[[i]]
  if (is.null(s[["fit"]]) || !is.null(s[["y"]]) || length(s[["x"]]) == 0) return(NULL)  # 1D smooths only
  v <- gam_final$smooth[[i]]$term
  data.frame(term = v, x = s$x, fit = as.numeric(s$fit),
             lwr = as.numeric(s$fit) - 2 * s$se, upr = as.numeric(s$fit) + 2 * s$se,
             label = sprintf("%s\nedf = %.2f%s", nice(v), edf[i],
                             as.character(cut(pv[i], c(-Inf, .001, .01, .05, Inf),
                                              c(", p < 0.001", ", p < 0.01",
                                                ", p < 0.05", ", n.s.")))))
}))

rug_df <- do.call(rbind, lapply(unique(smooth_df$term), function(v)
  data.frame(term = v, x = train[[v]])))

# order panels like the correlation figure
ord <- intersect(cor_tab$variable, nice(unique(smooth_df$term)))
lab_order <- unique(smooth_df$label)[match(ord, nice(unique(smooth_df$term)))]
smooth_df$label <- factor(smooth_df$label, levels = lab_order)
rug_df$label <- factor(smooth_df$label[match(rug_df$term, smooth_df$term)],
                       levels = lab_order)

plot_pe <- function(df, ref, ylab, log_y = FALSE) {
  p <- ggplot(df, aes(x, fit)) +
    geom_hline(yintercept = ref, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = col_band, alpha = 0.7) +
    geom_line(colour = col_main, linewidth = 0.7) +
    geom_rug(data = rug_df, aes(x = x), inherit.aes = FALSE, sides = "b",
             alpha = 0.25, linewidth = 0.2, length = unit(0.04, "npc")) +
    facet_wrap(~ label, scales = "free_x", ncol = 3) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
    labs(x = NULL, y = ylab) +
    theme_thesis()
  if (log_y) p <- p + scale_y_log10()
  p
}

p_pe <- plot_pe(smooth_df, 0, "Partial effect on PUD (link scale)")
save_fig(p_pe, "Fig_GAM_partial_effects", length(lab_order))

# Same figure as multiplicative effect (only meaningful with a log link):
# 2 = twice the PUD compared with an average cell, 0.5 = half.
if (family(gam_final)$link == "log") {
  smooth_exp <- transform(smooth_df, fit = exp(fit), lwr = exp(lwr), upr = exp(upr))
  p_pe_exp <- plot_pe(smooth_exp, 1, "Multiplicative effect on PUD", log_y = TRUE)
  save_fig(p_pe_exp, "Fig_GAM_partial_effects_multiplicative", length(lab_order))
}

# ---------------------------------------------------------------
# 9. THESIS FIGURE: CROSS-VALIDATION
# ---------------------------------------------------------------
cv_runs_csv <- paste0(dir_model, "PUD_spatialCV_all_runs.csv")   # from script 06

# --- 9.1 Spatial CV: train vs. test pseudo-R2 per model ------------
cv <- read.csv(cv_runs_csv)
model_levels <- c("LM", "GLM", "GAM", "RandomForest", "XGBoost")
model_labels <- c("LM", "GLM", "GAM", "Random\nForest", "XGBoost")
cv$Model <- factor(cv$Model, levels = model_levels, labels = model_labels)
cv$Set   <- factor(cv$Set, levels = c("Train", "Test"))

# 50 repetitions do not necessarily give 50 different fold layouts
# (k-means on the same coordinates often converges to the same solution)
n_distinct <- aggregate(Pseudo_R2 ~ Model, data = cv[cv$Set == "Test", ],
                        FUN = function(x) length(unique(round(x, 8))))
names(n_distinct)[2] <- "distinct_fold_layouts"
cat("\nDistinct spatial fold configurations (out of", max(cv$Repetition), "repetitions):\n")
print(n_distinct, row.names = FALSE)

cv_mean <- aggregate(Pseudo_R2 ~ Model + Set, data = cv, FUN = mean)

p_cv <- ggplot(cv, aes(Model, Pseudo_R2, colour = Set, fill = Set)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6,
               alpha = 0.25, linewidth = 0.4, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75,
                                             seed = 1),
             size = 0.9, alpha = 0.5, stroke = 0) +
  geom_point(data = cv_mean, position = position_dodge(width = 0.75),
             shape = 23, size = 2, colour = "black", stroke = 0.4) +
  scale_colour_manual(values = c(Train = "grey45", Test = col_main)) +
  scale_fill_manual(values   = c(Train = "grey45", Test = col_main)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 6)) +
  labs(x = NULL, y = expression("Pseudo-"*R^2), colour = NULL, fill = NULL) +
  theme_thesis(base_size = 9) +
  theme(legend.position = "top", legend.justification = "left",
        legend.margin = margin(0, 0, -4, 0),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"))

ggsave(paste0(dir_fig, "Fig_spatialCV_pseudoR2.png"), p_cv,
       width = 16, height = 8, units = "cm", dpi = 600, bg = "white")
ggsave(paste0(dir_fig, "Fig_spatialCV_pseudoR2.pdf"), p_cv,
       width = 16, height = 8, units = "cm")
cat(sprintf("Saved: %sFig_spatialCV_pseudoR2.png / .pdf\n", dir_fig))


# ---------------------------------------------------------------
# 10. SUMMARY STATISTICS OF THE FULL-AREA PREDICTION (for thesis text)
# ---------------------------------------------------------------
# Writes all values reported in the results section to CSV files:
#   fullarea_prediction_summary.csv        -- counts, range, median, quantiles
#   fullarea_prediction_by_distance.csv    -- median PUD by distance to Uhuru Peak

# Uhuru Peak (WGS84 lon/lat) -- replace with the point used in your QGIS map if different
uhuru_lonlat <- c(37.3533, -3.0758)

pred_sf <- grid %>%
  left_join(full_data %>% select(cell_id, PUD_predicted), by = "cell_id")

uhuru <- st_transform(st_sfc(st_point(uhuru_lonlat), crs = 4326), st_crs(pred_sf))
pred_sf$dist_summit_km <- as.numeric(st_distance(st_centroid(st_geometry(pred_sf)), uhuru)) / 1000

p <- pred_sf$PUD_predicted
i_max <- which.max(p)

# --- 10.1 Overall summary --------------------------------------
pred_summary <- data.frame(
  statistic = c("n_cells_total", "n_cells_predicted", "n_cells_not_predicted",
                "min", "q01", "q25", "median", "mean", "q75", "q95", "q99", "max",
                "share_below_0.30_percent",
                "max_cell_id", "max_dist_to_summit_km"),
  value = c(length(p), sum(!is.na(p)), sum(is.na(p)),
            min(p, na.rm = TRUE), quantile(p, c(0.01, 0.25), na.rm = TRUE),
            median(p, na.rm = TRUE), mean(p, na.rm = TRUE),
            quantile(p, c(0.75, 0.95, 0.99), na.rm = TRUE), max(p, na.rm = TRUE),
            100 * mean(p < 0.30, na.rm = TRUE),
            pred_sf$cell_id[i_max], pred_sf$dist_summit_km[i_max])
)
pred_summary$value <- round(pred_summary$value, 4)
write.csv(pred_summary, paste0(dir_fullarea, "fullarea_prediction_summary.csv"), row.names = FALSE)
print(pred_summary, row.names = FALSE)

# --- 10.2 Median PUD by distance to the summit -----------------
dist_breaks <- c(0, 3, 6, 10, 15, 20, 30, 40, 60, Inf)
pred_by_dist <- pred_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(PUD_predicted)) %>%
  mutate(distance_km = cut(dist_summit_km, dist_breaks, include.lowest = TRUE)) %>%
  group_by(distance_km) %>%
  summarise(n_cells = n(),
            median_PUD = round(median(PUD_predicted), 4),
            mean_PUD   = round(mean(PUD_predicted), 4),
            .groups = "drop")
write.csv(pred_by_dist, paste0(dir_fullarea, "fullarea_prediction_by_distance.csv"), row.names = FALSE)
print(as.data.frame(pred_by_dist), row.names = FALSE)

rho_dist <- cor(pred_sf$PUD_predicted, pred_sf$dist_summit_km,
                method = "spearman", use = "complete.obs")
cat(sprintf("Spearman correlation predicted PUD vs. distance to summit: %.3f\n", rho_dist))

cat("\nSaved summary tables to", dir_fullarea, "\n")
