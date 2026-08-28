install.packages("moments")
library(moments)
setwd("C:/Users/Lukas/masterthesis/thesis_2026/data/landscape_indicators")
setwd("C:/Users/Lukas/masterthesis/thesis_2026/data/PUD_Analysis/n_images_per_grid_cell.csv")
data <- read.csv("master_indicator_table_SCREENING_ONLY_4.csv")
data_2 <- read.csv("n_images_per_grid_cell.csv")

response_var <- c(
  "avg_annual_PUD"
)

response_var_2 <- c(
  "n_images"
)

predictor_vars <- c("dist_sd_near_m", "visible_area_middle_km2", "dist_sd_middle_m",
                     "ShapeIndex", "ForestEcotoneDensity", "ReliefDiversity", "NDSI",
                     "NDVI_var", "visible_naturalness_middle", "visible_water_middle",
                     "VerticalStructuralHeterogeneity", "LAI_var", "DistanceToWater_km",
                     "ContagionIndex", "ColorDiversity", "ColorDiversity_var",
                     "visible_ColorDiversity", "DistanceToWaterway_km")

png(paste0("histogram_", response_var, ".png"))
hist(data$avg_annual_PUD, main = response_var)
dev.off()


hist(data$avg_annual_PUD)
png(paste0("histogram_", response_var_2, ".png"))
hist(data_2$n_images, main = response_var_2)
dev.off()

first <- quantile(data$avg_annual_PUD, 0.25)
first
third <- quantile(data$avg_annual_PUD, 0.75)
third

IQR_value <- IQR(data$avg_annual_PUD)
IQR_value

lower <- max(0,first - 1.5 *IQR_value)
lower
upper <- third + 1.5 *IQR_value
upper

df_clean <- data[data$avg_annual_PUD>= lower& data$avg_annual_PUD<= upper,]
boxplot(df_clean$avg_annual_PUD, main = response_var)


png(paste0("boxplot_", response_var, ".png"))
boxplot(data$avg_annual_PUD, main = response_var)
dev.off()

png(paste0("boxplot_", response_var_2, ".png"))
boxplot(data_2$n_images, main = response_var_2)
dev.off()

distribution_result_table <- data.frame(
  Response  = response_var,
  Median    = median(data$avg_annual_PUD, na.rm = TRUE),
  Mean      = mean(data$avg_annual_PUD, na.rm = TRUE),
  Variance  = var(data$avg_annual_PUD, na.rm = TRUE),
  Min       = min(data$avg_annual_PUD, na.rm = TRUE),
  Max       = max(data$avg_annual_PUD, na.rm = TRUE),
  Number_0  = sum(data$avg_annual_PUD == 0, na.rm = TRUE),
  Share_0 = mean(data$avg_annual_PUD == 0, na.rm = TRUE),
  Number_NA = sum(is.na(data$avg_annual_PUD)),
  Share_NA = mean(is.na(data$avg_annual_PUD)),
  Skewness = moments::skewness(data$avg_annual_PUD, na.rm = TRUE)
)

distribution_result_table_2 <- data.frame(
  Response  = response_var_2,
  Median    = median(data_2$n_images, na.rm = TRUE),
  Mean      = mean(data_2$n_images, na.rm = TRUE),
  Variance  = var(data_2$n_images, na.rm = TRUE),
  Min       = min(data_2$n_images, na.rm = TRUE),
  Max       = max(data_2$n_images, na.rm = TRUE),
  Number_0  = sum(data_2$n_images == 0, na.rm = TRUE),
  Share_0 = mean(data_2$n_images == 0, na.rm = TRUE),
  Number_NA = sum(is.na(data_2$n_images)),
  Share_NA = mean(is.na(data_2$n_images)),
  Skewness = moments::skewness(data_2$n_images, na.rm = TRUE),
  Dispersion_Index = var(data_2$n_images, na.rm = TRUE) / mean(data_2$n_images, na.rm = TRUE)
)

distribution_result_table
write.csv(distribution_result_table, "distribution_result_table_PUD.csv", row.names = FALSE)

distribution_result_table_2
write.csv(distribution_result_table_2, "distribution_result_table_number_of_photos.csv", row.names = FALSE)
library(mgcv)

gam_explore <- gam(
  avg_annual_PUD ~ s(dist_sd_near_m) + s(visible_area_middle_km2) + s(dist_sd_middle_m) +
    s(ShapeIndex) + s(ForestEcotoneDensity) + s(ReliefDiversity) + s(NDSI) + s(NDVI_var) +
    s(visible_naturalness_middle) + s(visible_water_middle) + s(VerticalStructuralHeterogeneity) +
    s(LAI_var) + s(DistanceToWater_km) + s(ContagionIndex) + s(ColorDiversity) +
    s(ColorDiversity_var) + s(visible_ColorDiversity) + s(DistanceToWaterway_km),
  family = Gamma(link = "log"),
  data = data
)

summary(gam_explore)      
png("scatterplot.png", width = 1200, height = 800, res = 150)
plot(gam_explore, pages = 1, shade = TRUE)
dev.off()

s_table <- as.data.frame(summary(gam_explore)$s.table)
s_table$Term <- rownames(s_table)
write.csv(s_table, "gam_summary.csv", row.names = FALSE)

m_gamma <- glm(
  avg_annual_PUD ~ dist_sd_near_m + visible_area_middle_km2 + dist_sd_middle_m +
    ShapeIndex + ForestEcotoneDensity + ReliefDiversity + NDSI + NDVI_var +
    visible_naturalness_middle + visible_water_middle + VerticalStructuralHeterogeneity +
    LAI_var + DistanceToWater_km + ContagionIndex + ColorDiversity + ColorDiversity_var +
    visible_ColorDiversity + DistanceToWaterway_km,
  family = Gamma(link = "log"),
  data = data
)

m_lognorm <- glm(
  log(avg_annual_PUD) ~ dist_sd_near_m + visible_area_middle_km2 + dist_sd_middle_m +
    ShapeIndex + ForestEcotoneDensity + ReliefDiversity + NDSI + NDVI_var +
    visible_naturalness_middle + visible_water_middle + VerticalStructuralHeterogeneity +
    LAI_var + DistanceToWater_km + ContagionIndex + ColorDiversity + ColorDiversity_var +
    visible_ColorDiversity + DistanceToWaterway_km,
  family = gaussian(),
  data = data
)

sim_gamma <- simulateResiduals(m_gamma)
testDispersion(sim_gamma)
testOutliers(sim_gamma)

sim_lognorm <- simulateResiduals(m_lognorm)
testDispersion(sim_lognorm)
testOutliers(sim_lognorm)





####


library(DHARMa)

# ---- 1. Modelle fitten ----
m_gamma <- glm(
  avg_annual_PUD ~ dist_sd_near_m + visible_area_middle_km2 + dist_sd_middle_m +
    ShapeIndex + ForestEcotoneDensity + ReliefDiversity + NDSI + NDVI_var +
    visible_naturalness_middle + visible_water_middle + VerticalStructuralHeterogeneity +
    LAI_var + DistanceToWater_km + ContagionIndex + ColorDiversity + ColorDiversity_var +
    visible_ColorDiversity + DistanceToWaterway_km,
  family = Gamma(link = "log"),
  data = data
)

m_lognorm <- glm(
  log(avg_annual_PUD) ~ dist_sd_near_m + visible_area_middle_km2 + dist_sd_middle_m +
    ShapeIndex + ForestEcotoneDensity + ReliefDiversity + NDSI + NDVI_var +
    visible_naturalness_middle + visible_water_middle + VerticalStructuralHeterogeneity +
    LAI_var + DistanceToWater_km + ContagionIndex + ColorDiversity + ColorDiversity_var +
    visible_ColorDiversity + DistanceToWaterway_km,
  family = gaussian(),
  data = data
)

# ---- 2. DHARMa simulieren ----
sim_gamma   <- simulateResiduals(m_gamma)
sim_lognorm <- simulateResiduals(m_lognorm)

disp_gamma   <- testDispersion(sim_gamma, plot = FALSE)
out_gamma    <- testOutliers(sim_gamma, plot = FALSE)
disp_lognorm <- testDispersion(sim_lognorm, plot = FALSE)
out_lognorm  <- testOutliers(sim_lognorm, plot = FALSE)

# ---- 3. Diagnostik-Plots speichern ----
png("dharma_gamma.png", width = 1400, height = 700, res = 150)
plot(sim_gamma)
dev.off()

png("dharma_lognorm.png", width = 1400, height = 700, res = 150)
plot(sim_lognorm)
dev.off()

# ---- 4. Testergebnisse als Tabelle speichern ----
dharma_results <- data.frame(
  Model     = c("Gamma", "Gamma", "Lognormal", "Lognormal"),
  Test      = c("Dispersion", "Outliers", "Dispersion", "Outliers"),
  Statistic = c(as.numeric(disp_gamma$statistic), as.numeric(out_gamma$statistic),
                as.numeric(disp_lognorm$statistic), as.numeric(out_lognorm$statistic)),
  P_Value   = c(disp_gamma$p.value, out_gamma$p.value, disp_lognorm$p.value, out_lognorm$p.value)
)

dharma_results
write.csv(dharma_results, "dharma_results_gamma_vs_lognormal.csv", row.names = FALSE)


lm_test <- lm(log(avg_annual_PUD) ~  dist_sd_near_m + visible_area_middle_km2 + dist_sd_middle_m +
                ShapeIndex + ForestEcotoneDensity + ReliefDiversity + NDSI + NDVI_var +
                visible_naturalness_middle + visible_water_middle + VerticalStructuralHeterogeneity +
                LAI_var + DistanceToWater_km + ContagionIndex + ColorDiversity + ColorDiversity_var +
                visible_ColorDiversity + DistanceToWaterway_km, data = data)  

# Q-Q-Plot -- der eigentlich relevante Test
qqnorm(residuals(lm_test))
qqline(residuals(lm_test), col = "red")

# Histogramm der Residuen
hist(residuals(lm_test), breaks = 30)

# Formaler Test (optional)
shapiro.test(residuals(lm_test))


