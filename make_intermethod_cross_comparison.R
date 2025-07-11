# Created by Kellen Petersen, July 1, 2025

library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(here)
library(DescTools)
library(nopaco)

setwd("<SET WORKING DIRECTORY>")

which_dataset <- "ADNI" # "KADRC" or "ADNI"
which_dataset_title <- which_dataset

if(which_dataset == "ADNI"){
  df_adni <- read.csv(here("results","TrajPlotData_ADNI_dataset_TIRA_method_FINAL.csv")) %>% 
    arrange(ID, EXAMDATE)
  df_adni_SILA <- read.csv(here("results","TrajPlotData_ADNI_dataset_SILA_method_FINAL.csv")) %>% 
    arrange(ID, EXAMDATE)
  
  df <- df_adni
  df_SILA <- df_adni_SILA
  
} else {
  df_kadrc <- read.csv(here("results","TrajPlotData_KADRC_dataset_TIRA_method_FINAL.csv")) %>% 
    arrange(ID, EXAMDATE)
  df_kadrc_SILA <- read.csv(here("results","TrajPlotData_KADRC_dataset_SILA_method_FINAL.csv")) %>% 
    arrange(ID, EXAMDATE)
  
  df <- df_kadrc
  df_SILA <- df_kadrc_SILA
}

df_tira <- df %>% 
  group_by(ID) %>%
  slice(1) %>% 
  ungroup() %>%
  select(ID, est_onset_age) %>% 
  rename(est_onset_age_TIRA = est_onset_age)

df_sila <- df_SILA %>% 
  group_by(ID) %>%
  slice(1) %>% 
  ungroup() %>%
  select(ID, est_onset_age) %>% 
  rename(est_onset_age_SILA = est_onset_age)

combined_data <- df_tira %>%
  inner_join(df_sila, by = "ID")

lm_fit <- lm(est_onset_age_TIRA ~ est_onset_age_SILA, data = combined_data)

spearman_result <- cor.test(combined_data$est_onset_age_TIRA, 
                            combined_data$est_onset_age_SILA, 
                            method = "spearman")
pearson_result <- cor.test(combined_data$est_onset_age_TIRA, 
                           combined_data$est_onset_age_SILA, 
                           method = "pearson")

r2 <- summary(lm_fit)$r.squared
adj_r2 <- summary(lm_fit)$adj.r.squared

if (!requireNamespace("DescTools", quietly = TRUE)) install.packages("DescTools")
ccc_result <- DescTools::CCC(
  combined_data$est_onset_age_TIRA,
  combined_data$est_onset_age_SILA,
  na.rm = TRUE
)$rho.c
ccc <- ccc_result$est

if (!requireNamespace("nopaco", quietly = TRUE)) install.packages("nopaco")
mat <- cbind(combined_data$est_onset_age_TIRA, 
             combined_data$est_onset_age_SILA)
conc_test <- concordance.test(mat)
psi <- getPsi(conc_test)
np_ccc <- rfromPsi(psi)

cat("\n--- Cross-Comparison Metrics ---\n")
cat("R-squared:            ", round(r2, 4), "\n")
cat("Adjusted R-squared:   ", round(adj_r2, 4), "\n")
cat("Pearson r:            ", round(pearson_result$estimate, 4), "\n")
cat("Spearman r:           ", round(spearman_result$estimate, 4), "\n")
cat("CCC:                  ", round(ccc, 4), "\n")
cat("Non-parametric CCC:   ", round(np_ccc, 4), "\n")

comparison_plot <- ggplot() +
  geom_point(data = combined_data,
             aes(x = est_onset_age_SILA, y = est_onset_age_TIRA),
             size = 2, alpha = 0.8, color = "black") +
  annotate("text", x = min(combined_data$est_onset_age_SILA, combined_data$est_onset_age_TIRA, na.rm = TRUE), 
           y = max(combined_data$est_onset_age_SILA, combined_data$est_onset_age_TIRA, na.rm = TRUE) ,
           label = paste0(
             "Adjusted R² = ", sprintf("%.3f", summary(lm_fit)$adj.r.squared), "\n",
             "Pearson r = ", sprintf("%.3f", pearson_result$estimate)),
           color = "black", hjust = 0, size = 5) +
  labs(title = str_c(which_dataset_title),
       x = "SILA model\nEstimated age at %p-tau217 positivity (years)",
       y = "TIRA model\nEstimated age at %p-tau217 positivity (years)") +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", linewidth = 1) +
  coord_fixed(ratio = 1, xlim = c(min(combined_data$est_onset_age_SILA, combined_data$est_onset_age_TIRA) - 5,
                                  max(combined_data$est_onset_age_SILA, combined_data$est_onset_age_TIRA) + 5),
              ylim = c(min(combined_data$est_onset_age_SILA, combined_data$est_onset_age_TIRA) - 5,
                       max(combined_data$est_onset_age_SILA, combined_data$est_onset_age_TIRA) + 5))

comparison_plot

filename0 <- str_c("results/Inter-Method_comparison_", which_dataset, 
                   "_data_", "FINAL")
saveRDS(comparison_plot, str_c(filename0, ".rds"))
ggsave(str_c(filename0, ".png"), comparison_plot, width = 8, height = 8, dpi = 500, bg = "white")
