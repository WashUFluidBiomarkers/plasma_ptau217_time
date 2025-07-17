# Created by Kellen Petersen, July 1, 2025

library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(here)
library(DescTools)
library(nopaco)

setwd("<SET WORKING DIRECTORY>")

which_method <- "TIRA"

if (which_method == "TIRA") {
  clock_adni <- read.csv(here("results","<LOAD DATA>"))
  df_adni <- read.csv(here("results","<LOAD DATA>"))
  
  clock_kadrc <- read.csv(here("results","<LOAD DATA>"))
  df_kadrc <- read.csv(here("results","<LOAD DATA>"))
  
} else if (which_method == "SILA") {
  clock_adni <- read.csv(here("results_SILA","<LOAD DATA>")) %>% 
    mutate(clock_time = adtime)
  df_adni <- read.csv(here("results","<LOAD DATA>")) %>% 
    arrange(ID, EXAMDATE) %>% 
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>% 
    ungroup()
  
  clock_kadrc <- read.csv(here("results_SILA","<LOAD DATA>")) %>% 
    mutate(clock_time = adtime)
  df_kadrc <- read.csv(here("results","<LOAD DATA>"))  %>% 
    arrange(ID, EXAMDATE) %>% 
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>% 
    ungroup()
  
} else {
  stop("Invalid method selected. Please choose 'TIRA' or 'SILA'.")
}

df_adni <- df_adni %>% filter(plasma>1.06 & plasma<10.45)
df_kadrc <- df_kadrc %>% filter(plasma>1.06 & plasma<10.45)

find_closest_plasma <- function(plasma_value, reference_df) {
  if (is.na(plasma_value)) return(NA)
  idx <- which.min(abs(reference_df$plasma - plasma_value))
  return(as.numeric(reference_df$clock_time[idx])) 
}

DF_KADRC2 <- df_kadrc %>%
  group_by(ID) %>%
  mutate(
    plasma_time_KADRC = as.numeric(sapply(plasma, find_closest_plasma, clock_kadrc)),
    plasma_onset_each_KADRC = as.numeric(AGE) - plasma_time_KADRC,
    est_plasma_onset_KADRC = mean(plasma_onset_each_KADRC, na.rm = TRUE),
    plasma_time_ADNI = as.numeric(sapply(plasma, find_closest_plasma, clock_adni)),
    plasma_onset_each_ADNI = as.numeric(AGE) - plasma_time_ADNI,
    est_plasma_onset_ADNI = mean(plasma_onset_each_ADNI, na.rm = TRUE)
  ) %>%
  ungroup()

DF_ADNI2 <- df_adni %>%
  group_by(ID) %>%
  mutate(
    plasma_time_ADNI = as.numeric(sapply(plasma, find_closest_plasma, clock_adni)),
    plasma_onset_each_ADNI = as.numeric(AGE) - plasma_time_ADNI,
    est_plasma_onset_ADNI = mean(plasma_onset_each_ADNI, na.rm = TRUE),
    plasma_time_KADRC = as.numeric(sapply(plasma, find_closest_plasma, clock_kadrc)),
    plasma_onset_each_KADRC = as.numeric(AGE) - plasma_time_KADRC,
    est_plasma_onset_KADRC = mean(plasma_onset_each_KADRC, na.rm = TRUE)
  ) %>%
  ungroup()

combined_data <- bind_rows(
  DF_KADRC2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI),
  DF_ADNI2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI)
)

lm_fit <- lm(est_plasma_onset_KADRC ~ est_plasma_onset_ADNI, data = combined_data)

spearman_result <- cor.test(combined_data$est_plasma_onset_KADRC, 
                            combined_data$est_plasma_onset_ADNI, 
                            method = "spearman")
pearson_result <- cor.test(combined_data$est_plasma_onset_KADRC, 
                           combined_data$est_plasma_onset_ADNI, 
                           method = "pearson")

r2 <- summary(lm_fit)$r.squared
adj_r2 <- summary(lm_fit)$adj.r.squared

ccc_result <- CCC(combined_data$est_plasma_onset_KADRC, 
                  combined_data$est_plasma_onset_ADNI)$rho.c
ccc <- ccc_result$est

mat <- cbind(combined_data$est_plasma_onset_KADRC, 
             combined_data$est_plasma_onset_ADNI)
conc_test <- concordance.test(mat)
psi <- getPsi(conc_test)
np_ccc <- rfromPsi(psi)

cat("\n--- Cross-Comparison Metrics (", which_method, ") ---\n")
cat("R-squared:            ", round(r2, 4), "\n")
cat("Adjusted R-squared:   ", round(adj_r2, 4), "\n")
cat("Pearson r:            ", round(pearson_result$estimate, 4), "\n")
cat("Spearman r:           ", round(spearman_result$estimate, 4), "\n")
cat("CCC:                  ", round(ccc, 4), "\n")
cat("Non-parametric CCC:   ", round(np_ccc, 4), "\n")

comparison_plot <- ggplot() +
  geom_point(data = DF_KADRC2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI),
             aes(x = est_plasma_onset_ADNI, y = est_plasma_onset_KADRC, color = "Knight ADRC"),
             size = 2, alpha = 1) +
  geom_point(data = DF_ADNI2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI),
             aes(x = est_plasma_onset_ADNI, y = est_plasma_onset_KADRC, color = "ADNI"),
             size = 2, alpha = 0.6) +
  annotate("text", x = 40, y = 110,
           label = paste0(
             "Adjusted R² = ", sprintf("%.3f", summary(lm_fit)$adj.r.squared), "\n",
             "Pearson r = ", sprintf("%.3f", pearson_result$estimate)),
           color = "black", hjust = 0, size = 5) +
  scale_color_manual(values = c("ADNI" = "black", "Knight ADRC" = "darkgreen"),
                     breaks = c("Knight ADRC", "ADNI")) +
  labs(title = str_c("For all measured %p-tau217 values from both datasets using ", which_method),
       x = "ADNI model\nEstimated age at %p-tau217 positivity (years)",
       y = "Knight ADRC model\nEstimated age at %p-tau217 positivity (years)",
       color = "Dataset") +
  cowplot::theme_cowplot() +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(hjust = 0.5)) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype= "dashed", linewidth = 1) +
  xlim(25,115) +
  ylim(25,115) +
  theme(aspect.ratio = 1)
comparison_plot

filename0 <- str_c("results/intramethod_cross_comparison_plot_", which_method, 
                   "_", "FINAL")
saveRDS(comparison_plot, str_c(filename0, ".rds"))
ggsave(str_c(filename0,".png"), comparison_plot, width = 7.5, height = 7.5, dpi = 500, bg = "white")
