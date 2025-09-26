# Originally created by Kellen Petersen, July 1, 2025
# Updated by Kellen Petersen, September 12, 2025

# Load required libraries for data analysis and visualization
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(here)
library(DescTools)
library(nopaco)

setwd("")

# Set method for cross-dataset comparison analysis
which_method <- "TIRA"

# Load clock and trajectory data based on selected method
if (which_method == "TIRA") {
  clock_adni <- read.csv(here("results",""))
  df_adni <- read.csv(here("results",""))
  clock_kadrc <- read.csv(here("results",""))
  df_kadrc <- read.csv(here("results",""))
} else if (which_method == "SILA") {
  # Load SILA-specific data with appropriate transformations
  clock_adni <- read.csv(here("results_SILA","")) %>%
    mutate(clock_time = adtime)
  df_adni <- read.csv(here("results","")) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>%
    ungroup()
  
  clock_kadrc <- read.csv(here("results_SILA","")) %>%
    mutate(clock_time = adtime)
  df_kadrc <- read.csv(here("results","")) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>%
    ungroup()
} else {
  stop("Invalid method selected. Please choose 'TIRA' or 'SILA'.")
}

# Filter data to reliable plasma measurement range
df_adni <- df_adni %>% filter(plasma>1.06 & plasma<10.45)
df_kadrc <- df_kadrc %>% filter(plasma>1.06 & plasma<10.45)

# Function to find closest plasma time point in clock data
find_closest_plasma <- function(plasma_value, reference_df) {
  if (is.na(plasma_value)) return(NA)
  idx <- which.min(abs(reference_df$plasma - plasma_value))
  return(as.numeric(reference_df$clock_time[idx]))
}

# Apply both KADRC and ADNI clocks to KADRC data for cross-comparison
DF_KADRC2 <- df_kadrc %>%
  group_by(ID) %>%
  mutate(
    # Calculate onset age using KADRC-derived clock
    plasma_time_KADRC = as.numeric(sapply(plasma, find_closest_plasma, clock_kadrc)),
    plasma_onset_each_KADRC = as.numeric(AGE) - plasma_time_KADRC,
    est_plasma_onset_KADRC = mean(plasma_onset_each_KADRC, na.rm = TRUE),
    # Calculate onset age using ADNI-derived clock
    plasma_time_ADNI = as.numeric(sapply(plasma, find_closest_plasma, clock_adni)),
    plasma_onset_each_ADNI = as.numeric(AGE) - plasma_time_ADNI,
    est_plasma_onset_ADNI = mean(plasma_onset_each_ADNI, na.rm = TRUE)
  ) %>%
  ungroup()

# Apply both clocks to ADNI data for cross-comparison
DF_ADNI2 <- df_adni %>%
  group_by(ID) %>%
  mutate(
    # Calculate onset age using ADNI-derived clock
    plasma_time_ADNI = as.numeric(sapply(plasma, find_closest_plasma, clock_adni)),
    plasma_onset_each_ADNI = as.numeric(AGE) - plasma_time_ADNI,
    est_plasma_onset_ADNI = mean(plasma_onset_each_ADNI, na.rm = TRUE),
    # Calculate onset age using KADRC-derived clock
    plasma_time_KADRC = as.numeric(sapply(plasma, find_closest_plasma, clock_kadrc)),
    plasma_onset_each_KADRC = as.numeric(AGE) - plasma_time_KADRC,
    est_plasma_onset_KADRC = mean(plasma_onset_each_KADRC, na.rm = TRUE)
  ) %>%
  ungroup()

# Combine data from both datasets for comparison analysis
combined_data <- bind_rows(
  DF_KADRC2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI),
  DF_ADNI2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI)
)

# Perform statistical analyses comparing the two clock estimates
lm_fit <- lm(est_plasma_onset_KADRC ~ est_plasma_onset_ADNI, data = combined_data)

spearman_result <- cor.test(combined_data$est_plasma_onset_KADRC,
                            combined_data$est_plasma_onset_ADNI,
                            method = "spearman")

pearson_result <- cor.test(combined_data$est_plasma_onset_KADRC,
                           combined_data$est_plasma_onset_ADNI,
                           method = "pearson")

r2 <- summary(lm_fit)$r.squared
adj_r2 <- summary(lm_fit)$adj.r.squared

# Calculate concordance correlation coefficient for agreement assessment
ccc_result <- CCC(combined_data$est_plasma_onset_KADRC,
                  combined_data$est_plasma_onset_ADNI)$rho.c
ccc <- ccc_result$est

# Calculate non-parametric concordance correlation coefficient
mat <- cbind(combined_data$est_plasma_onset_KADRC,
             combined_data$est_plasma_onset_ADNI)
conc_test <- concordance.test(mat)
psi <- getPsi(conc_test)
np_ccc <- rfromPsi(psi)

# Display comparison metrics for cross-dataset validation
cat("\n--- Cross-Comparison Metrics (", which_method, ") ---\n")
cat("R-squared: ", round(r2, 4), "\n")
cat("Adjusted R-squared: ", round(adj_r2, 4), "\n")
cat("Pearson r: ", round(pearson_result$estimate, 4), "\n")
cat("Spearman r: ", round(spearman_result$estimate, 4), "\n")
cat("CCC: ", round(ccc, 4), "\n")
cat("Non-parametric CCC: ", round(np_ccc, 4), "\n")

# Create scatter plot comparing ADNI vs KADRC clock estimates
comparison_plot <- ggplot() +
  # Plot KADRC participants using both clocks
  geom_point(data = DF_KADRC2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI),
             aes(x = est_plasma_onset_ADNI, y = est_plasma_onset_KADRC, color = "Knight ADRC"),
             size = 2, alpha = 1) +
  # Plot ADNI participants using both clocks
  geom_point(data = DF_ADNI2 %>% distinct(ID, est_plasma_onset_KADRC, est_plasma_onset_ADNI),
             aes(x = est_plasma_onset_ADNI, y = est_plasma_onset_KADRC, color = "ADNI"),
             size = 2, alpha = 0.6) +
  # Add correlation statistics as annotation
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
  # Add identity line for perfect agreement reference
  geom_abline(intercept = 0, slope = 1, color = "red", linetype= "dashed", linewidth = 1) +
  xlim(25,115) +
  ylim(25,115) +
  theme(aspect.ratio = 1)

comparison_plot

# Save intra-method cross-dataset comparison plot
filename0 <- str_c("results/intramethod_cross_comparison_plot_", which_method,
                   "_", "FINAL")
saveRDS(comparison_plot, str_c(filename0, ".rds"))
ggsave(str_c(filename0,".png"), comparison_plot, width = 7.5, height = 7.5, dpi = 500, bg = "white")
