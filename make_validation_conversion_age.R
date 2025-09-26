# Originally created by Kellen Petersen, July 1, 2025
# Updated by Kellen Petersen, September 12, 2025

# Load required libraries for data manipulation and visualization
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(here)

setwd("")

# Resolve function conflicts with dplyr
conflicts_prefer(dplyr::lag)

# Analysis parameters
which_method <- "SILA"  # Which modeling approach to use
which_dataset <- "KADRC"  # Which dataset to analyze
avg_method <- "simple"  # Method for calculating conversion age
tip_point <- 4.06  # Plasma biomarker threshold for positivity

# Load clock and trajectory data based on selected method
if (which_method == "TIRA") {
  clock_adni <- read.csv(here("results","df_clock_ADNI_dataset_TIRA_method_FINAL.csv"))
  df_adni <- read.csv(here("results","TrajPlotData_ADNI_dataset_TIRA_method_FINAL.csv"))
  clock_kadrc <- read.csv(here("results","df_clock_KADRC_dataset_TIRA_method_FINAL.csv"))
  df_kadrc <- read.csv(here("results","TrajPlotData_KADRC_dataset_TIRA_method_FINAL.csv"))
} else if (which_method == "SILA") {
  clock_adni <- read.csv(here("results_SILA","clock_adni.csv")) %>%
    mutate(clock_time = adtime)
  df_adni <- read.csv(here("results_SILA","df_adni.csv")) %>%
    rename(est_onset_age = estaget0) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>%
    ungroup()
  
  # Create binary positivity indicator and identify converters
  df_adni <- df_adni %>%
    mutate(
      plasma_YN = if_else(plasma >= tip_point, 1, 0),
      plasma_YN = factor(plasma_YN, levels = c(0, 1), labels = c("Negative", "Positive"))
    ) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(Previous_plasma_YN = lag(plasma_YN)) %>%
    mutate(is_converter = any(Previous_plasma_YN == "Negative" & plasma_YN == "Positive", na.rm = TRUE),
           one_pos = any(plasma_YN == "Positive", na.rm = TRUE)) %>%
    select(-plasma_YN, -Previous_plasma_YN)
  
  clock_kadrc <- read.csv(here("results_SILA","clock_kadrc.csv")) %>%
    mutate(clock_time = adtime)
  df_kadrc <- read.csv(here("results_SILA","df_kadrc.csv")) %>%
    rename(est_onset_age = estaget0) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>%
    ungroup()
  
  # Apply same conversion logic to KADRC data
  df_kadrc <- df_kadrc %>%
    mutate(
      plasma_YN = if_else(plasma >= tip_point, 1, 0),
      plasma_YN = factor(plasma_YN, levels = c(0, 1), labels = c("Negative", "Positive"))
    ) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(Previous_plasma_YN = lag(plasma_YN)) %>%
    mutate(is_converter = any(Previous_plasma_YN == "Negative" & plasma_YN == "Positive", na.rm = TRUE),
           one_pos = any(plasma_YN == "Positive", na.rm = TRUE)) %>%
    select(-plasma_YN, -Previous_plasma_YN)
} else {
  stop("Invalid method selected. Please choose 'TIRA' or 'SILA'.")
}

# Select appropriate dataset for analysis
if (which_dataset == "ADNI") {
  df <- df_adni
  clock <- clock_adni
} else if (which_dataset == "KADRC") {
  df <- df_kadrc
  clock <- clock_kadrc
} else {
  stop("Invalid dataset selected. Please choose 'ADNI' or 'KADRC'.")
}

# Filter to participants who converted from negative to positive
converters_df <- df %>%
  filter(is_converter == TRUE) %>%
  arrange(ID, EXAMDATE)

# Function to determine age at which plasma biomarker crossed threshold
calculate_conversion_age <- function(id_data, tip_point, method = "simple") {
  id_data <- id_data %>% arrange(AGE)
  for (i in 1:(nrow(id_data)-1)) {
    if ((id_data$plasma[i] < tip_point && id_data$plasma[i+1] >= tip_point) ||
        (id_data$plasma[i] >= tip_point && id_data$plasma[i+1] < tip_point)) {
      age_before <- id_data$AGE[i]
      age_after <- id_data$AGE[i+1]
      plasma_before <- id_data$plasma[i]
      plasma_after <- id_data$plasma[i+1]
      if (method == "simple") {
        return((age_before + age_after) / 2)
      } else if (method == "weighted") {
        return((plasma_before * age_after + plasma_after * age_before) /
                 (plasma_before + plasma_after))
      }
    }
  }
  return(NA)
}

# Calculate observed conversion ages and compare to model estimates
conversion_ages <- converters_df %>%
  group_by(ID) %>%
  summarize(
    conversion_age = calculate_conversion_age(cur_data(), tip_point = 4.06, method = avg_method),
    est_onset_age = first(est_onset_age)
  ) %>%
  filter(!is.na(conversion_age))

# Fit linear regression model for validation
lm_fit <- lm(est_onset_age ~ conversion_age, data = conversion_ages)

# Calculate correlation statistics
pearson_r <- cor.test(conversion_ages$conversion_age,
                      conversion_ages$est_onset_age,
                      method = "pearson")
spearman_rho <- cor.test(conversion_ages$conversion_age,
                         conversion_ages$est_onset_age,
                         method = "spearman")

# Create scatter plot comparing observed vs estimated conversion ages
conversion_plot <- ggplot(conversion_ages, aes(x = conversion_age, y = est_onset_age)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fullrange = TRUE) +
  labs(
    title = "",
    x = "Observed conversion age (years)",
    y = "Estimated age at plasma %p-tau217 positivity (years)"
  ) +
  theme_cowplot() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 12)
  ) +
  xlim(55,95) +
  ylim(55,95)

# Add model statistics as plot annotation
plot_build <- ggplot_build(conversion_plot)
x_min <- plot_build$layout$panel_params[[1]]$x.range[1]
x_max <- plot_build$layout$panel_params[[1]]$x.range[2]
y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
y_max <- plot_build$layout$panel_params[[1]]$y.range[2]

conversion_plot <- conversion_plot +
  annotate("text",
           x = 55,
           y = 95,
           label = paste0(
             "\nAdjusted R² = ", sprintf("%.3f", summary(lm_fit)$adj.r.squared),
             "\nSpearman ρ = ", sprintf("%.3f", spearman_rho$estimate)
           ),
           hjust = 0,
           vjust = 0.5,
           size = 4)

print(conversion_plot)

# Save validation plot
filename0 <- str_c("results/","validation_conversion_age_plot_",
                   which_method, "_",
                   which_dataset, "_",
                   avg_method, "_avg_",
                   "FINAL")
saveRDS(conversion_plot, str_c(filename0, ".rds"))
ggsave(str_c(filename0, ".png"), conversion_plot, width = 8, height = 8, units = "in", dpi = 500, bg = "white")
