# Created by Kellen Petersen, July 1, 2025

library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(here)

setwd("<SET WORKING DIRECTORY>")

conflicts_prefer(dplyr::lag)

which_method <- "SILA"
which_dataset <- "KADRC"
avg_method <- "simple"
tip_point <- 4.06

if (which_method == "TIRA") {
  clock_adni <- read.csv(here("results","<LOAD DATA>"))
  df_adni <- read.csv(here("results","<LOAD DATA>"))
  
  clock_kadrc <- read.csv(here("results","<LOAD DATA>"))
  df_kadrc <- read.csv(here("results","<LOAD DATA>"))
  
} else if (which_method == "SILA") {
  clock_adni <- read.csv(here("results_SILA","<LOAD DATA>")) %>% 
    mutate(clock_time = adtime)
  df_adni <- read.csv(here("results_SILA","<LOAD DATA>")) %>% 
    rename(est_onset_age = estaget0) %>%
    arrange(ID, EXAMDATE) %>% 
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>% 
    ungroup()
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
  
  clock_kadrc <- read.csv(here("results_SILA","<LOAD DATA>")) %>% 
    mutate(clock_time = adtime)
  df_kadrc <- read.csv(here("results_SILA","<LOAD DATA>")) %>% 
    rename(est_onset_age = estaget0) %>%
    arrange(ID, EXAMDATE) %>% 
    group_by(ID) %>%
    mutate(years_since_onset = AGE - est_onset_age) %>% 
    ungroup()
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

if (which_dataset == "ADNI") {
  df <- df_adni
  clock <- clock_adni
} else if (which_dataset == "KADRC") {
  df <- df_kadrc
  clock <- clock_kadrc
} else {
  stop("Invalid dataset selected. Please choose 'ADNI' or 'KADRC'.")
}

converters_df <- df %>%
  filter(is_converter == TRUE) %>%
  arrange(ID, EXAMDATE)

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

conversion_ages <- converters_df %>%
  group_by(ID) %>%
  summarize(
    conversion_age = calculate_conversion_age(cur_data(), tip_point = 4.06, method = avg_method),
    est_onset_age = first(est_onset_age)
  ) %>%
  filter(!is.na(conversion_age))

lm_fit <- lm(est_onset_age ~ conversion_age, data = conversion_ages)

pearson_r <- cor.test(conversion_ages$conversion_age, 
                      conversion_ages$est_onset_age, 
                      method = "pearson")
spearman_rho <- cor.test(conversion_ages$conversion_age, 
                         conversion_ages$est_onset_age, 
                         method = "spearman")

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

filename0 <- str_c("results/","validation_conversion_age_plot_", 
                   which_method, "_",
                   which_dataset, "_", 
                   avg_method, "_avg_",
                   "FINAL")
saveRDS(conversion_plot, str_c(filename0, ".rds"))
ggsave(str_c(filename0, ".png"), conversion_plot, width = 8, height = 8, units = "in", dpi = 500, bg = "white")
