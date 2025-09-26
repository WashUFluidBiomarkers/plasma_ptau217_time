# Originally created by Kellen Petersen, July 1, 2025
# Updated by Kellen Petersen, September 12, 2025

# Load required libraries for data manipulation and visualization
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(here)

setwd("")

# Set analysis parameters for dataset, method, and outcome
which_dataset <- "ADNI"
which_method <- "SILA"
which_outcome <- "CDR_NEW"

filter_CDR_0_at_bl <- FALSE

# Load linear model coefficients for both SILA and TIRA methods
sila_coef <- read.csv(here("results",""))
tira_coef <- read.csv(here("results",""))

# Standardize dataset naming conventions
sila_coef$Dataset <- ifelse(sila_coef$dataset == "Knight ADRC", "KADRC", "ADNI")
tira_coef$Dataset <- ifelse(tira_coef$dataset == "Knight ADRC", "KADRC", "ADNI")

# Create structured data frames for model parameters
sila_df <- data.frame(
  Model = "SILA",
  Dataset = sila_coef$Dataset,
  intercept = sila_coef$intercept,
  slope = sila_coef$slope
)

tira_df <- data.frame(
  Model = "TIRA",
  Dataset = tira_coef$Dataset,
  intercept = tira_coef$intercept,
  slope = tira_coef$slope
)

# Combine coefficient data for both methods and datasets
df_linear <- rbind(
  tira_df[tira_df$Dataset == "KADRC", ],
  sila_df[sila_df$Dataset == "KADRC", ],
  tira_df[tira_df$Dataset == "ADNI", ],
  sila_df[sila_df$Dataset == "ADNI", ]
)

rownames(df_linear) <- NULL
print(df_linear)

# Load appropriate dataset and symptom files based on method and cohort
if (which_dataset == "ADNI" & which_method == "TIRA") {
  df <- read.csv(here("results",""))
  symps <- read.csv(here("data_final","")) %>%
    mutate(CDR_DX = CDR_DX_Imp)
  
  result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
  model_intercept <- result_tmp$intercept
  model_slope <- result_tmp$slope
  
} else if (which_dataset == "KADRC" & which_method == "TIRA") {
  df <- read.csv(here("results",""))
  symps <- read.csv(here("data_final",""))
  
  result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
  model_intercept <- result_tmp$intercept
  model_slope <- result_tmp$slope
  
} else if (which_dataset == "ADNI" & which_method == "SILA") {
  df <- read.csv(here("results",""))
  symps <- read.csv(here("data_final","")) %>%
    mutate(CDR_DX = CDR_DX_Imp)
  
  result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
  model_intercept <- result_tmp$intercept
  model_slope <- result_tmp$slope
  
} else if (which_dataset == "KADRC" & which_method == "SILA") {
  df <- read.csv(here("results",""))
  symps <- read.csv(here("data_final",""))
  
  result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
  model_intercept <- result_tmp$intercept
  model_slope <- result_tmp$slope
  
} else {
  stop("Invalid data type specified")
}

# Create composite outcome variable combining CDR and diagnosis
symps <- symps %>%
  mutate(
    CDR_NEW = case_when(
      CDR_ONLY == 0 & CDR_DX == 0 ~ 0,
      CDR_ONLY == 1 & CDR_DX == 0 ~ 1,
      CDR_ONLY == 1 & CDR_DX == 1 ~ 2,
      CDR_ONLY == 0 & CDR_DX == 1 ~ -99,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(outcome := .data[[which_outcome]])

# Extract onset age estimates (one per participant)
df <- df %>%
  select(ID, est_onset_age, one_pos) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

# Merge biomarker and symptom data
combo <- df %>%
  left_join(symps, by = "ID") %>%
  mutate(years_since_onset = AGE - est_onset_age) %>%
  arrange(est_onset_age, -ID) %>%
  mutate(est_age_symptom = model_slope * est_onset_age + model_intercept,
         line_values = est_age_symptom - est_onset_age) %>%
  filter(!is.na(CDR_NEW)) %>%
  filter(!is.na(years_since_onset))

table(combo$CDR_NEW)

# Recode outcome for cleaner visualization groupings
combo <- combo %>%
  mutate(CDR_NEW_0 = CDR_NEW) %>%
  mutate(CDR_NEW = ifelse(CDR_NEW==0,0,
                          ifelse(CDR_NEW==1,1,
                                 ifelse(((CDR_NEW==2)&(years_since_onset>=0)),3,2))))

# Add baseline CDR information for optional filtering
combo <- combo %>%
  group_by(ID) %>%
  mutate(first_CDR = first(CDR_ONLY)) %>%
  ungroup()

# Apply baseline filter if specified
if(filter_CDR_0_at_bl) {
  combo <- combo %>%
    filter(first_CDR == 0)
}

# Keep only participants with multiple observations
combo <- combo %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  ungroup()

n_distinct(combo$ID)

# Create age-based labels for y-axis using representative participants
target_ages <- c(40, 50, 60, 70, 80, 90, 100)
label_ids <- combo %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    interval_idx = findInterval(est_onset_age, target_ages),
    target_age = as.numeric(sapply(interval_idx, function(x) target_ages[x])),
    diff = abs(as.numeric(est_onset_age) - target_age)
  ) %>%
  group_by(target_age) %>%
  slice_min(diff, n = 1) %>%
  ungroup() %>%
  na.omit() %>%
  group_by(target_age) %>%
  slice(1) %>%
  ungroup()

title_label = ""
combo$outcome <- as.factor(combo$CDR_NEW)

# Create raster plot showing cognitive trajectories
# Each row represents a participant, ordered by estimated age at biomarker positivity
# X-axis shows years relative to biomarker positivity
# Brown line shows predicted symptom onset for each participant
dot_plot <- ggplot(combo,
                   aes(x = years_since_onset,
                       y = reorder(factor(ID), est_onset_age))) +
  theme_minimal() +
  geom_point(aes(color = factor(outcome)),
             size = 1,
             alpha = 1) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "black",
             alpha = 1,
             linewidth=1) +
  geom_line(aes(x = line_values,
                y = reorder(factor(ID), -est_onset_age)),
            color = "brown",
            size = 2.5) +
  geom_point(aes(x = line_values,
                 y = reorder(factor(ID), -est_onset_age)),
             color = "brown",
             size = .8) +
  scale_y_discrete(
    breaks = label_ids$ID,
    labels = label_ids$target_age
  ) +
  coord_cartesian(ylim = c(1, length(unique(combo$ID)) + 5)) +
  scale_x_continuous(
    limits = c(-35, 35),
    expand = c(0, 0)
  ) +
  scale_color_manual(values = c("#66B2FF", "darkorange","purple4", "red"),
                     labels = c("Cognitively Normal",
                                "Other Dementia Syndrome",
                                "Mild Cognitive Impairment",
                                "Typical AD Dementia Syndrome")) +
  labs(title = "",
       x = "Estimated years from %p-tau217 positivity",
       y = "Estimated age at %p-tau217 positivity") +
  cowplot::theme_cowplot() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0, "cm")
  )

dot_plot

# Generate output filename with analysis parameters
filename0 <- str_c("results/dot_plot_",
                   which_method, "_",
                   which_dataset, "_",
                   which_outcome, "_", "FINAL")

if (filter_CDR_0_at_bl) {
  filename0 <- str_c(filename0, "_filter_CDR_0_at_bl")
} else {
  filename0 <- str_c(filename0, "_no_filter_CDR_0_at_bl")
}

# Save plot in multiple formats
saveRDS(dot_plot, str_c(filename0, ".rds"))
ggsave(str_c(filename0, ".png"), dot_plot, width = 5, height = 10, units = "in", dpi = 500, bg = "white")
