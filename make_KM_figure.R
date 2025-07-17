# Created by Kellen Petersen, July 1, 2025

library(survival)
library(ggplot2)
library(survminer)
library(dplyr)
library(cowplot)
library(splines)
library(stringr)
library(here)
library(ggdist)
library(patchwork)
library(RColorBrewer)

setwd("<SET WORKING DIRECTORY>")

conflicts_prefer(dplyr::select)

which_dataset <- "ADNI"
which_method <- "SILA"
symptoms_after_positivity <- FALSE
use_linear_model <- TRUE
cox_df_ns <- 3
scale_max_density <- 2
set_max_density <- 0.5
palette_colors <- c("#1565C0", "#00B8D4", "#00B050")
y_label_text <- "Probability of cognitively unimpaired"

if (which_dataset == "KADRC") {
  which_outcome <- "CDR_DX"
} else if (which_dataset == "ADNI") {
  which_outcome <- "CDR_DX_Imp"
} else {
  stop("Invalid dataset specified")
}

if (which_dataset == "KADRC" & which_method == "TIRA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else if (which_dataset == "ADNI" & which_method == "TIRA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>")) %>% 
    mutate(CDR_DX = CDR_DX_Imp)
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else if (which_dataset == "KADRC" & which_method == "SILA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else if (which_dataset == "ADNI" & which_method == "SILA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>")) %>% 
    mutate(CDR_DX = CDR_DX_Imp)
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else {
  stop("Invalid dataset/method specified")
}

symps <- symps %>% 
  mutate(
    CDR_OUT = case_when(
      CDR_ONLY == 0 & CDR_DX == 0 ~ 0,
      CDR_ONLY == 1 & CDR_DX == 0 ~ 1,
      CDR_ONLY == 1 & CDR_DX == 1 ~ 2,
      CDR_ONLY == 0 & CDR_DX == 1 ~ -99,
      TRUE ~ NA_real_
    )
  )

symps <- symps %>%
  filter(ID %in% df$ID, !is.na(AGE)) %>%
  mutate(ID = as.factor(ID),
         EXAMDATE = as.Date(EXAMDATE),
         CDR_DX = as.numeric(CDR_DX),
         CDR_OUT = as.numeric(CDR_OUT)) %>%
  arrange(ID, EXAMDATE) %>%
  group_by(ID) %>%
  mutate(
    last_age = max(AGE),
    age_at_symptom_onset_AD = if (any(CDR_OUT == 2, na.rm = TRUE) & first(CDR_OUT == 0)) {
      first(AGE[CDR_OUT == 2 & !is.na(CDR_DX)])
    } else {
      NA
    },
    age_at_censoring_AD = ifelse(is.na(age_at_symptom_onset_AD), last_age, NA),
    
    age_at_symptom_onset_nonAD = if (any(CDR_OUT == 1, na.rm = TRUE) & first(CDR_OUT == 0)) {
      first(AGE[CDR_OUT == 1 & !is.na(CDR_DX)])
    } else {
      NA
    },
    age_at_censoring_nonAD = ifelse(is.na(age_at_symptom_onset_nonAD), last_age, NA)
    
  ) %>%
  ungroup() %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

n_distinct(symps$ID[!is.na(symps$age_at_symptom_onset_AD)])
n_distinct(symps$ID[!is.na(symps$age_at_symptom_onset_nonAD)])

df$ID <- as.factor(df$ID)
df_slice <- df %>%
  group_by(ID) %>%
  slice(1) %>%
  dplyr::select(ID, est_onset_age, PTGENDER, PTEDUCAT, AGE, APOE4_positive) %>%
  mutate(first_plasma_age = first(AGE)) %>%
  ungroup() %>%
  filter(!is.na(est_onset_age))

survival_data <- df_slice %>%
  left_join(
    symps %>% select(ID, 
                     age_at_symptom_onset_AD, age_at_censoring_AD, 
                     age_at_symptom_onset_nonAD, age_at_censoring_nonAD,
                     last_age,CDR_OUT),
    by = "ID"
  ) %>%
  mutate(
    converted_AD_pos = ((!is.na(age_at_symptom_onset_AD)) & ((age_at_symptom_onset_AD - est_onset_age>=0))),
    time_to_event_AD_pos = ifelse(((converted_AD_pos)&(age_at_symptom_onset_AD - est_onset_age>=0)), 
                                  age_at_symptom_onset_AD - est_onset_age, 
                                  age_at_censoring_AD - est_onset_age)
  ) %>%
  mutate(
    converted_AD_neg = ((!is.na(age_at_symptom_onset_AD)) & ((age_at_symptom_onset_AD - est_onset_age<0))),
    time_to_event_AD_neg = ifelse(((converted_AD_neg)&(age_at_symptom_onset_AD - est_onset_age<0)), 
                                  age_at_symptom_onset_AD - est_onset_age, 
                                  age_at_censoring_AD - est_onset_age)
  ) %>%
  mutate(
    converted_nonAD = !is.na(age_at_symptom_onset_nonAD),
    time_to_event_nonAD = ifelse(converted_nonAD, 
                                 age_at_symptom_onset_nonAD - est_onset_age, 
                                 age_at_censoring_nonAD - est_onset_age)
  )

table(survival_data$converted_AD_pos)
table(survival_data$converted_AD_neg)
table(survival_data$converted_nonAD)

ids_AD_pos <- survival_data %>% 
  filter(converted_AD_pos) %>% 
  arrange(time_to_event_AD_pos) %>% 
  pull(ID)
ids_AD_neg <- survival_data %>% 
  filter(converted_AD_neg) %>% 
  arrange(time_to_event_AD_neg) %>% 
  pull(ID)
ids_nonAD <- survival_data %>% 
  filter(converted_nonAD) %>% 
  arrange(time_to_event_nonAD) %>% 
  pull(ID)

survival_data_AD_pos <- survival_data
survival_data_AD_neg <- survival_data
survival_data_nonAD <- survival_data

shift_x <- 100

survival_data_AD_pos <- survival_data_AD_pos %>%
  mutate(time_to_event_shift = time_to_event_AD_pos + shift_x)
km_fit_AD_pos <- survfit(Surv(time_to_event_shift, converted_AD_pos) ~ 1, 
                         data = survival_data_AD_pos)

survival_data_AD_neg <- survival_data_AD_neg %>%
  mutate(time_to_event_shift = time_to_event_AD_neg + shift_x)
km_fit_AD_neg <- survfit(Surv(time_to_event_shift, converted_AD_neg) ~ 1,
                         data = survival_data_AD_neg)

survival_data_nonAD <- survival_data_nonAD %>%
  mutate(time_to_event_shift = time_to_event_nonAD + shift_x)
km_fit_nonAD <- survfit(Surv(time_to_event_shift, converted_nonAD) ~ 1,
                        data = survival_data_nonAD)

fit_list <- list(
  nonAD = km_fit_nonAD,
  AD_neg = km_fit_AD_neg,
  AD_pos = km_fit_AD_pos
)

all_times <- c(survival_data_nonAD$time_to_event_shift,
               survival_data_AD_neg$time_to_event_shift,
               survival_data_AD_pos$time_to_event_shift)

tick_breaks <- seq(
  floor(min(all_times, na.rm = TRUE)/5)*5,
  ceiling(max(all_times, na.rm = TRUE)/5)*5,
  by = 5
)
xlim_range <- c(min(all_times, na.rm = TRUE),
                max(all_times, na.rm = TRUE))

p1_three <- ggsurvplot_combine(
  fit_list,
  data = survival_data,
  palette = c("darkorange", "purple4", "red3"),
  legend.labs = c(
    "Non-AD syndrome", 
    "AD syndrome/\nbiomarker negative",
    "AD syndrome/\nbiomarker positive"
  ),
  legend.title = "none",
  xlab = "Estimated years from %p-tau217 positivity",
  ylab = y_label_text,
  conf.int = TRUE,
  ggtheme = theme_cowplot() +
    theme(
      legend.position = c(0.12, 0.15),
      legend.justification = c(0, 0),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.box.background = element_rect(colour = "black")
    ),
  risk.table = FALSE,
  xlim = xlim_range
)

p1_three$plot <- p1_three$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  theme(legend.position = "none")

calculate_median <- function(fit) {
  km_summary <- summary(fit)
  median_index <- which(km_summary$surv <= 0.5)[1]
  if(is.na(median_index)) return(c(time = NA, lower = NA, upper = NA))
  c(
    time = km_summary$time[median_index] - shift_x,
    lower = km_summary$time[which(km_summary$surv <= km_summary$upper[median_index])[1]] - shift_x,
    upper = km_summary$time[which(km_summary$surv <= km_summary$lower[median_index])[1]] - shift_x
  )
}

med_AD_pos <- calculate_median(km_fit_AD_pos)
med_AD_neg <- calculate_median(km_fit_AD_neg)
med_nonAD <- calculate_median(km_fit_nonAD)

survival_data_AD_pos <- survival_data_AD_pos %>%
  mutate(
    converter_IDs_AD_pos = ifelse(ID %in% ids_AD_pos, 1, 0),
    converter_IDs_AD_neg = ifelse(ID %in% ids_AD_neg, 1, 0),
    converter_IDs_nonAD = ifelse(ID %in% ids_nonAD, 1, 0)
  )

survival_data_AD_neg <- survival_data_AD_neg %>%
  mutate(
    converter_IDs_AD_pos = ifelse(ID %in% ids_AD_pos, 1, 0),
    converter_IDs_AD_neg = ifelse(ID %in% ids_AD_neg, 1, 0),
    converter_IDs_nonAD = ifelse(ID %in% ids_nonAD, 1, 0)
  )

survival_data_nonAD <- survival_data_nonAD %>%
  mutate(
    converter_IDs_AD_pos = ifelse(ID %in% ids_AD_pos, 1, 0),
    converter_IDs_AD_neg = ifelse(ID %in% ids_AD_neg, 1, 0),
    converter_IDs_nonAD = ifelse(ID %in% ids_nonAD, 1, 0)
  )

survival_data2_AD_pos_y <- survival_data_AD_pos %>% 
  filter(converter_IDs_AD_pos == 1) %>% 
  arrange(time_to_event_AD_pos)

survival_data2_AD_neg_y <- survival_data_AD_neg %>% 
  filter(converter_IDs_AD_neg == 1) %>% 
  arrange(time_to_event_AD_neg)

survival_data2_nonAD_y <- survival_data_nonAD %>% 
  filter(converter_IDs_nonAD == 1) %>% 
  arrange(time_to_event_nonAD)

max_density_AD_pos <- scale_max_density^2 * max(density(survival_data2_AD_pos_y$time_to_event_shift)$y)
max_density_AD_neg <- scale_max_density^2 * max(density(survival_data2_AD_neg_y$time_to_event_shift)$y)
max_density_nonAD <- scale_max_density^2 * max(density(survival_data2_nonAD_y$time_to_event_shift)$y)

max_density_AD_pos <- set_max_density
max_density_AD_neg <- set_max_density
max_density_nonAD <- set_max_density

p1_three$plot <- p1_three$plot +
  geom_density(
    data = survival_data2_nonAD_y,
    aes(x = time_to_event_shift, y = ..density../max_density_nonAD),
    color = "darkorange", fill = "#FFDEAD", alpha = 0.7
  ) +
  geom_density(
    data = survival_data2_AD_neg_y,
    aes(x = time_to_event_shift, y = ..density../max_density_AD_neg),
    color = "purple4", fill = "#E6CCFF", alpha = 0.7
  ) +
  geom_density(
    data = survival_data2_AD_pos_y,
    aes(x = time_to_event_shift, y = ..density../max_density_AD_pos),
    color = "red3", fill = "#FF9999", alpha = 0.7
  ) +
  geom_point(data = survival_data2_AD_pos_y,
             aes(x = time_to_event_shift, y = -0.01), color = "red3", size = 1) +
  geom_point(data = survival_data2_AD_neg_y,
             aes(x = time_to_event_shift, y = -0.02), color = "purple4", size = 1) +
  geom_point(data = survival_data2_nonAD_y,
             aes(x = time_to_event_shift, y = -0.03), color = "darkorange", size = 1)

p1_three$plot <- p1_three$plot +
  scale_y_continuous(
    name = y_label_text,
    sec.axis = sec_axis(
      ~ . * max(c(max_density_AD_pos, max_density_AD_neg, max_density_nonAD)),
      name = "Event density"
    )
  )

p1_three$plot <- p1_three$plot +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth = 1) 

p1_three$plot <- p1_three$plot +
  scale_fill_manual(values = c("darkorange","purple4","red3")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p1_three$plot <- p1_three$plot +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.05, 0.3),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white")
  )

p1_three$plot <- p1_three$plot +
  guides(
    color = guide_legend(reverse = TRUE),
    fill = guide_legend(reverse = TRUE)
  )

p1_three$plot <- p1_three$plot + theme(legend.key.spacing.y = unit(3, 'mm'))

print(p1_three)

sila_coef <- read.csv(here("results","<LOAD DATA>"))
tira_coef <- read.csv(here("results","<LOAD DATA>"))

sila_coef$Dataset <- ifelse(sila_coef$dataset == "Knight ADRC", "KADRC", "ADNI")
tira_coef$Dataset <- ifelse(tira_coef$dataset == "Knight ADRC", "KADRC", "ADNI")

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

df_linear <- rbind(
  tira_df[tira_df$Dataset == "KADRC", ],
  sila_df[sila_df$Dataset == "KADRC", ],
  tira_df[tira_df$Dataset == "ADNI", ],
  sila_df[sila_df$Dataset == "ADNI", ]
)
rownames(df_linear) <- NULL
print(df_linear)

result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
model_intercept <- result_tmp$intercept
model_slope <- result_tmp$slope

survival_data3 <- survival_data %>%
  mutate(
    predicted_symptom_onset_age = model_intercept + model_slope * est_onset_age
  )

df_last_plot <- survival_data3 %>% 
  mutate(
    event_AD_pos = converted_AD_pos,
    years_since_predicted_symptoms_AD_pos = ifelse(converted_AD_pos, 
                                                   age_at_symptom_onset_AD - predicted_symptom_onset_age, 
                                                   age_at_censoring_AD - predicted_symptom_onset_age),
    
    event_AD_neg = converted_AD_neg,
    years_since_predicted_symptoms_AD_neg = ifelse(converted_AD_neg, 
                                                   age_at_symptom_onset_AD - predicted_symptom_onset_age, 
                                                   age_at_censoring_AD - predicted_symptom_onset_age),
    
    event_nonAD = converted_nonAD,
    years_since_predicted_symptoms_nonAD = ifelse(converted_nonAD, 
                                                  age_at_symptom_onset_nonAD - predicted_symptom_onset_age, 
                                                  age_at_censoring_nonAD - predicted_symptom_onset_age)
  )

shift_x <- 100
survival_data_shift <- df_last_plot %>%
  mutate(
    years_since_predicted_symptoms_AD_pos_shift = years_since_predicted_symptoms_AD_pos + shift_x,
    years_since_predicted_symptoms_AD_neg_shift = years_since_predicted_symptoms_AD_neg + shift_x,
    years_since_predicted_symptoms_nonAD_shift = years_since_predicted_symptoms_nonAD + shift_x
  )

survival_data_AD_pos <- survival_data_shift
survival_data_AD_neg <- survival_data_shift
survival_data_nonAD <- survival_data_shift 

km_fit_AD_pos <- survfit(Surv(years_since_predicted_symptoms_AD_pos_shift, event_AD_pos) ~ 1, 
                         data = survival_data_AD_pos)
km_fit_AD_neg <- survfit(Surv(years_since_predicted_symptoms_AD_neg_shift, event_AD_neg) ~ 1, 
                         data = survival_data_AD_neg)
km_fit_nonAD <- survfit(Surv(years_since_predicted_symptoms_nonAD_shift, event_nonAD) ~ 1, 
                        data = survival_data_nonAD)

all_years <- c(
  survival_data_AD_pos$years_since_predicted_symptoms_AD_pos_shift,
  survival_data_AD_neg$years_since_predicted_symptoms_AD_neg_shift,
  survival_data_nonAD$years_since_predicted_symptoms_nonAD_shift
)
xlim_range <- c(min(all_years, na.rm=TRUE), max(all_years, na.rm=TRUE))

km_plot <- ggsurvplot(
  list(nonAD = km_fit_nonAD, AD_neg = km_fit_AD_neg, AD_pos = km_fit_AD_pos),
  combine = TRUE,
  data = survival_data_shift,
  palette = c("darkorange", "purple4", "red3"),
  legend.labs = c(
    "Non-AD syndrome", 
    "AD syndrome/\nbiomarker negative",
    "AD syndrome/\nbiomarker positive"
  ),
  xlab = "Estimated years from symptom onset",
  ylab = y_label_text,
  conf.int = TRUE,
  risk.table = TRUE,
  xlim = xlim_range,
  break.x.by = 5,
  ggtheme = theme_cowplot()
)

survival_data2_AD_pos_y <- survival_data_AD_pos %>% 
  filter(event_AD_pos == 1) %>% 
  arrange(years_since_predicted_symptoms_AD_pos_shift)

survival_data2_AD_neg_y <- survival_data_AD_neg %>% 
  filter(event_AD_neg == 1) %>% 
  arrange(years_since_predicted_symptoms_AD_neg_shift)

survival_data2_nonAD_y <- survival_data_nonAD %>% 
  filter(event_nonAD == 1) %>% 
  arrange(years_since_predicted_symptoms_nonAD_shift)

xlim_lower_limit <- floor(min(all_times, na.rm = TRUE)/5)*5
xlim_upper_limit <- ceiling(max(all_times, na.rm = TRUE)/5)*5

km_main_plot <- km_plot$plot +
  scale_x_continuous(
    breaks = seq(xlim_lower_limit, xlim_upper_limit, 5),
    labels = seq(xlim_lower_limit, xlim_upper_limit, 5) - shift_x
  )

km_main_plot <- km_plot$plot +
  scale_x_continuous(
    breaks = seq(xlim_lower_limit, xlim_upper_limit, 5),
    labels = seq(xlim_lower_limit, xlim_upper_limit, 5) - shift_x
  )

max_density_AD_pos2 <- scale_max_density^2 * max(density(survival_data2_AD_pos_y$years_since_predicted_symptoms_AD_pos_shift)$y)
max_density_AD_neg2 <- scale_max_density^2 * max(density(survival_data2_AD_neg_y$years_since_predicted_symptoms_AD_neg_shift)$y)
max_density_nonAD2 <- scale_max_density^2 * max(density(survival_data2_nonAD_y$years_since_predicted_symptoms_nonAD_shift)$y)
max_density2 <- max(c(max_density_AD_pos2, max_density_AD_neg2, max_density_nonAD2))

set_max_density2 <- 1.2
max_density2 <- set_max_density2
max_density_AD_pos2 <- set_max_density2
max_density_AD_neg2 <- set_max_density2
max_density_nonAD2 <- set_max_density2

km_main_plot <- km_main_plot +
  geom_density(
    data = survival_data2_nonAD_y,
    aes(x = years_since_predicted_symptoms_nonAD_shift, y = ..density../max_density2),
    color = "darkorange", fill = "#FFDEAD", alpha = 0.7
  ) +
  geom_density(
    data = survival_data2_AD_neg_y,
    aes(x = years_since_predicted_symptoms_AD_neg_shift, y = ..density../max_density2),
    color = "purple4", fill = "#E6CCFF", alpha = 0.7
  ) +
  geom_density(
    data = survival_data2_AD_pos_y,
    aes(x = years_since_predicted_symptoms_AD_pos_shift, y = ..density../max_density2),
    color = "red3", fill = "#FF9999", alpha = 0.7
  )

km_main_plot <- km_main_plot +
  scale_y_continuous(
    name = y_label_text,
    sec.axis = sec_axis(~ . * max_density2, name = "Event density"),
    limits = c(-0.03, 1)
  )

km_main_plot <- km_main_plot +
  geom_point(data = survival_data2_AD_pos_y,
             aes(x = years_since_predicted_symptoms_AD_pos_shift, y = -0.01), 
             color = "red3", size = 1) +
  geom_point(data = survival_data2_AD_neg_y,
             aes(x = years_since_predicted_symptoms_AD_neg_shift, y = -0.02), 
             color = "purple4", size = 1) +
  geom_point(data = survival_data2_nonAD_y,
             aes(x = years_since_predicted_symptoms_nonAD_shift, y = -0.03), 
             color = "darkorange", size = 1)

km_main_plot <- km_main_plot +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth = 1)

km_main_plot <- km_main_plot +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.05, 0.3),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white"),
    legend.key.spacing.y = unit(3, 'mm')
  ) +
  guides(
    color = guide_legend(reverse = TRUE),
    fill = guide_legend(reverse = TRUE)
  )

survival_data <- survival_data %>%
  mutate(
    age_bin = case_when(
      est_onset_age < 70 ~ "< 70",
      est_onset_age >= 70 & est_onset_age < 80 ~ "70-80",
      est_onset_age >= 80 ~ "≥ 80",
      TRUE ~ NA_character_
    ),
    age_bin = factor(age_bin, levels = c("< 70", "70-80", "≥ 80"))
  )

survival_data_AD_pos_bin <- survival_data %>%
  mutate(time_to_event_shift = time_to_event_AD_pos + shift_x)

km_fit_by_bin <- survfit(Surv(time_to_event_shift, converted_AD_pos) ~ age_bin, 
                         data = survival_data_AD_pos_bin)

p_age_bins <- ggsurvplot(
  km_fit_by_bin,
  data = survival_data_AD_pos_bin,
  palette = palette_colors,
  xlab = "Estimated years from %p-tau217 positivity",
  ylab = y_label_text,
  title = "",
  conf.int = TRUE,
  legend.labs = c("< 70", "70-80", "≥ 80"),
  legend.title = "Age at \n%p-tau217 positivity",
  ggtheme = theme_cowplot(),
  xlim = c(100,130)
)

p_age_bins$plot <- p_age_bins$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  theme(legend.position = "right")

p_age_bins$plot <- p_age_bins$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth=1) +
  theme(
    legend.position = c(0.01, 0.1),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white")
  )

print(p_age_bins)

survival_data_AD_pos <- survival_data_shift %>%
  mutate(
    age_bin = case_when(
      est_onset_age < 70 ~ "< 70",
      est_onset_age >= 70 & est_onset_age < 80 ~ "70-80",
      est_onset_age >= 80 ~ "≥ 80",
      TRUE ~ NA_character_
    ),
    age_bin = factor(age_bin, levels = c("< 70", "70-80", "≥ 80"))
  )

km_fit_AD_pos_bins <- survfit(Surv(years_since_predicted_symptoms_AD_pos_shift, converted_AD_pos) ~ age_bin, 
                              data = survival_data_AD_pos)

fit_list_binned <- list(
  "AD_pos < 70" = km_fit_AD_pos_bins[1],
  "AD_pos 70-80" = km_fit_AD_pos_bins[2],
  "AD_pos ≥ 80" = km_fit_AD_pos_bins[3]
)

greens <- brewer.pal(3, "Greens")

km_main_binned <- ggsurvplot(
  fit_list_binned,
  combine = TRUE,
  data = survival_data_shift,
  palette = palette_colors,
  xlab = "Estimated years from symptom onset",
  ylab = y_label_text,
  conf.int = TRUE,
  ggtheme = theme_cowplot(),
  xlim = c(90,115),
  break.x.by = 5,
  legend.title = "Age at \n%p-tau217 positivity",
  legend.labs = c("< 70", "70-80", "≥ 80"),
)

km_main_binned$plot <- km_main_binned$plot +
  scale_x_continuous(
    breaks = seq(xlim_lower_limit, xlim_upper_limit, 5),
    labels = seq(xlim_lower_limit, xlim_upper_limit, 5) - shift_x
  ) +
  theme(legend.position = "right")

km_main_binned$plot <- km_main_binned$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth = 1) +
  theme(
    legend.position = c(0.01, 0.1),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white")
  )

print(km_main_binned)

if (which_dataset == "KADRC") {
  which_outcome <- "CDR_DX"
} else if (which_dataset == "ADNI") {
  which_outcome <- "CDR_DX_Imp"
} else {
  stop("Invalid dataset specified")
}

if (which_dataset == "KADRC" & which_method == "TIRA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else if (which_dataset == "ADNI" & which_method == "TIRA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>")) %>% 
    mutate(CDR_DX = CDR_DX_Imp)
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else if (which_dataset == "KADRC" & which_method == "SILA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else if (which_dataset == "ADNI" & which_method == "SILA") {
  df <- read.csv(here::here("results","<LOAD DATA>"))
  symps <- read.csv(here::here("data_final","<LOAD DATA>")) %>% 
    mutate(CDR_DX = CDR_DX_Imp)
  ids <- read.csv(here::here("data_final","<LOAD DATA>"))
  ids_sens <- read.csv(here::here("data_final","<LOAD DATA>"))
} else {
  stop("Invalid dataset/method specified")
}

symps <- symps %>% 
  mutate(
    CDR_OUT = case_when(
      CDR_ONLY == 0 & CDR_DX == 0 ~ 0,
      CDR_ONLY == 1 & CDR_DX == 0 ~ 1,
      CDR_ONLY == 1 & CDR_DX == 1 ~ 2,
      CDR_ONLY == 0 & CDR_DX == 1 ~ -99,
      TRUE ~ NA_real_
    )
  )

symps <- symps %>%
  filter(ID %in% df$ID, !is.na(AGE)) %>%
  mutate(ID = as.factor(ID),
         EXAMDATE = as.Date(EXAMDATE),
         CDR_DX = as.numeric(CDR_DX),
         CDR_OUT = as.numeric(CDR_OUT)) %>%
  arrange(ID, EXAMDATE) %>%
  group_by(ID) %>%
  mutate(
    last_age = max(AGE),
    age_at_symptom_onset_AD = if (any(CDR_OUT == 2, na.rm = TRUE) & first(CDR_OUT == 0)) {
      first(AGE[CDR_OUT == 2 & !is.na(CDR_DX)])
    } else {
      NA
    },
    age_at_censoring_AD = ifelse(is.na(age_at_symptom_onset_AD), last_age, NA),
    
    age_at_symptom_onset_nonAD = if (any(CDR_OUT == 1, na.rm = TRUE) & first(CDR_OUT == 0)) {
      first(AGE[CDR_OUT == 1 & !is.na(CDR_DX)])
    } else {
      NA
    },
    age_at_censoring_nonAD = ifelse(is.na(age_at_symptom_onset_nonAD), last_age, NA)
    
  ) %>%
  ungroup() %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

n_distinct(symps$ID[!is.na(symps$age_at_symptom_onset_AD)])
n_distinct(symps$ID[!is.na(symps$age_at_symptom_onset_nonAD)])

df$ID <- as.factor(df$ID)
df_slice <- df %>%
  group_by(ID) %>%
  slice(1) %>%
  dplyr::select(ID, plasma, est_onset_age, PTGENDER, PTEDUCAT, AGE, APOE4_positive) %>%
  mutate(first_plasma_age = first(AGE)) %>%
  ungroup() %>%
  filter(!is.na(est_onset_age))

survival_data <- df_slice %>%
  left_join(
    symps %>% select(ID, 
                     age_at_symptom_onset_AD, age_at_censoring_AD, 
                     age_at_symptom_onset_nonAD, age_at_censoring_nonAD,
                     last_age,CDR_OUT),
    by = "ID"
  ) %>%
  mutate(
    converted_AD_pos = ((!is.na(age_at_symptom_onset_AD)) & ((age_at_symptom_onset_AD - est_onset_age>=0))),
    time_to_event_AD_pos = ifelse(((converted_AD_pos)&(age_at_symptom_onset_AD - est_onset_age>=0)), 
                                  age_at_symptom_onset_AD - first_plasma_age, 
                                  age_at_censoring_AD - first_plasma_age)
  ) %>%
  mutate(
    converted_AD_neg = ((!is.na(age_at_symptom_onset_AD)) & ((age_at_symptom_onset_AD - est_onset_age<0))),
    time_to_event_AD_neg = ifelse(((converted_AD_neg)&(age_at_symptom_onset_AD - est_onset_age<0)), 
                                  age_at_symptom_onset_AD - first_plasma_age, 
                                  age_at_censoring_AD - first_plasma_age)
  ) %>%
  mutate(
    converted_nonAD = !is.na(age_at_symptom_onset_nonAD),
    time_to_event_nonAD = ifelse(converted_nonAD, 
                                 age_at_symptom_onset_nonAD - first_plasma_age, 
                                 age_at_censoring_nonAD - first_plasma_age)
  )

survival_data <- survival_data %>% 
  mutate(diff = age_at_symptom_onset_AD - est_onset_age)

filter_bl_CDRge0 <- TRUE
survival_data0 <- survival_data
if (filter_bl_CDRge0) {
  survival_data <- survival_data0 %>%
    filter(plasma >= 4.06)
}

ids_AD_pos <- survival_data %>% 
  filter(converted_AD_pos) %>% 
  arrange(time_to_event_AD_pos) %>% 
  pull(ID)
ids_AD_neg <- survival_data %>% 
  filter(converted_AD_neg) %>% 
  arrange(time_to_event_AD_neg) %>% 
  pull(ID)
ids_nonAD <- survival_data %>% 
  filter(converted_nonAD) %>% 
  arrange(time_to_event_nonAD) %>% 
  pull(ID)

survival_data_AD_pos <- survival_data
survival_data_AD_neg <- survival_data
survival_data_nonAD <- survival_data

shift_x <- 100

survival_data_AD_pos <- survival_data_AD_pos %>%
  mutate(time_to_event_shift = time_to_event_AD_pos + shift_x)
km_fit_AD_pos <- survfit(Surv(time_to_event_shift, converted_AD_pos) ~ 1, 
                         data = survival_data_AD_pos)

survival_data_AD_neg <- survival_data_AD_neg %>%
  mutate(time_to_event_shift = time_to_event_AD_neg + shift_x)

survival_data_nonAD <- survival_data_nonAD %>%
  mutate(time_to_event_shift = time_to_event_nonAD + shift_x)
km_fit_nonAD <- survfit(Surv(time_to_event_shift, converted_nonAD) ~ 1,
                        data = survival_data_nonAD)

fit_list <- list(
  nonAD = km_fit_nonAD,
  AD_pos = km_fit_AD_pos
)

all_times <- c(survival_data_nonAD$time_to_event_shift,
               survival_data_AD_pos$time_to_event_shift)

tick_breaks <- seq(
  floor(min(all_times, na.rm = TRUE)/5)*5,
  ceiling(max(all_times, na.rm = TRUE)/5)*5,
  by = 5
)
xlim_range <- c(min(all_times, na.rm = TRUE),
                max(all_times, na.rm = TRUE))

p1_baseline <- ggsurvplot_combine(
  fit_list,
  data = survival_data,
  palette = c("darkorange", "red3"),
  legend.labs = c(
    "AD syndrome/\nbiomarker negative",
    "AD syndrome/\nbiomarker positive"
  ),
  xlab = "Time from baseline %p-tau217 collection (years)",
  ylab = y_label_text,
  conf.int = TRUE,
  ggtheme = theme_cowplot(),
  risk.table = FALSE,
  xlim = xlim_range
)

p1_baseline$plot <- p1_baseline$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  theme(legend.position = "none")

calculate_median <- function(fit) {
  km_summary <- summary(fit)
  median_index <- which(km_summary$surv <= 0.5)[1]
  if(is.na(median_index)) return(c(time = NA, lower = NA, upper = NA))
  c(
    time = km_summary$time[median_index] - shift_x,
    lower = km_summary$time[which(km_summary$surv <= km_summary$upper[median_index])[1]] - shift_x,
    upper = km_summary$time[which(km_summary$surv <= km_summary$lower[median_index])[1]] - shift_x
  )
}

med_AD_pos <- calculate_median(km_fit_AD_pos)
med_nonAD <- calculate_median(km_fit_nonAD)

survival_data_AD_pos <- survival_data_AD_pos %>%
  mutate(
    converter_IDs_AD_pos = ifelse(ID %in% ids_AD_pos, 1, 0),
    converter_IDs_nonAD = ifelse(ID %in% ids_nonAD, 1, 0)
  )

survival_data_AD_neg <- survival_data_AD_neg %>%
  mutate(
    converter_IDs_AD_pos = ifelse(ID %in% ids_AD_pos, 1, 0),
    converter_IDs_nonAD = ifelse(ID %in% ids_nonAD, 1, 0)
  )

survival_data_nonAD <- survival_data_nonAD %>%
  mutate(
    converter_IDs_AD_pos = ifelse(ID %in% ids_AD_pos, 1, 0),
    converter_IDs_nonAD = ifelse(ID %in% ids_nonAD, 1, 0)
  )

survival_data2_AD_pos_y <- survival_data_AD_pos %>% 
  filter(converter_IDs_AD_pos == 1) %>% 
  arrange(time_to_event_AD_pos)

survival_data2_nonAD_y <- survival_data_nonAD %>% 
  filter(converter_IDs_nonAD == 1) %>% 
  arrange(time_to_event_nonAD)

max_density_AD_pos <- scale_max_density^2 * max(density(survival_data2_AD_pos_y$time_to_event_shift)$y)
max_density_nonAD <- scale_max_density^2 * max(density(survival_data2_nonAD_y$time_to_event_shift)$y)

max_density_AD_pos <- set_max_density
max_density_AD_neg <- set_max_density

p1_baseline$plot <- p1_baseline$plot +
  geom_density(
    data = survival_data2_nonAD_y,
    aes(x = time_to_event_shift, y = ..density../max_density_nonAD),
    color = "darkorange", fill = "#FFDEAD", alpha = 0.7
  ) +
  geom_density(
    data = survival_data2_AD_pos_y,
    aes(x = time_to_event_shift, y = ..density../max_density_AD_pos),
    color = "red3", fill = "#FF9999", alpha = 0.7
  ) +
  geom_point(data = survival_data2_AD_pos_y,
             aes(x = time_to_event_shift, y = -0.01), color = "red3", size = 1) +
  geom_point(data = survival_data2_nonAD_y,
             aes(x = time_to_event_shift, y = -0.03), color = "darkorange", size = 1)

p1_baseline$plot <- p1_baseline$plot +
  scale_y_continuous(
    name = y_label_text,
    sec.axis = sec_axis(
      ~ . * max(c(max_density_AD_pos, max_density_nonAD)),
      name = "Event density"
    )
  )

p1_baseline$plot <- p1_baseline$plot +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth = 1) 

p1_baseline$plot <- p1_baseline$plot +
  scale_fill_manual(values = c("darkorange","red3")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p1_baseline <- p1_baseline$plot 

p1_baseline <- p1_baseline +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.05, 0.3),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white"),
    legend.key.spacing.y = unit(3, 'mm')
  ) +
  guides(
    color = guide_legend(reverse = TRUE),
    fill = guide_legend(reverse = TRUE)
  )

print(p1_baseline)

survival_data <- survival_data %>%
  mutate(
    age_bin = case_when(
      est_onset_age < 70 ~ "< 70",
      est_onset_age >= 70 & est_onset_age < 80 ~ "70-80",
      est_onset_age >= 80 ~ "≥ 80",
      TRUE ~ NA_character_
    ),
    age_bin = factor(age_bin, levels = c("< 70", "70-80", "≥ 80"))
  )

survival_data_AD_pos_bin <- survival_data %>%
  mutate(time_to_event_shift = time_to_event_AD_pos + shift_x)

km_fit_by_bin <- survfit(Surv(time_to_event_shift, converted_AD_pos) ~ age_bin, 
                         data = survival_data_AD_pos_bin)

p_age_bins_baseline <- ggsurvplot(
  km_fit_by_bin,
  data = survival_data_AD_pos_bin,
  palette = palette_colors,
  xlab = "Time from baseline %p-tau217 collection (years)",
  ylab = y_label_text,
  title = "",
  conf.int = TRUE,
  legend.labs = c("< 70", "70-80", "≥ 80"),
  legend.title = "Age at \n%p-tau217 positivity",
  ggtheme = theme_cowplot(),
  xlim = xlim_range
)

p_age_bins_baseline$plot <- p_age_bins_baseline$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth = 1) +
  theme(
    legend.position = c(0.01, 0.1),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white")
  )

print(p_age_bins_baseline)

p_age_bins$plot <- p_age_bins$plot +
  scale_x_continuous(
    breaks = tick_breaks,
    labels = tick_breaks - shift_x
  ) +
  geom_vline(xintercept = shift_x, linetype = "dashed", color = "black", linewidth=1) +
  theme(
    legend.position = c(0.01, 0.1),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "white", linewidth = 0.3),
    legend.box.background = element_rect(color = "white")
  )

pp <- (p1_baseline + p1_three$plot + km_main_plot) / (p_age_bins_baseline$plot + p_age_bins$plot + km_main_binned$plot)
pp <- pp +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold")) +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(hjust = 0, vjust = 1))

filename0 <- str_c("results/figure_KM_",
                   which_method, "_",
                   which_dataset, "_", 
                   format(Sys.Date(), "%Y-%m-%d"), "_",
                   format(Sys.time(), "%H-%M-%S"))
ggsave(
  filename = str_c(filename0,".png"),
  plot = pp,
  width = 20, height = 12, units = "in",
  dpi = 500
)
