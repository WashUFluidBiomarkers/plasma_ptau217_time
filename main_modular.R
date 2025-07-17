# Created by Kellen Petersen, July 1, 2025

library(tidyverse)  # Data manipulation and visualization
library(cowplot)    # Publication-ready plots
library(here)       # Project-relative file paths
library(nlme)       # Linear mixed-effects models
library(lme4)       # Mixed-effects models
library(mgcv)       # Generalized additive models

setwd("<SET WORKING DIRECTORY>")
i_am("<SET I_AM>")
here()
source("main_functions.R")

which_data <- "ADNI"
choose_plasma <- "C2N_plasma_ptau217_ratio"

year <- 365.25
dx <- 0.001
tip_point <- 4.06
tip_point2 <- tip_point + dx/2
which_model <- 1 # 1 = LME, 2 = LinReg

clock_parameters = list(year = year,
                        dx = dx,
                        tip_point = tip_point,
                        tip_point2 = tip_point2,
                        which_model = which_model)

if (which_data == "KADRC") {
  print("Loading KADRC data")
  df0 <- read.csv(here("data_final","<LOAD DATA>"), header = TRUE)
  dataset_name <- which_data
} else if (which_data == "ADNI") {
  print("Loading ADNI data")
  df0 <- read.csv(here("data_final","<LOAD DATA>"), header = TRUE)
  dataset_name <- which_data
} else {
  print("Error: Choose correct dataset")
}
n_distinct(df0$ID)

pre_process_data_out <- pre_process_data(df0, choose_plasma)
df <- pre_process_data_out$df

clock_out <- make_clock(df, clock_parameters)
df<- clock_out$df
plasma_midpoint <- clock_out$plasma_midpoint
model_time <- clock_out$model_time

df <- df %>% 
  mutate(plasma_time2 = predict(model_time, newdata = ., type = "response"))

df <- df %>% 
  mutate(est_onset_age_i = AGE - plasma_time) %>% 
  group_by(ID) %>%
  mutate( est_onset_age = ifelse(n() == 1, est_onset_age_i, mean(est_onset_age_i)),
          years_since_onset = AGE - est_onset_age) %>% 
  ungroup() %>%
  arrange(ID, EXAMDATE)

df <- df %>% select(ID, EXAMDATE, time, fu_time, 
                    AGE, PTGENDER, PTEDUCAT,
                    plasma, plasma_time, plasma_time2, 
                    est_onset_age_i,est_onset_age, years_since_onset, 
                    one_pos, is_converter,
                    CDR, CDR_SOB, CDR_10)

filename_df <- str_c("results/", "df_", dataset_name, "_", 
                     format(Sys.Date(), "%Y-%m-%d"), "_", 
                     format(Sys.time(), "%H-%M-%S"), ".csv")
write.csv(df, filename_df, row.names = FALSE)

df_clock <- plasma_midpoint %>% 
  select(plasma_midpoint, plasma_time) %>% 
  rename(plasma = plasma_midpoint, 
         clock_time = plasma_time)

filename_df_clock <- str_c("results/", "df_clock_", dataset_name, "_", 
                           format(Sys.Date(), "%Y-%m-%d"), "_", 
                           format(Sys.time(), "%H-%M-%S"), ".csv")
write.csv(df_clock, filename_df_clock, row.names = FALSE)

y_limits <- c(0,12)
y_label <- "Plasma %ptau217 (pg/mL)"

p1 <- make_spaghetti_plot(
  data = df,
  x_var = "AGE",
  y_var = "plasma",
  group_var = "ID",
  color_var = "one_pos",
  y_lim = y_limits,
  tip_point = tip_point,
  title = "",
  x_label = "Age (years)",
  y_label = y_label
)

p2 <- make_spaghetti_plot(
  data = df,
  x_var = "years_since_onset",
  y_var = "plasma",
  group_var = "ID",
  color_var = "one_pos",
  y_lim = y_limits,
  tip_point = tip_point,
  title = "",
  x_label = "TIRA time (years)",
  y_label = y_label,
  add_line = TRUE,
  line_data = plasma_midpoint,
  line_x = "plasma_time",
  line_y = "plasma_midpoint"
)

combined_spaghetti_plot <- p1 + p2
combined_spaghetti_plot

intervals_table <- create_clock_intervals(
  midpoint_data = plasma_midpoint,
  biomarker_col = "plasma_midpoint",
  time_col = "plasma_time",
  start_value = 2.06,
  end_value = 19.06,
  by_value = 1,
  reference_value = 4.06)
intervals_table

filename_clock_table <- str_c("results/","clock_intervals_", dataset_name, "_", 
                              format(Sys.Date(), "%Y-%m-%d"), "_", 
                              format(Sys.time(), "%H-%M-%S"), ".csv")
write.csv(intervals_table, filename_clock_table, row.names = FALSE)
