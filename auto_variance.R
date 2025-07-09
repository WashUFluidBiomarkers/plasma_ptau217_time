# Created by Kellen Petersen, July 1, 2025

auto_variance <- function(which_biomarker, 
                          df0,
                          choose_plasma,
                          variance_cut = 0.95, 
                          tip_point) {
  
  library(tidyverse)
  library(here)
  library(nlme)
  library(mgcv)
  library(cowplot)

  setwd("<SET WORKING DIRECTORY>")
  i_am("<SET I_AM>")
  here()
  source("main_functions.R")
  
  year <- 365.25
  dx <- 0.001
  tip_point2 <- tip_point + dx/2
  time_since_last_cut <- 1
  which_model <- 1 # 1 = LME, 2 = LinReg
  
  clock_parameters = list(tip_point = tip_point)
  
  df0$ID <- as.factor(df0$ID)
  df0 <- df0 %>% 
    select(ID, EXAMDATE, C2N_plasma_ptau217_ratio, AGE) %>% 
    na.omit()
  df0$EXAMDATE <- as.Date(df0$EXAMDATE, format = "%Y-%m-%d")
  df <- df0 %>% 
    relocate(EXAMDATE, .after = ID)
  
  df <- df %>% 
    group_by(ID) %>% 
    mutate(
      time_days = EXAMDATE - first(EXAMDATE),
      time = as.numeric(time_days, units = "days") / year,
      fu_time = as.numeric(max(time_days)) / year
    ) %>% 
    relocate(time, .after = EXAMDATE)
  
  df <- df %>% 
    group_by(ID) %>% 
    mutate(
      time_since_last = time - dplyr::lag(time, default = NA)
    ) %>% 
    ungroup() %>% 
    relocate(time_since_last, .after = time)
  
  df <- df %>% 
    filter(is.na(time_since_last) | time_since_last > time_since_last_cut)
  
  df <- df %>% 
    group_by(ID) %>% 
    filter(n() > 1) %>% 
    ungroup()
  
  n_distinct(df$ID)
  
  df$plasma <- df$C2N_plasma_ptau217_ratio
  df$plasma_YN <- ifelse(df$plasma >= tip_point, 1, 0)
  df$plasma_YN <- factor(df$plasma_YN, levels=c(0,1), labels = c("Negative", "Positive"))
  
  df <- df %>%
    arrange(ID,EXAMDATE) %>% 
    group_by(ID) %>% 
    mutate(Previous_plasma_YN = dplyr::lag(plasma_YN))
  conversion_plasma_YN <- df %>%
    filter(Previous_plasma_YN  == "Negative", plasma_YN == "Positive") %>%
    group_by(ID) %>% 
    distinct(ID)
  converters_subset <- df %>%
    semi_join(conversion_plasma_YN, by = "ID")
  
  midpoint_ages <- converters_subset %>%
    group_by(ID) %>%
    summarize(
      Last_0_Age = dplyr::last(AGE[plasma_YN == "Negative"]),
      First_1_Age = first(AGE[plasma_YN == "Positive"]),
      Midpoint_Age = (Last_0_Age + First_1_Age) / 2
    )
  df <- merge(df, midpoint_ages, by = "ID", all.x = TRUE)
  df$plasma_time_conv <- (df$AGE - df$Midpoint_Age)
  df <- df %>% 
    arrange(ID, EXAMDATE)
  
  fit_SUVRrate <- lme(plasma ~ time,
                      random = ~1 + time | ID,
                      data = df,
                      na.action = na.omit,
                      method = "ML",
                      control = lmeControl(opt = "optim"))
  
  slopes <- coef(fit_SUVRrate)$time
  slopes_df <- data.frame(ID = unique(df$ID), slopes = slopes)
  df <- left_join(df, slopes_df, by = "ID")  
  
  df <- df %>%
    group_by(ID) %>%
    mutate(
      Baseline_SUVR = first(plasma),
      SUVR_midpoint = ((fu_time / 2) * slopes) + Baseline_SUVR
    ) %>%
    mutate(slopes_over_midpoint = slopes / SUVR_midpoint) %>% 
    ungroup()
  
  df2 <- df %>% distinct(ID, .keep_all = TRUE)
  
  rate_plot <- ggplot(df2, aes(x = SUVR_midpoint, y = slopes)) + 
    geom_point(color = "grey20", alpha = 0.6) + 
    geom_smooth(method = "gam", color = "black") + 
    theme_cowplot() +
    labs(
      x = "Estimated biomarker at midpoint", 
      y = "Biomarker rate of change"
    )
  df2_later <- df2
  
  gam_model <- gam(slopes ~ s(SUVR_midpoint), data = df2)
  new_data <- data.frame(
    SUVR_midpoint = seq(min(df2$SUVR_midpoint), max(df2$SUVR_midpoint), length = 1000)
  )
  new_data$predicted <- predict(gam_model, new_data)
  zero_crossing <- approxfun(new_data$SUVR_midpoint, new_data$predicted)
  
  fit_SUVRrate <- gam(slopes ~ s(SUVR_midpoint, bs = "cr"), data = df2)
  
  df2 <- df2[order(df2$SUVR_midpoint), ]
  predictions <- predict(fit_SUVRrate, newdata = df2, se.fit = TRUE)
  df2$variance <- predictions$se.fit^2
  
  variance <- df2$variance
  mean_variance <- mean(variance, na.rm = TRUE)
  sd_variance <- sd(variance, na.rm = TRUE)
  
  cutpoint <- quantile(variance, variance_cut) 
  low_variance_points <- df2$SUVR_midpoint[variance <= cutpoint]
  high_variance_points <- df2$SUVR_midpoint[variance > cutpoint]
  options(scipen = 999)
  
  df2$variance <- as.numeric(df2$variance)
  p_variance <- ggplot() +
    geom_line(data = data.frame(x = df2$SUVR_midpoint, y = variance), 
              aes(x = x, y = y), linewidth = 1) +
    geom_hline(yintercept = cutpoint, color = "red", linetype = "dashed") +
    geom_point(data = data.frame(x = low_variance_points, 
                                 y = variance[variance <= cutpoint]),
               aes(x = x, y = y), color = "blue", size = 3) +
    geom_point(data = data.frame(x = high_variance_points, 
                                 y = variance[variance > cutpoint]),
               aes(x = x, y = y), color = "red", size = 3) +
    geom_vline(xintercept = c(min(low_variance_points), max(low_variance_points)), 
               color = "red", linetype = "dashed") +
    theme_cowplot() +
    labs(x = str_c("Estimated biomarker at midpoint"),
         y = "Variance in rate of change")
  p_variance
  
  min_lv <- min(low_variance_points)
  max_lv <- max(low_variance_points)
  min_lv <-floor(min_lv*100)/100
  max_lv <-ceiling(max_lv*100)/100
  cut_points <- c(min_lv,max_lv)
  cut_points
  
  rate_plot <- rate_plot + 
    geom_vline(xintercept = cut_points, linetype = "dashed", color = "red")
  
  pp <- rate_plot | p_variance
  
  return(list(cut_points = cut_points, rate_plot = rate_plot, p_variance = p_variance))
}
