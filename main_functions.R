# Created by Kellen Petersen, July 1, 2025

pre_process_data <- function(df0, choose_plasma) {
  year <- 365.25
  
  df <- df0 %>% 
    select(ID, EXAMDATE, all_of(choose_plasma), 
           AGE, PTGENDER, PTEDUCAT, 
           CDR, CDR_SOB, CDR_10) %>% 
    rename(plasma = all_of(choose_plasma)) %>%
    arrange(ID, EXAMDATE)
  head(df)
  
  df <- df %>% 
    mutate(
      ID = as.factor(ID),
      EXAMDATE = as.Date(EXAMDATE, format = "%Y-%m-%d"),
      PTGENDER = as.factor(PTGENDER),
      PTEDUCAT = as.numeric(PTEDUCAT),
      CDR_10 = as.factor(CDR_10)
    )
  
  df <- df %>%
    group_by(ID) %>%
    mutate(
      time_days = EXAMDATE - first(EXAMDATE),
      time = as.numeric(time_days, units = "days") / year,
      fu_time = as.numeric(max(time_days)) / year
    ) %>%
    relocate(time, .after = EXAMDATE) %>%
    group_by(ID) %>%
    mutate(time_since_last = time - dplyr::lag(time, default = NA)) %>%
    ungroup() %>%
    relocate(time_since_last, .after = time) %>%
    relocate(fu_time, .after = time_since_last)
  
  df_all <- df
  df <- df %>% group_by(ID) %>% filter(n() > 1) %>% ungroup()
  
  return(list(df = df, df_all = df_all))
}

make_clock <- function(df, clock_parameters) {
  
  year <- clock_parameters$year
  dx <- clock_parameters$dx
  tip_point <- clock_parameters$tip_point
  tip_point2 <- clock_parameters$tip_point2
  which_model <- clock_parameters$which_model
  
  df <- df %>%
    mutate(
      plasma_YN = if_else(plasma >= tip_point, 1, 0),
      plasma_YN = factor(plasma_YN, levels = c(0, 1), labels = c("Negative", "Positive"))
    ) %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(Previous_plasma_YN = dplyr::lag(plasma_YN)) %>%
    mutate(is_converter = any(Previous_plasma_YN == "Negative" & plasma_YN == "Positive", na.rm = TRUE),
           one_pos = any(plasma_YN == "Positive", na.rm = TRUE))
  
  safe_last_neg <- function(ages, plasma_yn) {
    neg_ages <- ages[plasma_yn == "Negative"]
    if(length(neg_ages) > 0) {
      return(neg_ages[length(neg_ages)])
    } else {
      return(NA_real_)
    }
  }
  
  safe_first_pos <- function(ages, plasma_yn) {
    pos_ages <- ages[plasma_yn == "Positive"]
    if(length(pos_ages) > 0) {
      return(pos_ages[1])
    } else {
      return(NA_real_)
    }
  }
  
  df <- df %>%
    group_by(ID) %>%
    mutate(
      Last_0_Age = if(first(is_converter)) safe_last_neg(AGE, plasma_YN) else NA_real_,
      First_1_Age = if(first(is_converter)) safe_first_pos(AGE, plasma_YN) else NA_real_,
      Midpoint_Age = (Last_0_Age + First_1_Age) / 2,
      plasma_time_conv = AGE - Midpoint_Age
    ) %>%
    ungroup() %>%
    arrange(ID, EXAMDATE)
  
  if (which_model == 1) {
    # Linear Mixed-Effects Model
    fit_plasma_rate <- lme(plasma ~ time,  
                           random = ~1 + time | ID, 
                           data = df, 
                           na.action = na.omit, 
                           method = "ML",
                           control = lmeControl(opt = "optim"))
    
    df <- df %>%
      left_join(
        tibble(
          ID = rownames(coef(fit_plasma_rate)),
          intercepts = coef(fit_plasma_rate)[,"(Intercept)"],
          slopes = coef(fit_plasma_rate)$time
        ),
        by = "ID"
      )
    
  } else if (which_model == 2) {
    # Linear Regression for each ID
    slopes_df <- df %>%
      group_by(ID) %>%
      summarize(
        slopes = if(n() >= 2) {
          model <- lm(plasma ~ time, data = cur_data())
          coef(model)[2]
        } else {
          NA_real_
        }
      )
    
    df <- df %>%
      left_join(slopes_df, by = "ID")
  }
  
  df <- df %>%
    group_by(ID) %>%
    mutate(
      Baseline_plasma = first(plasma),
      Baseline_plasma = intercepts,
      plasma_midpoint = ((fu_time / 2) * slopes) + Baseline_plasma
    ) %>%
    ungroup()
  
  df %>% 
    distinct(ID, .keep_all = TRUE) %>%
    ggplot(aes(x = plasma_midpoint, y = slopes)) + 
    geom_point() + 
    geom_smooth(method = "gam", color = "darkgreen") + 
    theme_cowplot() +
    labs(
      x = "Estimated %p-tau-217 at midpoint", 
      y = "%p-tau-217 rate of change",
      title = "Relationship between p-tau-217 levels and rate of change"
    )
  
  gam_model <- df %>%
    distinct(ID, .keep_all = TRUE) %>%
    gam(slopes ~ s(plasma_midpoint, bs = "cr"), data = .)
  
  seq_min <- min(df$plasma)
  seq_max <- max(df$plasma)
  plasma_midpoint <- tibble(plasma_midpoint = seq(from = seq_min, to = seq_max, by = dx) + dx/2) %>%
    mutate(
      estim_rate_i = predict(gam_model, newdata = ., type = "response"),
      recip_rate_i = 1 / estim_rate_i,
      plasma_TIME_INT_i = dx * recip_rate_i,
      TimeSum_i = cumsum(plasma_TIME_INT_i),
      plasma_time = TimeSum_i - TimeSum_i[which.min(abs(plasma_midpoint - tip_point2))]
    ) %>% 
    mutate(
      estim_rate_i = as.numeric(estim_rate_i),
      recip_rate_i = as.numeric(recip_rate_i),
      plasma_TIME_INT_i = as.numeric(plasma_TIME_INT_i)
    )
  
  find_nearest <- function(x, y) {
    y[findInterval(x, y, all.inside = TRUE)]
  }
  
  all_nearest <- find_nearest(df$plasma, plasma_midpoint$plasma_midpoint)
  
  df <- df %>%
    mutate(Nearest_Value = all_nearest) %>%
    left_join(plasma_midpoint, by = c("Nearest_Value" = "plasma_midpoint")) %>%
    arrange(ID, EXAMDATE)
  
  model_time <- gam(plasma_time ~ s(plasma, bs = "cr"), data = df)
  summary(model_time)
  
  return(list(df = df,
              plasma_midpoint = plasma_midpoint,
              model_time = model_time))
}

make_spaghetti_plot <- function(data, x_var, y_var, group_var, color_var, 
                                x_lim = NULL, y_lim = NULL, tip_point = NULL,
                                title, x_label, y_label,
                                grey_color = "#999999", red_color = "#E41A1C",
                                add_line = FALSE, line_data = NULL, 
                                line_x = NULL, line_y = NULL,
                                line_color = "black", line_width = 1.5) {
  
  x_var_sym <- sym(x_var)
  y_var_sym <- sym(y_var)
  group_var_sym <- sym(group_var)
  color_var_sym <- sym(color_var)
  
  p <- ggplot(data, aes(x = !!x_var_sym, y = !!y_var_sym, 
                        group = !!group_var_sym, color = as.factor(!!color_var_sym))) +
    geom_line(linewidth = .5) +
    geom_point(size = .75) +
    scale_color_manual(values = c("FALSE" = grey_color, "TRUE" = red_color))
  
  if (!is.null(tip_point)) {
    p <- p + geom_hline(yintercept = tip_point, linetype = "dashed", color = "black")
  }
  
  if (add_line && !is.null(line_data) && !is.null(line_x) && !is.null(line_y)) {
    line_x_sym <- sym(line_x)
    line_y_sym <- sym(line_y)
    
    p <- p + geom_line(data = line_data, 
                       aes(x = !!line_x_sym, y = !!line_y_sym, group = NULL, color = NULL), 
                       color = line_color, linewidth = line_width)
  }
  
  p <- p + theme_cowplot() +
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme(legend.position = "none")
  
  if (!is.null(y_lim)) {
    p <- p + ylim(y_lim[1], y_lim[2])
  }
  
  if (!is.null(x_lim)) {
    p <- p + xlim(x_lim[1], x_lim[2])
  }
  
  return(p)
}

create_clock_intervals <- function(midpoint_data, 
                                   biomarker_col = NULL,
                                   time_col = NULL,
                                   start_value = NULL, 
                                   end_value = NULL,
                                   by_value = 1,
                                   reference_value = NULL) {
  
  clock_intervals <- data.frame(Beginning = seq(start_value, end_value, by = by_value)) %>%
    mutate(
      End = Beginning + by_value,
      Time_Start = sapply(Beginning, function(x) {
        midpoint_data[[time_col]][which.min(abs(midpoint_data[[biomarker_col]] - x))]
      }),
      Time_End = sapply(End, function(x) {
        midpoint_data[[time_col]][which.min(abs(midpoint_data[[biomarker_col]] - x))]
      }),
      Interval = round(Time_End - Time_Start, 1)
    ) %>%
    mutate(
      Cumulative_raw = cumsum(Interval),
      ref_row = which.min(abs(Beginning - reference_value)),
      Cumulative = round(Cumulative_raw - Cumulative_raw[ref_row - 1], 1)
    ) %>%
    select(-Cumulative_raw, -ref_row)
  
  clock_intervals_export <- clock_intervals[, c("Beginning", "End", "Interval", "Cumulative")]
  colnames(clock_intervals_export) <- c("Plasma beginning", "Plasma end", "Interval", "Cumulative (from 4.06%)")
  
  return(clock_intervals_export)
}
