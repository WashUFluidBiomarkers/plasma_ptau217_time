# Created by Kellen Petersen, July 1, 2025

library(tidyverse)
library(ggdist)
library(cowplot)
library(ggpubr)
library(here)
library(rstatix)
library(conflicted)
library(patchwork)

setwd("<SET WORKING DIRECTORY>")

conflict_prefer("filter", "dplyr")
conflict_prefer("first", "dplyr")

df_long_tira_0 <- read.csv(here("results","<LOAD DATA>"))
df_long_sila_0 <- read.csv(here("results","<LOAD DATA>"))

df_long_tira <- df_long_tira_0 %>% 
  dplyr::select(ID, plasma, AGE,
                est_onset_age, years_since_onset,
                years_since_predicted_symptoms) %>%
  rename(plasma_TIRA = plasma) %>%
  rename(est_onset_age_TIRA = est_onset_age) %>%
  rename(years_since_onset_TIRA = years_since_onset) %>%
  rename(years_since_predicted_symptoms_TIRA = years_since_predicted_symptoms) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

df_long_sila <- df_long_sila_0 %>% 
  dplyr::select(ID, plasma, AGE, 
                est_onset_age, years_since_onset,
                years_since_predicted_symptoms) %>%
  rename(plasma_SILA = plasma) %>%
  rename(est_onset_age_SILA = est_onset_age) %>%
  rename(years_since_onset_SILA = years_since_onset) %>%
  rename(years_since_predicted_symptoms_SILA = years_since_predicted_symptoms) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup()

df_long <- df_long_tira %>%
  left_join(df_long_sila, by = "ID") %>% 
  dplyr::select(-AGE.y) %>% 
  rename(AGE=AGE.x)

df0 <- read.csv(here("data_raw","<LOAD DATA>"))

df <- df0 %>%
  dplyr::select(RID, EXAMDATE,
                C2N_plasma_ptau217_ratio,
                CENTILOIDS, CENTILOIDS_10,
                MesialTemporal, TAU_MesialTemporal_10,
                TemporoParietal, TAU_TemporoParietal_10) %>%
  rename(ID = RID) 

df <- df %>%
  filter(!is.na(C2N_plasma_ptau217_ratio))

df$ID <- as.factor(df$ID)
df$EXAMDATE <- as.Date(df$EXAMDATE, format = "%Y-%m-%d")

df <- df %>%
  mutate(
    AT_stage = case_when(
      CENTILOIDS_10 == 0 & TAU_MesialTemporal_10 == 0 & TAU_TemporoParietal_10 == 0 ~ 0,
      CENTILOIDS_10 == 1 & TAU_MesialTemporal_10 == 0 & TAU_TemporoParietal_10 == 0 ~ 1,
      CENTILOIDS_10 == 1 & TAU_MesialTemporal_10 == 1 & TAU_TemporoParietal_10 == 0 ~ 2,
      CENTILOIDS_10 == 1 & TAU_MesialTemporal_10 == 1 & TAU_TemporoParietal_10 == 1 ~ 3,
      TRUE ~ NA_integer_
    )
  )

df_stages <- df

create_dataframe <- function(which_model, df = df_stages) {
  if (which_model == "tira") {
    model <- "TIRA"
    dataset <- "ADNI"
    clock <- read.csv(here("results","<LOAD DATA>"))
    df <- read.csv(here("data_raw","<LOAD DATA>"))
  } else if (which_model == "sila") {
    model <- "SILA"
    dataset <- "ADNI"
    clock <- read.csv(here("results_SILA","<LOAD DATA>")) %>%
      mutate(clock_time = adtime)
    df <- read.csv(here("data_raw","<LOAD DATA>"))
  }
  
  df  <- df %>%
    rename(ID = RID,
           plasma = C2N_plasma_ptau217_ratio) %>%
    mutate(ID = as.factor(ID),
           EXAMDATE = as.Date(EXAMDATE, format = "%Y-%m-%d"),
           AGE = as.numeric(AGE),
           plasma = as.numeric(plasma))
  
  clock_plasma_min <- min(clock$plasma)
  clock_plasma_max <- max(clock$plasma)
  df$plasma_time <- sapply(df$plasma, function(x) {
    if (is.na(x)) return(NA)
    if (x < clock_plasma_min || x > clock_plasma_max) {
      return(NA)
    }
    idx <- which.min(abs(clock$plasma - x))
    clock$clock_time[idx]
  })
  df$plasma_time <- as.numeric(df$plasma_time)
  
  df <- df %>%
    arrange(ID, EXAMDATE) %>%
    group_by(ID) %>%
    mutate(time = (difftime(EXAMDATE, first(EXAMDATE), units = "days") / 365.25),
           time = as.numeric(time),
           fu_time = max(time)) %>%
    ungroup() %>%
    relocate(time, .after = EXAMDATE) %>%
    relocate(fu_time, .after = time)
  
  df <- df %>%
    mutate(est_onset_age_i = AGE - plasma_time) %>%
    group_by(ID) %>%
    mutate(est_onset_age = ifelse(n() == 1, est_onset_age_i, mean(est_onset_age_i, na.rm = TRUE)),
           years_since_onset = AGE - est_onset_age) %>%
    ungroup() %>%
    arrange(ID, EXAMDATE)
  
  df2 <- df %>%
    dplyr::select(ID, EXAMDATE, AGE,
                  plasma_time, est_onset_age_i, est_onset_age, years_since_onset)
  df2$ID <- as.factor(df2$ID)
  df2$EXAMDATE <- as.Date(df2$EXAMDATE, format = "%Y-%m-%d")
  
  DF <- df_stages %>%
    left_join(df2, by = c("ID" = "ID", "EXAMDATE" = "EXAMDATE")) %>%
    mutate(
      AT_stage = factor(
        AT_stage,
        levels = c(0, 1, 2, 3),
        labels = c(
          "A-T_early-T_late-",
          "A+T_early-T_late-",
          "A+T_early+T_late-",
          "A+T_early+T_late+"
        )
      ),
      AT_stage = as.character(AT_stage)
    ) %>%
    arrange(ID, EXAMDATE)
  
  return(DF)
}

DF_tira <- create_dataframe("tira")
DF_sila <- create_dataframe("sila")

at_stage_colors <- c(
  "A-T_early-T_late-" = "#B0B0B0",
  "A+T_early-T_late-" = "#4682B4",
  "A+T_early+T_late-" = "#90EE90",
  "A+T_early+T_late+" = "#006400"
)

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

which_dataset <- "ADNI"
which_method <- "TIRA"
result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
model_intercept <- result_tmp$intercept
model_slope <- result_tmp$slope
DF_tira <- DF_tira %>% 
  mutate(est_symtoms_age = est_onset_age*model_slope + model_intercept) %>% 
  mutate(years_since_symptoms = AGE - est_symtoms_age)

which_dataset <- "ADNI"
which_method <- "SILA"
result_tmp <- df_linear[df_linear$Dataset == which_dataset & df_linear$Model == which_method, c("intercept", "slope")]
model_intercept <- result_tmp$intercept
model_slope <- result_tmp$slope
DF_sila <- DF_sila %>% 
  mutate(est_symtoms_age = est_onset_age*model_slope + model_intercept) %>% 
  mutate(years_since_symptoms = AGE - est_symtoms_age)

DF_tira_2 <- DF_tira %>% filter(!is.na(AT_stage), !is.na(years_since_onset))
DF_sila_2 <- DF_sila %>% filter(!is.na(AT_stage), !is.na(years_since_onset))

DF_tira_p <- DF_tira %>% filter(!is.na(AT_stage), !is.na(C2N_plasma_ptau217_ratio))
DF_sila_p <- DF_tira %>% filter(!is.na(AT_stage), !is.na(C2N_plasma_ptau217_ratio))

table(DF_tira$AT_stage)
table(DF_tira_2$AT_stage)
table(DF_tira_p$AT_stage)

table(DF_sila$AT_stage)
table(DF_sila_2$AT_stage)
table(DF_sila_p$AT_stage)

makeATStageRaincloud_0 <- function(data,
                                   value_var,
                                   y_label = "Years since %p-tau217 positivity",
                                   titlename = "TIRA",
                                   normalization = "none",
                                   group_var = AT_stage,
                                   show_n = TRUE) {
  value_var <- ensym(value_var)
  group_var <- ensym(group_var)
  
  df_plot <- data %>%
    filter(!is.na(!!group_var), !is.na(!!value_var))
  
  df_plot[[quo_name(group_var)]] <- factor(
    df_plot[[quo_name(group_var)]],
    levels = c("A-T_early-T_late-",
               "A+T_early-T_late-",
               "A+T_early+T_late-",
               "A+T_early+T_late+")
  )
  
  n_counts <- df_plot %>%
    group_by(AT_stage) %>%
    summarise(n = n()) %>%
    deframe()
  
  if (show_n) {
    at_stage_labels_with_n <- c(
      "A-T_early-T_late-" = bquote(atop("A-T"["early"]*"-T"["late"]*"-", .(paste0("n=", n_counts["A-T_early-T_late-"])))),
      "A+T_early-T_late-" = bquote(atop("A+T"["early"]*"-T"["late"]*"-", .(paste0("n=", n_counts["A+T_early-T_late-"])))),
      "A+T_early+T_late-" = bquote(atop("A+T"["early"]*"+"*"T"["late"]*"-", .(paste0("n=", n_counts["A+T_early+T_late-"])))),
      "A+T_early+T_late+" = bquote(atop("A+T"["early"]*"+"*"T"["late"]*"+", .(paste0("n=", n_counts["A+T_early+T_late+"]))))
    )
  } else {
    at_stage_labels_with_n <- c(
      "A-T_early-T_late-" = expression("A-T"["early"]*"-T"["late"]*"-"),
      "A+T_early-T_late-" = expression("A+T"["early"]*"-T"["late"]*"-"),
      "A+T_early+T_late-" = expression("A+T"["early"]*"+"*"T"["late"]*"-"),
      "A+T_early+T_late+" = expression("A+T"["early"]*"+"*"T"["late"]*"+")
    )
  }
  
  plot <- ggplot(df_plot, aes(x = !!group_var, y = !!value_var, fill = !!group_var)) +
    ggdist::geom_dots(
      side = "bottom",
      dotsize = 0.8,
      stackratio = 1.05,
      alpha = 0.7
    ) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA
    ) +
    ggdist::stat_halfeye(
      side = "top",
      alpha = 0.7,
      normalize = normalization
    ) +
    scale_fill_manual(values = at_stage_colors) +
    scale_x_discrete(labels = at_stage_labels_with_n) +
    labs(
      title = titlename,
      y = y_label,
      x = ""
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    ) +
    ggpubr::stat_compare_means(
      comparisons = list(
        c("A-T_early-T_late-", "A+T_early-T_late-"),
        c("A-T_early-T_late-", "A+T_early+T_late-"),
        c("A-T_early-T_late-", "A+T_early+T_late+"),
        c("A+T_early-T_late-", "A+T_early+T_late-"),
        c("A+T_early-T_late-", "A+T_early+T_late+"),
        c("A+T_early+T_late-", "A+T_early+T_late+")
      ),
      method = "wilcox.test",
      p.adjust.method = "holm",
      label = "p.signif",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      ),
      hide.ns = TRUE,
      step.increase = 0.03,
      vjust = 1,
      tip.length = 0.01
    )
  
  return(plot)
}

makeATStageRaincloud_1 <- function(data,
                                   value_var,
                                   y_label = "Years since %p-tau217 positivity",
                                   titlename = "TIRA",
                                   normalization = "none",
                                   group_var = AT_stage,
                                   show_n = TRUE) {
  value_var <- ensym(value_var)
  group_var <- ensym(group_var)
  
  df_plot <- data %>%
    filter(!is.na(!!group_var), !is.na(!!value_var))
  
  df_plot[[quo_name(group_var)]] <- factor(
    df_plot[[quo_name(group_var)]],
    levels = c("A-T_early-T_late-",
               "A+T_early-T_late-",
               "A+T_early+T_late-",
               "A+T_early+T_late+")
  )
  
  n_counts <- df_plot %>%
    group_by(AT_stage) %>%
    summarise(n = n()) %>%
    deframe()
  
  if (show_n) {
    at_stage_labels_with_n <- c(
      "A-T_early-T_late-" = bquote(atop("A-T_early-T_late-", .(paste0("n=", n_counts["A-T_early-T_late-"])))),
      "A+T_early-T_late-" = bquote(atop("A+T"["early"]*"-T"["late"]*"-", .(paste0("n=", n_counts["A+T_early-T_late-"])))),
      "A+T_early+T_late-" = bquote(atop("A+T"["early"]*"+"*"T"["late"]*"-", .(paste0("n=", n_counts["A+T_early+T_late-"])))),
      "A+T_early+T_late+" = bquote(atop("A+T"["early"]*"+"*"T"["late"]*"+", .(paste0("n=", n_counts["A+T_early+T_late+"]))))
    )
  } else {
    at_stage_labels_with_n <- c(
      "A-T_early-T_late-" = expression("A-T_early-T_late-"),
      "A+T_early-T_late-" = expression("A+T"["early"]*"-T"["late"]*"-"),
      "A+T_early+T_late-" = expression("A+T"["early"]*"+"*"T"["late"]*"-"),
      "A+T_early+T_late+" = expression("A+T"["early"]*"+"*"T"["late"]*"+")
    )
  }
  
  kw_test <- df_plot %>%
    kruskal_test(as.formula(paste(quo_name(value_var), "~", quo_name(group_var))))
  
  plot <- ggplot(df_plot, aes(x = !!group_var, y = !!value_var, fill = !!group_var)) +
    ggdist::geom_dots(
      side = "bottom",
      dotsize = 0.8,
      stackratio = 1.05,
      alpha = 0.7
    ) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA
    ) +
    ggdist::stat_halfeye(
      side = "top",
      alpha = 0.7,
      normalize = normalization
    ) +
    scale_fill_manual(values = at_stage_colors) +
    scale_x_discrete(labels = at_stage_labels_with_n) +
    labs(
      title = titlename,
      y = y_label,
      x = ""
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
  
  if (kw_test$p < 0.05) {
    dunn_test_results <- df_plot %>%
      dunn_test(as.formula(paste(quo_name(value_var), "~", quo_name(group_var))),
                p.adjust.method = "holm")
    
    dunn_test_results <- dunn_test_results %>%
      add_xy_position(x = quo_name(group_var))
    
    plot <- plot +
      stat_pvalue_manual(
        dunn_test_results,
        hide.ns = TRUE,
        label = "p.adj.signif",
        tip.length = 0.01,
        step.increase = 0.03,
        vjust = 1
      ) +
      labs(subtitle = paste("Kruskal-Wallis, p =", signif(kw_test$p, 3)))
  } else {
    plot <- plot +
      labs(subtitle = paste("Kruskal-Wallis, p =", signif(kw_test$p, 3), "(ns)"))
  }
  
  return(plot)
}

makeATStageRaincloud <- makeATStageRaincloud_0

DF_tira_plot <- DF_tira %>% 
  dplyr::filter((C2N_plasma_ptau217_ratio>1.06) & (C2N_plasma_ptau217_ratio<10.45)) %>% 
  dplyr::filter(!is.na(CENTILOIDS)) %>% 
  dplyr::filter(!is.na(MesialTemporal)) %>% 
  dplyr::filter(!is.na(TemporoParietal))

DF_sila_plot <- DF_sila %>% 
  dplyr::filter((C2N_plasma_ptau217_ratio>1.06) & (C2N_plasma_ptau217_ratio<10.45)) %>% 
  dplyr::filter(!is.na(CENTILOIDS_10)) %>% 
  dplyr::filter(!is.na(MesialTemporal)) %>% 
  dplyr::filter(!is.na(TemporoParietal))

DF_plasma_plot <- DF_tira %>% 
  dplyr::filter((C2N_plasma_ptau217_ratio>1.06) & (C2N_plasma_ptau217_ratio<10.45)) %>% 
  dplyr::filter(!is.na(CENTILOIDS_10)) %>% 
  dplyr::filter(!is.na(MesialTemporal)) %>% 
  dplyr::filter(!is.na(TemporoParietal))

DF_tira_plot2 <- DF_tira_plot %>% 
  group_by(ID) %>% 
  slice(n()) %>% 
  ungroup()

DF_sila_plot2 <- DF_sila_plot %>%
  group_by(ID) %>% 
  slice(n()) %>% 
  ungroup()

DF_plasma_plot2 <- DF_plasma_plot %>%
  group_by(ID) %>% 
  slice(n()) %>% 
  ungroup()

y_label <- "Years since estimated %p-tau217 positivity"
at_stage_plot_tira <- makeATStageRaincloud(DF_tira_plot2, years_since_onset, y_label, "TIRA", "groups", show_n = TRUE)
at_stage_plot_sila <- makeATStageRaincloud(DF_sila_plot2, years_since_onset, y_label, "SILA", "groups", show_n = TRUE)

y_label <- "Plasma %p-tau217"
at_stage_plot_plasma <- makeATStageRaincloud(DF_plasma_plot2, C2N_plasma_ptau217_ratio, y_label, "Plasma", "groups", show_n = TRUE)

y_label <- "Estimated years since symptoms onset"
at_stage_plot_tira_symps <- makeATStageRaincloud(DF_tira_plot2, years_since_symptoms, y_label, "TIRA", "groups", show_n = TRUE)
at_stage_plot_sila_symps <- makeATStageRaincloud(DF_sila_plot2, years_since_symptoms, y_label, "SILA", "groups", show_n = TRUE)

pp <- (at_stage_plot_plasma) / 
  ((at_stage_plot_tira) + (at_stage_plot_sila)) / 
  ((at_stage_plot_tira_symps) + (at_stage_plot_sila_symps))
pp <- pp +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold")) +
  theme(plot.tag.position = c(0, 1), plot.tag = element_text(hjust = 0, vjust = 1))
pp
ggsave("results/at_stage_raincloud_combined_FINAL.png", pp, width = 14, height = 16, dpi = 500, bg = "white")