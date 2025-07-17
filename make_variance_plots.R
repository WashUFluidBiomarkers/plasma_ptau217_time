# Created by Kellen Petersen, July 1, 2025

library(tidyverse)  # Data manipulation and visualization
library(data.table)
library(cowplot)    # Publication-ready plots
library(patchwork)

library(here)       # Project-relative file paths

library(nlme)       # Linear mixed-effects models
library(lme4)       # Mixed-effects models
library(mgcv)       # Generalized additive models
library(effects)    # Effect displays
library(performance) # Model diagnostics
library(rmcorr)     # Repeated measures correlation

library(tictoc)     # Timing
library(demography) # Demographic analysis

setwd("<SET WORKING DIRECTORY>")
i_am("<SET I_AM>")
here()
source("main_functions.R")
source("auto_variance.R")

which_dataset <- "KADRC"
if (which_dataset == "ADNI") {
  df0 <- read.csv(here("data_raw","<LOAD DATA>"), header = TRUE)
} else if (which_dataset == "KADRC") {
  df0 <- read.csv(here("data_raw","<LOAD DATA>"), header = TRUE)
} else {
  stop("Invalid method selected. Please choose 'TIRA' or 'SILA'.")
}
df0 <- df0 %>% rename(ID = RID)
variance_cut <- 0.90
which_biomarker <- "BLANK_ENTRY"
tip_point <- 4.06
choose_plasma <- "C2N_plasma_ptau217_ratio"

df0 <- df0 %>% filter(!is.na(df0[[choose_plasma]]))

results_auto_variance <- auto_variance(which_biomarker, 
                                       df0,
                                       choose_plasma,
                                       variance_cut,
                                       tip_point)
rate_plot <- results_auto_variance$rate_plot + ylim(0, 2) + xlim(0, 22)
plot_variance <- results_auto_variance$p_variance + ylim(0, 0.005) + xlim(0, 22)
restrict_range_0 <- results_auto_variance$cut_points

rate_plot | plot_variance
restrict_range_0
