# Originally created by Kellen Petersen, July 1, 2025
# Updated by Kellen Petersen, September 12, 2025

# Load required libraries for data analysis and visualization
library(tidyverse) # Data manipulation and visualization
library(data.table)
library(cowplot) # Publication-ready plots
library(patchwork)
library(here) # Project-relative file paths
library(nlme) # Linear mixed-effects models
library(lme4) # Mixed-effects models
library(mgcv) # Generalized additive models
library(effects) # Effect displays
library(performance) # Model diagnostics
library(rmcorr) # Repeated measures correlation
library(tictoc) # Timing
library(demography) # Demographic analysis

# Set working directory and source functions
setwd("")
i_am("")
here()
source("main_functions.R")
source("auto_variance.R")

# Set analysis parameters
which_dataset <- "KADRC"

# Load dataset based on selection
if (which_dataset == "ADNI") {
  df0 <- read.csv(here("data_raw",""), header = TRUE)
} else if (which_dataset == "KADRC") {
  df0 <- read.csv(here("data_raw",""), header = TRUE)
} else {
  stop("Invalid method selected. Please choose 'TIRA' or 'SILA'.")
}

# Standardize ID column name across datasets
df0 <- df0 %>% rename(ID = RID)

# Define analysis parameters
variance_cut <- 0.90
which_biomarker <- "BLANK_ENTRY"
tip_point <- 4.06
choose_plasma <- "C2N_plasma_ptau217_ratio"

# Filter out missing plasma values
df0 <- df0 %>% filter(!is.na(df0[[choose_plasma]]))

# Run automatic variance analysis to identify stable biomarker regions
results_auto_variance <- auto_variance(which_biomarker,
                                       df0,
                                       choose_plasma,
                                       variance_cut,
                                       tip_point)

# Create customized plots with specified axis limits
rate_plot <- results_auto_variance$rate_plot + ylim(0, 2) + xlim(0, 22)
plot_variance <- results_auto_variance$p_variance + ylim(0, 0.005) + xlim(0, 22)

# Extract cut points for biomarker range restriction
restrict_range_0 <- results_auto_variance$cut_points

# Display combined plots
rate_plot | plot_variance

# Output the calculated restriction range
restrict_range_0
