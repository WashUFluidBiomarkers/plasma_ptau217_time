# Originally created by Kellen Petersen, July 1, 2025
# Updated by Kellen Petersen, September 12, 2025

# Load required packages for plotting and statistical analysis
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)
library(here)
library(DescTools)
library(nopaco)

setwd("")

# Analysis parameters - specify which method and dataset scope to use
which_method <- "TIRA"
AD_or_ALL <- "ALL"

# Basic preprocessing function to prepare symptom conversion data
preprocess_data <- function(df, symps, outcome, converter_ids) {
  combined <- df %>%
    mutate(symp_conversion = age_at_symptom_onset)
  return(combined)
}

# Main function to create combined scatter plot comparing multiple datasets
create_combined_plot <- function(lines_config, include_legend = TRUE, AD_or_ALL) {
  all_data <- data.frame()
  
  # Process each dataset configuration
  for (i in seq_along(lines_config)) {
    config <- lines_config[[i]]
    cohort <- case_when(
      config$annotation_title == "Knight ADRC" ~ "KADRC",
      config$annotation_title == "ADNI" ~ "ADNI",
      TRUE ~ NA_character_
    )
    
    # Load appropriate converter IDs based on analysis scope
    if (AD_or_ALL == "AD") {
      id_file <- paste0(
        here("data_final"),
        "/"
      )
    } else if (AD_or_ALL == "ALL") {
      id_file <- paste0(
        here("data_final"),
        "/"
      )
    } else {
      stop("Error: Choose correct dataset")
    }
    
    converter_ids <- read.csv(id_file)$ID
    print(length(converter_ids))
    
    # Get data from global environment and preprocess
    df <- get(config$df_name)
    symps <- get(config$symps_name)
    data <- preprocess_data(df, symps, config$outcome, converter_ids)
    data$dataset <- config$annotation_title
    all_data <- rbind(all_data, data)
  }
  
  # Ensure consistent factor ordering for datasets
  all_data$dataset <- factor(all_data$dataset,
                             levels = sapply(lines_config, function(config) config$annotation_title))
  
  title_label <- ""
  
  # Create base scatter plot with consistent theme
  combined_plot <- ggplot(all_data, aes(x = est_onset_age, y = symp_conversion, color = dataset)) +
    labs(
      title = title_label,
      x = "Estimated age at plasma %p-tau217 positivity (years)",
      y = "Age at symptom onset (years)"
    ) +
    theme_cowplot() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = if(include_legend) c(0.95, 0.05) else "none",
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0.5, linetype = "solid"),
      legend.margin = margin(6, 6, 6, 6)
    )
  
  # Add regression lines and points for non-black colors first
  for (i in seq_along(lines_config)) {
    config <- lines_config[[i]]
    if (config$line_color == "black") {
      next
    }
    
    dataset_data <- all_data[all_data$dataset == config$annotation_title, ]
    combined_plot <- combined_plot +
      geom_smooth(data = dataset_data,
                  method = "lm",
                  se = config$ci_include,
                  linetype = config$line_type,
                  size = config$line_size,
                  fill = config$ci_color) +
      geom_point(data = dataset_data,
                 shape = config$point_shape,
                 size = config$point_size)
  }
  
  # Add regression lines and points for black color last (plotting order)
  for (i in seq_along(lines_config)) {
    config <- lines_config[[i]]
    if (config$line_color != "black") {
      next
    }
    
    dataset_data <- all_data[all_data$dataset == config$annotation_title, ]
    combined_plot <- combined_plot +
      geom_smooth(data = dataset_data,
                  method = "lm",
                  se = config$ci_include,
                  linetype = config$line_type,
                  size = config$line_size,
                  fill = config$ci_color) +
      geom_point(data = dataset_data,
                 shape = config$point_shape,
                 size = config$point_size)
  }
  
  # Apply color scheme and legend
  combined_plot <- combined_plot +
    scale_color_manual(values = sapply(lines_config, function(config) config$line_color)) +
    guides(color=guide_legend(title="Dataset")) +
    theme(legend.title.align = 0.5)
  
  return(combined_plot)
}

# Function to add statistical annotations to the plot
add_annotations <- function(plot, lines_config, AD_or_ALL) {
  # Get plot dimensions for annotation positioning
  plot_build <- ggplot_build(plot)
  x_min <- plot_build$layout$panel_params[[1]]$x.range[1]
  x_max <- plot_build$layout$panel_params[[1]]$x.range[2]
  y_min <- plot_build$layout$panel_params[[1]]$y.range[1]
  y_max <- plot_build$layout$panel_params[[1]]$y.range[2]
  annotation_x <- x_min + (x_max - x_min) * 0.05
  annotation_y <- y_max - (y_max - y_min) * 0.05
  annotation_spacing <- (y_max - y_min) * 0.2
  
  results_list <- list()
  
  # Process each dataset for statistical calculations
  for (i in seq_along(lines_config)) {
    config <- lines_config[[i]]
    cohort <- case_when(
      config$annotation_title == "Knight ADRC" ~ "KADRC",
      config$annotation_title == "ADNI" ~ "ADNI",
      TRUE ~ NA_character_
    )
    
    # Load converter IDs based on analysis scope
    if (AD_or_ALL == "AD") {
      id_file <- paste0(
        here("data_final"),
        "/"
      )
    } else if (AD_or_ALL == "ALL") {
      id_file <- paste0(
        here("data_final"),
        "/"
      )
    } else {
      stop("Error: Choose correct dataset")
    }
    
    converter_ids <- read.csv(id_file)$ID
    df <- get(config$df_name)
    symps <- get(config$symps_name)
    data <- preprocess_data(df, symps, config$outcome, converter_ids)
    
    # Fit linear model for regression statistics
    lm_fit <- lm(symp_conversion ~ est_onset_age, data = data)
    
    # Ensure consistent factor levels for demographic variables
    data$APOE4_positive <- factor(data$APOE4_positive, levels = c("0", "1"))
    
    if (i==1) {
      data$PTGENDER <- factor(data$PTGENDER, levels = c("M","F"))
    } else {
      data$PTGENDER <- factor(data$PTGENDER, levels = c("1","2"))
    }
    
    lm_fit2 <- lm(symp_conversion ~ est_onset_age, data = data)
    n_points <- nrow(data)
    
    # Extract key regression statistics
    adj_r2 <- summary(lm_fit)$adj.r.squared
    intercept <- coef(lm_fit)[1]
    slope <- coef(lm_fit)[2]
    spearman_r <- cor.test(data$symp_conversion, data$est_onset_age, method = "spearman", exact=FALSE)$estimate
    r2 <- summary(lm_fit)$r.squared
    
    # Calculate concordance correlation coefficients
    ccc_result <- DescTools::CCC(data$est_onset_age, data$symp_conversion)$rho.c
    print(ccc_result)
    ccc <- ccc_result$est
    print(ccc)
    
    # Non-parametric concordance correlation coefficient
    mat <- cbind(data$est_onset_age, data$symp_conversion)
    conc_test <- nopaco::concordance.test(mat)
    psi <- nopaco::getPsi(conc_test)
    np_ccc <- nopaco::rfromPsi(psi)
    print(np_ccc)
    
    # Print statistical summary to console
    cat("\n--- Metrics for", config$annotation_title, "---\n")
    cat("R-squared: ", round(r2, 4), "\n")
    cat("Adjusted R-squared: ", round(adj_r2, 4), "\n")
    cat("Spearman r: ", round(spearman_r, 4), "\n")
    cat("CCC: ", round(ccc, 4), "\n")
    cat("Non-parametric CCC: ", round(np_ccc, 4), "\n")
    
    # Store results for export
    results_list[[i]] <- data.frame(
      dataset = config$annotation_title,
      r2 = r2,
      adj_r2 = adj_r2,
      spearman_r = spearman_r,
      ccc = ccc,
      np_ccc = np_ccc,
      intercept = intercept,
      slope = slope,
      n = n_points,
      stringsAsFactors = FALSE
    )
    
    # Add regression equation and statistics as text annotation
    plot <- plot +
      annotate("text",
               x = annotation_x,
               y = annotation_y - (i - 1) * annotation_spacing,
               label = paste0(
                 "\ny = ", sprintf("%.2f", slope), "x + ",
                 sprintf("%.2f", intercept),
                 "\nAdjusted R² = ", sprintf("%.3f", adj_r2),
                 "\nSpearman ρ = ", sprintf("%.3f", spearman_r),
                 "\nN = ", n_points),
               hjust = 0,
               vjust = 1,
               color = config$line_color)
  }
  
  return(list(
    plot = plot,
    results = do.call(rbind, results_list)
  ))
}

# Load datasets based on analysis scope (AD converters only vs all participants)
if (AD_or_ALL == "AD") {
  df_ADNI_TIRA <- read.csv(here("data_final",""))
  df_ADNI_TIRA$ID <- as.factor(df_ADNI_TIRA$ID)
  df_KADRC_TIRA <- read.csv(here("data_final",""))
  df_KADRC_TIRA$ID <- as.factor(df_KADRC_TIRA$ID)
  df_ADNI_SILA <- read.csv(here("data_final",""))
  df_ADNI_SILA$ID <- as.factor(df_ADNI_SILA$ID)
  df_KADRC_SILA <- read.csv(here("data_final",""))
  df_KADRC_SILA$ID <- as.factor(df_KADRC_SILA$ID)
} else if (AD_or_ALL == "ALL") {
  df_ADNI_TIRA <- read.csv(here("data_final",""))
  df_ADNI_TIRA$ID <- as.factor(df_ADNI_TIRA$ID)
  df_KADRC_TIRA <- read.csv(here("data_final",""))
  df_KADRC_TIRA$ID <- as.factor(df_KADRC_TIRA$ID)
  df_ADNI_SILA <- read.csv(here("data_final",""))
  df_ADNI_SILA$ID <- as.factor(df_ADNI_SILA$ID)
  df_KADRC_SILA <- read.csv(here("data_final",""))
  df_KADRC_SILA$ID <- as.factor(df_KADRC_SILA$ID)
} else {
  stop("Error: Choose correct dataset")
}

# Load symptom data for each cohort and method
symps_ADNI_TIRA <- read.csv(here("data_final",""))
symps_ADNI_TIRA$ID <- as.factor(symps_ADNI_TIRA$ID)
symps_KADRC_TIRA <- read.csv(here("data_final",""))
symps_KADRC_TIRA$ID <- as.factor(symps_KADRC_TIRA$ID)
symps_ADNI_SILA <- read.csv(here("data_final","")) %>%
  filter(!is.na(EXAMDATE))
symps_ADNI_SILA$ID <- as.factor(symps_ADNI_SILA$ID)
symps_KADRC_SILA <- read.csv(here("data_final",""))
symps_KADRC_SILA$ID <- as.factor(symps_KADRC_SILA$ID)

# Configure plot appearance and data sources for each dataset
lines_config <- list(
  list(df_name=str_c("df_KADRC_",which_method), symps_name=str_c("symps_KADRC_",which_method),
       outcome="CDR_DX", method = which_method,
       line_color="darkgreen", line_type="solid", line_size=1.5,
       ci_include=TRUE, ci_color="lightgreen",
       point_shape=16, point_color="darkgreen", point_size=3,
       annotation_title="Knight ADRC"),
  list(df_name=str_c("df_ADNI_",which_method), symps_name=str_c("symps_ADNI_",which_method),
       outcome="CDR_DX_Imp", method = which_method,
       line_color="black", line_type="solid", line_size=1.5,
       ci_include=TRUE, ci_color="darkgray",
       point_shape=16, point_color="black", point_size=3,
       annotation_title="ADNI")
)

# Extract method and outcome information for file naming
tmp <- lines_config[[2]]
fn_method <- tmp$method
fn_outcome <- tmp$outcome

# Create base plot
base_plot <- create_combined_plot(lines_config=lines_config,
                                  include_legend=TRUE,
                                  AD_or_ALL)

# Add statistical annotations to the plot
add_annotations <- add_annotations(base_plot,
                                   lines_config,
                                   AD_or_ALL)

final_plot <- add_annotations$plot
results_df <- add_annotations$results

print(final_plot)

# Generate output filenames with timestamp and method information
filename0 <- str_c("results/simple_symptoms_MULTIPLE_byID_",
                   fn_method, "_",
                   fn_outcome, "_",
                   format(Sys.Date(), "%Y-%m-%d"), "_",
                   format(Sys.time(), "%H-%M-%S"))

filename0 <- str_c("results/simple_symptoms_MULTIPLE_byID_",
                   fn_method, "_",
                   fn_outcome, "_",
                   "FINAL")

# Add sensitivity analysis suffix if analyzing all participants
if (AD_or_ALL == "ALL") {
  filename0 <- str_c(filename0,"_sens")
}

# Generate filename for statistical results
filename_linear_models <- str_c("results/linear_models_",
                                fn_method, "_",
                                fn_outcome, "_",
                                "FINAL")

if (AD_or_ALL == "ALL") {
  filename_linear_models <- str_c(filename_linear_models,"_sens")
}

# Export statistical results to CSV
write.csv(results_df, str_c(filename_linear_models, ".csv"), row.names = FALSE)
