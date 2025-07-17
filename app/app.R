# Created by Kellen Petersen, July 1, 2025

library(shiny)
library(ggplot2)
library(dplyr)
library(icenReg)
library(cowplot)
library(scales)

parameter_sets <- list(
  "C2N Diagnostics PrecivityAD2 %p-tau217_KADRC_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-20, 20), y_lim = c(0, 13),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "C2N Diagnostics PrecivityAD2 %p-tau217_KADRC_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-15, 15), y_lim = c(0, 11),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "C2N Diagnostics PrecivityAD2 %p-tau217_ADNI_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "C2N Diagnostics PrecivityAD2 %p-tau217_ADNI_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "C2N Diagnostics PrecivityAD2 p-tau217_ADNI_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "C2N Diagnostics PrecivityAD2 p-tau217_ADNI_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42_ADNI_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42_ADNI_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "Fujirebio Diagnostics Lumipulse p-tau217_ADNI_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "Fujirebio Diagnostics Lumipulse p-tau217_ADNI_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "Janssen LucentAD Quanterix p-tau217_ADNI_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "Janssen LucentAD Quanterix p-tau217_ADNI_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "ALZpath Quanterix p-tau217_ADNI_TIRA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  ),
  "ALZpath Quanterix p-tau217_ADNI_SILA" = list(
    plasma_min = NULL, plasma_max = NULL, x_lim = c(-18, 18), y_lim = c(0, 12),
    clock_file = "<LOAD DATA>",
    cox_file = "<LOAD DATA>"
  )
)

assay_thresholds <- list(
  "C2N Diagnostics PrecivityAD2 %p-tau217" = 4.06,
  "C2N Diagnostics PrecivityAD2 p-tau217" = 2.34,
  "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42" = 0.006312,
  "Fujirebio Diagnostics Lumipulse p-tau217" = 0.158,
  "Janssen LucentAD Quanterix p-tau217" = 0.0615,
  "ALZpath Quanterix p-tau217" = 0.444
)

assay_plasma_params <- list(
  "C2N Diagnostics PrecivityAD2 %p-tau217" = list(min = 1.1, max = 10.3, value = 7.0, step = 0.1),
  "C2N Diagnostics PrecivityAD2 p-tau217" = list(min = 1.3, max = 7.9, value = 4, step = 0.1),
  "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42" = list(min = 0.001, max = 0.018, value = 0.01, step = 0.001),
  "Fujirebio Diagnostics Lumipulse p-tau217" = list(min = 0.03, max = 0.49, value = 0.3, step = 0.001),
  "Janssen LucentAD Quanterix p-tau217" = list(min = 0.012, max = 0.156, value = 0.1, step = 0.001),
  "ALZpath Quanterix p-tau217" = list(min = 0.11, max = 1.10, value = 0.7, step = 0.01)
)

get_param_key <- function(assay, cohort, method) {
  paste(assay, cohort, method, sep = "_")
}

get_available_assays <- function(cohort) {
  if (cohort == "KADRC") {
    return(c("C2N Diagnostics PrecivityAD2 %p-tau217" = "C2N Diagnostics PrecivityAD2 %p-tau217"))
  } else if (cohort == "ADNI") {
    return(c(
      "C2N Diagnostics PrecivityAD2 %p-tau217" = "C2N Diagnostics PrecivityAD2 %p-tau217",
      "C2N Diagnostics PrecivityAD2 p-tau217" = "C2N Diagnostics PrecivityAD2 p-tau217",
      "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42" = "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42",
      "Fujirebio Diagnostics Lumipulse p-tau217" = "Fujirebio Diagnostics Lumipulse p-tau217",
      "Janssen LucentAD Quanterix p-tau217" = "Janssen LucentAD Quanterix p-tau217",
      "ALZpath Quanterix p-tau217" = "ALZpath Quanterix p-tau217"
    ))
  }
  return(c())
}

get_assay_threshold <- function(assay) {
  threshold <- assay_thresholds[[assay]]
  if (is.null(threshold)) return(4.06)
  return(threshold)
}

create_smooth_clock_function <- function(clock_data) {
  clock_data <- clock_data[!duplicated(clock_data$plasma), ]
  clock_data <- clock_data[order(clock_data$plasma), ]
  plasma_range <- range(clock_data$plasma, na.rm = TRUE)
  smooth_fit <- smooth.spline(clock_data$plasma, clock_data$clock_time, spar = 0.5)
  return(function(plasma_values) {
    result <- predict(smooth_fit, plasma_values)$y
    outside_range <- plasma_values < plasma_range[1] | plasma_values > plasma_range[2]
    result[outside_range] <- NA
    return(result)
  })
}

load_model_data <- function(assay, cohort, method) {
  param_key <- get_param_key(assay, cohort, method)
  params <- parameter_sets[[param_key]]
  if (is.null(params)) {
    stop("Invalid assay-cohort-method combination: ", assay, "-", cohort, "-", method)
  }
  
  tryCatch({
    clock_data <- read.csv(params$clock_file)
    if (method == "SILA") {
      if ("adtime" %in% names(clock_data) && !"clock_time" %in% names(clock_data)) {
        clock_data <- clock_data %>% mutate(clock_time = adtime)
      }
    }
    if (!"plasma" %in% names(clock_data) || !"clock_time" %in% names(clock_data)) {
      stop("Clock data file must contain 'plasma' and 'clock_time' columns")
    }
    clock_data <- clock_data[complete.cases(clock_data[c("plasma", "clock_time")]), ]
    smooth_clock_fn <- create_smooth_clock_function(clock_data)
  }, error = function(e) {
    stop("Failed to load clock file: ", params$clock_file, ". Error: ", e$message)
  })
  
  tryCatch({
    cox_model <- readRDS(params$cox_file)
  }, error = function(e) {
    stop("Failed to load Cox model file: ", params$cox_file, ". Error: ", e$message)
  })
  
  params$plasma_min <- min(clock_data$plasma, na.rm = TRUE)
  params$plasma_max <- max(clock_data$plasma, na.rm = TRUE)
  
  list(
    clock_data = clock_data,
    cox_model = cox_model,
    params = params,
    smooth_clock_fn = smooth_clock_fn
  )
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
      
      :root {
        --primary-red: #BA0C2F;
        --highlight-red: #dc2626;
        --light-red: #fecaca;
        --dark-red: #8b0000;
        --accent-orange: #f97316;
        --success-green: black;
        --neutral-gray: #6b7280;
        --light-gray: #f9fafb;
      }
      
      body {
        font-family: 'Inter', Arial, sans-serif;
        background: linear-gradient(135deg, #fef2f2 0%, #f9fafb 100%);
        color: #1f2937;
        margin: 0;
        padding: 0;
      }
      
      .custom-header {
        background: linear-gradient(135deg, var(--primary-red) 0%, var(--dark-red) 100%);
        padding: 25px 0;
        box-shadow: 0 6px 20px rgba(186, 12, 47, 0.3);
        border-bottom: 4px solid var(--dark-red);
      }
      
      .header-content {
        max-width: 1400px;
        margin: 0 auto;
        display: flex;
        align-items: center;
        justify-content: space-between;
        padding: 0 30px;
      }
      
      .title-section {
        display: flex;
        flex-direction: column;
      }
      
      .app-title {
        color: white;
        font-size: 28px;
        font-weight: 700;
        margin: 0;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.4);
        line-height: 1.2;
      }
      
      .subtitle {
        color: rgba(255, 255, 255, 0.95);
        font-size: 16px;
        font-weight: 500;
        margin: 8px 0 0 0;
        text-shadow: 1px 1px 2px rgba(0, 0, 0, 0.3);
      }
      
      .header-nav {
        display: flex;
        gap: 20px;
      }
      
      .nav-button {
        background: rgba(255, 255, 255, 0.25);
        color: white;
        border: 2px solid rgba(255, 255, 255, 0.6);
        padding: 12px 24px;
        border-radius: 12px;
        text-decoration: none;
        font-weight: 600;
        font-size: 16px;
        transition: all 0.3s ease;
        cursor: pointer;
        backdrop-filter: blur(10px);
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.15);
        text-shadow: 1px 1px 2px rgba(0, 0, 0, 0.3);
      }
      
      .nav-button:hover {
        background: rgba(255, 255, 255, 0.4);
        border-color: rgba(255, 255, 255, 0.8);
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0, 0, 0, 0.25);
      }
      
      .nav-button.active {
        background: white;
        color: var(--primary-red);
        border-color: white;
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2);
        text-shadow: none;
        font-weight: 700;
      }
      
      .container-fluid {
        max-width: 1400px;
        margin: 0 auto;
        padding: 30px 20px;
      }
      
      .card {
        background: white;
        border-radius: 12px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
        border: 1px solid #e5e7eb;
        margin-bottom: 25px;
        transition: all 0.3s ease;
      }
      
      .card:hover {
        box-shadow: 0 8px 25px rgba(186, 12, 47, 0.15);
        transform: translateY(-2px);
      }
      
      .card-header {
        background: linear-gradient(90deg, var(--primary-red) 0%, var(--highlight-red) 100%);
        color: white;
        padding: 20px;
        border-radius: 12px 12px 0 0;
        font-weight: 600;
        font-size: 18px;
        border: none;
      }
      
      .card-body {
        padding: 25px;
      }
      
      .reset-btn {
        background: linear-gradient(90deg, var(--primary-red) 0%, var(--dark-red) 100%);
        color: white;
        border: none;
        width: 100%;
        padding: 12px 20px;
        margin-bottom: 25px;
        font-weight: 600;
        border-radius: 8px;
        cursor: pointer;
        font-size: 14px;
        transition: all 0.3s ease;
        box-shadow: 0 2px 4px rgba(186, 12, 47, 0.3);
      }
      
      .reset-btn:hover {
        background: linear-gradient(90deg, var(--dark-red) 0%, #660000 100%);
        transform: translateY(-1px);
        box-shadow: 0 4px 8px rgba(186, 12, 47, 0.4);
      }
      
      .section-header {
        font-weight: 600;
        font-size: 16px;
        margin-bottom: 15px;
        margin-top: 25px;
        color: var(--primary-red);
        border-bottom: 2px solid var(--light-red);
        padding-bottom: 8px;
        display: flex;
        align-items: center;
      }
      
      .section-header:first-of-type {
        margin-top: 0;
      }
      
      .section-icon {
        margin-right: 8px;
        font-size: 18px;
      }
      
      .form-group {
        margin-bottom: 20px;
        width: 100%;
      }
      
      .shiny-input-select {
        width: 100%;
        padding: 12px 16px;
        border: 2px solid #e5e7eb;
        border-radius: 8px;
        font-size: 14px;
        background: #fafafa;
        color: #374151;
        transition: all 0.3s ease;
        appearance: none;
        background-image: url('data:image/svg+xml;utf8,<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 24 24\" fill=\"none\" stroke=\"%236b7280\" stroke-width=\"2\" stroke-linecap=\"round\" stroke-linejoin=\"round\"><polyline points=\"6,9 12,15 18,9\"></polyline></svg>');
        background-repeat: no-repeat;
        background-position: right 12px center;
        background-size: 16px;
        padding-right: 45px;
      }
      
      .shiny-input-select:focus {
        border-color: var(--primary-red);
        box-shadow: 0 0 0 3px rgba(186, 12, 47, 0.1);
        background: white;
        outline: none;
      }
      
      .form-control {
        border: 2px solid #e5e7eb;
        border-radius: 8px;
        padding: 12px 16px;
        font-size: 14px;
        transition: all 0.3s ease;
        background: #fafafa;
      }
      
      .form-control:focus {
        border-color: var(--primary-red);
        box-shadow: 0 0 0 3px rgba(186, 12, 47, 0.1);
        background: white;
      }
      
      .prediction-results {
        background: linear-gradient(135deg, #fef2f2 0%, #fef7f7 100%);
        padding: 20px;
        border-radius: 12px;
        margin: 20px 0;
        border: 2px solid var(--light-red);
        box-shadow: 0 2px 4px rgba(186, 12, 47, 0.1);
      }
      
      .prediction-item {
        margin-bottom: 12px;
        font-size: 14px;
        display: flex;
        justify-content: space-between;
        align-items: center;
        width: 100%;
      }
      
      .prediction-label {
        font-weight: 500;
        color: var(--neutral-gray);
        flex: 1;
        text-align: left;
      }
      
      .prediction-value {
        font-weight: 700;
        color: var(--primary-red);
        font-size: 16px;
        text-align: right;
        min-width: 80px;
      }
      
      .explanatory-text {
        margin-top: 25px;
        padding: 20px;
        background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%);
        border-left: 4px solid #0ea5e9;
        color: #0c4a6e;
        font-size: 14px;
        border-radius: 0 8px 8px 0;
        line-height: 1.6;
      }
      
      .plot-section-title {
        font-size: 20px;
        font-weight: 700;
        margin: 30px 0 20px 0;
        color: var(--primary-red);
        border-bottom: 3px solid var(--primary-red);
        padding-bottom: 10px;
        display: flex;
        align-items: center;
      }
      
      .plot-container {
        margin-bottom: 30px;
        border: 2px solid #e5e7eb;
        border-radius: 12px;
        overflow: hidden;
        background: white;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
      }
      
      .comparison-table {
        margin: 40px auto;
        max-width: 1200px;
        background: white;
        border-radius: 12px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
        border: 2px solid var(--light-red);
        overflow: hidden;
      }
      
      .comparison-table .plot-section-title {
        text-align: center;
        margin: 20px 0;
        padding: 0 20px;
      }
      
      .table-container {
        display: flex;
        justify-content: center;
        padding: 0 20px 20px 20px;
      }
      
      .table {
        width: auto;
        max-width: 100%;
        margin: 0 auto;
        color: #374151;
      }
      
      .table th {
        background: linear-gradient(90deg, var(--primary-red) 0%, var(--highlight-red) 100%);
        color: white;
        font-weight: 600;
        padding: 15px;
        border: none;
        text-align: center;
      }
      
      .table td {
        padding: 12px 15px;
        border-top: 1px solid #e5e7eb;
        text-align: center;
        font-weight: 500;
      }
      
      .table tbody tr:hover {
        background-color: #fef2f2;
      }
      
      .content-section {
        display: none;
      }
      
      .content-section.active {
        display: block;
      }
      
      .animation-section {
        background: white;
        border-radius: 12px;
        padding: 25px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.05);
        border: 1px solid #e5e7eb;
        height: fit-content;
        margin-top: 0;
      }
      
      .animation-section h3 {
        font-weight: 600;
        margin-bottom: 20px;
      }
      
      .animation-section img {
        transition: transform 0.3s ease;
      }
      
      .animation-section img:hover {
        transform: scale(1.02);
      }
      
      .podcast-link {
        color: var(--primary-red);
        text-decoration: none;
        font-weight: 600;
      }
      
      .podcast-link:hover {
        text-decoration: underline;
      }
      
      .footer {
        background: linear-gradient(135deg, var(--primary-red) 0%, var(--dark-red) 100%);
        color: white;
        padding: 40px 20px;
        text-align: center;
        margin-top: 50px;
        box-shadow: 0 -4px 20px rgba(186, 12, 47, 0.3);
      }

      .footer-logo {
        margin-bottom: 25px;
        display: flex;
        justify-content: center;
        align-items: center;
        gap: 40px;
        padding: 20px;
        background: rgba(255, 255, 255, 0.1);
        border-radius: 12px;
        backdrop-filter: blur(10px);
        max-width: 600px;
        margin: 0 auto 25px auto;
      }

      .footer-logo a {
        display: block;
        transition: all 0.3s ease;
      }

      .footer-logo a:hover {
        transform: scale(1.1);
      }

      .footer-logo img {
        height: 70px;
        filter: brightness(0) invert(1) drop-shadow(2px 2px 4px rgba(0, 0, 0, 0.3));
        transition: all 0.3s ease;
      }

      .footer-text {
        font-family: 'Inter', Arial, sans-serif;
        font-size: 18px;
        font-weight: 600;
        line-height: 1.4;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.3);
        letter-spacing: 0.5px;
      }

      .footer-institution {
        font-size: 20px;
        font-weight: 700;
        margin-bottom: 8px;
      }

      .footer-department {
        font-size: 16px;
        font-weight: 500;
        opacity: 0.95;
      }
    "))
  ),
  
  div(class = "custom-header",
      div(class = "header-content",
          div(class = "title-section",
              h1("Plasma p-tau217 Biomarker Clocks", class = "app-title"),
              p("Time since biomarker positivity and risk for symptomatic Alzheimer's disease", class = "subtitle")
          ),
          div(class = "header-nav",
              tags$a("Analysis", class = "nav-button active", href = "#", onclick = "showSection('analysis')"),
              tags$a("References", class = "nav-button", href = "#", onclick = "showSection('references')")
          )
      )
  ),
  
  tags$script(HTML("
    function showSection(sectionName) {
      document.querySelectorAll('.content-section').forEach(function(section) {
        section.classList.remove('active');
      });
      document.getElementById(sectionName + '-section').classList.add('active');
      document.querySelectorAll('.nav-button').forEach(function(button) {
        button.classList.remove('active');
      });
      event.target.classList.add('active');
    }
  ")),
  
  div(id = "analysis-section", class = "content-section active",
      div(class = "container-fluid",
          fluidRow(
            column(5,
                   div(class = "card",
                       div(class = "card-header",
                           HTML("🩸 Biomarker Analysis Tool")
                       ),
                       div(class = "card-body",
                           actionButton("reset", "🔄 Reset All Inputs", class = "reset-btn"),
                           
                           div(class = "section-header",
                               span(class = "section-icon", "🏥"),
                               "Cohort Selection"
                           ),
                           div(class = "form-group",
                               selectInput("cohort", NULL,
                                           choices = c("Knight ADRC" = "KADRC", "ADNI" = "ADNI"),
                                           selected = "ADNI")
                           ),
                           
                           div(class = "section-header",
                               span(class = "section-icon", "🧬"),
                               "Assay selection"
                           ),
                           
                           div(class = "form-group",
                               selectInput("assay", NULL,
                                           choices = get_available_assays("ADNI"),
                                           selected = "C2N Diagnostics PrecivityAD2 %p-tau217")
                           ),
                           
                           div(class = "section-header",
                               span(class = "section-icon", "⚙️"),
                               "Method Selection"
                           ),
                           div(class = "form-group",
                               selectInput("method", NULL,
                                           choices = c("TIRA" = "TIRA", "SILA" = "SILA"),
                                           selected = "TIRA")
                           ),
                           
                           div(class = "section-header",
                               span(class = "section-icon", "👤"),
                               "Individual Information"
                           ),
                           div(class = "form-group",
                               numericInput("plasma",
                                            HTML("🩸 Plasma level"),
                                            value = 7.0, min = 1.1, max = 10, step = 0.1)
                           ),
                           div(class = "form-group",
                               numericInput("age_plasma",
                                            HTML("📅 Age at plasma collection"),
                                            value = 80, min = 60, max = 100, step = 1)
                           ),
                           div(class = "form-group",
                               numericInput("time_since_plasma",
                                            HTML("🕐 Time since plasma collection (years)"),
                                            value = 5, min = 0, max = 40, step = 1)
                           ),
                           
                           div(class = "section-header",
                               span(class = "section-icon", "📊"),
                               "Plasma clock estimates"
                           ),
                           div(class = "prediction-results",
                               div(class = "prediction-item",
                                   span(class = "prediction-label", "Age at biomarker positivity:"),
                                   span(class = "prediction-value", textOutput("onset_age_display", inline = TRUE))
                               ),
                               div(class = "prediction-item",
                                   span(class = "prediction-label", "Time since biomarker positivity at plasma collection:"),
                                   span(class = "prediction-value", textOutput("clock_time_display", inline = TRUE))
                               ),
                               div(class = "prediction-item",
                                   span(class = "prediction-label", "Time since biomarker positivity at follow-up:"),
                                   span(class = "prediction-value", textOutput("years_since_onset_display", inline = TRUE))
                               )
                           ),
                           
                           div(class = "explanatory-text",
                               HTML("ℹ️ This tool uses plasma biomarker values to predict time since biomarker positivity and probability of developing symptomatic Alzheimer's disease")
                           )
                       )
                   )
            ),
            
            column(7,
                   div(class = "card",
                       div(class = "card-body",
                           div(class = "plot-section-title",
                               HTML("⏰ Clock Visualization")),
                           div(class = "plot-container",
                               plotOutput("clock_plot", height = "500px")),
                           div(class = "plot-section-title",
                               HTML("📈 Probability of developing symptomatic Alzheimer's disease")),
                           div(class = "plot-container",
                               plotOutput("survival_plot", height = "500px"))
                       )
                   )
            )
          ),
          
          fluidRow(
            column(12,
                   div(class = "comparison-table",
                       div(class = "plot-section-title",
                           HTML("🔍 Results for selected biomarker analysis")),
                       div(class = "table-container",
                           tableOutput("comparison_table")
                       )
                   )
            )
          )
      )
  ),
  
  div(id = "references-section", class = "content-section",
      div(class = "container-fluid",
          fluidRow(
            column(8,
                   div(class = "references-section",
                       h2("References", style = "color: var(--primary-red); margin-bottom: 30px;"),
                       div(class = "reference-container",
                           p(class = "hangingindent",
                             strong("1. "), "Schindler SE, Li Y, Buckles VD, Gordon BA, Benzinger TLS, Wang G, Coble D, Klunk WE, Fagan AM, Holtzman DM, Bateman RJ, Morris JC, Xiong C. Predicting Symptom Onset in Sporadic Alzheimer Disease With Amyloid PET. Neurology. 2021 Nov 2;97(18):e1823-e1834. doi: 10.1212/WNL.0000000000012775. PMID: 34504028. ",
                             tags$a("https://pubmed.ncbi.nlm.nih.gov/34504028/", href = "https://pubmed.ncbi.nlm.nih.gov/34504028/", target = "_blank", class = "podcast-link")),
                           p(class = "hangingindent",
                             strong("2. "), "Betthauser TJ, Bilgel M, Koscik RL, Jedynak BM, An Y, Kellett KA, Moghekar A, Jonaitis EM, Stone CK, Engelman CD, Asthana S, Christian BT, Wong DF, Albert M, Resnick SM, Johnson SC; Alzheimer's Disease Neuroimaging Initiative. Multi-method investigation of factors influencing amyloid onset and impairment in three cohorts. Brain. 2022 Nov 21;145(11):4065-4079. doi: 10.1093/brain/awac213. PMID: 35856240. ",
                             tags$a("https://pubmed.ncbi.nlm.nih.gov/35856240/", href = "https://pubmed.ncbi.nlm.nih.gov/35856240/", target = "_blank", class = "podcast-link")),
                           p(class = "hangingindent",
                             strong("3. "), "Li Y, Yen D, Hendrix RD, Gordon BA, Dlamini S, Barthélemy NR, Aschenbrenner AJ, Henson RL, Herries EM, Volluz K, Kirmess K, Eastwood S, Meyer M, Heller M, Jarrett L, McDade E, Holtzman DM, Benzinger TLS, Morris JC, Bateman RJ, Xiong C, Schindler SE. Timing of Biomarker Changes in Sporadic Alzheimer's Disease in Estimated Years from Symptom Onset. Ann Neurol. 2024 May;95(5):951-965. doi: 10.1002/ana.26891. PMID: 38400792. ",
                             tags$a("https://pubmed.ncbi.nlm.nih.gov/38400792/", href = "https://pubmed.ncbi.nlm.nih.gov/38400792/", target = "_blank", class = "podcast-link")),
                           p(class = "hangingindent",
                             strong("4. "), "Milà-Alomà M, Tosun D, Schindler SE, Hausle I, Petersen KK, Li Y, Dage JL, Du-Cuny L, Saad ZS, Saef B, Triana-Baltzer G, Raunig DL, Coomaraswamy J, Baratta M, Meyers EA, Mordashova Y, Rubel CE, Ferber K, Kolb H, Ashton NJ, Zetterberg H, Rosenbaugh EG, Sabandal M, Shaw LM, Bannon AW, Potter WZ; Alzheimer's Disease Neuroimaging Initiative (ADNI); Foundation for the National Institutes of Health (FNIH) Biomarkers Consortium Plasma Aβ and Phosphorylated Tau as Predictors of Amyloid and Tau Positivity in Alzheimer's Disease Project Team. Timing of Changes in Alzheimer's Disease Plasma Biomarkers as Assessed by Amyloid and Tau PET Clocks. Ann Neurol. 2025 Jun 20. doi: 10.1002/ana.27285. Epub ahead of print. PMID: 40539416. ",
                             tags$a("https://pubmed.ncbi.nlm.nih.gov/40539416/", href = "https://pubmed.ncbi.nlm.nih.gov/40539416/", target = "_blank", class = "podcast-link")),
                           p(class = "hangingindent",
                             strong("5. "), "Schindler SE, Petersen KK, Saef B, Tosun D, Shaw LM, Zetterberg H, Dage JL, Ferber K, Triana-Baltzer G, Du-Cuny L, Li Y, Coomaraswamy J, Baratta M, Mordashova Y, Saad ZS, Raunig DL, Ashton NJ, Meyers EA, Rubel CE, Rosenbaugh EG, Bannon AW, Potter WZ; et al. Head-to-head comparison of leading blood tests for Alzheimer's disease pathology. Alzheimers Dement. 2024 Nov;20(11):8074-8096. doi: 10.1002/alz.14315. PMID: 39394841. ",
                             tags$a("https://pubmed.ncbi.nlm.nih.gov/39394841/", href = "https://pubmed.ncbi.nlm.nih.gov/39394841/", target = "_blank", class = "podcast-link")),
                           p(class = "hangingindent",
                             strong("6. "), "Petersen KK, Milà-Alomà ML, Yan; Du, Lianlian; Xiong, Chengjie; Tosun, Duygu; Saef, Benjamin; Saad, Ziad S.; Du-Cuny, Lei; Coomaraswamy, Janaky; Mordashova, Yulia; Rubel, Carrie E.; Meyers, Emily A.; Shaw, Leslie M.; Dage, Jeffrey L.; Ashton, Nicholas J.; Zetterberg, Henrik; Ferber, Kyle; Triana-Baltzer, Gallen; Baratta, Michael; Rosenbaugh, Erin G.; Cruchaga, Carlos; McDade, Eric; Holtzman, David M.; Morris, John C.; Sabandal, J. Martin; Bateman, Randall J.; Bannon, Anthony W.; Potter, William Z.; Schindler, Suzanne E.; Alzheimer's Disease Neuroimaging Initiative (ADNI); Foundation for the National Institutes of Health (FNIH) Biomarkers Consortium Plasma Aβ and Phosphorylated Tau as Predictors of Amyloid and Tau Positivity in Alzheimer's Disease Project Team. Predicting onset of symptomatic Alzheimer disease with a plasma %p-tau217 clock. Research Square. 09 July 2025 2025;[Preprint]. Version 1doi:10.21203/rs.3.rs-7059258/v1. ",
                             tags$a("https://www.researchsquare.com/article/rs-7059258/v1", href = "https://www.researchsquare.com/article/rs-7059258/v1", target = "_blank", class = "podcast-link"))
                       ),
                       h3("Methodology", style = "color: var(--primary-red); margin-top: 40px; margin-bottom: 20px;"),
                       p("This biomarker clock tool implements advanced statistical modeling techniques to predict the onset of Alzheimer's disease symptoms based on plasma biomarker levels."),
                       h3("Data Sources", style = "color: var(--primary-red); margin-top: 30px; margin-bottom: 20px;"),
                       p("Models are trained on data from the Knight Alzheimer Disease Research Center (Knight ADRC) and the Alzheimer's Disease Neuroimaging Initiative (ADNI).")
                   )
            ),
            
            column(4,
                   div(class = "animation-section",
                       h3("", style = "color: var(--primary-red); margin-bottom: 20px; text-align: center;"),
                       div(style = "text-align: center; padding: 20px;",
                           img(src = "trajectory_animation_hd_450x300.gif",
                               style = "max-width: 100%; height: auto; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);",
                               alt = "")
                       ),
                       p(style = "text-align: center; font-size: 14px; color: var(--neutral-gray); margin-top: 15px;",
                         "Visualization showing plasma trajectories aligning when plotted against time from biomarker positivity.")
                   )
            )
          )
      )
  ),
  
  div(class = "footer",
      div(class = "footer-logo",
          tags$a(href = "https://medicine.washu.edu/", target = "_blank",
                 img(src = "washu_logo.png", height = "70px", alt = "Washington University")),
          tags$a(href = "https://knightadrc.wustl.edu/", target = "_blank",
                 img(src = "knight_adrc_logo.png", height = "70px", alt = "Knight ADRC"))
      )
  )
)

server <- function(input, output, session) {
  
  assay_plasma_params <- list(
    "C2N Diagnostics PrecivityAD2 %p-tau217" = list(min = 1.1, max = 10.3, value = 7.0, step = 0.1),
    "C2N Diagnostics PrecivityAD2 p-tau217" = list(min = 1.3, max = 7.9, value = 4, step = 0.1),
    "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42" = list(min = 0.001, max = 0.018, value = 0.01, step = 0.001),
    "Fujirebio Diagnostics Lumipulse p-tau217" = list(min = 0.03, max = 0.49, value = 0.3, step = 0.001),
    "Janssen LucentAD Quanterix p-tau217" = list(min = 0.012, max = 0.156, value = 0.1, step = 0.001),
    "ALZpath Quanterix p-tau217" = list(min = 0.11, max = 1.10, value = 0.7, step = 0.01)
  )
  
  get_linear_model_filename <- function(cohort, assay, method) {
    assay_file_map <- list(
      "C2N Diagnostics PrecivityAD2 %p-tau217" = "C2N_ptau217ratio",
      "C2N Diagnostics PrecivityAD2 p-tau217" = "C2N_ptau217",
      "Fujirebio Diagnostics Lumipulse p-tau217/Aβ42" = "Fuji_ptau217_Abeta42",
      "Fujirebio Diagnostics Lumipulse p-tau217" = "Fuji_ptau217",
      "Janssen LucentAD Quanterix p-tau217" = "Janssen_ptau217",
      "ALZpath Quanterix p-tau217" = "AlzPath_ptau217"
    )
    
    assay_suffix <- assay_file_map[[assay]]
    if (is.null(assay_suffix)) return(NULL)
    
    filename <- paste0("linear_models_", cohort, "_", method, "_CDR_DX_Imp_FINAL_", assay_suffix, ".csv")
    return(file.path("data", filename))
  }
  
  load_linear_model_params <- function(cohort, assay, method) {
    filename <- get_linear_model_filename(cohort, assay, method)
    
    if (is.null(filename) || !file.exists(filename)) {
      return(list(intercept = NA, slope = NA))
    }
    
    tryCatch({
      linear_data <- read.csv(filename)
      list(
        intercept = linear_data$intercept[1],
        slope = linear_data$slope[1]
      )
    }, error = function(e) {
      cat("Error loading linear model file:", filename, ". Error:", e$message, "\n")
      list(intercept = NA, slope = NA)
    })
  }
  
  observe({
    available_assays <- get_available_assays("ADNI")
    updateSelectInput(session, "assay",
                      choices = available_assays,
                      selected = "C2N Diagnostics PrecivityAD2 %p-tau217")
  })
  
  observeEvent(input$cohort, {
    available_assays <- get_available_assays(input$cohort)
    updateSelectInput(session, "assay",
                      choices = available_assays,
                      selected = available_assays[[1]])
  })
  
  observeEvent(input$assay, {
    req(input$assay)
    
    params <- assay_plasma_params[[input$assay]]
    
    if (!is.null(params)) {
      updateNumericInput(session, "plasma",
                         label = HTML(paste0("🩸 Plasma level (", input$assay, ")")),
                         min = params$min,
                         max = params$max,
                         value = params$value,
                         step = params$step
      )
    }
  }, ignoreInit = FALSE)
  
  model_data <- reactive({
    req(input$assay, input$cohort, input$method)
    tryCatch({
      load_model_data(input$assay, input$cohort, input$method)
    }, error = function(e) {
      cat("Error loading model data:", e$message, "\n")
      return(NULL)
    })
  })
  
  params <- reactive({
    req(model_data())
    model_data()$params
  })
  
  observeEvent(input$reset, {
    updateSelectInput(session, "cohort", selected = "ADNI")
    updateSelectInput(session, "assay", selected = "C2N Diagnostics PrecivityAD2 %p-tau217")
    updateSelectInput(session, "method", selected = "TIRA")
    updateNumericInput(session, "age_plasma", value = 80)
    updateNumericInput(session, "time_since_plasma", value = 5)
    updateNumericInput(session, "plasma", value = 7.0, step = 0.1, min = 1.1, max = 10.3)
  })
  
  get_clock_time <- reactive({
    req(model_data())
    function(plasma_level) {
      smooth_fn <- model_data()$smooth_clock_fn
      result <- smooth_fn(plasma_level)
      if (is.na(result)) return(0)
      return(result)
    }
  })
  
  create_survival_curve <- reactive({
    req(model_data())
    function(est_onset_age, age_range = c(50, 100)) {
      time_points <- seq(age_range[1], age_range[2], by = 0.5)
      newdata <- data.frame(est_onset_age = est_onset_age)
      fit <- model_data()$cox_model
      surv_probs <- sapply(time_points, function(t) {
        tryCatch({
          1 - getFitEsts(fit, newdata = newdata, q = t)
        }, error = function(e) {
          exp(-0.1 * (t - est_onset_age))
        })
      })
      surv_probs <- pmax(0, pmin(1, surv_probs))
      data.frame(time = time_points, surv = surv_probs)
    }
  })
  
  calculations <- reactive({
    req(model_data(), params())
    req(input$plasma, input$age_plasma, input$time_since_plasma)
    
    current_age <- input$age_plasma + input$time_since_plasma
    
    assay_params <- assay_plasma_params[[input$assay]]
    
    shiny::validate(
      need(input$time_since_plasma >= 0, "Time since plasma collection must be ≥ 0"),
      need(!is.null(assay_params), "Invalid assay selected"),
      need(input$plasma >= assay_params$min && input$plasma <= assay_params$max,
           paste("Plasma level must be between", assay_params$min, "and", assay_params$max, "for", input$assay))
    )
    
    clock_time <- get_clock_time()(input$plasma)
    est_onset_age <- input$age_plasma - clock_time
    years_since_onset_current <- current_age - est_onset_age
    list(
      clock_time = clock_time,
      est_onset_age = est_onset_age,
      years_since_onset_current = years_since_onset_current,
      current_age = current_age
    )
  })
  
  output$onset_age_display <- renderText({
    calc <- calculations()
    sprintf("%.1f years", calc$est_onset_age)
  })
  
  output$clock_time_display <- renderText({
    calc <- calculations()
    sprintf("%.1f years", calc$clock_time)
  })
  
  output$years_since_onset_display <- renderText({
    calc <- calculations()
    sprintf("%.1f years", calc$years_since_onset_current)
  })
  
  output$clock_plot <- renderPlot({
    req(params(), calculations(), model_data())
    p <- params()
    calc <- calculations()
    
    plasma_seq <- seq(p$plasma_min, p$plasma_max, length = 1000)
    clock_times <- model_data()$smooth_clock_fn(plasma_seq)
    clock_plot_data <- data.frame(plasma = plasma_seq, clock_time = clock_times)
    clock_plot_data <- clock_plot_data[complete.cases(clock_plot_data), ]
    clock_plot_data <- clock_plot_data[
      clock_plot_data$clock_time >= p$x_lim[1] &
        clock_plot_data$clock_time <= p$x_lim[2] &
        clock_plot_data$plasma >= p$plasma_min &
        clock_plot_data$plasma <= p$plasma_max, ]
    
    current_time_plasma <- NA
    if (!is.na(calc$years_since_onset_current)) {
      time_diffs <- abs(clock_times - calc$years_since_onset_current)
      if (min(time_diffs, na.rm = TRUE) < Inf) {
        closest_idx <- which.min(time_diffs)
        current_time_plasma <- plasma_seq[closest_idx]
      }
    }
    
    plasma_label_y <- input$plasma + (p$plasma_max - p$plasma_min) * 0.08
    
    if (!is.na(current_time_plasma)) {
      current_label_y <- current_time_plasma - (p$plasma_max - p$plasma_min) * 0.08
      current_label_y <- max(current_label_y, p$plasma_min + (p$plasma_max - p$plasma_min) * 0.05)
    } else {
      current_label_y <- input$plasma - (p$plasma_max - p$plasma_min) * 0.08
    }
    
    plasma_label_y <- min(plasma_label_y, p$plasma_max * 0.95)
    
    plot <- ggplot(clock_plot_data, aes(x = clock_time, y = plasma)) +
      geom_line(color = "#BA0C2F", linewidth = 2, alpha = 0.8) +
      
      annotate("point", x = calc$clock_time, y = input$plasma,
               color = "black", size = 5, shape = 21, stroke = 3, fill = "#dcfce7") +
      
      geom_vline(xintercept = 0, color = "#374151", linetype = "dashed", alpha = 0.7, linewidth = 1) +
      geom_vline(xintercept = calc$clock_time, color = "black", linetype = "dashed", alpha = 0.7, linewidth = 1) +
      
      geom_hline(yintercept = get_assay_threshold(input$assay), color = "#374151", linetype = "dashed", alpha = 0.7, linewidth = 1) +
      
      annotate("text", x = calc$clock_time, y = plasma_label_y,
               label = "Plasma\nCollection", color = "black", fontface = "bold", hjust = 0.5, size = 4) +
      
      labs(
        title = paste0(ifelse(input$cohort == "KADRC", "Knight ADRC", "ADNI"), " - ", input$method),
        subtitle = paste("Assay:", input$assay),
        x = "Estimated time since biomarker positivity (years)",
        y = "Plasma level"
      ) +
      coord_cartesian(xlim = p$x_lim, ylim = c(p$plasma_min, p$plasma_max)) +
      theme_cowplot() +
      theme(
        plot.title = element_text(size = 18, face = "bold", color = "#BA0C2F"),
        plot.subtitle = element_text(size = 14, color = "#6b7280"),
        axis.title = element_text(size = 12, face = "bold", color = "#374151"),
        axis.text = element_text(size = 11, color = "#6b7280"),
        panel.grid.major = element_line(color = "#f3f4f6", linewidth = 0.5),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    if (!is.na(current_time_plasma)) {
      plot <- plot +
        annotate("point", x = calc$years_since_onset_current, y = current_time_plasma,
                 color = "#BA0C2F", size = 5, shape = 21, stroke = 3, fill = "#fecaca") +
        geom_vline(xintercept = calc$years_since_onset_current, color = "#BA0C2F", linetype = "dashed", alpha = 0.7, linewidth = 1) +
        annotate("text", x = calc$years_since_onset_current, y = current_label_y,
                 label = "Follow-up", color = "#BA0C2F", fontface = "bold", hjust = 0.5, size = 4)
    }
    
    plot
  })
  
  output$survival_plot <- renderPlot({
    req(calculations(), create_survival_curve(), model_data())
    calc <- calculations()
    surv_data <- create_survival_curve()(calc$est_onset_age)
    surv_data <- surv_data[complete.cases(surv_data), ]
    
    plasma_collection_idx <- which.min(abs(surv_data$time - input$age_plasma))
    current_age_idx <- which.min(abs(surv_data$time - calc$current_age))
    plasma_collection_surv <- surv_data$surv[plasma_collection_idx]
    current_age_surv <- surv_data$surv[current_age_idx]
    
    plasma_label_y <- plasma_collection_surv + 0.1
    current_label_y <- current_age_surv - 0.1
    
    plasma_label_y <- min(plasma_label_y, 0.95)
    current_label_y <- max(current_label_y, 0.05)
    
    plot <- ggplot(surv_data, aes(x = time, y = surv)) +
      geom_line(color = "#BA0C2F", linewidth = 2) +
      geom_ribbon(aes(ymin = pmax(0, surv - 0.1), ymax = pmin(1, surv + 0.1)),
                  alpha = 0.2, fill = "#fecaca") +
      annotate("point", x = input$age_plasma, y = plasma_collection_surv,
               color = "black", size = 5, shape = 21, stroke = 3, fill = "#dcfce7") +
      annotate("point", x = calc$current_age, y = current_age_surv,
               color = "#BA0C2F", size = 5, shape = 21, stroke = 3, fill = "#fecaca") +
      geom_vline(xintercept = input$age_plasma, color = "black", linetype = "dashed", alpha = 0.7, linewidth = 1) +
      geom_vline(xintercept = calc$current_age, color = "#BA0C2F", linetype = "dashed", alpha = 0.7, linewidth = 1) +
      
      annotate("text", x = input$age_plasma, y = plasma_label_y,
               label = "Plasma\nCollection", color = "black", fontface = "bold", hjust = 0.5, size = 4) +
      annotate("text", x = calc$current_age, y = current_label_y,
               label = "Follow-up", color = "#BA0C2F", fontface = "bold", hjust = 0.5, size = 4) +
      
      scale_x_continuous(limits = c(50, 100), breaks = seq(50, 100, 10)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent) +
      labs(
        title = "",
        subtitle = sprintf("Based on estimated age at biomarker positivity of %.1f years", calc$est_onset_age),
        x = "Age (years)",
        y = "Likelihood of remaining cognitively unimpaired from AD"
      ) +
      theme_cowplot() +
      theme(
        plot.title = element_text(size = 18, face = "bold", color = "#BA0C2F"),
        plot.subtitle = element_text(size = 14, color = "#6b7280"),
        axis.title = element_text(size = 12, face = "bold", color = "#374151"),
        axis.text = element_text(size = 11, color = "#6b7280"),
        panel.grid.major = element_line(color = "#f3f4f6", linewidth = 0.5),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
    plot
  })
  
  output$comparison_table <- renderTable({
    req(input$plasma, input$age_plasma, input$time_since_plasma, input$cohort, input$assay, input$method)
    req(model_data())
    
    current_data <- model_data()
    results <- data.frame(
      Cohort = character(),
      Assay = character(),
      Method = character(),
      `Estimated age at AD symptom onset` = character(),
      `Probability of symptomatic AD at plasma collection` = character(),
      `Probability of symptomatic AD at follow-up` = character(),
      stringsAsFactors = FALSE
    )
    
    tryCatch({
      clock_time <- current_data$smooth_clock_fn(input$plasma)
      if (is.na(clock_time)) clock_time <- 0
      est_onset_age <- input$age_plasma - clock_time
      current_age <- input$age_plasma + input$time_since_plasma
      
      linear_params <- load_linear_model_params(input$cohort, input$assay, input$method)
      if (!is.na(linear_params$slope) && !is.na(linear_params$intercept)) {
        est_symptom_onset_age <- est_onset_age * linear_params$slope + linear_params$intercept
      } else {
        est_symptom_onset_age <- NA
      }
      
      time_points <- seq(50, 100, by = 0.5)
      newdata <- data.frame(est_onset_age = est_onset_age)
      surv_probs <- sapply(time_points, function(t) {
        tryCatch({
          1 - getFitEsts(current_data$cox_model, newdata = newdata, q = t)
        }, error = function(e) {
          exp(-0.1 * (t - est_onset_age))
        })
      })
      surv_probs <- pmax(0, pmin(1, surv_probs))
      surv_data <- data.frame(time = time_points, surv = surv_probs)
      
      plasma_idx <- which.min(abs(surv_data$time - input$age_plasma))
      current_idx <- which.min(abs(surv_data$time - current_age))
      prob_at_collection <- 1 - surv_data$surv[plasma_idx]
      current_prob <- 1 - surv_data$surv[current_idx]
      
      results <- rbind(results, data.frame(
        Cohort = ifelse(input$cohort == "KADRC", "Knight ADRC", "ADNI"),
        Assay = input$assay,
        Method = input$method,
        `Estimated age at AD symptom onset` = sprintf("%.1f years", est_symptom_onset_age),
        `Probability of symptomatic AD at plasma collection` = sprintf("%.1f%%", prob_at_collection * 100),
        `Probability of symptomatic AD at follow-up` = sprintf("%.1f%%", current_prob * 100),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ))
    }, error = function(e) {
      return(results)
    })
    
    results
  }, striped = TRUE, bordered = TRUE, spacing = "s", align = "c", width = "100%")
}

shinyApp(ui = ui, server = server)
