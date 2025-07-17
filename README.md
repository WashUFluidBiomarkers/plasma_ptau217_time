# Overview
This repository contains R code for analyzing plasma p-tau217 biomarker trajectories to estimate the timing of Alzheimer's disease biomarker positivity and predicting symptomatic AD using the modeling approach TIRA (Temporal Integration of Rate Accumulation). TIRA estimates individual plasma %p-tau217 rates of change using linear mixed-effects modeling with random slopes and intercepts. The rates of change are used in generalized additive models (GAM) with cubic splines to characterize non-linear relationships between the rates of change and plasma %p-tau217 levels at the estimated midpoint of follow-up. The inverse of the modeled rate of change is integrated to derive the time between plasma %p-tau217 values.


# Interactive Tool
Explore the plasma p-tau217 biomarker clocks using our interactive Shiny application:
https://amyloid.shinyapps.io/plasma_ptau217_time/
