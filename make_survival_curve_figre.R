# Created by Kellen Petersen, July 1, 2025

library(tidyverse)
library(survival)
library(here)
library(icenReg)
library(doParallel)

setwd("<SET WORKING DIRECTORY>")

conflicts_prefer(dplyr::last)
conflicts_prefer(dplyr::lag)

all_or_healthy <- "HEALTHY"
which_dataset <-  "KADRC"
which_method <- "TIRA"
n_boot <- 5000
n_cores <- 8

if (which_dataset == "KADRC" & which_method == "TIRA") {
  df0 <- read.csv(here("results","<LOAD DATA>"))
} else if (which_dataset == "KADRC" & which_method == "SILA") {
  df0 <- read.csv(here("results","<LOAD DATA>"))
} else if (which_dataset == "ADNI" & which_method == "TIRA") {
  df0 <- read.csv(here("results","<LOAD DATA>"))
} else if (which_dataset == "ADNI" & which_method == "SILA") {
  df0 <- read.csv(here("results","<LOAD DATA>"))
}

if (which_dataset =="ADNI") {
  cdr0 <- read.csv(here("data_raw","<LOAD DATA>"))
  dx0 <- read.csv(here("data_raw","<LOAD DATA>"))
  demo0 <- read.csv(here("data_raw","<LOAD DATA>"))
  
  cdr <- cdr0 %>% 
    select(RID, VISCODE, VISCODE2, VISDATE, CDGLOBAL) %>% 
    rename(ID = RID, CDR = CDGLOBAL, TESTDATE_CDR = VISDATE) %>%
    mutate(VISCODE3 = if_else(VISCODE2 %in% c("f","sc","bl"), "bl0", VISCODE2)) %>% 
    filter(!is.na(CDR) & CDR != "") %>%
    select(-c(VISCODE, VISCODE2)) %>% 
    relocate(ID, VISCODE3, TESTDATE_CDR, CDR)
  
  dx <- dx0 %>% 
    select(RID, VISCODE, VISCODE2, EXAMDATE, DIAGNOSIS, DXMDUE, DXDDUE) %>% 
    rename(ID = RID) %>%
    mutate(VISCODE3 = if_else(VISCODE2 %in% c("f","sc","bl"), "bl0", VISCODE2)) %>% 
    filter(!is.na(DIAGNOSIS) & DIAGNOSIS != "") %>%
    select(-c(VISCODE, VISCODE2, EXAMDATE)) %>% 
    relocate(ID, VISCODE3, DIAGNOSIS, DXMDUE, DXDDUE)
  
  CDR_DX <- cdr %>% 
    full_join(dx, by = c("ID", "VISCODE3"), relationship = "many-to-many") %>% 
    select(-VISCODE3) %>% 
    relocate(ID, TESTDATE_CDR, CDR, DIAGNOSIS, DXMDUE, DXDDUE) %>%
    mutate(
      CDR_01 = ifelse(!is.na(CDR) & CDR > 0, 1, 0),
      CDGLOBAL = factor(CDR, levels = c(-1, 0, 0.5, 1.0, 2.0, 3.0)),
      DIAGNOSIS = factor(DIAGNOSIS, levels = c(1, 2, 3)),
      EXAMDATE = TESTDATE_CDR
    ) %>%
    group_by(ID) %>%
    arrange(EXAMDATE, .by_group = TRUE) %>%
    mutate(
      last_DIAGNOSIS = data.table::last(DIAGNOSIS, na.rm = TRUE),
      last_DXMDUE = data.table::last(DXMDUE, na.rm = TRUE),
      last_DXDDUE = data.table::last(DXDDUE, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      CDR_DX_Imp = ifelse(
        CDR_01 == 1 &
          (last_DIAGNOSIS %in% c(2, 3)) &
          (last_DXMDUE == 1 | last_DXDDUE == 1),
        1, 0
      ),
      CDR_DX_MCI = ifelse(
        CDR_01 == 1 &
          (last_DIAGNOSIS %in% c(2)) &
          (last_DXMDUE == 1),
        1, 0
      ),
      CDR_DX_AD = ifelse(
        CDR_01 == 1 &
          (last_DIAGNOSIS %in% c(3)) &
          (last_DXDDUE == 1),
        1, 0
      ),
      CDR_ONLY = ifelse(CDR>0, 1,0)
    )
  
  demo <- demo0 %>% 
    select(RID,PTDOB) %>% 
    rename(ID = RID) %>%
    arrange(ID)
  demo$PTDOB <- as.Date(demo$PTDOB, format = "%m/%d/%Y")
  demo <- demo %>% 
    na.omit() %>% 
    group_by(ID) %>%
    slice(1) %>%
    ungroup()
  
  CDR_DX$EXAMDATE <- as.Date(CDR_DX$EXAMDATE)
  CDR_DX_demo <- CDR_DX %>% 
    left_join(demo, by = "ID") %>% 
    mutate(AGE = as.numeric(difftime(EXAMDATE, PTDOB, units = "days"))/365.25) %>% 
    arrange(ID, EXAMDATE)
  
  dk <- CDR_DX_demo %>%
    select(ID, EXAMDATE, AGE, CDR_ONLY, CDR_DX_Imp) %>% 
    mutate(CDR_DX = CDR_DX_Imp,
           age = AGE) %>% 
    filter(!is.na(CDR_DX) & CDR_DX != "")
} else if (which_dataset =="KADRC") {
  dk0 <- read.csv(here("data_raw","<LOAD DATA>"), header = TRUE)
  
  dk0$TESTDATE <- as.Date(dk0$TESTDATE, format = "%d-%b-%y")
  
  dk <- dk0 %>% 
    filter(dx1 != "" & dx1 != ".") %>%  
    rename(EXAMDATE = TESTDATE)     
  
  dx1_AD <- c(
    "AD Dementia", "AD dem distrubed social, after", "AD dem distrubed social, with",
    "AD dem Language dysf after", "AD dem Language dysf prior", "AD dem Language dysf with",
    "AD dem w/CVD not contrib", "AD dem w/depresss, contribut", "AD dem w/depresss, not contribut",
    "AD dem w/Frontal lobe/demt at onset", "AD dem w/oth (list B) contribut",
    "AD dem w/oth (list B) not contrib", "AD dem w/oth unusual feat/subs demt",
    "AD dem w/PDI after AD dem contribut", "AD dem w/PDI after AD dem not contrib"
  )
  
  dk <- dk %>%
    mutate(
      dx_01 = ifelse(dx1 %in% dx1_AD, "AD", "Non-AD"),
      CDR_ONLY = ifelse(cdr > 0, 1, 0)
    ) %>%
    group_by(ID) %>%
    arrange(EXAMDATE) %>%
    mutate(dx_01_last_visit = last(dx_01)) %>%
    ungroup() %>%
    mutate(CDR_DX = ifelse(cdr > 0 & dx_01_last_visit == "AD", 1, 0)) %>% 
    mutate(ID = as.factor(ID)) %>%  
    arrange(ID, EXAMDATE)
  
  m0 <- read.csv(here("data_raw","<LOAD DATA>"), header = TRUE)
  mm <- m0 %>% 
    select(MAP_ID, BIRTH) %>%
    rename(ID = MAP_ID) %>%
    distinct(ID, .keep_all = TRUE) %>%  
    mutate(BIRTH = as.Date(BIRTH, format = "%d%b%Y")) %>%  
    mutate(ID = as.factor(ID)) %>%  
    arrange(ID)
  
  dk <- dk %>%
    left_join(mm, by = "ID") %>%
    mutate(
      age = difftime(EXAMDATE, BIRTH, units = "days"),  
      age = round(as.numeric(age/365.25), 2)
    ) %>% 
    filter(!is.na(BIRTH)) %>% 
    select(ID, EXAMDATE, BIRTH, age, dx1, dx_01, cdr, CDR_ONLY, CDR_DX)
}

if (all_or_healthy == "ALL") {
  surv_data <- dk %>%
    group_by(ID) %>%
    arrange(EXAMDATE) %>%
    mutate(
      event = CDR_DX,
      prev_age = lag(age, default = NA)
    ) %>%
    summarise(
      first_event = if(any(event == 1, na.rm = TRUE)) which.max(event == 1) else NA_integer_,
      first_age = first(age),
      last_age = last(age),
      event_age = if(any(event == 1, na.rm = TRUE)) age[first_event] else NA_real_,
      prev_event_age = if(any(event == 1, na.rm = TRUE)) prev_age[first_event] else NA_real_,
      left = case_when(
        any(event == 1, na.rm = TRUE) && first_event == 1 ~ 0,
        any(event == 1, na.rm = TRUE) && first_event > 1 ~ prev_event_age,
        TRUE ~ last_age
      ),
      right = case_when(
        any(event == 1, na.rm = TRUE) && first_event == 1 ~ first_age,
        any(event == 1, na.rm = TRUE) && first_event > 1 ~ event_age,
        TRUE ~ Inf
      ),
      status = case_when(
        any(event == 1, na.rm = TRUE) && first_event == 1 ~ 2L,
        any(event == 1, na.rm = TRUE) && first_event > 1 ~ 1L,
        TRUE ~ 0L
      )
    ) %>%
    select(ID, left, right, status) %>%
    ungroup()
} else if (all_or_healthy == "HEALTHY") {
  healthy_ids <- dk %>%
    group_by(ID) %>%
    arrange(EXAMDATE) %>%
    filter(row_number() == 1) %>%
    filter(CDR_DX == 0) %>%
    pull(ID) %>%
    unique()
  
  surv_data <- dk %>%
    filter(ID %in% healthy_ids) %>%
    group_by(ID) %>%
    arrange(EXAMDATE) %>%
    mutate(
      event = CDR_DX,
      prev_age = lag(age, default = NA)
    ) %>%
    summarise(
      first_event = if(any(event == 1, na.rm = TRUE)) which.max(event == 1) else NA_integer_,
      first_age = first(age),
      last_age = last(age),
      event_age = if(any(event == 1, na.rm = TRUE)) age[first_event] else NA_real_,
      prev_event_age = if(any(event == 1, na.rm = TRUE)) prev_age[first_event] else NA_real_,
      left = case_when(
        any(event == 1, na.rm = TRUE) && first_event > 1 ~ prev_event_age,
        TRUE ~ last_age
      ),
      right = case_when(
        any(event == 1, na.rm = TRUE) && first_event > 1 ~ event_age,
        TRUE ~ Inf
      ),
      status = case_when(
        any(event == 1, na.rm = TRUE) && first_event > 1 ~ 1L,
        TRUE ~ 0L
      )
    ) %>%
    select(ID, left, right, status) %>%
    ungroup()
}
surv_data$ID <- as.factor(surv_data$ID)

df <- df0 %>%
  select(ID, est_onset_age) %>%
  group_by(ID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(ID = as.factor(ID))
df$ID <- as.factor(df$ID)

surv_data <- surv_data %>%
  left_join(df, by = "ID") %>%
  na.omit()

cl <- makeCluster(n_cores)
registerDoParallel(cl)
fit <- ic_sp(Surv(left, right, type = "interval2") ~ est_onset_age,
             data = surv_data,
             bs_samples = n_boot,
             useMCores = TRUE)
stopCluster(cl)

summary(fit)

quantiles <- c(60, 70, 80, 90)

median_surv_direct <- sapply(quantiles, function(q) {
  newdata <- data.frame(est_onset_age = q)
  getFitEsts(fit, newdata = newdata)
})

time_points <- seq(50, 100, by = 0.5)
surv_list <- lapply(1:length(quantiles), function(i) {
  newdata <- data.frame(est_onset_age = quantiles[i])
  probs <- getFitEsts(fit, newdata = newdata, q = time_points)
  
  data.frame(
    time = time_points,
    surv = 1 - probs,
    Onset_Age = factor(paste("Positivity age of", quantiles[i]),
                       levels = paste("Positivity age of", quantiles))
  )
})

surv_df <- bind_rows(surv_list)

median_df <- data.frame(
  Onset_Age = paste("Positivity Age", quantiles),
  Median_Survival = median_surv_direct,
  Quantile = quantiles
) %>%
  mutate(Median_Survival_diff = Median_Survival - Quantile)

print(median_df)

segment_data <- data.frame(
  x_start = median_df$Quantile,
  x_end = median_df$Median_Survival,
  y_pos = c(-0.02, -0.04, -0.06, -0.08),
  colors = c("red", "blue", "darkgreen", "purple"),
  onset_age = median_df$Quantile
)

vline_data <- data.frame(
  x_pos = quantiles,
  colors = c("red", "blue", "darkgreen", "purple")
)

p_survival <- ggplot(surv_df, aes(x = time, y = surv, color = Onset_Age)) +
  geom_step(linewidth = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.7) +
  geom_vline(data = vline_data, 
             aes(xintercept = x_pos), 
             color = vline_data$colors,
             linetype = "dashed", linewidth = 0.5, alpha = 0.7, 
             inherit.aes = FALSE) +
  geom_segment(data = segment_data,
               aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos),
               color = segment_data$colors,
               linewidth = 3,
               inherit.aes = FALSE) +
  scale_color_manual(values = c("red", "blue", "darkgreen", "purple")) +
  scale_x_continuous(limits = c(50, 100), breaks = seq(50, 100, 10)) +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Age (years)",
    y = "Probability of cognitively unimpaired",
    color = "Estimated age of \n%p-tau217 positivity"
  ) +
  cowplot::theme_cowplot() +
  theme(
    legend.position = c(0.025, 0.15),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "gray90", size = 0.5),
    legend.margin = margin(6, 6, 6, 6),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p_survival)
ggsave(
  filename = here("results", paste0("survival_plot_", which_dataset, "_", which_method, ".png")),
  plot = p_survival,
  width = 8,
  height = 6,
  dpi = 500,
  bg = "white"
)