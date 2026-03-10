#setup memory
rm(list = ls())
gc()
gc()

# load library-----
library(pacman)
p_load(tidyverse, here, lubridate, skimr, naniar, mice, tictoc,
       doParallel, foreach, benchmarkme)

get_ram()

# load data-----
dt <- readRDS("./data/ROX/preIM_cohort.rds")

glimpse(dt)

# Check missing pattern by variable
miss_summary <- miss_var_summary(dt)
print(miss_summary, n = 49)

# Visualize: missing rate by variable and time (time >= 0) -----
plot_vars <- c(
  "heart_rate", "resp_rate", "sbp", "dbp", "mbp", "spo2", "temperature", "fio2",
  "glucose", "gcs", "ph", "po2", "pco2", "sofa_24hours",
  "vasopressor", "crrt", "invasive", "noninvasive", "highflow",
  "sepsis3", "antibiotic_use", "lactate"
)
plot_vars <- intersect(plot_vars, names(dt))

# Variable grouping
var_groups <- tibble(variable = plot_vars) %>%
  mutate(group = case_when(
    variable %in% c("heart_rate","resp_rate","sbp","dbp","mbp","spo2","temperature","fio2") ~ "Vital Signs",
    variable %in% c("glucose","gcs","ph","po2","pco2","sofa_24hours") ~ "Lab / Score",
    TRUE ~ "Treatment / Intervention"
  ))

# Calculate missing rate by time point
miss_by_time <- dt %>%
  filter(time >= 0) %>%
  group_by(time) %>%
  summarise(across(all_of(plot_vars), ~mean(is.na(.)) * 100), n = n(), .groups = "drop") %>%
  pivot_longer(-c(time, n), names_to = "variable", values_to = "pct_missing") %>%
  left_join(var_groups, by = "variable")

# Color palettes (distinct colors per group)
cols_vital <- c(
  "heart_rate"   = "#E41A1C",
  "resp_rate"    = "#377EB8",
  "sbp"          = "#4DAF4A",
  "dbp"          = "#984EA3",
  "mbp"          = "#FF7F00",
  "spo2"         = "#A65628",
  "temperature"  = "#F781BF",
  "fio2"         = "#999999"
)

cols_lab <- c(
  "glucose"      = "#E41A1C",
  "gcs"          = "#377EB8",
  "ph"           = "#4DAF4A",
  "po2"          = "#984EA3",
  "pco2"         = "#FF7F00",
  "sofa_24hours" = "#A65628"
)

cols_tx <- c(
  "vasopressor"    = "#E41A1C",
  "crrt"           = "#377EB8",
  "invasive"       = "#4DAF4A",
  "noninvasive"    = "#984EA3",
  "highflow"       = "#FF7F00",
  "sepsis3"        = "#A65628",
  "antibiotic_use" = "#F781BF",
  "lactate"        = "#999999"
)

# Common theme
theme_miss <- theme_bw() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"))

# Plot 1: Vital Signs
print(
  ggplot(miss_by_time %>% filter(group == "Vital Signs"),
         aes(x = time, y = pct_missing, color = variable)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1) +
    scale_color_manual(values = cols_vital) +
    labs(title = "Missing Rate: Vital Signs",
         x = "Time (hours from HFNC start)", y = "Missing (%)", color = NULL) +
    theme_miss
)

# Plot 2: Lab / Score
print(
  ggplot(miss_by_time %>% filter(group == "Lab / Score"),
         aes(x = time, y = pct_missing, color = variable)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1) +
    scale_color_manual(values = cols_lab) +
    labs(title = "Missing Rate: Lab / Score",
         x = "Time (hours from HFNC start)", y = "Missing (%)", color = NULL) +
    theme_miss
)

# Plot 3: Treatment / Intervention
print(
  ggplot(miss_by_time %>% filter(group == "Treatment / Intervention"),
         aes(x = time, y = pct_missing, color = variable)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1) +
    scale_color_manual(values = cols_tx) +
    labs(title = "Missing Rate: Treatment / Intervention",
         x = "Time (hours from HFNC start)", y = "Missing (%)", color = NULL) +
    theme_miss
)

# Step 1: Split data (time <= 0: MICE, time > 0: LOCF) -----
cat("=== Splitting data into time <= 0 and time > 0 ===\n")

dt_pre  <- dt %>% filter(time <= 0)  # MICE target
dt_post <- dt %>% filter(time > 0)   # LOCF target

cat(sprintf("time <= 0: %s rows\n", format(nrow(dt_pre), big.mark = ",")))
cat(sprintf("time > 0 : %s rows\n", format(nrow(dt_post), big.mark = ",")))
cat(sprintf("Total    : %s rows (original: %s rows)\n",
            format(nrow(dt_pre) + nrow(dt_post), big.mark = ","),
            format(nrow(dt), big.mark = ",")))

# Step 2: MICE imputation (time <= 0 only, m=5) -----

# Variables to exclude from imputation
exclude_vars <- c(
  "subject_id",
  "elig_1",
  "elig_2",
  "death_datetime",      # Death datetime (no imputation needed)
  "death_location",      # Death location (no imputation needed)
  "intime",             # ICU admission time
  "outtime",            # ICU discharge time
  "admittime",          # Hospital admission time
  "dischtime"            # Hospital discharge time
)

dt_pre_for_impute <- dt_pre %>%
  dplyr::select(-all_of(exclude_vars))

# Check data size
cat(sprintf("\nMICE target data size: %s\n", format(object.size(dt_pre_for_impute), units = "MB")))

# Run mice (m=5: 5 imputed datasets)
cat("=== Starting mice-pmm imputation (m=5, time <= 0 only) ===\n")
cat(sprintf("Data: %s rows x %d columns\n",
            format(nrow(dt_pre_for_impute), big.mark = ","),
            ncol(dt_pre_for_impute)))
cat("Generating 5 imputed datasets\n\n")

set.seed(123)
start_time <- Sys.time()

tic()

imp <- mice(
  dt_pre_for_impute,
  m = 5,
  method = 'pmm',
  maxit = 50,
  printFlag = TRUE
)

total_time <- toc()
end_time <- Sys.time()

# Check MICE results
cat("\n=== MICE imputation results (time <= 0) ===\n")
cat(sprintf("Number of imputed datasets: %d\n", imp$m))

for(i in 1:imp$m) {
  dt_imputed_i <- complete(imp, i)
  cat(sprintf("Imputed dataset %d missing values: %d\n", i, sum(is.na(dt_imputed_i))))
}

# Step 3: Apply LOCF to each imputed dataset and save -----
cat("\n=== Restoring excluded variables + applying LOCF ===\n")

# LOCF target variables
locf_vars <- c(
  # Vital signs
  "heart_rate", "resp_rate", "sbp", "dbp", "mbp",
  "spo2", "temperature", "fio2", "glucose", "gcs",
  "ph", "po2", "pco2",
  # Treatment / status variables
  "vasopressor", "crrt", "sofa_24hours",
  "invasive", "noninvasive", "highflow",
  "sepsis3", "antibiotic_use",
  # Individual vasopressors
  "norepinephrine", "epinephrine", "dopamine",
  "phenylephrine", "dobutamine",
  # Code status
  "fullcode", "dnr", "dni", "dnr_dni", "cmo",
  # Demographics (constant within patient, but missing at time > 0)
  "race_grouped", "Married", "medicare_medicaid",
  "admission_location_grouped", "language_grouped",
  # MICE-imputed variables (may vary over time)
  "lactate", "baseline_sofa"
)

# Use only variables that exist in the data
locf_vars <- intersect(locf_vars, names(dt))
cat(sprintf("Number of LOCF target variables: %d\n", length(locf_vars)))

dt_final_list <- list()

for(i in 1:imp$m) {
  cat(sprintf("\n--- Imputed dataset %d ---\n", i))

  # Step 3a: Get imputed time <= 0 data + restore excluded variables
  dt_pre_imputed <- complete(imp, i) %>%
    mutate(
      subject_id = dt_pre$subject_id,
      elig_1 = dt_pre$elig_1,
      elig_2 = dt_pre$elig_2,
      death_datetime = dt_pre$death_datetime,
      death_location = dt_pre$death_location,
      intime = dt_pre$intime,
      outtime = dt_pre$outtime,
      admittime = dt_pre$admittime,
      dischtime = dt_pre$dischtime
    ) %>%
    dplyr::select(all_of(names(dt)))  # Restore original column order

  # Step 3b: Combine with time > 0
  dt_combined <- bind_rows(dt_pre_imputed, dt_post) %>%
    arrange(subject_id, time)

  cat(sprintf("After combining: %s rows x %d columns\n",
              format(nrow(dt_combined), big.mark = ","), ncol(dt_combined)))

  # Step 3c: Apply LOCF (forward fill within each subject_id)
  dt_final <- dt_combined %>%
    group_by(subject_id) %>%
    fill(all_of(locf_vars), .direction = "down") %>%
    ungroup()

  dt_final_list[[i]] <- dt_final

  cat(sprintf("Total missing values after LOCF: %d\n", sum(is.na(dt_final))))

  # Check missing at time=0
  miss_at_t0 <- dt_final %>%
    filter(time == 0) %>%
    summarise(across(everything(), ~sum(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "n_miss") %>%
    filter(n_miss > 0)

  if(nrow(miss_at_t0) > 0) {
    cat("  Variables with remaining missing at time=0:\n")
    print(miss_at_t0)
  } else {
    cat("  Missing at time=0: none (LOCF target variables)\n")
  }
}

# Final check -----
cat("\n=== Final data ===\n")
for(i in 1:imp$m) {
  cat(sprintf("\nImputed dataset %d:\n", i))
  cat(sprintf("  Size: %s rows x %d columns\n",
              format(nrow(dt_final_list[[i]]), big.mark = ","),
              ncol(dt_final_list[[i]])))
  cat(sprintf("  Total missing values: %d\n", sum(is.na(dt_final_list[[i]]))))

  # Variables with remaining missing
  final_miss_i <- miss_var_summary(dt_final_list[[i]]) %>%
    filter(n_miss > 0)

  if(nrow(final_miss_i) > 0) {
    cat("  Variables with remaining missing:\n")
    print(final_miss_i)
  }
}

# Step 4: Save data -----
cat("\n=== Saving data ===\n")

# 1. Save mids object (imputation results for time <= 0 only)
output_mids <- "./data/ROX/imputed_cohort_mice_m5_mids.rds"
saveRDS(imp, file = output_mids)
cat(sprintf("mids object saved: %s\n", output_mids))

# 2. Save all imputed datasets as a list
output_list <- "./data/ROX/imputed_cohort_mice_m5_list.rds"
saveRDS(dt_final_list, file = output_list)
cat(sprintf("Imputed dataset list saved: %s\n", output_list))

# 3. Save each imputed dataset as individual files
for(i in 1:imp$m) {
  output_i <- sprintf("./data/ROX/imputed_cohort_mice_m5_imp%d.rds", i)
  saveRDS(dt_final_list[[i]], file = output_i)
  cat(sprintf("Imputed dataset %d saved: %s\n", i, output_i))
}

cat("\n=== Processing complete ===\n")
cat(sprintf("MICE: imputed m=5 for time <= 0 (%s rows)\n",
            format(nrow(dt_pre), big.mark = ",")))
cat(sprintf("LOCF: forward-filled %d variables for time > 0 (%s rows)\n",
            length(locf_vars), format(nrow(dt_post), big.mark = ",")))
