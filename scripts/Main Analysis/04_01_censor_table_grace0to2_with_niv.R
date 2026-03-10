#setup memory
rm(list = ls())
gc()
gc()

# load library -----
library(pacman)

p_load(tidyverse, here, data.table, boot, lubridate, skimr, naniar,gtsummary,officer, flextable,
       patchwork, ggsci, parallel,speedglm, doParallel, foreach)



# load data
#dt <- readRDS("./data/ROX/imputed_data_for_analysis.rds")
dt <-readRDS("./data/ROX/imputed_data_for_analysis_m5_imp1.rds")
dt<-as.data.table(dt)
class(dt)

# filter eligible patients
dt<-dt |> filter(elig_1==1)

# Disable printing results in scientific notation
options(scipen=999)

# Set time point for estimation of risks (time 0:720 hours; 30 days)
K <- 721

# Cloning ----
# Clone each eligible individual twice
arm0 <- copy(dt)[,arm:=0]
arm1 <- copy(dt)[,arm:=1]

# Check the number of clones in each arm
length(arm0$time[which(arm0$time==0)])
length(arm1$time[which(arm1$time==0)])

# Check the number of events in each arm prior to censoring
table(arm0$death)
table(arm1$death)


# Censoring ----

### Censor clones in arm 0 (control: Ususal care) when non-adherent
### In this arm no non-adherent
length(arm0)

# Match structure with arm1
arm0.c <- copy(arm0)[, my.censor := 0]

### Censor clones in arm 1 when non-adherent
## When ROX < 3.85, treat 1
### Non adherent
#### - if and when they were intubated rox > 3.85
#### - if and when they were not intubated rox < 3.85
### Censoring function for "intubate when ROX < 3.85" strategy


# Censor function -----
### Arm 1: ROX-guided strategy WITH GRACE PERIOD
### Grace period: Allow specified hours after ROX first drops below cutoff

arm1.censor.grace <- function(rox, treat_cum, treat, highflow_lag, cutoff = 3.85, grace_hours = 1){
  n <- length(rox)
  my.censor <- rep(0, n)
  first.cens <- NA

  # Create lagged treat_cum (to check if pre-intubation)
  treat_cum_lag <- c(0, treat_cum[1:(n-1)])

  # Record the first time point when ROX < cutoff
  first_below_cutoff <- NA

  for(i in 1:n){
    # Check only pre-intubation (treat_cum == 0 at previous time) and only while on HFNC (aligning with lagged ROX)
    if(!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
       !is.na(highflow_lag[i]) && highflow_lag[i] == 1){

      # === Non-adherence 1: ROX < cutoff response (WITH GRACE PERIOD) ===
      if(!is.na(rox[i]) && rox[i] < cutoff){

        # Record the first time point when ROX < cutoff
        if(is.na(first_below_cutoff)){
          first_below_cutoff <- i
        }
      }

      # Censor if not intubated after grace period expires since ROX dropped below cutoff
      # (counter is not reset even if ROX recovers)
      if(!is.na(first_below_cutoff) && i >= first_below_cutoff + grace_hours &&
         !is.na(treat[i]) && treat[i] == 0){
        first.cens <- i
        break
      }

      # === Non-adherence 2: Intubated despite ROX >= cutoff ===
      if(!is.na(rox[i]) && !is.na(treat[i]) &&
         rox[i] >= cutoff && treat[i] == 1){
        first.cens <- i
        break
      }
    }
  }

  if(!is.na(first.cens)){
    my.censor[first.cens:n] <- 1
  }

  return(my.censor)
}

## Apply censoring with grace period --------
# Default is 1 hour; sensitivity analysis possible by changing grace_hours
arm1.c <- copy(arm1)[
  , my.censor := arm1.censor.grace(rox, treat_cum, treat, highflow_lag, cutoff = 3.85, grace_hours = 1), by = id][
    , good := cumsum(my.censor), by = id][
      good <= 1][
        , good := NULL]

# Censor example
arm1.c |>filter(id==5) |>
  dplyr::select(id,time,rox,treat,treat_cum,death,my.censor) |> head(50)


# Check example of first censored patient
censored_ids <- arm1[, .(censored = max(arm1.censor.grace(rox, treat_cum, treat, highflow_lag, 3.85))), by = id][censored == 1, id]
if(length(censored_ids) > 0){
  cat("Example of censored patient (first case):\n")
  print(arm1[id == censored_ids[1], .(id, time, rox, treat, treat_cum, death)])
}

# Check censoring status
table(arm1.c$my.censor)

# Censoring table any cutoff  ----------

# Set cutoff values
cutoff_values <- seq(2, 13, by = 1.0)
grace_hours_values <- c(0, 1, 2)

# Store all results
all_censoring_summary <- data.frame()
cat("Calculating...\n")

for(grace_hours in grace_hours_values) {

  cat(sprintf("\n=== Grace Period: %d hour(s) ===\n", grace_hours))

  for(cutoff in cutoff_values) {

    # Apply censoring for each patient
    arm1_temp <- copy(arm1)

    arm1_temp[, my.censor := arm1.censor.grace(rox, treat_cum, treat, highflow_lag,
                                               cutoff = cutoff,
                                               grace_hours = grace_hours),
              by = id]

    # Get information at censoring time point (for reason classification)
    censor_info <- arm1_temp[my.censor == 1, .SD[1], by = id]
    censor_info[, censor_reason := ifelse(rox < cutoff & treat == 0,
                                          "Not Intubated",
                                          "Intubated")]

    # Aggregate censoring and intubation status per patient
    patient_censor <- arm1_temp[, .(
      censored = any(my.censor == 1),
      intubated = any(treat == 1),
      rox_at_intubation = {
        idx <- which(treat == 1 & highflow_lag == 1)[1]
        if(!is.na(idx) && length(idx) > 0) {
          as.numeric(rox[idx])
        } else {
          NA_real_
        }
      },
      time_to_intubation = {
        idx <- which(treat == 1 & highflow_lag == 1)[1]
        if(!is.na(idx) && length(idx) > 0) {
          as.numeric(time[idx])
        } else {
          NA_real_
        }
      }
    ), by = id]

    # Merge censoring reasons
    patient_censor <- merge(patient_censor,
                            censor_info[, .(id, censor_reason)],
                            by = "id", all.x = TRUE)

    # Overall statistics
    total_patients <- uniqueN(arm1$id)
    censored_patients <- sum(patient_censor$censored)
    censored_pct <- (censored_patients / total_patients) * 100

    # Aggregation by censoring reason
    n_censor_not_intubated <- sum(patient_censor$censor_reason == "Not Intubated", na.rm = TRUE)
    n_censor_intubated <- sum(patient_censor$censor_reason == "Intubated", na.rm = TRUE)

    # Breakdown of protocol-adherent patients
    following_protocol <- patient_censor[censored == FALSE]
    n_following <- nrow(following_protocol)
    n_intubated <- sum(following_protocol$intubated)
    n_not_intubated <- n_following - n_intubated

    # ROX and time for protocol-adherent intubated patients
    following_intubated <- following_protocol[intubated == TRUE]
    rox_median <- median(following_intubated$rox_at_intubation, na.rm = TRUE)
    rox_q25 <- quantile(following_intubated$rox_at_intubation, 0.25, na.rm = TRUE)
    rox_q75 <- quantile(following_intubated$rox_at_intubation, 0.75, na.rm = TRUE)

    time_median <- median(following_intubated$time_to_intubation, na.rm = TRUE)
    time_q25 <- quantile(following_intubated$time_to_intubation, 0.25, na.rm = TRUE)
    time_q75 <- quantile(following_intubated$time_to_intubation, 0.75, na.rm = TRUE)

    # Append results
    all_censoring_summary <- rbind(all_censoring_summary, data.frame(
      Grace_Hours = grace_hours,
      Cutoff = cutoff,
      N_Censored = censored_patients,
      Percentage = censored_pct,
      N_Censor_Not_Intubated = n_censor_not_intubated,
      N_Censor_Intubated = n_censor_intubated,
      N_Following = n_following,
      N_Intubated = n_intubated,
      N_Not_Intubated = n_not_intubated,
      ROX_Median = rox_median,
      ROX_Q25 = rox_q25,
      ROX_Q75 = rox_q75,
      Time_Median = time_median,
      Time_Q25 = time_q25,
      Time_Q75 = time_q75
    ))

    cat(sprintf("Grace %dhr, Cutoff %.1f: Censored %d (Not intubated: %d, Intubated: %d), Adherent: %d (Intubated: %d, Not intubated: %d), ROX: %.2f, Time: %.1fh\n",
                grace_hours, cutoff, censored_patients, n_censor_not_intubated, n_censor_intubated,
                n_following, n_intubated, n_not_intubated, rox_median, time_median))
  }
}

# Create table
total_patients <- uniqueN(arm1$id)

censoring_table <- all_censoring_summary |>
  mutate(
    `Grace Period (hours)` = sprintf("%d", Grace_Hours),
    `ROX Cutoff` = sprintf("%.1f", Cutoff),
    `Total N` = total_patients,
    `N Censored` = N_Censored,
    `% Censored` = sprintf("%.1f%%", Percentage),
    `Censored (Not Intubated)` = N_Censor_Not_Intubated,
    `Censored (Intubated)` = N_Censor_Intubated,
    `N Following Protocol` = N_Following,
    `% Following Protocol` = sprintf("%.1f%%", (N_Following / total_patients) * 100),
    `N Intubated` = N_Intubated,
    `N Not Intubated` = N_Not_Intubated,
    `ROX at Intubation` = sprintf("%.2f (%.2f-%.2f)", ROX_Median, ROX_Q25, ROX_Q75),
    `Time to Intubation, hours` = sprintf("%.1f (%.1f-%.1f)", Time_Median, Time_Q25, Time_Q75)
  ) |>
  dplyr::select(`Grace Period (hours)`, `ROX Cutoff`, `Total N`,
                `N Censored`, `% Censored`,
                `Censored (Not Intubated)`, `Censored (Intubated)`,
                `N Following Protocol`, `% Following Protocol`,
                `N Intubated`, `N Not Intubated`,
                `ROX at Intubation`, `Time to Intubation, hours`)

# Create flextable
ft <- flextable(censoring_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 2), align = "center", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  add_header_row(
    values = c("", "", "", "Censored (Non-adherence)", "Following Protocol (Adherence)", "Among Adherent Intubated"),
    colwidths = c(1, 1, 1, 4, 4, 2)
  ) %>%
  align(part = "header", i = 1, align = "center") %>%
  bold(part = "header") %>%
  add_footer_lines(
    values = c(
      sprintf("Total patients in Arm 1 (ROX-guided strategy): %d", total_patients),
      "Grace period: Time allowed after ROX first drops below cutoff before censoring",
      "Results shown for grace periods of 0, 1, and 2 hours",
      "Censoring reasons:",
      "  - Not Intubated: ROX < cutoff but patient was not intubated",
      "  - Intubated: ROX ≥ cutoff but patient was intubated",
      "Following Protocol:",
      "  - N Intubated: ROX < cutoff and patient was appropriately intubated",
      "  - N Not Intubated: ROX remained ≥ cutoff and patient was not intubated",
      "ROX at Intubation and Time to Intubation: Median (IQR) among adherent intubated patients"
    )
  ) %>%
  fontsize(part = "footer", size = 9) %>%
  italic(part = "footer")

# Display
print(ft)

ft1<-ft

# Save as Word document
doc <- read_docx() %>%
  body_add_par("Supplementary Table. Censoring Status by ROX Cutoff Value and Grace Period",
               style = "heading 1") %>%
  body_add_flextable(ft1)

print(doc, target = "outputs/ROX/S.Table_censoring_by_cutoff_grace0to2hr.docx")

cat("\n✓ Table saved: outputs/ROX/S.Table_censoring_by_cutoff_grace0to2hr.docx\n")



# Ususal care, 3.85, 4.88, time depended table summary (WITH NIV CENSOR INFO) -----
library(data.table)
library(dplyr)
library(flextable)
library(officer)

# Set cutoff values (0 = Usual Care, "time_dependent" = Time-varying ROX)
cutoff_values <- c(0, 3.85, 4.88, "time_dependent")
grace_hours_values <- c(0, 1, 2)

# Censoring function for fixed threshold
arm1.censor.grace.fixed <- function(rox, treat_cum, treat, highflow_lag, cutoff, grace_hours = 0) {
  n <- length(rox)
  my.censor <- rep(0, n)
  first.cens <- NA

  treat_cum_lag <- c(0, treat_cum[1:(n-1)])
  first_below_cutoff <- NA

  for(i in 1:n) {
    if(!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
       !is.na(highflow_lag[i]) && highflow_lag[i] == 1) {

      if(!is.na(rox[i]) && rox[i] < cutoff) {
        if(is.na(first_below_cutoff)) {
          first_below_cutoff <- i
        }
      }

      if(!is.na(first_below_cutoff) && i >= first_below_cutoff + grace_hours &&
         !is.na(treat[i]) && treat[i] == 0) {
        first.cens <- i
        break
      }

      if(!is.na(rox[i]) && !is.na(treat[i]) &&
         rox[i] >= cutoff && treat[i] == 1) {
        first.cens <- i
        break
      }
    }
  }

  if(!is.na(first.cens)) {
    my.censor[first.cens:n] <- 1
  }

  return(my.censor)
}

# Censoring function for time-varying ROX
arm1.censor.grace.timedep <- function(rox, treat_cum, treat, highflow_lag, cal_time, grace_hours = 0) {
  n <- length(rox)
  my.censor <- rep(0, n)
  first.cens <- NA

  treat_cum_lag <- c(0, treat_cum[1:(n-1)])

  get_cutoff <- function(time_hr) {
    if (time_hr >= 1 && time_hr <= 5) {
      return(2.85)
    } else if (time_hr >= 6 && time_hr <= 11) {
      return(3.47)
    } else {
      return(3.85)
    }
  }

  first_below_cutoff <- NA

  for (i in 1:n) {
    if (!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
        !is.na(highflow_lag[i]) && highflow_lag[i] == 1) {

      current_time <- cal_time[i]
      cutoff <- get_cutoff(current_time)
      current_rox <- rox[i]
      current_treat <- treat[i]

      if (is.na(current_rox) || is.na(current_treat)) {
        next
      }

      if (current_rox < cutoff) {
        if (grace_hours > 0) {
          if (is.na(first_below_cutoff)) {
            first_below_cutoff <- i
          }
          if (i >= first_below_cutoff + grace_hours && current_treat == 0) {
            first.cens <- i
            break
          }
        } else {
          if (current_treat == 0) {
            first.cens <- i
            break
          }
        }
      } else {
        if (current_treat == 1) {
          if (grace_hours > 0 && !is.na(first_below_cutoff) &&
              i < first_below_cutoff + grace_hours) {
            # Adherent within grace period
          } else if (grace_hours > 0 && !is.na(first_below_cutoff)) {
            first.cens <- i
            break
          } else {
            first.cens <- i
            break
          }
        }
      }
    }
  }

  if (!is.na(first.cens)) {
    my.censor[first.cens:n] <- 1
  }

  return(my.censor)
}

# List to store results
all_summary_results <- list()

for(grace_hours in grace_hours_values) {

  cat(sprintf("\n\n========== Grace Period: %d hour(s) ==========\n", grace_hours))

  for(cutoff in cutoff_values) {

    if(cutoff == 0) {
      cat("\n=== Calculating Usual Care (no cutoff) ===\n")

      arm_temp <- copy(arm0)
      arm_temp$my.censor <- 0

      strategy_name <- "Usual Care"

      n_censor_not_intubated <- NA
      n_censor_intubated <- NA

    } else if(cutoff == "time_dependent") {
      cat("\n=== Calculating Time-varying ROX ===\n")

      arm_temp <- copy(arm1)
      arm_temp[, my.censor := arm1.censor.grace.timedep(rox, treat_cum, treat, highflow_lag,
                                                        cal_time,
                                                        grace_hours = grace_hours), by = id]

      # Get information at censoring time point (for reason classification)
      # For time-varying, use the threshold at that time point
      censor_info <- arm_temp[my.censor == 1, .SD[1], by = id]
      censor_info[, cutoff_at_censor := fifelse(cal_time >= 1 & cal_time <= 5, 2.85,
                                                fifelse(cal_time >= 6 & cal_time <= 11, 3.47, 3.85))]
      censor_info[, censor_reason := ifelse(rox < cutoff_at_censor & treat == 0,
                                            "Not Intubated",
                                            "Intubated")]

      n_censor_not_intubated <- sum(censor_info$censor_reason == "Not Intubated", na.rm = TRUE)
      n_censor_intubated <- sum(censor_info$censor_reason == "Intubated", na.rm = TRUE)

      arm_temp[, good := cumsum(my.censor), by = id]
      arm_temp <- arm_temp[good <= 1]
      arm_temp[, good := NULL]

      strategy_name <- "Time-varying ROX"

    } else {
      cutoff_num <- as.numeric(cutoff)
      cat(sprintf("\n=== Calculating ROX-guided (cutoff %.2f) ===\n", cutoff_num))

      arm_temp <- copy(arm1)
      arm_temp[, my.censor := arm1.censor.grace.fixed(rox, treat_cum, treat, highflow_lag,
                                                      cutoff = cutoff_num,
                                                      grace_hours = grace_hours), by = id]

      # Get information at censoring time point
      censor_info <- arm_temp[my.censor == 1, .SD[1], by = id]
      censor_info[, censor_reason := ifelse(rox < cutoff_num & treat == 0,
                                            "Not Intubated",
                                            "Intubated")]

      n_censor_not_intubated <- sum(censor_info$censor_reason == "Not Intubated", na.rm = TRUE)
      n_censor_intubated <- sum(censor_info$censor_reason == "Intubated", na.rm = TRUE)

      arm_temp[, good := cumsum(my.censor), by = id]
      arm_temp <- arm_temp[good <= 1]
      arm_temp[, good := NULL]

      strategy_name <- sprintf("ROX < %.2f", cutoff_num)
    }

  # Patient-level aggregation (with NIV censor info)
  patient_stats <- arm_temp[, .(
    intubated = any(treat == 1),
    rox_at_intubation = {
      idx <- which(treat == 1 & highflow_lag == 1)[1]
      if(!is.na(idx) && length(idx) > 0) {
        as.numeric(rox[idx])
      } else {
        NA_real_
      }
    },
    spo2_at_intubation = {
      idx <- which(treat == 1 & highflow_lag == 1)[1]
      if(!is.na(idx) && length(idx) > 0) {
        as.numeric(spo2[idx])
      } else {
        NA_real_
      }
    },
    fio2_at_intubation = {
      idx <- which(treat == 1 & highflow_lag == 1)[1]
      if(!is.na(idx) && length(idx) > 0) {
        as.numeric(fio2[idx])
      } else {
        NA_real_
      }
    },
    resp_rate_at_intubation = {
      idx <- which(treat == 1 & highflow_lag == 1)[1]
      if(!is.na(idx) && length(idx) > 0) {
        as.numeric(resp_rate[idx])
      } else {
        NA_real_
      }
    },
    time_to_intubation = {
      idx <- which(treat == 1 & highflow_lag == 1)[1]
      if(!is.na(idx) && length(idx) > 0) {
        as.numeric(time[idx])
      } else {
        NA_real_
      }
    },
    censored = any(my.censor == 1),
    # Add NIV censor
    censored_niv = any(censor == 1, na.rm = TRUE),
    died = any(death == 1, na.rm = TRUE)
  ), by = id]

  # Calculate statistics
  total_n <- nrow(patient_stats)
  n_censored <- sum(patient_stats$censored)
  pct_censored <- (n_censored / total_n) * 100

  # NIV censor statistics
  n_censored_niv <- sum(patient_stats$censored_niv)
  pct_censored_niv <- (n_censored_niv / total_n) * 100

  patient_stats_not_censored <- patient_stats[censored == FALSE]
  n_not_censored <- nrow(patient_stats_not_censored)

  n_following_intubated <- sum(patient_stats_not_censored$intubated)
  n_following_not_intubated <- n_not_censored - n_following_intubated

  n_intubated <- sum(patient_stats_not_censored$intubated)
  pct_intubated <- if(n_not_censored > 0) (n_intubated / n_not_censored) * 100 else NA

  rox_intubated <- patient_stats_not_censored$rox_at_intubation[patient_stats_not_censored$intubated]
  rox_median <- median(rox_intubated, na.rm = TRUE)
  rox_q25 <- quantile(rox_intubated, 0.25, na.rm = TRUE)
  rox_q75 <- quantile(rox_intubated, 0.75, na.rm = TRUE)

  time_intubated <- patient_stats_not_censored$time_to_intubation[patient_stats_not_censored$intubated]
  time_median <- median(time_intubated, na.rm = TRUE)
  time_q25 <- quantile(time_intubated, 0.25, na.rm = TRUE)
  time_q75 <- quantile(time_intubated, 0.75, na.rm = TRUE)

  # SpO2 at intubation
  spo2_intubated <- patient_stats_not_censored$spo2_at_intubation[patient_stats_not_censored$intubated]
  spo2_median <- median(spo2_intubated, na.rm = TRUE)
  spo2_q25 <- quantile(spo2_intubated, 0.25, na.rm = TRUE)
  spo2_q75 <- quantile(spo2_intubated, 0.75, na.rm = TRUE)

  # FiO2 at intubation
  fio2_intubated <- patient_stats_not_censored$fio2_at_intubation[patient_stats_not_censored$intubated]
  fio2_median <- median(fio2_intubated, na.rm = TRUE)
  fio2_q25 <- quantile(fio2_intubated, 0.25, na.rm = TRUE)
  fio2_q75 <- quantile(fio2_intubated, 0.75, na.rm = TRUE)

  # Respiratory rate at intubation
  resp_rate_intubated <- patient_stats_not_censored$resp_rate_at_intubation[patient_stats_not_censored$intubated]
  resp_rate_median <- median(resp_rate_intubated, na.rm = TRUE)
  resp_rate_q25 <- quantile(resp_rate_intubated, 0.25, na.rm = TRUE)
  resp_rate_q75 <- quantile(resp_rate_intubated, 0.75, na.rm = TRUE)

  n_died <- sum(patient_stats$died)
  pct_died <- (n_died / total_n) * 100

    # Save results (with NIV censor info)
    list_key <- paste0("grace", grace_hours, "_", cutoff)
    all_summary_results[[list_key]] <- data.frame(
      Grace_Hours = grace_hours,
      Strategy = strategy_name,
      ROX_Cutoff = as.character(cutoff),
      Total_N = total_n,
      N_Censored = n_censored,
      Pct_Censored = pct_censored,
      N_Censor_Not_Intubated = n_censor_not_intubated,
      N_Censor_Intubated = n_censor_intubated,
      N_Censored_NIV = n_censored_niv,
      Pct_Censored_NIV = pct_censored_niv,
      N_Not_Censored = n_not_censored,
      N_Following_Intubated = n_following_intubated,
      N_Following_Not_Intubated = n_following_not_intubated,
      N_Intubated = n_intubated,
      Pct_Intubated = pct_intubated,
      ROX_at_Intub_Median = rox_median,
      ROX_at_Intub_Q25 = rox_q25,
      ROX_at_Intub_Q75 = rox_q75,
      SpO2_at_Intub_Median = spo2_median,
      SpO2_at_Intub_Q25 = spo2_q25,
      SpO2_at_Intub_Q75 = spo2_q75,
      FiO2_at_Intub_Median = fio2_median,
      FiO2_at_Intub_Q25 = fio2_q25,
      FiO2_at_Intub_Q75 = fio2_q75,
      RR_at_Intub_Median = resp_rate_median,
      RR_at_Intub_Q25 = resp_rate_q25,
      RR_at_Intub_Q75 = resp_rate_q75,
      Time_to_Intub_Median = time_median,
      Time_to_Intub_Q25 = time_q25,
      Time_to_Intub_Q75 = time_q75,
      N_Died = n_died,
      Pct_Died = pct_died
    )

    # Progress display
    cat(sprintf("Total patients: %d\n", total_n))
    cat(sprintf("Censored for protocol violation: %d (%.1f%%)\n", n_censored, pct_censored))
    if(cutoff != 0) {
      cat(sprintf("  - Censored (not intubated): %d\n", n_censor_not_intubated))
      cat(sprintf("  - Censored (intubated): %d\n", n_censor_intubated))
    }
    cat(sprintf("Censored for NIV transition: %d (%.1f%%)\n", n_censored_niv, pct_censored_niv))
    cat(sprintf("Adherent patients: %d\n", n_not_censored))
    cat(sprintf("  - Adherent and intubated: %d\n", n_following_intubated))
    cat(sprintf("  - Adherent and not intubated: %d\n", n_following_not_intubated))
    cat(sprintf("ROX at intubation (median): %.2f (IQR: %.2f-%.2f)\n",
                rox_median, rox_q25, rox_q75))
    cat(sprintf("Time to intubation (median): %.1f hours (IQR: %.1f-%.1f)\n",
                time_median, time_q25, time_q75))
    cat(sprintf("Deaths: %d (%.1f%%)\n", n_died, pct_died))
  }
}

# Combine into data frame
summary_df <- do.call(rbind, all_summary_results)

# Create table (with NIV censor column)
comparison_table <- summary_df |>
  mutate(
    `Grace Period (hours)` = sprintf("%d", Grace_Hours),
    `Strategy` = Strategy,
    `Total N` = Total_N,
    `Censored (Protocol violation), n (%)` = sprintf("%d (%.1f%%)", N_Censored, Pct_Censored),
    `Censored (Not Intubated)` = ifelse(is.na(N_Censor_Not_Intubated), "-", as.character(N_Censor_Not_Intubated)),
    `Censored (Intubated)` = ifelse(is.na(N_Censor_Intubated), "-", as.character(N_Censor_Intubated)),
    `Censored (NIV), n (%)` = sprintf("%d (%.1f%%)", N_Censored_NIV, Pct_Censored_NIV),
    `Following Protocol, n (%)` = sprintf("%d (%.1f%%)", N_Not_Censored, (N_Not_Censored / Total_N) * 100),
    `Following (Intubated), n (%)` = sprintf("%d (%.1f%%)", N_Following_Intubated, (N_Following_Intubated / N_Not_Censored) * 100),
    `Following (Not Intubated), n (%)` = sprintf("%d (%.1f%%)", N_Following_Not_Intubated, (N_Following_Not_Intubated / N_Not_Censored) * 100),
    `ROX at Intubation` = sprintf("%.2f (%.2f-%.2f)",
                                  ROX_at_Intub_Median,
                                  ROX_at_Intub_Q25,
                                  ROX_at_Intub_Q75),
    `Time to Intubation, hours` = sprintf("%.1f (%.1f-%.1f)",
                                          Time_to_Intub_Median,
                                          Time_to_Intub_Q25,
                                          Time_to_Intub_Q75),
    `Deaths, n (%)` = sprintf("%d (%.1f%%)", N_Died, Pct_Died)
  ) |>
  dplyr::select(`Grace Period (hours)`, `Strategy`, `Total N`,
                `Censored (Protocol violation), n (%)`, `Censored (Not Intubated)`, `Censored (Intubated)`,
                `Censored (NIV), n (%)`,
                `Following Protocol, n (%)`, `Following (Intubated), n (%)`, `Following (Not Intubated), n (%)`,
                `ROX at Intubation`, `Time to Intubation, hours`, `Deaths, n (%)`)

# Create flextable
ft <- flextable(comparison_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 2), align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  add_header_row(
    values = c("", "", "", "Censored (Non-adherence)", "", "Following Protocol (Adherence)", "", "", ""),
    colwidths = c(1, 1, 1, 3, 1, 3, 1, 1, 1)
  ) %>%
  align(part = "header", i = 1, align = "center") %>%
  bold(part = "header") %>%
  add_footer_lines(
    values = c(
      "Grace period: Time allowed after ROX first drops below cutoff before censoring",
      "Results shown for grace periods of 0, 1, and 2 hours",
      "Strategies:",
      "  - Usual Care: No ROX-based intubation criteria (grace period not applicable)",
      "  - ROX < 3.85: Intubate when ROX falls below 3.85",
      "  - ROX < 4.88: Intubate when ROX falls below 4.88",
      "  - Time-varying ROX: Intubate when ROX < 2.85 (1-5h), < 3.47 (6-11h), or < 3.85 (≥12h)",
      "Censored (Protocol violation):",
      "  - Not Intubated: ROX < cutoff but patient was not intubated",
      "  - Intubated: ROX ≥ cutoff but patient was intubated",
      "Censored (NIV): Patients who switched from HFNC to non-invasive ventilation",
      "Following Protocol (Adherence):",
      "  - Percentage is relative to Total N",
      "  - Intubated: ROX < cutoff and patient was appropriately intubated (% of Following Protocol)",
      "  - Not Intubated: ROX remained ≥ cutoff and patient was not intubated (% of Following Protocol)",
      "ROX at intubation and Time to intubation: Median (IQR) among adherent intubated patients",
      "'-' indicates not applicable for Usual Care strategy"
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

# Display
print(ft)

# Save as Word document
doc <- read_docx() %>%
  body_add_par("Table. Comparison of Treatment Strategies by Grace Period (with NIV Censoring)",
               style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "outputs/ROX/Table_strategy_comparison_3arms_grace0to2hr_with_niv.docx")

cat("\n✓ Table saved: outputs/ROX/Table_strategy_comparison_3arms_grace0to2hr_with_niv.docx\n")

# Save raw data (for additional analysis)
saveRDS(summary_df, "outputs/ROX/strategy_comparison_3arms_data_grace0to2hr_with_niv.rds")
cat("✓ Raw data saved: outputs/ROX/strategy_comparison_3arms_data_grace0to2hr_with_niv.rds\n")

# Display summary
cat("\n=== Strategy Comparison Summary (Grace Periods: 0, 1, 2 hours) ===\n")
print(summary_df |> dplyr::select(Grace_Hours, Strategy, Total_N, N_Not_Censored,
                                  N_Censored_NIV, Pct_Censored_NIV,
                                  Pct_Intubated,
                                  ROX_at_Intub_Median, Time_to_Intub_Median,
                                  Pct_Censored, Pct_Died))

# Create simplified comparison table with NIV censor ----
cat("\n\n=== Creating simplified comparison table with NIV censor ===\n")

simple_comparison_table <- summary_df |>
  mutate(
    `Grace Period (hours)` = sprintf("%d", Grace_Hours),
    `Strategy` = Strategy,
    `Total N` = Total_N,
    `Censored (Protocol), n (%)` = sprintf("%d (%.1f%%)", N_Censored, Pct_Censored),
    `Censored (NIV), n (%)` = sprintf("%d (%.1f%%)", N_Censored_NIV, Pct_Censored_NIV),
    `ROX at Intubation*` = sprintf("%.2f (%.2f-%.2f)",
                                  ROX_at_Intub_Median,
                                  ROX_at_Intub_Q25,
                                  ROX_at_Intub_Q75),
    `SpO2 at Intubation, %*` = sprintf("%.1f (%.1f-%.1f)",
                                       SpO2_at_Intub_Median,
                                       SpO2_at_Intub_Q25,
                                       SpO2_at_Intub_Q75),
    `FiO2 at Intubation, %*` = sprintf("%.1f (%.1f-%.1f)",
                                       FiO2_at_Intub_Median,
                                       FiO2_at_Intub_Q25,
                                       FiO2_at_Intub_Q75),
    `RR at Intubation, /min*` = sprintf("%.1f (%.1f-%.1f)",
                                        RR_at_Intub_Median,
                                        RR_at_Intub_Q25,
                                        RR_at_Intub_Q75),
    `Time to Intubation, hours*` = sprintf("%.1f (%.1f-%.1f)",
                                          Time_to_Intub_Median,
                                          Time_to_Intub_Q25,
                                          Time_to_Intub_Q75)
  ) |>
  dplyr::select(`Grace Period (hours)`, `Strategy`, `Total N`,
                `Censored (Protocol), n (%)`, `Censored (NIV), n (%)`,
                `ROX at Intubation*`,
                `SpO2 at Intubation, %*`, `FiO2 at Intubation, %*`, `RR at Intubation, /min*`,
                `Time to Intubation, hours*`)

# Create flextable
ft_simple <- flextable(simple_comparison_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 2), align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(
    values = c(
      "Grace period: Time allowed after ROX first drops below cutoff before censoring",
      "Results shown for grace periods of 0, 1, and 2 hours",
      "Strategies:",
      "  - Usual Care: No ROX-based intubation criteria (grace period not applicable)",
      "  - ROX < 3.85: Intubate when ROX falls below 3.85",
      "  - ROX < 4.88: Intubate when ROX falls below 4.88",
      "  - Time-varying ROX: Intubate when ROX < 2.85 (1-5h), < 3.47 (6-11h), or < 3.85 (≥12h)",
      "Censored (Protocol): Patients who deviated from the assigned treatment strategy",
      "Censored (NIV): Patients who switched from HFNC to non-invasive ventilation",
      "Values are median (IQR)",
      "* Among patients who initiate intubation without being censored."
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

# Display
print(ft_simple)

# Save as Word document
doc_simple <- read_docx() %>%
  body_add_par("Table. Simplified Comparison of Treatment Strategies by Grace Period (with NIV)",
               style = "heading 1") %>%
  body_add_flextable(ft_simple)

print(doc_simple, target = "outputs/ROX/Table_strategy_comparison_simplified_grace0to2hr_with_niv.docx")

cat("\n✓ Simplified table saved: outputs/ROX/Table_strategy_comparison_simplified_grace0to2hr_with_niv.docx\n")

# Create Grace period 1 only table ----
cat("\n\n=== Creating Grace period 1 only table with NIV censor ===\n")

grace1_only_table <- summary_df |>
  filter(Grace_Hours == 1) |>
  mutate(
    `Strategy` = Strategy,
    `Total N` = Total_N,
    `Censored (Protocol), n (%)` = sprintf("%d (%.1f%%)", N_Censored, Pct_Censored),
    `Censored (NIV), n (%)` = sprintf("%d (%.1f%%)", N_Censored_NIV, Pct_Censored_NIV),
    `ROX at Intubation*` = sprintf("%.2f (%.2f-%.2f)",
                                  ROX_at_Intub_Median,
                                  ROX_at_Intub_Q25,
                                  ROX_at_Intub_Q75),
    `SpO2 at Intubation, %*` = sprintf("%.1f (%.1f-%.1f)",
                                       SpO2_at_Intub_Median,
                                       SpO2_at_Intub_Q25,
                                       SpO2_at_Intub_Q75),
    `FiO2 at Intubation, %*` = sprintf("%.1f (%.1f-%.1f)",
                                       FiO2_at_Intub_Median,
                                       FiO2_at_Intub_Q25,
                                       FiO2_at_Intub_Q75),
    `RR at Intubation, /min*` = sprintf("%.1f (%.1f-%.1f)",
                                        RR_at_Intub_Median,
                                        RR_at_Intub_Q25,
                                        RR_at_Intub_Q75),
    `Time to Intubation, hours*` = sprintf("%.1f (%.1f-%.1f)",
                                          Time_to_Intub_Median,
                                          Time_to_Intub_Q25,
                                          Time_to_Intub_Q75)
  ) |>
  dplyr::select(`Strategy`, `Total N`,
                `Censored (Protocol), n (%)`, `Censored (NIV), n (%)`,
                `ROX at Intubation*`,
                `SpO2 at Intubation, %*`, `FiO2 at Intubation, %*`, `RR at Intubation, /min*`,
                `Time to Intubation, hours*`)

# Create flextable
ft_grace1 <- flextable(grace1_only_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(
    values = c("Median (IQR)",
      "Usual Care: No ROX-based intubation criteria, ROX < 3.85: Intubate when ROX falls below 3.85, ROX < 4.88: Intubate when ROX falls below 4.88, Time-varying ROX: Intubate when ROX < 2.85 (1-5h), < 3.47 (6-11h), or < 3.85 (≥12h), Censored (Protocol): Patients who deviated from the assigned treatment strategy, Censored (NIV): Patients who switched from HFNC to non-invasive ventilation",
      "* Among patients who were intubated from HFNC without being censored."
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

# Display
print(ft_grace1)

# Save as Word document
# Output paths
path1 <- "outputs/ROX/S.Table 5_strategy_comparison_censor_rox_time_with_niv.docx"
path2 <- "C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Table 5_strategy_comparison_censor_rox_time_with_niv.docx"

# Optimize table for Word document
ft_word <- ft_grace1 %>%
  set_table_properties(layout = "autofit", width = 1) %>%
  fontsize(size = 9, part = "all")

save_as_docx(
  "S.Table. Comparison of Treatment Strategies with Grace Period of 1 Hour (with NIV Censoring)" = ft_word,
  path = path1
)

save_as_docx(
  "S.Table. Comparison of Treatment Strategies with Grace Period of 1 Hour (with NIV Censoring)" = ft_word,
  path = path2
)

cat("\n✓ Saved to:", path1)
cat("\n✓ Saved to:", path2, "\n")

# ROX at intubation by anchor_year_group (Usual Care, Grace 1hr) ----
cat("\n\n=== ROX at intubation by anchor_year_group (Usual Care, Grace 1hr) ===\n")

# Usual Care uses arm0 (no censoring by protocol)
arm0_temp <- copy(arm0)
arm0_temp[, my.censor := 0]

# Patient-level stats for Usual Care
patient_stats_uc <- arm0_temp[, .(
  intubated = any(treat == 1),
  rox_at_intubation = {
    idx <- which(treat == 1 & highflow_lag == 1)[1]
    if(!is.na(idx) && length(idx) > 0) as.numeric(rox[idx]) else NA_real_
  },
  spo2_at_intubation = {
    idx <- which(treat == 1 & highflow_lag == 1)[1]
    if(!is.na(idx) && length(idx) > 0) as.numeric(spo2[idx]) else NA_real_
  },
  fio2_at_intubation = {
    idx <- which(treat == 1 & highflow_lag == 1)[1]
    if(!is.na(idx) && length(idx) > 0) as.numeric(fio2[idx]) else NA_real_
  },
  resp_rate_at_intubation = {
    idx <- which(treat == 1 & highflow_lag == 1)[1]
    if(!is.na(idx) && length(idx) > 0) as.numeric(resp_rate[idx]) else NA_real_
  },
  time_to_intubation = {
    idx <- which(treat == 1 & highflow_lag == 1)[1]
    if(!is.na(idx) && length(idx) > 0) as.numeric(time[idx]) else NA_real_
  },
  anchor_year_group = anchor_year_group[1]
), by = id]

# Filter to intubated patients from HFNC
intubated_uc <- patient_stats_uc[intubated == TRUE & !is.na(rox_at_intubation)]

# Total patients and intubated patients by anchor_year_group
total_by_year <- patient_stats_uc[, .N, by = anchor_year_group][order(anchor_year_group)]
setnames(total_by_year, "N", "Total_N")

# Summarise by anchor_year_group
stats_by_year <- intubated_uc[, .(
  N_intubated = .N,
  ROX_med = median(rox_at_intubation, na.rm = TRUE),
  ROX_q25 = quantile(rox_at_intubation, 0.25, na.rm = TRUE),
  ROX_q75 = quantile(rox_at_intubation, 0.75, na.rm = TRUE),
  SpO2_med = median(spo2_at_intubation, na.rm = TRUE),
  SpO2_q25 = quantile(spo2_at_intubation, 0.25, na.rm = TRUE),
  SpO2_q75 = quantile(spo2_at_intubation, 0.75, na.rm = TRUE),
  FiO2_med = median(fio2_at_intubation, na.rm = TRUE),
  FiO2_q25 = quantile(fio2_at_intubation, 0.25, na.rm = TRUE),
  FiO2_q75 = quantile(fio2_at_intubation, 0.75, na.rm = TRUE),
  RR_med = median(resp_rate_at_intubation, na.rm = TRUE),
  RR_q25 = quantile(resp_rate_at_intubation, 0.25, na.rm = TRUE),
  RR_q75 = quantile(resp_rate_at_intubation, 0.75, na.rm = TRUE),
  Time_med = median(time_to_intubation, na.rm = TRUE),
  Time_q25 = quantile(time_to_intubation, 0.25, na.rm = TRUE),
  Time_q75 = quantile(time_to_intubation, 0.75, na.rm = TRUE)
), by = anchor_year_group][order(anchor_year_group)]

# Merge total N
stats_by_year <- merge(stats_by_year, total_by_year, by = "anchor_year_group")

cat("\nIntubation stats from HFNC (Usual Care) by anchor_year_group:\n")
print(stats_by_year)

# Create flextable
year_table <- stats_by_year[, .(
  `Year Group` = anchor_year_group,
  `Total N` = Total_N,
  `N Intubated from HFNC` = N_intubated,
  `ROX at Intubation` = sprintf("%.2f (%.2f-%.2f)", ROX_med, ROX_q25, ROX_q75),
  `SpO2 at Intubation, %` = sprintf("%.1f (%.1f-%.1f)", SpO2_med, SpO2_q25, SpO2_q75),
  `FiO2 at Intubation, %` = sprintf("%.1f (%.1f-%.1f)", FiO2_med, FiO2_q25, FiO2_q75),
  `RR at Intubation, /min` = sprintf("%.1f (%.1f-%.1f)", RR_med, RR_q25, RR_q75),
  `Time to Intubation, hours` = sprintf("%.1f (%.1f-%.1f)", Time_med, Time_q25, Time_q75)
)]

ft_rox_year <- flextable(year_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  add_footer_lines(
    values = c(
      "Usual Care strategy",
      "Values are median (IQR) among patients intubated from HFNC"
    )
  ) %>%
  fontsize(part = "footer", size = 9) %>%
  italic(part = "footer")

print(ft_rox_year)

# Save as Word document
ft_rox_year_word <- ft_rox_year %>%
  set_table_properties(layout = "autofit", width = 1) %>%
  fontsize(size = 9, part = "all")

path1_year <- "outputs/ROX/S.Table 4_intubation_stats_by_anchor_year_group.docx"
path2_year <- "C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Table 4_intubation_stats_by_anchor_year_group.docx"

save_as_docx(
  #"S.Table 4. Intubation Stats from HFNC by Year Group (Usual Care)" = 
  ft_rox_year_word,
  path = path1_year
)

save_as_docx(
  #"S.Table 4. Intubation Stats from HFNC by Year Group (Usual Care)" = 
  ft_rox_year_word,
  path = path2_year
)

cat("\n✓ Saved to:", path1_year)
cat("\n✓ Saved to:", path2_year, "\n")

# Supplementary Figure 1. Histogram of ROX index at intubation in HFNC patients ----

# Extract ROX index at intubation for patients intubated from HFNC
rox_at_intubation <- dt[, .(
  rox_at_intub = {
    idx <- which(treat == 1 & highflow_lag == 1)[1]
    if (!is.na(idx)) as.numeric(rox[idx]) else NA_real_
  }
), by = id][!is.na(rox_at_intub)]

summary(rox_at_intubation$rox_at_intub)

cat(sprintf("\nHFNC patients intubated: n = %d\n", nrow(rox_at_intubation)))
cat(sprintf("Median (IQR): %.2f (%.2f-%.2f)\n",
            median(rox_at_intubation$rox_at_intub),
            quantile(rox_at_intubation$rox_at_intub, 0.25),
            quantile(rox_at_intubation$rox_at_intub, 0.75)))

sfig1 <- ggplot(rox_at_intubation, aes(x = rox_at_intub)) +
  geom_histogram(binwidth = 1, fill = "#2E9FDF", color = "white", alpha = 0.8) +
  labs(
    #title = "Supplementary Figure 1.",
    #subtitle = "Histogram of ROX index at intubation in HFNC patients",
    x = "ROX index at intubation",
    y = "Number of patients"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

print(sfig1)

ggsave("outputs/ROX/S.Fig1_histogram_rox_at_intubation_hfnc.png",
       sfig1, width = 8, height = 6, dpi = 300)

ggsave("C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Fig1_histogram_rox_at_intubation_hfnc.png",
       sfig1, width = 8, height = 6, dpi = 300)

cat("✓ Supplementary Figure 1 saved\n")
