# ============================================================================
# Bootstrap Analysis for ROX < 3.85 Strategy with Multiple Imputation (m=5)
# ============================================================================
# Purpose: Run bootstrap analysis across m=5 multiply imputed datasets
# Bootstrap: 100 iterations per imputation (500 total across all imputations)
# Grace period: 1 hour
# Strategy: ROX-guided intubation with cutoff ROX < 3.85
# ============================================================================

# Setup memory
rm(list = ls())
gc()
gc()

# Load libraries -----
library(pacman)
p_load(tidyverse, here, data.table, boot, lubridate, skimr, naniar, gtsummary,
       patchwork, ggsci, parallel, speedglm, doParallel, foreach,
       officer, flextable)

# Disable printing results in scientific notation
options(scipen = 999)

# Set time point for estimation of risks (time 0:720 hours; 30 days)
K <- 721

# Strategy parameters -----
index <- 3.85  # ROX cutoff
GRACE_HOURS <- 1  # Grace period in hours
N_IMPUTATIONS <- 5  # Number of imputed datasets
N_BOOTSTRAP <- 200  # Bootstrap iterations per imputation

cat("\n============================================================================\n")
cat("  Bootstrap Analysis: ROX < 3.85 Strategy (Multiple Imputation m=5)\n")
cat("============================================================================\n")
cat(sprintf("ROX Cutoff: %.2f\n", index))
cat(sprintf("Grace Period: %d hour(s)\n", GRACE_HOURS))
cat(sprintf("Number of Imputations: %d\n", N_IMPUTATIONS))
cat(sprintf("Bootstrap per Imputation: %d\n", N_BOOTSTRAP))
cat(sprintf("Total Bootstrap Samples: %d\n", N_IMPUTATIONS * N_BOOTSTRAP))
cat("============================================================================\n\n")

# Censoring function for "intubate when ROX < cutoff" strategy WITH GRACE PERIOD -----
arm1.censor.grace <- function(rox, treat_cum, treat, highflow_lag, cutoff = index,
                              grace_hours = 1){
  n <- length(rox)
  my.censor <- rep(0, n)
  first.cens <- NA

  treat_cum_lag <- c(0, treat_cum[1:(n-1)])
  first_below_cutoff <- NA

  for(i in 1:n){
    # Check only before intubation AND when HFNC is in use (lag1 to align with ROX timing)
    if(!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
       !is.na(highflow_lag[i]) && highflow_lag[i] == 1){

      # === Non-adherence 1: Response to ROX < cutoff (WITH GRACE PERIOD) ===
      if(!is.na(rox[i]) && rox[i] < cutoff){
        if(is.na(first_below_cutoff)){
          first_below_cutoff <- i
        }
      }

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
    # When highflow_lag=0, loop continues but adherence check is skipped
    # → No censoring occurs, follow-up continues
  }

  if(!is.na(first.cens)){
    my.censor[first.cens:n] <- 1
  }

  return(my.censor)
}

### Grace period flag function
arm1.grace.flag <- function(rox, treat_cum, highflow_lag, cutoff = index,
                            grace_hours = 1){
  n <- length(rox)
  grace <- rep(0, n)
  treat_cum_lag <- c(0, treat_cum[1:(n-1)])
  first_below_cutoff <- NA

  for(i in 1:n){
    if(!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
       !is.na(highflow_lag[i]) && highflow_lag[i] == 1){
      if(!is.na(rox[i]) && rox[i] < cutoff && is.na(first_below_cutoff)){
        first_below_cutoff <- i
      }
      if(!is.na(first_below_cutoff) &&
         i >= first_below_cutoff && i < first_below_cutoff + grace_hours){
        grace[i] <- 1
      }
    }
  }
  return(grace)
}

# Bootstrap function for cumulative incidence curves -----
clone.boot.graph <- function(data, indices, dt_full, grace_hours = 1){

  # Select individuals into each bootstrapped sample
  boot.ids <- data.frame(id = data$id[indices])
  boot.ids$bid <- 1:nrow(boot.ids)

  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- merge(boot.ids, dt_full, by = "id")
  d <- as.data.table(d)

  # Cloning
  arm0.d <- copy(d)[, arm := 0]
  arm1.d <- copy(d)[, arm := 1]

  # Censoring (arm1 only) WITH GRACE PERIOD
  arm0.d.c <- copy(arm0.d)
  arm0.d.c$my.censor <- 0

  arm1.d.c <- copy(arm1.d)
  arm1.d.c[, my.censor := arm1.censor.grace(rox, treat_cum, treat, highflow_lag, cutoff = index, grace_hours = grace_hours), by = bid]
  arm1.d.c[, good := cumsum(my.censor), by = bid]
  arm1.d.c <- arm1.d.c[good <= 1]
  arm1.d.c[, good := NULL]

  both.arms.d <- rbind(arm0.d.c, arm1.d.c)
  both.arms.d[, death := ifelse(my.censor == 1, NA, death)]

  # === IPTW estimation ===
  # Treatment weights - denominator
  psc.denom.a <- suppressWarnings(glm(
    treat ~ age + female + medicare_medicaid + language_grouped + Married +
      race_grouped + antibiotic_use + elixhauser_vanwalraven + hr_at_hfnc_start +
      heart_rate + mbp + temperature + resp_rate + spo2 + fio2 +
      gcs_category + ph + po2 + pco2 + sofa_24hours + crrt + vasopressor +
      treatment_limitation +
      cal_time + cal_timesqr,
    family = binomial(link = "logit"),
    data = d[d$treat_cum_lag1 == 0, ]
  ))

  d$psc.denom.a <- predict(psc.denom.a, d, type = "response")
  d$psc.denom.a.t <- ifelse(d$treat == 1, d$psc.denom.a, 1 - d$psc.denom.a)
  d$psc.denom.a.t <- ifelse(d$treat_cum_lag1 == 0, d$psc.denom.a.t, 1)

  # Treatment weights - numerator
  psc.numer.a <- suppressWarnings(glm(
    treat ~ cal_time + cal_timesqr,
    family = binomial(link = "logit"),
    data = d[d$treat_cum_lag1 == 0, ]
  ))

  d$psc.numer.a <- predict(psc.numer.a, d, type = "response")
  d$psc.numer.a.t <- ifelse(d$treat == 1, d$psc.numer.a, 1 - d$psc.numer.a)
  d$psc.numer.a.t <- ifelse(d$treat_cum_lag1 == 0, d$psc.numer.a.t, 1)

  # Censoring weights - denominator
  psc.denom.c <- suppressWarnings(glm(
    censor == 0 ~ age + female + medicare_medicaid + language_grouped + Married +
      race_grouped + antibiotic_use + elixhauser_vanwalraven + hr_at_hfnc_start +
      heart_rate + mbp + temperature + resp_rate + spo2 + fio2 +
      gcs_category + ph + po2 + pco2 + sofa_24hours + crrt + vasopressor +
      treatment_limitation +
      cal_time + cal_timesqr,
    family = binomial(link = "logit"),
    data = d[d$treat_cum_lag1 == 0, ]
  ))

  d$psc.denom.c <- predict(psc.denom.c, d, type = "response")
  d$psc.denom.c.t <- ifelse(d$censor == 0, d$psc.denom.c, 1 - d$psc.denom.c)
  d$psc.denom.c.t <- ifelse(d$treat_cum_lag1 == 0, d$psc.denom.c.t, 1)

  # Censoring weights - numerator
  psc.numer.c <- suppressWarnings(glm(
    censor == 0 ~ cal_time + cal_timesqr,
    family = binomial(link = "logit"),
    data = d[d$treat_cum_lag1 == 0, ]
  ))

  d$psc.numer.c <- predict(psc.numer.c, d, type = "response")
  d$psc.numer.c.t <- ifelse(d$censor == 0, d$psc.numer.c, 1 - d$psc.numer.c)
  d$psc.numer.c.t <- ifelse(d$treat_cum_lag1 == 0, d$psc.numer.c.t, 1)

  # Censoring weight: cumulative product on pre-cloning data (stabilized)
  d <- d %>%
    arrange(bid, cal_time) %>%
    group_by(bid) %>%
    mutate(sw_c_boot = cumprod(psc.numer.c.t) / cumprod(psc.denom.c.t)) %>%
    ungroup()

  # Merge per-hour treatment weight contributions to cloned data
  d_w <- d[, c("bid", "cal_time", "psc.numer.a.t", "psc.denom.a.t", "sw_c_boot")]
  both.arms.sw.d <- merge(d_w, both.arms.d, by = c("bid", "cal_time"))

  # Grace period: set treatment weight contribution = 1
  both.arms.sw.d <- as.data.table(both.arms.sw.d)
  both.arms.sw.d[arm == 1, grace_flag := arm1.grace.flag(rox, treat_cum, highflow_lag,
                 cutoff = index, grace_hours = grace_hours), by = bid]
  both.arms.sw.d[is.na(grace_flag), grace_flag := 0]
  both.arms.sw.d[grace_flag == 1, psc.numer.a.t := 1]
  both.arms.sw.d[grace_flag == 1, psc.denom.a.t := 1]

  # Treatment weight: cumulative product on cloned data by (bid, arm) (stabilized)
  both.arms.sw.d <- both.arms.sw.d %>%
    arrange(bid, arm, cal_time) %>%
    group_by(bid, arm) %>%
    mutate(sw_a_boot = cumprod(psc.numer.a.t) / cumprod(psc.denom.a.t)) %>%
    ungroup()

  # Arm-specific final weights (control arm is not artificially censored)
  both.arms.sw.d <- both.arms.sw.d %>%
    mutate(sw_ac_boot = ifelse(arm == 1, sw_a_boot * sw_c_boot, sw_c_boot))

  # Truncate weights at 99th percentile
  threshold_99_boot <- quantile(both.arms.sw.d$sw_ac_boot, 0.99, na.rm = TRUE)
  both.arms.sw.d$sw_ac_99_boot <- pmin(both.arms.sw.d$sw_ac_boot, threshold_99_boot)

  # Fit weighted pooled logistic regression
  fit.pool.boot <- suppressWarnings(glm(
    death == 1 ~ arm + cal_time + cal_timesqr +
      I(arm * cal_time) + I(arm * cal_timesqr),
    family = binomial(link = 'logit'),
    data = both.arms.sw.d,
    weights = sw_ac_99_boot
  ))

  # === Create dataset for prediction (all timepoints) ===
  K_boot <- 721  # 0-720 hours

  # Arm 0 (Usual Care)
  arm0 <- data.frame(
    cal_time = seq(0, K_boot - 1),
    arm = 0,
    cal_timesqr = (seq(0, K_boot - 1))^2
  )

  # Arm 1 (ROX Strategy)
  arm1 <- data.frame(
    cal_time = seq(0, K_boot - 1),
    arm = 1,
    cal_timesqr = (seq(0, K_boot - 1))^2
  )

  # Extract predicted values (discrete-time hazards)
  arm0$p.event0 <- predict(fit.pool.boot, arm0, type = "response")
  arm1$p.event1 <- predict(fit.pool.boot, arm1, type = "response")

  # Estimate survival probabilities from hazards
  arm0$surv0 <- cumprod(1 - arm0$p.event0)
  arm1$surv1 <- cumprod(1 - arm1$p.event1)

  # Estimate risks from survival probabilities
  arm0$risk0 <- 1 - arm0$surv0
  arm1$risk1 <- 1 - arm1$surv1

  # Merge data from two treatment groups
  graph.pred <- merge(
    arm0[, c("cal_time", "cal_timesqr", "p.event0", "surv0", "risk0")],
    arm1[, c("cal_time", "cal_timesqr", "p.event1", "surv1", "risk1")],
    by = c("cal_time", "cal_timesqr")
  )

  # Edit data frame to reflect that risks are estimated at the END of each interval
  graph.pred$time_0 <- graph.pred$cal_time + 1

  # Add zero at time 0
  zero <- data.frame(
    cal_time = 0,
    cal_timesqr = 0,
    p.event0 = 0,
    surv0 = 1,
    risk0 = 0,
    p.event1 = 0,
    surv1 = 1,
    risk1 = 0,
    time_0 = 0
  )

  graph <- rbind(zero, graph.pred)

  # Calculate risk difference and risk ratio
  graph$rd <- graph$risk1 - graph$risk0
  graph$rr <- graph$risk1 / graph$risk0

  # Return the complete graph (all timepoints)
  return(graph[, c("time_0", "risk0", "risk1", "rd", "rr")])
}

# ============================================================================
# Main Analysis Loop: Process Each Imputed Dataset
# ============================================================================

# Record script start time
script_start_time <- Sys.time()

# Check available CPU cores
n_cores <- detectCores()
cat("Available CPU cores:", n_cores, "\n")

# Use 80% of cores for parallel processing
use_cores <- max(1, floor(n_cores * 0.6))
cat("Using CPU cores:", use_cores, "\n\n")

# Loop through each imputation
for(imp in 1:N_IMPUTATIONS) {

  cat("\n============================================================================\n")
  cat(sprintf("  Processing Imputation %d/%d\n", imp, N_IMPUTATIONS))
  cat("============================================================================\n\n")

  # Load imputed data
  input_file <- sprintf("./data/ROX/imputed_data_for_analysis_m5_imp%d.rds", imp)
  cat("Loading data:", input_file, "\n")
  dt <- readRDS(input_file)

  # Convert to data.table and filter eligible patients
  dt <- as.data.table(dt)
  dt <- dt %>% filter(elig_1 == 1)

  cat(sprintf("Sample size: %s patients\n", format(nrow(dt %>% distinct(id)), big.mark = ",")))

  # Create patient ID list for bootstrapping
  dt_ids <- data.frame(id = unique(dt$id))

  # Set seed for reproducibility (different seed per imputation)
  seed_value <- 100 + imp * 11
  set.seed(seed_value)
  cat(sprintf("Random seed set to: %d\n\n", seed_value))

  # Setup parallel processing
  cl <- makeCluster(use_cores)
  registerDoParallel(cl)

  clusterExport(cl, c("dt", "dt_ids", "arm1.censor.grace", "arm1.grace.flag", "clone.boot.graph",
                      "GRACE_HOURS", "index"))
  clusterEvalQ(cl, {
    library(data.table)
    library(dplyr)
  })

  # Run bootstrap
  cat(sprintf("Starting bootstrap (R=%d, Grace Period=%d hour)\n", N_BOOTSTRAP, GRACE_HOURS))
  cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  start_time <- Sys.time()

  boot_graphs_list <- foreach(
    i = 1:N_BOOTSTRAP,
    .packages = c("data.table", "dplyr"),
    .errorhandling = "pass"
  ) %dopar% {
    set.seed(seed_value + i)
    indices <- sample(1:nrow(dt_ids), replace = TRUE)
    clone.boot.graph(dt_ids, indices, dt, grace_hours = GRACE_HOURS)
  }

  end_time <- Sys.time()
  elapsed_time <- end_time - start_time

  cat("\nEnd time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Elapsed time:", round(as.numeric(elapsed_time, units = "mins"), 1), "minutes\n")

  # Stop cluster
  stopCluster(cl)

  # Save results
  output_file <- sprintf("outputs/ROX/bootdata/bootstrap_results_3.85_imp%d.rds", imp)
  saveRDS(boot_graphs_list, output_file)
  cat("\nResults saved:", output_file, "\n")

  # Clean up
  rm(dt, dt_ids, boot_graphs_list, cl)
  gc()
}

cat("\n============================================================================\n")
cat("  All Imputations Completed Successfully\n")
cat("============================================================================\n")
cat(sprintf("\nTotal bootstrap samples generated: %d\n", N_IMPUTATIONS * N_BOOTSTRAP))
cat(sprintf("Results saved in: outputs/ROX/bootdata/bootstrap_results_3.85_imp*.rds\n\n"))

# Calculate and print total execution time
script_end_time <- Sys.time()
total_elapsed_time <- script_end_time - script_start_time

cat("============================================================================\n")
cat("  Total Execution Time\n")
cat("============================================================================\n")
cat("Script start time:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Script end time:  ", format(script_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat(sprintf("Total elapsed time: %.1f minutes (%.2f hours)\n",
    as.numeric(total_elapsed_time, units = "mins"),
    as.numeric(total_elapsed_time, units = "hours")))
cat("============================================================================\n\n")

cat("Next step: Run pooling analysis using Rubin's Rule (05_02_Bootstrap script)\n\n")
