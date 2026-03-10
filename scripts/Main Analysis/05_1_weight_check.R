# =============================================================================
# Weight Distribution Check for 05_1 (ROX < 3.85, v4)
# Runs a single iteration (no bootstrap) to inspect w_boot / w_boot_99
# =============================================================================

rm(list = ls())
gc()

library(pacman)
p_load(tidyverse, here, data.table, lubridate, skimr)

# --- Load data ---------------------------------------------------------------
dt <- readRDS("./data/ROX/imputed_data_for_analysis_m5_imp1.rds")
dt <- as.data.table(dt)
dt <- dt[elig_1 == 1]

options(scipen = 999)
K <- 721
index <- 3.85
GRACE_HOURS <- 1

# --- Functions (same as 05_1_Bootstrap_CCW_3.85.R) ---------------------------

arm1.censor.grace <- function(rox, treat_cum, treat, highflow_lag, cutoff = index,
                              grace_hours = 1){
  n <- length(rox)
  my.censor <- rep(0, n)
  first.cens <- NA
  treat_cum_lag <- c(0, treat_cum[1:(n-1)])
  first_below_cutoff <- NA

  for(i in 1:n){
    if(!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
       !is.na(highflow_lag[i]) && highflow_lag[i] == 1){
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

# --- Run on original sample (no resampling) ----------------------------------
cat("=== Weight Distribution Check (ROX <", index, ", Grace =", GRACE_HOURS, "hr) ===\n\n")

dt_ids <- data.frame(id = unique(dt$id))
indices <- 1:nrow(dt_ids)   # original sample, no resampling

boot.ids <- data.frame(id = dt_ids$id[indices])
boot.ids$bid <- 1:nrow(boot.ids)

d <- merge(boot.ids, dt, by = "id")
d <- as.data.table(d)

# Cloning
arm0.d <- copy(d)[, arm := 0]
arm1.d <- copy(d)[, arm := 1]

# Censoring (arm1 only)
arm0.d.c <- copy(arm0.d)
arm0.d.c$my.censor <- 0

arm1.d.c <- copy(arm1.d)
arm1.d.c[, my.censor := arm1.censor.grace(rox, treat_cum, treat, highflow_lag,
                                           cutoff = index, grace_hours = GRACE_HOURS), by = bid]
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
    #gcs_category + ph + po2 + pco2 + sofa_24hours + 
    #crrt + 
    #vasopressor +
    #treatment_limitation +
    cal_time + cal_timesqr,
  family = binomial(link = "logit"),
  data = d[d$treat_cum_lag1 == 0, ]
))

d$psc.denom.a <- predict(psc.denom.a, d, type = "response")
d$psc.denom.a.t <- ifelse(d$treat == 1, d$psc.denom.a, 1 - d$psc.denom.a)
d$psc.denom.a.t <- ifelse(d$treat_cum_lag1 == 0, d$psc.denom.a.t, 1)

# Censoring weights - denominator
psc.denom.c <- suppressWarnings(glm(
  censor == 0 ~ age + female + medicare_medicaid + language_grouped + Married +
    race_grouped + antibiotic_use +
    elixhauser_vanwalraven + hr_at_hfnc_start +
    heart_rate + mbp + temperature + resp_rate + spo2 + fio2 +
    #gcs_category + ph + po2 + pco2 + sofa_24hours + 
    #crrt + 
    #vasopressor +
    #treatment_limitation +
    cal_time + cal_timesqr,
  family = binomial(link = "logit"),
  data = d[d$treat_cum_lag1 == 0, ]
))

d$psc.denom.c <- predict(psc.denom.c, d, type = "response")
d$psc.denom.c.t <- ifelse(d$censor == 0, d$psc.denom.c, 1 - d$psc.denom.c)
d$psc.denom.c.t <- ifelse(d$treat_cum_lag1 == 0, d$psc.denom.c.t, 1)

# Censoring weight: cumulative product on pre-cloning data
d <- d %>%
  arrange(bid, cal_time) %>%
  group_by(bid) %>%
  mutate(w_c_boot = 1 / cumprod(psc.denom.c.t)) %>%
  ungroup()

# Merge per-hour treatment weight contributions to cloned data
d_w <- d[, c("bid", "cal_time", "psc.denom.a.t", "w_c_boot")]
both.arms.sw.d <- merge(d_w, both.arms.d, by = c("bid", "cal_time"))

# Grace period: set treatment weight contribution = 1
both.arms.sw.d <- as.data.table(both.arms.sw.d)
both.arms.sw.d[arm == 1, grace_flag := arm1.grace.flag(rox, treat_cum, highflow_lag,
               cutoff = index, grace_hours = GRACE_HOURS), by = bid]
both.arms.sw.d[is.na(grace_flag), grace_flag := 0]
both.arms.sw.d[grace_flag == 1, psc.denom.a.t := 1]

# Treatment weight: cumulative product on cloned data by (bid, arm)
both.arms.sw.d <- both.arms.sw.d %>%
  arrange(bid, arm, cal_time) %>%
  group_by(bid, arm) %>%
  mutate(w_a_boot = 1 / cumprod(psc.denom.a.t)) %>%
  ungroup()

# Arm-specific final weights
both.arms.sw.d <- both.arms.sw.d %>%
  mutate(w_boot = ifelse(arm == 1, w_a_boot * w_c_boot, w_c_boot))

# Truncate at 99th percentile
threshold_99 <- quantile(both.arms.sw.d$w_boot, 0.99, na.rm = TRUE)
both.arms.sw.d$w_boot_99 <- pmin(both.arms.sw.d$w_boot, threshold_99)

# =============================================================================
# Weight diagnostics
# =============================================================================

cat("\n--- 1. Overall w_boot ---\n")
cat("N person-hours:", nrow(both.arms.sw.d), "\n")
print(summary(both.arms.sw.d$w_boot))
cat("SD:", sd(both.arms.sw.d$w_boot, na.rm = TRUE), "\n")

cat("\n--- 2. Overall w_boot_99 (truncated at 99th pctl) ---\n")
cat("99th percentile threshold:", threshold_99, "\n")
cat("N truncated:", sum(both.arms.sw.d$w_boot > threshold_99, na.rm = TRUE), "\n")
print(summary(both.arms.sw.d$w_boot_99))
cat("SD:", sd(both.arms.sw.d$w_boot_99, na.rm = TRUE), "\n")

cat("\n--- 3. Percentiles of w_boot ---\n")
probs <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975,
           0.98,0.99, 0.995, 0.999, 1.00)
print(quantile(both.arms.sw.d$w_boot, probs, na.rm = TRUE))

cat("\n--- 4. By arm ---\n")
both.arms.sw.d %>%
  group_by(arm) %>%
  summarise(
    n = n(),
    mean_w = mean(w_boot, na.rm = TRUE),
    median_w = median(w_boot, na.rm = TRUE),
    sd_w = sd(w_boot, na.rm = TRUE),
    max_w = max(w_boot, na.rm = TRUE),
    p99_w = quantile(w_boot, 0.99, na.rm = TRUE),
    mean_w99 = mean(w_boot_99, na.rm = TRUE),
    median_w99 = median(w_boot_99, na.rm = TRUE),
    sd_w99 = sd(w_boot_99, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

cat("\n--- 5. Component weights ---\n")
cat("w_a_boot (treatment weight):\n")
print(summary(both.arms.sw.d$w_a_boot))
cat("w_c_boot (censoring weight):\n")
print(summary(both.arms.sw.d$w_c_boot))

cat("\n--- 6. Grace period flag ---\n")
cat("Total grace_flag == 1:", sum(both.arms.sw.d$grace_flag == 1, na.rm = TRUE), "\n")
cat("Unique patients with grace_flag:",
    length(unique(both.arms.sw.d$bid[both.arms.sw.d$grace_flag == 1])), "\n")

cat("\n--- 7. w_boot by arm and time period ---\n")
both.arms.sw.d %>%
  mutate(time_group = case_when(
    cal_time <= 24  ~ "0-24h",
    cal_time <= 72  ~ "25-72h",
    cal_time <= 168 ~ "73-168h (3-7d)",
    TRUE            ~ "169h+ (>7d)"
  )) %>%
  group_by(arm, time_group) %>%
  summarise(
    n = n(),
    mean_w = mean(w_boot, na.rm = TRUE),
    p99_w = quantile(w_boot, 0.99, na.rm = TRUE),
    max_w = max(w_boot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(arm, time_group) %>%
  print(n = Inf)

both.arms.sw.d %>%
  group_by(arm) %>%
  summarise(
    mean_w99 = mean(w_boot_99, na.rm = TRUE),
    median_w99 = median(w_boot_99, na.rm = TRUE),
    max_w99 = max(w_boot_99, na.rm = TRUE),
    pct_truncated = mean(w_boot > w_boot_99, na.rm = TRUE) * 100,
    .groups = "drop"
  )

both.arms.sw.d %>%
  filter(arm == 1) %>%
  summarise(
    p95 = quantile(w_boot, 0.95, na.rm = TRUE),
    p97.5 = quantile(w_boot, 0.975, na.rm = TRUE),
    p99 = quantile(w_boot, 0.99, na.rm = TRUE),
    pct_above_20 = mean(w_boot > 20, na.rm = TRUE) * 100
  )

cat("\n=== Done ===\n")
