#setup memory
rm(list = ls())
gc()
gc()

# load library -----
library(pacman)

p_load(tidyverse, here, data.table, boot, lubridate, skimr, naniar, gtsummary, officer, flextable,
       patchwork, ggsci, parallel, speedglm, doParallel, foreach)


# load data
#dt <- readRDS("./data/ROX/imputed_data_for_analysis.rds")
dt <-readRDS("./data/ROX/imputed_data_for_analysis_m5_imp1.rds")

dt <- as.data.table(dt)
class(dt)

# filter eligible patients
dt <- dt |> filter(elig_1 == 1)

# Disable printing results in scientific notation
options(scipen = 999)

# Set time point for estimation of risks (time 0:720 hours; 30 days)
K <- 721

# Cloning ----
# Clone each eligible individual twice
arm0 <- copy(dt)[, arm := 0]
arm1 <- copy(dt)[, arm := 1]

# Censoring functions ----
# 固定閾値用のセンサー関数
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

# Time-Dependent ROX用のセンサー関数
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
            # Grace期間内で遵守
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

# Function to create cumulative status data ----
create_cumulative_status <- function(data, max_time = 720) {
  # Get patient-level event times
  patient_events <- data[, .(
    death_time = if(any(!is.na(death) & death == 1)) as.numeric(min(time[!is.na(death) & death == 1])) else NA_real_,
    niv_censor_time = if(any(censor == 1 & my.censor == 0)) as.numeric(min(time[censor == 1 & my.censor == 0])) else NA_real_,
    rox_censor_time = if(any(my.censor == 1)) as.numeric(min(time[my.censor == 1])) else NA_real_
  ), by = id]

  # Create all time points for all patients (using data.table's CJ for cross join)
  all_times <- CJ(
    id = unique(patient_events$id),
    time = 0:max_time
  )

  # Merge with event times
  status_data <- merge(all_times, patient_events, by = "id", all.x = TRUE)
  setDT(status_data)  # Ensure it's a data.table

  # Assign status for each time point
  # Priority: Death > ROX censor > NIV censor > Following
  status_data[, status := fcase(
    !is.na(death_time) & time >= death_time, "Death",
    !is.na(rox_censor_time) & time >= rox_censor_time, "Censored: Protocol deviation (ROX index)",
    !is.na(niv_censor_time) & time >= niv_censor_time, "Censored: NIV use",
    default = "Following intubation strategy"
  )]

  return(status_data)
}

# Cumulative status plots for all strategies ----
grace_hours_values <- c(0, 1, 2)

# Store plots for final display
all_cumulative_plots <- list()
all_cumulative_plots_72h <- list()

for(grace_hours in grace_hours_values) {

  cat(sprintf("\n\n========== Creating cumulative status plot for Grace Period: %d hour(s) ==========\n", grace_hours))

  # List to store cumulative data for all strategies
  all_strategies_cumulative <- list()

  # Strategy 1: Usual Care ----
  cat("\n=== Processing Usual Care ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  both.arms.uc <- arm0.c[, death := ifelse(my.censor == 1, NA, death)]

  cumulative_uc <- create_cumulative_status(both.arms.uc)
  cumulative_uc[, strategy := "Usual Care"]
  all_strategies_cumulative[["Usual Care"]] <- cumulative_uc

  # Strategy 2: ROX < 3.85 ----
  cat("\n=== Processing ROX < 3.85 ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  arm1.c <- copy(arm1)[
    , my.censor := arm1.censor.grace.fixed(rox, treat_cum, treat, highflow_lag, cutoff = 3.85, grace_hours = grace_hours), by = id][
      , good := cumsum(my.censor), by = id][
        good <= 1][
          , good := NULL]

  both.arms.385 <- rbind(arm0.c, arm1.c)[, death := ifelse(my.censor == 1, NA, death)]

  # Only use intervention arm (arm = 1) data
  cumulative_385 <- create_cumulative_status(both.arms.385[arm == 1])
  cumulative_385[, strategy := "ROX < 3.85"]
  all_strategies_cumulative[["ROX < 3.85"]] <- cumulative_385

  # Strategy 3: ROX < 4.88 ----
  cat("\n=== Processing ROX < 4.88 ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  arm1.c <- copy(arm1)[
    , my.censor := arm1.censor.grace.fixed(rox, treat_cum, treat, highflow_lag, cutoff = 4.88, grace_hours = grace_hours), by = id][
      , good := cumsum(my.censor), by = id][
        good <= 1][
          , good := NULL]

  both.arms.488 <- rbind(arm0.c, arm1.c)[, death := ifelse(my.censor == 1, NA, death)]

  cumulative_488 <- create_cumulative_status(both.arms.488[arm == 1])
  cumulative_488[, strategy := "ROX < 4.88"]
  all_strategies_cumulative[["ROX < 4.88"]] <- cumulative_488

  # Strategy 4: Time-Dependent ROX ----
  cat("\n=== Processing Time-Dependent ROX ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  arm1.c <- copy(arm1)[
    , my.censor := arm1.censor.grace.timedep(rox, treat_cum, treat, highflow_lag, cal_time, grace_hours = grace_hours), by = id][
      , good := cumsum(my.censor), by = id][
        good <= 1][
          , good := NULL]

  both.arms.timedep <- rbind(arm0.c, arm1.c)[, death := ifelse(my.censor == 1, NA, death)]

  cumulative_timedep <- create_cumulative_status(both.arms.timedep[arm == 1])
  cumulative_timedep[, strategy := "Time-Dependent ROX"]
  all_strategies_cumulative[["Time-Dependent ROX"]] <- cumulative_timedep

  # Combine all strategies ----
  cumulative_combined <- rbindlist(all_strategies_cumulative)

  # Count patients in each status at each time point
  status_counts <- cumulative_combined[, .N, by = .(strategy, time, status)]

  # Set factor levels for proper ordering
  status_counts$strategy <- factor(status_counts$strategy,
                                    levels = c("Usual Care", "ROX < 3.85", "ROX < 4.88", "Time-Dependent ROX"))

  status_counts$status <- factor(status_counts$status,
                                  levels = c("Censored: Protocol deviation (ROX index)",
                                             "Censored: NIV use",
                                             "Death",
                                             "Following intubation strategy"))

  # Create plot ----
  cumulative_plot <- ggplot(status_counts, aes(x = time, y = N, fill = status)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    xlab("Time since baseline (hours)") +
    scale_x_continuous(breaks = c(0, 120, 240, 360, 480, 600, 720)) +
    ylab("Number of patients") +
    scale_fill_manual(
      values = c(
        "Censored: Protocol deviation (ROX index)" = "#878787",
        "Censored: NIV use" = "#14B5D0",
        "Death" = "#000000",
        "Following intubation strategy" = "#800000"
      ),
      breaks = c("Censored: Protocol deviation (ROX index)",
                 "Censored: NIV use",
                 "Death",
                 "Following intubation strategy")
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    facet_wrap(~strategy, scales = "free_y", ncol = 2)
    # Note: Title removed for cleaner presentation

  # Store plot for later display
  all_cumulative_plots[[paste0("grace_", grace_hours, "hr")]] <- cumulative_plot

  # Save full time range plot
  ggsave(sprintf("outputs/ROX/S.Fig_cumulative_status_4strategies_grace%dhr.png", grace_hours),
         cumulative_plot,
         dpi = 1200,
         width = 14,
         height = 10)

  cat(sprintf("✓ Saved: outputs/ROX/S.Fig_cumulative_status_4strategies_grace%dhr.png\n", grace_hours))

  # Create 0-72 hours plot ----
  status_counts_72h <- status_counts[time <= 72]

  cumulative_plot_72h <- ggplot(status_counts_72h,
                                 aes(x = time, y = N, fill = status)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    xlab("Time since baseline (hours)") +
    scale_x_continuous(
      breaks = c(0, 12, 24, 36, 48, 60, 72),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(xlim = c(0, 72)) +
    ylab("Number of patients") +
    scale_fill_manual(
      values = c(
        "Censored: Protocol deviation (ROX index)" = "#878787",
        "Censored: NIV use" = "#14B5D0",
        "Death" = "#000000",
        "Following intubation strategy" = "#800000"
      ),
      breaks = c("Censored: Protocol deviation (ROX index)",
                 "Censored: NIV use",
                 "Death",
                 "Following intubation strategy")
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      axis.line = element_line(colour = "black"),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    facet_wrap(~strategy, scales = "free_y", ncol = 2)
    # Note: Title removed for cleaner presentation

  # Store 72h plot for later display
  all_cumulative_plots_72h[[paste0("grace_", grace_hours, "hr_72h")]] <- cumulative_plot_72h

  # Save 72h plot
  ggsave(sprintf("outputs/ROX/S.Fig_cumulative_status_4strategies_72h_grace%dhr.png", grace_hours),
         cumulative_plot_72h,
         dpi = 1200,
         width = 14,
         height = 10)

  cat(sprintf("✓ Saved: outputs/ROX/S.Fig_cumulative_status_4strategies_72h_grace%dhr.png\n", grace_hours))
}

# Display all plots ----
cat("\n\n========== Displaying all cumulative status plots ==========\n")

cat("\n=== Full time range plots ===\n")
for(i in seq_along(all_cumulative_plots)) {
  cat(sprintf("Displaying: %s\n", names(all_cumulative_plots)[i]))
  print(all_cumulative_plots[[i]])
}

cat("\n=== 0-72 hours plots ===\n")
for(i in seq_along(all_cumulative_plots_72h)) {
  cat(sprintf("Displaying: %s\n", names(all_cumulative_plots_72h)[i]))
  print(all_cumulative_plots_72h[[i]])
}

# Create combined comparison plots using patchwork ----
cat("\n\n=== Creating combined comparison plots ===\n")

# Add grace period labels to plots for combined display
plot_0hr <- all_cumulative_plots[[1]] +
  labs(subtitle = "Grace Period: 0 hour", tag = "A") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_1hr <- all_cumulative_plots[[2]] +
  labs(subtitle = "Grace Period: 1 hour", tag = "B") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_2hr <- all_cumulative_plots[[3]] +
  labs(subtitle = "Grace Period: 2 hours", tag = "C") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

# Combined plot: All grace periods (full time range)
combined_full <- plot_0hr / plot_1hr / plot_2hr

print(combined_full)

ggsave("outputs/ROX/S.Fig_cumulative_status_comparison_all_grace_periods.png",
       combined_full,
       dpi = 1200,
       width = 14,
       height = 24)

cat("✓ Saved: outputs/ROX/S.Fig_cumulative_status_comparison_all_grace_periods.png\n")

# Add grace period labels to 72h plots
plot_0hr_72h <- all_cumulative_plots_72h[[1]] +
  labs(subtitle = "Grace Period: 0 hour", tag = "A") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_1hr_72h <- all_cumulative_plots_72h[[2]] +
  labs(subtitle = "Grace Period: 1 hour", tag = "B") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_2hr_72h <- all_cumulative_plots_72h[[3]] +
  labs(subtitle = "Grace Period: 2 hours", tag = "C") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

# Combined plot: All grace periods (0-72h)
combined_72h <- plot_0hr_72h / plot_1hr_72h / plot_2hr_72h

print(combined_72h)

ggsave("outputs/ROX/S.Fig_cumulative_status_comparison_all_grace_periods_72h.png",
       combined_72h,
       dpi = 1200,
       width = 14,
       height = 24)

cat("✓ Saved: outputs/ROX/S.Fig_cumulative_status_comparison_all_grace_periods_72h.png\n")

cat("\n\n========== All cumulative status plots completed! ==========\n")
cat("\nGenerated files:\n")
cat("\nIndividual grace period plots:\n")
for(grace_hours in grace_hours_values) {
  cat(sprintf("  - S.Fig_cumulative_status_4strategies_grace%dhr.png\n", grace_hours))
  cat(sprintf("  - S.Fig_cumulative_status_4strategies_72h_grace%dhr.png\n", grace_hours))
}
cat("\nCombined comparison plots:\n")
cat("  - S.Fig_cumulative_status_comparison_all_grace_periods.png\n")
cat("  - S.Fig_cumulative_status_comparison_all_grace_periods_72h.png\n")

cat("\n✓ All cumulative status plots have been saved to outputs/ROX/ directory\n")
cat("✓ All plots are now displayed in the Plots panel for review\n")
