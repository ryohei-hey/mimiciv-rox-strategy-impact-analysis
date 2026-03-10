#setup memory
rm(list = ls())
gc()
gc()

# load library -----
library(pacman)

p_load(tidyverse, here, data.table, boot, lubridate, skimr, naniar, gtsummary, officer, flextable,
       patchwork, ggsci, parallel, speedglm, doParallel, foreach)


# load data
dt <- readRDS("./data/ROX/imputed_data_for_analysis.rds")

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

# Time-varying ROX用のセンサー関数
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

# Censor plots for all strategies ----
grace_hours_values <- c(0, 1, 2)

# Store plots for final display
all_plots <- list()
all_plots_72h <- list()

for(grace_hours in grace_hours_values) {

  cat(sprintf("\n\n========== Creating censoring plot for Grace Period: %d hour(s) ==========\n", grace_hours))

  # List to store data for all strategies
  all_strategies_data <- list()

  # Strategy 1: Usual Care ----
  cat("\n=== Processing Usual Care ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]

  both.arms.uc <- arm0.c[, death := ifelse(my.censor == 1, NA, death)]

  # Prepare for plotting
  both.arms.new.cens <- both.arms.uc[which(both.arms.uc$censor == 1 & both.arms.uc$my.censor == 0), ]
  both.arms.new.cens$time <- both.arms.new.cens$time + 1
  both.arms.no.cens <- both.arms.uc
  both.arms.no.cens$censor <- 0
  both.arms.plot <- rbind(both.arms.no.cens, both.arms.new.cens)

  both.arms.plot$status <- ifelse(
    !is.na(both.arms.plot$death) & both.arms.plot$death == 1, 3,
    ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 0, 2,
           ifelse(both.arms.plot$censor == 1 & both.arms.plot$my.censor == 0, 1,
                  ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 1, 0, NA))))

  plot.counts.uc <- both.arms.plot %>%
    group_by(time, status) %>%
    dplyr::summarize(total_count = n(), .groups = 'drop') %>%
    mutate(strategy = "Usual Care")

  all_strategies_data[["Usual Care"]] <- plot.counts.uc

  # Strategy 2: ROX < 3.85 ----
  cat("\n=== Processing ROX < 3.85 ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  arm1.c <- copy(arm1)[
    , my.censor := arm1.censor.grace.fixed(rox, treat_cum, treat, highflow_lag, cutoff = 3.85, grace_hours = grace_hours), by = id][
      , good := cumsum(my.censor), by = id][
        good <= 1][
          , good := NULL]

  both.arms.385 <- rbind(arm0.c, arm1.c)[, death := ifelse(my.censor == 1, NA, death)]

  # Prepare for plotting
  both.arms.new.cens <- both.arms.385[which(both.arms.385$censor == 1 & both.arms.385$my.censor == 0), ]
  both.arms.new.cens$time <- both.arms.new.cens$time + 1
  both.arms.no.cens <- both.arms.385
  both.arms.no.cens$censor <- 0
  both.arms.plot <- rbind(both.arms.no.cens, both.arms.new.cens)

  both.arms.plot$status <- ifelse(
    !is.na(both.arms.plot$death) & both.arms.plot$death == 1, 3,
    ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 0, 2,
           ifelse(both.arms.plot$censor == 1 & both.arms.plot$my.censor == 0, 1,
                  ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 1, 0, NA))))

  plot.counts.385 <- both.arms.plot %>%
    group_by(time, status, arm) %>%
    dplyr::summarize(total_count = n(), .groups = 'drop') %>%
    mutate(strategy = ifelse(arm == 0, "Usual Care", "ROX < 3.85"))

  all_strategies_data[["ROX < 3.85"]] <- plot.counts.385

  # Strategy 3: ROX < 4.88 ----
  cat("\n=== Processing ROX < 4.88 ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  arm1.c <- copy(arm1)[
    , my.censor := arm1.censor.grace.fixed(rox, treat_cum, treat, highflow_lag, cutoff = 4.88, grace_hours = grace_hours), by = id][
      , good := cumsum(my.censor), by = id][
        good <= 1][
          , good := NULL]

  both.arms.488 <- rbind(arm0.c, arm1.c)[, death := ifelse(my.censor == 1, NA, death)]

  # Prepare for plotting
  both.arms.new.cens <- both.arms.488[which(both.arms.488$censor == 1 & both.arms.488$my.censor == 0), ]
  both.arms.new.cens$time <- both.arms.new.cens$time + 1
  both.arms.no.cens <- both.arms.488
  both.arms.no.cens$censor <- 0
  both.arms.plot <- rbind(both.arms.no.cens, both.arms.new.cens)

  both.arms.plot$status <- ifelse(
    !is.na(both.arms.plot$death) & both.arms.plot$death == 1, 3,
    ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 0, 2,
           ifelse(both.arms.plot$censor == 1 & both.arms.plot$my.censor == 0, 1,
                  ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 1, 0, NA))))

  plot.counts.488 <- both.arms.plot %>%
    group_by(time, status, arm) %>%
    dplyr::summarize(total_count = n(), .groups = 'drop') %>%
    mutate(strategy = ifelse(arm == 0, "Usual Care", "ROX < 4.88"))

  all_strategies_data[["ROX < 4.88"]] <- plot.counts.488

  # Strategy 4: Time-varying ROX ----
  cat("\n=== Processing Time-varying ROX ===\n")
  arm0.c <- copy(arm0)[, my.censor := 0]
  arm1.c <- copy(arm1)[
    , my.censor := arm1.censor.grace.timedep(rox, treat_cum, treat, highflow_lag, cal_time, grace_hours = grace_hours), by = id][
      , good := cumsum(my.censor), by = id][
        good <= 1][
          , good := NULL]

  both.arms.timedep <- rbind(arm0.c, arm1.c)[, death := ifelse(my.censor == 1, NA, death)]

  # Prepare for plotting
  both.arms.new.cens <- both.arms.timedep[which(both.arms.timedep$censor == 1 & both.arms.timedep$my.censor == 0), ]
  both.arms.new.cens$time <- both.arms.new.cens$time + 1
  both.arms.no.cens <- both.arms.timedep
  both.arms.no.cens$censor <- 0
  both.arms.plot <- rbind(both.arms.no.cens, both.arms.new.cens)

  both.arms.plot$status <- ifelse(
    !is.na(both.arms.plot$death) & both.arms.plot$death == 1, 3,
    ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 0, 2,
           ifelse(both.arms.plot$censor == 1 & both.arms.plot$my.censor == 0, 1,
                  ifelse(both.arms.plot$censor == 0 & both.arms.plot$my.censor == 1, 0, NA))))

  plot.counts.timedep <- both.arms.plot %>%
    group_by(time, status, arm) %>%
    dplyr::summarize(total_count = n(), .groups = 'drop') %>%
    mutate(strategy = ifelse(arm == 0, "Usual Care", "Time-varying ROX"))

  all_strategies_data[["Time-varying ROX"]] <- plot.counts.timedep

  # Combine all strategies ----
  # For strategies with 2 arms, keep only the intervention arm for comparison
  plot.counts.385.arm1 <- plot.counts.385 %>% filter(strategy == "ROX < 3.85")
  plot.counts.488.arm1 <- plot.counts.488 %>% filter(strategy == "ROX < 4.88")
  plot.counts.timedep.arm1 <- plot.counts.timedep %>% filter(strategy == "Time-varying ROX")

  # Combine Usual Care + 3 intervention strategies
  plot_combined <- bind_rows(
    plot.counts.uc,
    plot.counts.385.arm1,
    plot.counts.488.arm1,
    plot.counts.timedep.arm1
  )

  # Set strategy as factor with specific order
  plot_combined$strategy <- factor(plot_combined$strategy,
                                    levels = c("Usual Care", "ROX < 3.85", "ROX < 4.88", "Time-varying ROX"))

  # Set status as factor with specific order (Death before Following strategy)
  plot_combined$status <- factor(plot_combined$status,
                                  levels = c("0", "1", "3", "2"))

  # Create plot ----
  cens_plot <- ggplot(plot_combined, aes(x = time, y = total_count, fill = status)) +
    geom_bar(position = "stack", stat = "identity") +
    xlab("Time since baseline (hours)") +
    scale_x_continuous(breaks = c(0, 120, 240, 360, 480, 600, 720)) +
    ylab("Number of observations") +
    scale_fill_manual(
      values = c("0" = "#878787", "1" = "#14B5D0", "3" = "#000000", "2" = "#800000"),
      labels = c(
        "0" = "Censored: Protocol deviation (ROX index)",
        "1" = "Censored: NIV use",
        "3" = "Death",
        "2" = "Following intubation strategy"
      ),
      breaks = c("0", "1", "3", "2")
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
      panel.grid.major.y = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    facet_wrap(~strategy, scales = "free_y", ncol = 2)
    # Note: Title removed for cleaner presentation

  # Store plot for later display
  all_plots[[paste0("grace_", grace_hours, "hr")]] <- cens_plot

  # Save full time range plot
  ggsave(sprintf("outputs/ROX/S.Fig_censoring_4strategies_grace%dhr.png", grace_hours),
         cens_plot,
         dpi = 1200,
         width = 14,
         height = 10)

  cat(sprintf("✓ Saved: outputs/ROX/S.Fig_censoring_4strategies_grace%dhr.png\n", grace_hours))

  # Save to manuscript folder if grace_hours == 1
  if(grace_hours == 1) {
    ggsave("C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Fig_censoring_4strategies_grace1hr.png",
           cens_plot,
           dpi = 1200,
           width = 14,
           height = 10)
    cat("✓ Also saved to manuscript folder: S.Fig_censoring_4strategies_grace1hr.png\n")
  }

  # Create 0-72 hours plot ----
  plot_combined_72h <- plot_combined %>%
    filter(time <= 72)

  # Status is already a factor from the full time range plot
  cens_plot_72h <- ggplot(plot_combined_72h,
                          aes(x = time, y = total_count, fill = status)) +
    geom_bar(position = "stack", stat = "identity", width = 0.95) +
    xlab("Time since baseline (hours)") +
    scale_x_continuous(
      breaks = c(0, 12, 24, 36, 48, 60, 72),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(xlim = c(-1, 73)) +
    ylab("Number of observations") +
    scale_fill_manual(
      values = c("0" = "#878787", "1" = "#14B5D0", "3" = "#000000", "2" = "#800000"),
      labels = c(
        "0" = "Censored: Protocol deviation (ROX index)",
        "1" = "Censored: NIV use",
        "3" = "Death",
        "2" = "Following intubation strategy"
      ),
      breaks = c("0", "1", "3", "2")
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
      panel.grid.major.y = element_blank(),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    facet_wrap(~strategy, scales = "free_y", ncol = 2)
    # Note: Title removed for cleaner presentation

  # Store 72h plot for later display
  all_plots_72h[[paste0("grace_", grace_hours, "hr_72h")]] <- cens_plot_72h

  # Save 72h plot
  ggsave(sprintf("outputs/ROX/S.Fig_censoring_4strategies_72h_grace%dhr.png", grace_hours),
         cens_plot_72h,
         dpi = 1200,
         width = 14,
         height = 10)

  cat(sprintf("✓ Saved: outputs/ROX/S.Fig_censoring_4strategies_72h_grace%dhr.png\n", grace_hours))

  # Save to manuscript folder if grace_hours == 1
  if(grace_hours == 1) {
    ggsave("C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Fig_censoring_4strategies_72h_grace1hr.png",
           cens_plot_72h,
           dpi = 1200,
           width = 14,
           height = 10)
    cat("✓ Also saved to manuscript folder: S.Fig_censoring_4strategies_72h_grace1hr.png\n")
  }
}

# Display all plots ----
cat("\n\n========== Displaying all plots ==========\n")

cat("\n=== Full time range plots ===\n")
for(i in seq_along(all_plots)) {
  cat(sprintf("Displaying: %s\n", names(all_plots)[i]))
  print(all_plots[[i]])
}

cat("\n=== 0-72 hours plots ===\n")
for(i in seq_along(all_plots_72h)) {
  cat(sprintf("Displaying: %s\n", names(all_plots_72h)[i]))
  print(all_plots_72h[[i]])
}

# Create combined comparison plots using patchwork ----
cat("\n\n=== Creating combined comparison plots ===\n")

# Add grace period labels to plots for combined display
plot_0hr <- all_plots[[1]] +
  labs(subtitle = "Grace Period: 0 hour", tag = "A") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_1hr <- all_plots[[2]] +
  labs(subtitle = "Grace Period: 1 hour", tag = "B") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_2hr <- all_plots[[3]] +
  labs(subtitle = "Grace Period: 2 hours", tag = "C") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

# Combined plot: All grace periods (full time range)
combined_full <- plot_0hr / plot_1hr / plot_2hr

print(combined_full)

ggsave("outputs/ROX/S.Fig_censoring_comparison_all_grace_periods.png",
       combined_full,
       dpi = 1200,
       width = 14,
       height = 24)

cat("✓ Saved: outputs/ROX/S.Fig_censoring_comparison_all_grace_periods.png\n")

# Add grace period labels to 72h plots
plot_0hr_72h <- all_plots_72h[[1]] +
  labs(subtitle = "Grace Period: 0 hour", tag = "A") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_1hr_72h <- all_plots_72h[[2]] +
  labs(subtitle = "Grace Period: 1 hour", tag = "B") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

plot_2hr_72h <- all_plots_72h[[3]] +
  labs(subtitle = "Grace Period: 2 hours", tag = "C") +
  theme(plot.subtitle = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.tag = element_text(size = 14, face = "bold"))

# Combined plot: All grace periods (0-72h)
combined_72h <- plot_0hr_72h / plot_1hr_72h / plot_2hr_72h

print(combined_72h)

ggsave("outputs/ROX/S.Fig_censoring_comparison_all_grace_periods_72h.png",
       combined_72h,
       dpi = 1200,
       width = 14,
       height = 24)

cat("✓ Saved: outputs/ROX/S.Fig_censoring_comparison_all_grace_periods_72h.png\n")

cat("\n\n========== All censoring plots completed! ==========\n")
cat("\nGenerated files:\n")
cat("\nIndividual grace period plots:\n")
for(grace_hours in grace_hours_values) {
  cat(sprintf("  - S.Fig_censoring_4strategies_grace%dhr.png\n", grace_hours))
  cat(sprintf("  - S.Fig_censoring_4strategies_72h_grace%dhr.png\n", grace_hours))
}
cat("\nCombined comparison plots:\n")
cat("  - S.Fig_censoring_comparison_all_grace_periods.png\n")
cat("  - S.Fig_censoring_comparison_all_grace_periods_72h.png\n")

cat("\n✓ All plots have been saved to outputs/ROX/ directory\n")
cat("✓ All plots are now displayed in the Plots panel for review\n")
