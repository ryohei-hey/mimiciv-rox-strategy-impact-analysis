#setup memory
rm(list = ls())
gc()
gc()


# load library -----
library(pacman)

p_load(tidyverse, here, data.table,boot,lubridate, skimr, naniar,gtsummary,EValue,
       patchwork, ggsci, parallel,speedglm, doParallel, foreach,
       officer,flextable)

# ===================================================
# Sensitivity Analysis: Excluding Respiratory Variables (SpO2, FiO2, Respiratory Rate)
# ===================================================
GRACE_HOURS <- 1  # Grace period setting

cat("=== Sensitivity Analysis Settings ===\n")
cat("Analysis type: Excluding respiratory variables from outcome model\n")
cat("Grace period:", GRACE_HOURS, "hour(s)\n\n")

# ===================================================
# Load bootstrap data for no_resp_spo2_fio2 sensitivity analysis
# ===================================================

# File paths for sensitivity analysis (note: slightly different naming conventions)
file_3.85 <- paste0("outputs/ROX/bootdata/bootstrap_graphs_3.85_no_resp_vars_r500_grace", GRACE_HOURS, "hr.rds")
file_4.88 <- paste0("outputs/ROX/bootdata/bootstrap_graphs_4.88_r500_grace", GRACE_HOURS, "hr_sensitivity_no_resp_spo2_fio2.rds")
file_complex <- paste0("outputs/ROX/bootdata/bootstrap_graphs_complex_r500_grace", GRACE_HOURS, "hr_sensitivity_no_resp_spo2_fio2.rds")

cat("Loading files:\n")
cat("  ", file_3.85, "\n")
cat("  ", file_4.88, "\n")
cat("  ", file_complex, "\n\n")

# Read bootstrap graph data
boot_data_3.85 <- readRDS(file_3.85)
boot_data_4.88 <- readRDS(file_4.88)
boot_data_complex <- readRDS(file_complex)

cat("Data loaded successfully\n\n")

length(boot_data_3.85)
length(boot_data_4.88)
length(boot_data_complex)

## Step 2: calc 95% CI for 3.85 -----

n_boot_3.85 <- length(boot_data_3.85)
n_time_3.85 <- nrow(boot_data_3.85[[1]])

risk0_boot_3.85 <- matrix(NA, nrow = n_boot_3.85, ncol = n_time_3.85)
risk1_boot_3.85 <- matrix(NA, nrow = n_boot_3.85, ncol = n_time_3.85)

for (i in 1:n_boot_3.85) {
  if (is.data.frame(boot_data_3.85[[i]])) {
    risk0_boot_3.85[i, ] <- boot_data_3.85[[i]]$risk0
    risk1_boot_3.85[i, ] <- boot_data_3.85[[i]]$risk1
  }
}

time_0_3.85 <- boot_data_3.85[[1]]$time_0

graph_3.85 <- data.frame(
  time_0 = time_0_3.85,
  time_days = time_0_3.85 / 24,

  risk0 = colMeans(risk0_boot_3.85, na.rm = TRUE),
  risk1 = colMeans(risk1_boot_3.85, na.rm = TRUE),

  risk0_lower = apply(risk0_boot_3.85, 2, quantile, 0.025, na.rm = TRUE),
  risk0_upper = apply(risk0_boot_3.85, 2, quantile, 0.975, na.rm = TRUE),

  risk1_lower = apply(risk1_boot_3.85, 2, quantile, 0.025, na.rm = TRUE),
  risk1_upper = apply(risk1_boot_3.85, 2, quantile, 0.975, na.rm = TRUE)
)


## Step 2: calc 95% CI for 4.88 -----

n_boot_4.88 <- length(boot_data_4.88)
n_time_4.88 <- nrow(boot_data_4.88[[1]])

risk0_boot_4.88 <- matrix(NA, nrow = n_boot_4.88, ncol = n_time_4.88)
risk1_boot_4.88 <- matrix(NA, nrow = n_boot_4.88, ncol = n_time_4.88)

for (i in 1:n_boot_4.88) {
  if (is.data.frame(boot_data_4.88[[i]])) {
    risk0_boot_4.88[i, ] <- boot_data_4.88[[i]]$risk0
    risk1_boot_4.88[i, ] <- boot_data_4.88[[i]]$risk1
  }
}

time_0_4.88 <- boot_data_4.88[[1]]$time_0

graph_4.88 <- data.frame(
  time_0 = time_0_4.88,
  time_days = time_0_4.88 / 24,

  risk0 = colMeans(risk0_boot_4.88, na.rm = TRUE),
  risk1 = colMeans(risk1_boot_4.88, na.rm = TRUE),

  risk0_lower = apply(risk0_boot_4.88, 2, quantile, 0.025, na.rm = TRUE),
  risk0_upper = apply(risk0_boot_4.88, 2, quantile, 0.975, na.rm = TRUE),

  risk1_lower = apply(risk1_boot_4.88, 2, quantile, 0.025, na.rm = TRUE),
  risk1_upper = apply(risk1_boot_4.88, 2, quantile, 0.975, na.rm = TRUE)
)


## Step 2: calc 95% CI for Complex -----

n_boot_complex <- length(boot_data_complex)
n_time_complex <- nrow(boot_data_complex[[1]])

risk0_boot_complex <- matrix(NA, nrow = n_boot_complex, ncol = n_time_complex)
risk1_boot_complex <- matrix(NA, nrow = n_boot_complex, ncol = n_time_complex)

for (i in 1:n_boot_complex) {
  if (is.data.frame(boot_data_complex[[i]])) {
    risk0_boot_complex[i, ] <- boot_data_complex[[i]]$risk0
    risk1_boot_complex[i, ] <- boot_data_complex[[i]]$risk1
  }
}

time_0_complex <- boot_data_complex[[1]]$time_0

graph_complex <- data.frame(
  time_0 = time_0_complex,
  time_days = time_0_complex / 24,

  risk0 = colMeans(risk0_boot_complex, na.rm = TRUE),
  risk1 = colMeans(risk1_boot_complex, na.rm = TRUE),

  risk0_lower = apply(risk0_boot_complex, 2, quantile, 0.025, na.rm = TRUE),
  risk0_upper = apply(risk0_boot_complex, 2, quantile, 0.975, na.rm = TRUE),

  risk1_lower = apply(risk1_boot_complex, 2, quantile, 0.025, na.rm = TRUE),
  risk1_upper = apply(risk1_boot_complex, 2, quantile, 0.975, na.rm = TRUE)
)


# Check the graph data structures
cat("Graph data for 3.85:\n")
cat("  Time points:", nrow(graph_3.85), "\n")
cat("  Max time (days):", max(graph_3.85$time_days), "\n\n")

cat("Graph data for 4.88:\n")
cat("  Time points:", nrow(graph_4.88), "\n")
cat("  Max time (days):", max(graph_4.88$time_days), "\n\n")

cat("Graph data for Complex:\n")
cat("  Time points:", nrow(graph_complex), "\n")
cat("  Max time (days):", max(graph_complex$time_days), "\n\n")


## Step 3a: Create combined graph -----

# Grace period label for legend
grace_label <- ifelse(GRACE_HOURS == 0, "", paste0(", ", GRACE_HOURS, "hr grace"))

# Color palette
col_black <- "#000000"
col_navy <- "#0c385c"
col_light_blue <- "#1770b8"
col_dark_mauve <- "#6b8ac5"

# Combine data for plotting
plot_combined <- ggplot() +

  # Usual Care - ribbon and main line
  geom_ribbon(data = graph_3.85,
              aes(x = time_days, ymin = risk0_lower, ymax = risk0_upper),
              fill = col_black, alpha = 0.15) +
  geom_line(data = graph_3.85,
            aes(x = time_days, y = risk0, color = "Usual Care", linetype = "Usual Care"),
            linewidth = 1.5) +

  # ROX Strategy 3.85 - ribbon and main line
  geom_ribbon(data = graph_3.85,
              aes(x = time_days, ymin = risk1_lower, ymax = risk1_upper),
              fill = col_navy, alpha = 0.2) +
  geom_line(data = graph_3.85,
            aes(x = time_days, y = risk1, color = "ROX < 3.85", linetype = "ROX < 3.85"),
            linewidth = 1.5) +

  # ROX Strategy 4.88 - ribbon and main line
  geom_ribbon(data = graph_4.88,
              aes(x = time_days, ymin = risk1_lower, ymax = risk1_upper),
              fill = col_light_blue, alpha = 0.2) +
  geom_line(data = graph_4.88,
            aes(x = time_days, y = risk1, color = "ROX < 4.88", linetype = "ROX < 4.88"),
            linewidth = 1.5) +

  # Complex Strategy - ribbon and main line
  geom_ribbon(data = graph_complex,
              aes(x = time_days, ymin = risk1_lower, ymax = risk1_upper),
              fill = col_dark_mauve, alpha = 0.25) +
  geom_line(data = graph_complex,
            aes(x = time_days, y = risk1, color = "Time-varying", linetype = "Time-varying"),
            linewidth = 1.5) +

  # Labels and scales
  xlab("Days since HFNC initiation") +
  ylab("Cumulative Incidence of 30-Day Mortality") +

  scale_x_continuous(limits = c(0, 30.5),
                     breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(
    limits = c(0, 0.40),
    breaks = seq(0, 0.40, by = 0.05),
    labels = scales::percent_format(accuracy = 1)
  ) +

  # Theme with Helvetica font
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.30, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +

  # Color scale
  scale_color_manual(
    values = c(
      "Usual Care" = col_black,
      "ROX < 3.85" = col_navy,
      "ROX < 4.88" = col_light_blue,
      "Time-varying" = col_dark_mauve
    ),
    breaks = c("Usual Care", "ROX < 3.85", "ROX < 4.88", "Time-varying"),
    labels = c(
      "Usual Care" = "Usual Care",
      "ROX < 3.85" = paste0("ROX < 3.85", grace_label),
      "ROX < 4.88" = paste0("ROX < 4.88", grace_label),
      "Time-varying" = paste0("Time-varying ROX", grace_label)
    )
  ) +

  # Linetype scale (all solid, legend hidden)
  scale_linetype_manual(
    values = c(
      "Usual Care" = 1,
      "ROX < 3.85" = 1,
      "ROX < 4.88" = 1,
      "Time-varying" = 1
    ),
    guide = "none"
  )

print(plot_combined)

# Save the combined plot
plot_filename <- paste0("outputs/ROX/sensitivity/Fig.cumulative_incidence_all_strategies_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png")

# Create directory if it doesn't exist
if (!dir.exists("outputs/ROX/sensitivity")) {
  dir.create("outputs/ROX/sensitivity", recursive = TRUE)
}

ggsave(plot_filename, plot_combined, width = 12, height = 7, dpi = 300)

cat("\nCombined plot saved:", plot_filename, "\n")


## Step 3b: Create combined graph without CI -----

# Combine data for plotting (without confidence interval ribbons)
plot_combined_noCI <- ggplot() +

  # Usual Care - line only
  geom_line(data = graph_3.85,
            aes(x = time_days, y = risk0, color = "Usual Care", linetype = "Usual Care"),
            linewidth = 1.5) +

  # ROX Strategy 3.85 - line only
  geom_line(data = graph_3.85,
            aes(x = time_days, y = risk1, color = "ROX < 3.85", linetype = "ROX < 3.85"),
            linewidth = 1.5) +

  # ROX Strategy 4.88 - line only
  geom_line(data = graph_4.88,
            aes(x = time_days, y = risk1, color = "ROX < 4.88", linetype = "ROX < 4.88"),
            linewidth = 1.5) +

  # Complex Strategy - line only
  geom_line(data = graph_complex,
            aes(x = time_days, y = risk1, color = "Time-varying", linetype = "Time-varying"),
            linewidth = 1.5) +

  # Labels and scales
  xlab("Days since HFNC initiation") +
  ylab("Cumulative Incidence of 30-Day Mortality") +

  scale_x_continuous(limits = c(0, 30.5),
                     breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(
    limits = c(0, 0.40),
    breaks = seq(0, 0.40, by = 0.05),
    labels = scales::percent_format(accuracy = 1)
  ) +

  # Theme with Helvetica font
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.30, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +

  # Color scale
  scale_color_manual(
    values = c(
      "Usual Care" = col_black,
      "ROX < 3.85" = col_navy,
      "ROX < 4.88" = col_light_blue,
      "Time-varying" = col_dark_mauve
    ),
    breaks = c("Usual Care", "ROX < 3.85", "ROX < 4.88", "Time-varying"),
    labels = c("Usual Care", "ROX < 3.85", "ROX < 4.88", "Time-varying ROX")
  ) +

  # Linetype scale (all solid, legend hidden)
  scale_linetype_manual(
    values = c(
      "Usual Care" = 1,
      "ROX < 3.85" = 1,
      "ROX < 4.88" = 1,
      "Time-varying" = 1
    ),
    guide = "none"
  )

print(plot_combined_noCI)

# Save the combined plot without CI
plot_filename_noCI <- paste0("outputs/ROX/sensitivity/Fig.cumulative_incidence_all_strategies_no_resp_spo2_fio2_noCI_grace", GRACE_HOURS, "hr.png")
ggsave(plot_filename_noCI, plot_combined_noCI, width = 12, height = 7, dpi = 300)

cat("\nCombined plot (no CI) saved:", plot_filename_noCI, "\n")


## Step 3c: Create individual plots for each emulated target trial -----

# Common theme for individual plots
theme_trial <- theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.35, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# Emulated Target Trial 1: ROX < 3.85 vs Usual Care
plot_trial1 <- ggplot() +
  # Usual Care - ribbon and line
  geom_ribbon(data = graph_3.85,
              aes(x = time_days, ymin = risk0_lower, ymax = risk0_upper),
              fill = col_black, alpha = 0.15) +
  geom_line(data = graph_3.85,
            aes(x = time_days, y = risk0, color = "Usual Care"),
            linewidth = 1.2) +
  # ROX < 3.85 - ribbon and line
  geom_ribbon(data = graph_3.85,
              aes(x = time_days, ymin = risk1_lower, ymax = risk1_upper),
              fill = col_navy, alpha = 0.2) +
  geom_line(data = graph_3.85,
            aes(x = time_days, y = risk1, color = "ROX < 3.85"),
            linewidth = 1.2) +
  # Scales
  xlab("Days since HFNC initiation") +
  ylab("Cumulative Incidence of 30-Day Mortality") +
  scale_x_continuous(limits = c(0, 30.5), breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits = c(0, 0.40), breaks = seq(0, 0.40, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Usual Care" = col_black, "ROX < 3.85" = col_navy)) +
  theme_trial

print(plot_trial1)

# Emulated Target Trial 2: ROX < 4.88 vs Usual Care
plot_trial2 <- ggplot() +
  # Usual Care - ribbon and line
  geom_ribbon(data = graph_4.88,
              aes(x = time_days, ymin = risk0_lower, ymax = risk0_upper),
              fill = col_black, alpha = 0.15) +
  geom_line(data = graph_4.88,
            aes(x = time_days, y = risk0, color = "Usual Care"),
            linewidth = 1.2) +
  # ROX < 4.88 - ribbon and line
  geom_ribbon(data = graph_4.88,
              aes(x = time_days, ymin = risk1_lower, ymax = risk1_upper),
              fill = col_light_blue, alpha = 0.2) +
  geom_line(data = graph_4.88,
            aes(x = time_days, y = risk1, color = "ROX < 4.88"),
            linewidth = 1.2) +
  # Scales
  xlab("Days since HFNC initiation") +
  ylab("Cumulative Incidence of 30-Day Mortality") +
  scale_x_continuous(limits = c(0, 30.5), breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits = c(0, 0.40), breaks = seq(0, 0.40, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Usual Care" = col_black, "ROX < 4.88" = col_light_blue)) +
  theme_trial

print(plot_trial2)

# Emulated Target Trial 3: Time-varying ROX vs Usual Care
plot_trial3 <- ggplot() +
  # Usual Care - ribbon and line
  geom_ribbon(data = graph_complex,
              aes(x = time_days, ymin = risk0_lower, ymax = risk0_upper),
              fill = col_black, alpha = 0.15) +
  geom_line(data = graph_complex,
            aes(x = time_days, y = risk0, color = "Usual Care"),
            linewidth = 1.2) +
  # Time-varying ROX - ribbon and line
  geom_ribbon(data = graph_complex,
              aes(x = time_days, ymin = risk1_lower, ymax = risk1_upper),
              fill = col_dark_mauve, alpha = 0.25) +
  geom_line(data = graph_complex,
            aes(x = time_days, y = risk1, color = "Time-varying ROX"),
            linewidth = 1.2) +
  # Scales
  xlab("Days since HFNC initiation") +
  ylab("Cumulative Incidence of 30-Day Mortality") +
  scale_x_continuous(limits = c(0, 30.5), breaks = seq(0, 30, by = 5)) +
  scale_y_continuous(limits = c(0, 0.40), breaks = seq(0, 0.40, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = c("Usual Care" = col_black, "Time-varying ROX" = col_dark_mauve)) +
  theme_trial

print(plot_trial3)

# Save individual plots
ggsave(paste0("outputs/ROX/sensitivity/Fig_trial1_ROX3.85_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png"),
       plot_trial1, width = 8, height = 6, dpi = 300)
ggsave(paste0("outputs/ROX/sensitivity/Fig_trial2_ROX4.88_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png"),
       plot_trial2, width = 8, height = 6, dpi = 300)
ggsave(paste0("outputs/ROX/sensitivity/Fig_trial3_TimeVarying_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png"),
       plot_trial3, width = 8, height = 6, dpi = 300)

cat("\nIndividual trial plots saved\n")

# Combine all three plots into one figure using patchwork
plot_combined_trials <- plot_trial1 + plot_trial2 + plot_trial3 +
  plot_layout(nrow = 3) +
  plot_annotation(tag_levels = 'A')

print(plot_combined_trials)

# Save combined figure
ggsave(paste0("outputs/ROX/sensitivity/Fig_all_trials_combined_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png"),
       plot_combined_trials, width = 8, height = 16, dpi = 300)

cat("\nCombined trial figure saved:", paste0("outputs/ROX/sensitivity/Fig_all_trials_combined_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png"), "\n")


## Step 4: Create comparison table for all cutoffs -----

# Helper function to calculate RD with bootstrap percentile CI
calc_rd_bootstrap <- function(time_h, risk0_boot, risk1_boot, graph_data) {
  idx <- which(graph_data$time_0 == time_h)
  rd_boot <- (risk1_boot[, idx] - risk0_boot[, idx]) * 100

  list(
    est = mean(rd_boot, na.rm = TRUE),
    lower = quantile(rd_boot, 0.025, na.rm = TRUE),
    upper = quantile(rd_boot, 0.975, na.rm = TRUE)
  )
}

# Calculate RD for each strategy and time point using bootstrap percentile method
rd_3.85_7d <- calc_rd_bootstrap(168, risk0_boot_3.85, risk1_boot_3.85, graph_3.85)
rd_3.85_14d <- calc_rd_bootstrap(336, risk0_boot_3.85, risk1_boot_3.85, graph_3.85)
rd_3.85_30d <- calc_rd_bootstrap(720, risk0_boot_3.85, risk1_boot_3.85, graph_3.85)

rd_4.88_7d <- calc_rd_bootstrap(168, risk0_boot_4.88, risk1_boot_4.88, graph_4.88)
rd_4.88_14d <- calc_rd_bootstrap(336, risk0_boot_4.88, risk1_boot_4.88, graph_4.88)
rd_4.88_30d <- calc_rd_bootstrap(720, risk0_boot_4.88, risk1_boot_4.88, graph_4.88)

rd_complex_7d <- calc_rd_bootstrap(168, risk0_boot_complex, risk1_boot_complex, graph_complex)
rd_complex_14d <- calc_rd_bootstrap(336, risk0_boot_complex, risk1_boot_complex, graph_complex)
rd_complex_30d <- calc_rd_bootstrap(720, risk0_boot_complex, risk1_boot_complex, graph_complex)

# Complete results table comparing all strategies
comparison_results <- data.frame(
  Time_Point = c("7-Day", "14-Day", "30-Day"),

  # Usual Care
  Usual_Care = sprintf("%.1f (%.1f-%.1f)",
                       c(graph_3.85$risk0[graph_3.85$time_0 == 168] * 100,
                         graph_3.85$risk0[graph_3.85$time_0 == 336] * 100,
                         graph_3.85$risk0[graph_3.85$time_0 == 720] * 100),
                       c(graph_3.85$risk0_lower[graph_3.85$time_0 == 168] * 100,
                         graph_3.85$risk0_lower[graph_3.85$time_0 == 336] * 100,
                         graph_3.85$risk0_lower[graph_3.85$time_0 == 720] * 100),
                       c(graph_3.85$risk0_upper[graph_3.85$time_0 == 168] * 100,
                         graph_3.85$risk0_upper[graph_3.85$time_0 == 336] * 100,
                         graph_3.85$risk0_upper[graph_3.85$time_0 == 720] * 100)),

  # ROX Strategy 3.85
  ROX_3.85 = sprintf("%.1f (%.1f-%.1f)",
                     c(graph_3.85$risk1[graph_3.85$time_0 == 168] * 100,
                       graph_3.85$risk1[graph_3.85$time_0 == 336] * 100,
                       graph_3.85$risk1[graph_3.85$time_0 == 720] * 100),
                     c(graph_3.85$risk1_lower[graph_3.85$time_0 == 168] * 100,
                       graph_3.85$risk1_lower[graph_3.85$time_0 == 336] * 100,
                       graph_3.85$risk1_lower[graph_3.85$time_0 == 720] * 100),
                     c(graph_3.85$risk1_upper[graph_3.85$time_0 == 168] * 100,
                       graph_3.85$risk1_upper[graph_3.85$time_0 == 336] * 100,
                       graph_3.85$risk1_upper[graph_3.85$time_0 == 720] * 100)),

  # ROX Strategy 4.88
  ROX_4.88 = sprintf("%.1f (%.1f-%.1f)",
                     c(graph_4.88$risk1[graph_4.88$time_0 == 168] * 100,
                       graph_4.88$risk1[graph_4.88$time_0 == 336] * 100,
                       graph_4.88$risk1[graph_4.88$time_0 == 720] * 100),
                     c(graph_4.88$risk1_lower[graph_4.88$time_0 == 168] * 100,
                       graph_4.88$risk1_lower[graph_4.88$time_0 == 336] * 100,
                       graph_4.88$risk1_lower[graph_4.88$time_0 == 720] * 100),
                     c(graph_4.88$risk1_upper[graph_4.88$time_0 == 168] * 100,
                       graph_4.88$risk1_upper[graph_4.88$time_0 == 336] * 100,
                       graph_4.88$risk1_upper[graph_4.88$time_0 == 720] * 100)),

  # Complex Strategy
  Complex = sprintf("%.1f (%.1f-%.1f)",
                    c(graph_complex$risk1[graph_complex$time_0 == 168] * 100,
                      graph_complex$risk1[graph_complex$time_0 == 336] * 100,
                      graph_complex$risk1[graph_complex$time_0 == 720] * 100),
                    c(graph_complex$risk1_lower[graph_complex$time_0 == 168] * 100,
                      graph_complex$risk1_lower[graph_complex$time_0 == 336] * 100,
                      graph_complex$risk1_lower[graph_complex$time_0 == 720] * 100),
                    c(graph_complex$risk1_upper[graph_complex$time_0 == 168] * 100,
                      graph_complex$risk1_upper[graph_complex$time_0 == 336] * 100,
                      graph_complex$risk1_upper[graph_complex$time_0 == 720] * 100)),

  # Risk Difference (3.85 vs Usual Care) - Bootstrap percentile method
  RD_3.85 = sprintf("%.1f (%.1f to %.1f)",
                    c(rd_3.85_7d$est, rd_3.85_14d$est, rd_3.85_30d$est),
                    c(rd_3.85_7d$lower, rd_3.85_14d$lower, rd_3.85_30d$lower),
                    c(rd_3.85_7d$upper, rd_3.85_14d$upper, rd_3.85_30d$upper)),

  # Risk Difference (4.88 vs Usual Care) - Bootstrap percentile method
  RD_4.88 = sprintf("%.1f (%.1f to %.1f)",
                    c(rd_4.88_7d$est, rd_4.88_14d$est, rd_4.88_30d$est),
                    c(rd_4.88_7d$lower, rd_4.88_14d$lower, rd_4.88_30d$lower),
                    c(rd_4.88_7d$upper, rd_4.88_14d$upper, rd_4.88_30d$upper)),

  # Risk Difference (Complex vs Usual Care) - Bootstrap percentile method
  RD_Complex = sprintf("%.1f (%.1f to %.1f)",
                       c(rd_complex_7d$est, rd_complex_14d$est, rd_complex_30d$est),
                       c(rd_complex_7d$lower, rd_complex_14d$lower, rd_complex_30d$lower),
                       c(rd_complex_7d$upper, rd_complex_14d$upper, rd_complex_30d$upper)),

  # Risk Ratio (3.85 vs Usual Care)
  RR_3.85 = sprintf("%.2f (%.2f-%.2f)",
                    sapply(c(168, 336, 720), function(h) {
                      idx <- which(graph_3.85$time_0 == h)
                      mean(risk1_boot_3.85[, idx] / risk0_boot_3.85[, idx], na.rm = TRUE)
                    }),
                    sapply(c(168, 336, 720), function(h) {
                      idx <- which(graph_3.85$time_0 == h)
                      quantile(risk1_boot_3.85[, idx] / risk0_boot_3.85[, idx], 0.025, na.rm = TRUE)
                    }),
                    sapply(c(168, 336, 720), function(h) {
                      idx <- which(graph_3.85$time_0 == h)
                      quantile(risk1_boot_3.85[, idx] / risk0_boot_3.85[, idx], 0.975, na.rm = TRUE)
                    })),

  # Risk Ratio (4.88 vs Usual Care)
  RR_4.88 = sprintf("%.2f (%.2f-%.2f)",
                    sapply(c(168, 336, 720), function(h) {
                      idx <- which(graph_4.88$time_0 == h)
                      mean(risk1_boot_4.88[, idx] / risk0_boot_4.88[, idx], na.rm = TRUE)
                    }),
                    sapply(c(168, 336, 720), function(h) {
                      idx <- which(graph_4.88$time_0 == h)
                      quantile(risk1_boot_4.88[, idx] / risk0_boot_4.88[, idx], 0.025, na.rm = TRUE)
                    }),
                    sapply(c(168, 336, 720), function(h) {
                      idx <- which(graph_4.88$time_0 == h)
                      quantile(risk1_boot_4.88[, idx] / risk0_boot_4.88[, idx], 0.975, na.rm = TRUE)
                    })),

  # Risk Ratio (Complex vs Usual Care)
  RR_Complex = sprintf("%.2f (%.2f-%.2f)",
                       sapply(c(168, 336, 720), function(h) {
                         idx <- which(graph_complex$time_0 == h)
                         mean(risk1_boot_complex[, idx] / risk0_boot_complex[, idx], na.rm = TRUE)
                       }),
                       sapply(c(168, 336, 720), function(h) {
                         idx <- which(graph_complex$time_0 == h)
                         quantile(risk1_boot_complex[, idx] / risk0_boot_complex[, idx], 0.025, na.rm = TRUE)
                       }),
                       sapply(c(168, 336, 720), function(h) {
                         idx <- which(graph_complex$time_0 == h)
                         quantile(risk1_boot_complex[, idx] / risk0_boot_complex[, idx], 0.975, na.rm = TRUE)
                       }))
)

# Set column names
colnames(comparison_results) <- c(
  "Time Point",
  "Usual Care",
  "ROX < 3.85",
  "ROX < 4.88",
  "Time-varying",
  "RD (3.85)",
  "RD (4.88)",
  "RD (Time-varying)",
  "RR (3.85)",
  "RR (4.88)",
  "RR (Time-varying)"
)

# Create flextable
ft_comparison <- flextable(comparison_results) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "center", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  add_header_row(
    values = c("", "Cumulative Mortality, % (95% CI)", "Risk Difference (95% CI)", "Risk Ratio (95% CI)"),
    colwidths = c(1, 4, 3, 3)
  ) %>%
  align(part = "header", i = 1, align = "center") %>%
  bold(part = "header") %>%
  add_footer_lines(
    values = c(
      "Sensitivity Analysis: Excluding respiratory variables (SpO2, FiO2, respiratory rate) from outcome model",
      "CI indicates confidence interval; RD, risk difference; RR, risk ratio.",
      "Time-varying strategy: ROX < 2.85 (hours 1-5), ROX < 3.47 (hours 6-11), ROX < 3.85 (hours 12+).",
      "Values are percentages with 95% confidence intervals.",
      "Risk differences are in percentage points.",
      paste0("Grace period: ", GRACE_HOURS, " hour(s)."),
      paste0("Bootstrap resamples: ", n_boot_3.85)
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

ft_comparison

# Save to Word document
doc_comparison <- read_docx()

doc_comparison <- doc_comparison %>%
  body_add_par(paste0("Table. Comparison of Cumulative Mortality Outcomes Across ROX-Guided Intubation Strategies - Sensitivity Analysis (Excluding Respiratory Variables, Grace Period: ", GRACE_HOURS, " hour)"),
               style = "heading 1") %>%
  body_add_flextable(ft_comparison)

doc_filename <- paste0("outputs/ROX/sensitivity/mortality_comparison_all_strategies_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.docx")
print(doc_comparison, target = doc_filename)

cat("\nComparison table saved:", doc_filename, "\n")


## Step 5: Table for sensitivity analysis results ------
# Effect of ROX-Guided Strategy vs Usual Care on 30-Day Risk of Mortality

# Extract 30-day results (time_0 == 720 hours = 30 days)
time_30day <- 720

# Calculate 30-day Risk Difference and Risk Ratio with 95% CI from bootstrap
# ROX 3.85
idx_30d_3.85 <- which(graph_3.85$time_0 == time_30day)
rd_boot_3.85 <- (risk1_boot_3.85[, idx_30d_3.85] - risk0_boot_3.85[, idx_30d_3.85]) * 100
rr_boot_3.85 <- risk1_boot_3.85[, idx_30d_3.85] / risk0_boot_3.85[, idx_30d_3.85]

rd_3.85_est <- mean(rd_boot_3.85, na.rm = TRUE)
rd_3.85_ci <- quantile(rd_boot_3.85, c(0.025, 0.975), na.rm = TRUE)
rr_3.85_est <- mean(rr_boot_3.85, na.rm = TRUE)
rr_3.85_ci <- quantile(rr_boot_3.85, c(0.025, 0.975), na.rm = TRUE)

# ROX 4.88
idx_30d_4.88 <- which(graph_4.88$time_0 == time_30day)
rd_boot_4.88 <- (risk1_boot_4.88[, idx_30d_4.88] - risk0_boot_4.88[, idx_30d_4.88]) * 100
rr_boot_4.88 <- risk1_boot_4.88[, idx_30d_4.88] / risk0_boot_4.88[, idx_30d_4.88]

rd_4.88_est <- mean(rd_boot_4.88, na.rm = TRUE)
rd_4.88_ci <- quantile(rd_boot_4.88, c(0.025, 0.975), na.rm = TRUE)
rr_4.88_est <- mean(rr_boot_4.88, na.rm = TRUE)
rr_4.88_ci <- quantile(rr_boot_4.88, c(0.025, 0.975), na.rm = TRUE)

# Complex (Time-varying)
idx_30d_complex <- which(graph_complex$time_0 == time_30day)
rd_boot_complex <- (risk1_boot_complex[, idx_30d_complex] - risk0_boot_complex[, idx_30d_complex]) * 100
rr_boot_complex <- risk1_boot_complex[, idx_30d_complex] / risk0_boot_complex[, idx_30d_complex]

rd_complex_est <- mean(rd_boot_complex, na.rm = TRUE)
rd_complex_ci <- quantile(rd_boot_complex, c(0.025, 0.975), na.rm = TRUE)
rr_complex_est <- mean(rr_boot_complex, na.rm = TRUE)
rr_complex_ci <- quantile(rr_boot_complex, c(0.025, 0.975), na.rm = TRUE)

# Extract absolute risks at 30 days with 95% CI
# Usual care (risk0)
risk0_3.85_30d <- graph_3.85$risk0[idx_30d_3.85] * 100
risk0_3.85_30d_ci <- c(graph_3.85$risk0_lower[idx_30d_3.85], graph_3.85$risk0_upper[idx_30d_3.85]) * 100

risk0_4.88_30d <- graph_4.88$risk0[idx_30d_4.88] * 100
risk0_4.88_30d_ci <- c(graph_4.88$risk0_lower[idx_30d_4.88], graph_4.88$risk0_upper[idx_30d_4.88]) * 100

risk0_complex_30d <- graph_complex$risk0[idx_30d_complex] * 100
risk0_complex_30d_ci <- c(graph_complex$risk0_lower[idx_30d_complex], graph_complex$risk0_upper[idx_30d_complex]) * 100

# ROX strategy (risk1)
risk1_3.85_30d <- graph_3.85$risk1[idx_30d_3.85] * 100
risk1_3.85_30d_ci <- c(graph_3.85$risk1_lower[idx_30d_3.85], graph_3.85$risk1_upper[idx_30d_3.85]) * 100

risk1_4.88_30d <- graph_4.88$risk1[idx_30d_4.88] * 100
risk1_4.88_30d_ci <- c(graph_4.88$risk1_lower[idx_30d_4.88], graph_4.88$risk1_upper[idx_30d_4.88]) * 100

risk1_complex_30d <- graph_complex$risk1[idx_30d_complex] * 100
risk1_complex_30d_ci <- c(graph_complex$risk1_lower[idx_30d_complex], graph_complex$risk1_upper[idx_30d_complex]) * 100

# Create Table data frame
table_sensitivity <- data.frame(
  Trial = c("Emulated Trial 1", "Emulated Trial 2", "Emulated Trial 3"),

  ROX_strategy = c("ROX < 3.85", "ROX < 4.88", "Time-varying threshold*"),

  # Usual care: Absolute risk (95% CI)
  Usual_care = c(
    sprintf("%.1f (%.1f to %.1f)", risk0_3.85_30d, risk0_3.85_30d_ci[1], risk0_3.85_30d_ci[2]),
    sprintf("%.1f (%.1f to %.1f)", risk0_4.88_30d, risk0_4.88_30d_ci[1], risk0_4.88_30d_ci[2]),
    sprintf("%.1f (%.1f to %.1f)", risk0_complex_30d, risk0_complex_30d_ci[1], risk0_complex_30d_ci[2])
  ),

  # ROX: Absolute risk (95% CI)
  ROX_risk = c(
    sprintf("%.1f (%.1f to %.1f)", risk1_3.85_30d, risk1_3.85_30d_ci[1], risk1_3.85_30d_ci[2]),
    sprintf("%.1f (%.1f to %.1f)", risk1_4.88_30d, risk1_4.88_30d_ci[1], risk1_4.88_30d_ci[2]),
    sprintf("%.1f (%.1f to %.1f)", risk1_complex_30d, risk1_complex_30d_ci[1], risk1_complex_30d_ci[2])
  ),

  # Risk difference, % (95% CI)
  Risk_difference = c(
    sprintf("%.1f (%.1f to %.1f)", rd_3.85_est, rd_3.85_ci[1], rd_3.85_ci[2]),
    sprintf("%.1f (%.1f to %.1f)", rd_4.88_est, rd_4.88_ci[1], rd_4.88_ci[2]),
    sprintf("%.1f (%.1f to %.1f)", rd_complex_est, rd_complex_ci[1], rd_complex_ci[2])
  ),

  # Risk ratio (95% CI)
  Risk_ratio = c(
    sprintf("%.2f (%.2f to %.2f)", rr_3.85_est, rr_3.85_ci[1], rr_3.85_ci[2]),
    sprintf("%.2f (%.2f to %.2f)", rr_4.88_est, rr_4.88_ci[1], rr_4.88_ci[2]),
    sprintf("%.2f (%.2f to %.2f)", rr_complex_est, rr_complex_ci[1], rr_complex_ci[2])
  )
)

# Set column names for the table
colnames(table_sensitivity) <- c(
  "Trial",
  "ROX strategy",
  "Usual care: Absolute risk, % (95% CI)",
  "ROX: Absolute risk, % (95% CI)",
  "Risk difference, % (95% CI)",
  "Risk ratio (95% CI)"
)

# Create flextable
ft_table_sensitivity <- flextable(table_sensitivity) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  add_footer_lines(
    values = c(
      "Sensitivity Analysis: Excluding respiratory variables (SpO2, FiO2, respiratory rate) from outcome model",
      "*Time-varying threshold: ROX < 2.85 (hours 1-5), ROX < 3.47 (hours 6-11), ROX < 3.85 (hours 12+).",
      paste0("Grace period: ", GRACE_HOURS, " hour(s). Bootstrap resamples: ", n_boot_3.85, ".")
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

ft_table_sensitivity

# Save Table to Word document
path1 <- paste0("outputs/ROX/sensitivity/S.Tab_ROX_strategy_effect_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.docx")

ft_table_sensitivity |> autofit() |>
  set_table_properties(layout = "autofit") |>
  save_as_docx(path = path1)

cat("\nSensitivity analysis table saved:", path1, "\n")


## Step 6: Supplementary Table - Multiple time points -----
# Table with 7-day, 14-day, and 30-day mortality for each strategy

# Time points
time_points <- c(168, 336, 720)  # 7, 14, 30 days in hours
time_labels <- c("Day 7", "Day 14", "Day 30*")

# Function to extract results for a given time point
extract_results <- function(time_h, graph_data, risk0_boot, risk1_boot) {
  idx <- which(graph_data$time_0 == time_h)

  # Absolute risks
  risk0 <- graph_data$risk0[idx] * 100
  risk0_ci <- c(graph_data$risk0_lower[idx], graph_data$risk0_upper[idx]) * 100
  risk1 <- graph_data$risk1[idx] * 100
  risk1_ci <- c(graph_data$risk1_lower[idx], graph_data$risk1_upper[idx]) * 100

  # Risk difference from bootstrap
  rd_boot <- (risk1_boot[, idx] - risk0_boot[, idx]) * 100
  rd_est <- mean(rd_boot, na.rm = TRUE)
  rd_ci <- quantile(rd_boot, c(0.025, 0.975), na.rm = TRUE)

  # Risk ratio from bootstrap
  rr_boot <- risk1_boot[, idx] / risk0_boot[, idx]
  rr_est <- mean(rr_boot, na.rm = TRUE)
  rr_ci <- quantile(rr_boot, c(0.025, 0.975), na.rm = TRUE)

  list(
    risk0 = sprintf("%.1f (%.1f to %.1f)", risk0, risk0_ci[1], risk0_ci[2]),
    risk1 = sprintf("%.1f (%.1f to %.1f)", risk1, risk1_ci[1], risk1_ci[2]),
    rd = sprintf("%.1f (%.1f to %.1f)", rd_est, rd_ci[1], rd_ci[2]),
    rr = sprintf("%.2f (%.2f to %.2f)", rr_est, rr_ci[1], rr_ci[2])
  )
}

# Extract results for all strategies and time points
# ROX 3.85
results_3.85 <- lapply(time_points, function(t) {
  extract_results(t, graph_3.85, risk0_boot_3.85, risk1_boot_3.85)
})

# ROX 4.88
results_4.88 <- lapply(time_points, function(t) {
  extract_results(t, graph_4.88, risk0_boot_4.88, risk1_boot_4.88)
})

# Complex
results_complex <- lapply(time_points, function(t) {
  extract_results(t, graph_complex, risk0_boot_complex, risk1_boot_complex)
})

# Create table data frame
table_step6 <- data.frame(
  Trial = c(
    "Emulated Target Trial 1", "", "",
    "Emulated Target Trial 2", "", "",
    "Emulated Target Trial 3", "", ""
  ),

  ROX_strategy = c(
    "ROX < 3.85", "", "",
    "ROX < 4.88", "", "",
    "Time-varying ROX strategy", "", ""
  ),

  Outcome = rep(time_labels, 3),

  Usual_care = c(
    results_3.85[[1]]$risk0, results_3.85[[2]]$risk0, results_3.85[[3]]$risk0,
    results_4.88[[1]]$risk0, results_4.88[[2]]$risk0, results_4.88[[3]]$risk0,
    results_complex[[1]]$risk0, results_complex[[2]]$risk0, results_complex[[3]]$risk0
  ),

  ROX_risk = c(
    results_3.85[[1]]$risk1, results_3.85[[2]]$risk1, results_3.85[[3]]$risk1,
    results_4.88[[1]]$risk1, results_4.88[[2]]$risk1, results_4.88[[3]]$risk1,
    results_complex[[1]]$risk1, results_complex[[2]]$risk1, results_complex[[3]]$risk1
  ),

  Risk_difference = c(
    results_3.85[[1]]$rd, results_3.85[[2]]$rd, results_3.85[[3]]$rd,
    results_4.88[[1]]$rd, results_4.88[[2]]$rd, results_4.88[[3]]$rd,
    results_complex[[1]]$rd, results_complex[[2]]$rd, results_complex[[3]]$rd
  ),

  Risk_ratio = c(
    results_3.85[[1]]$rr, results_3.85[[2]]$rr, results_3.85[[3]]$rr,
    results_4.88[[1]]$rr, results_4.88[[2]]$rr, results_4.88[[3]]$rr,
    results_complex[[1]]$rr, results_complex[[2]]$rr, results_complex[[3]]$rr
  )
)

# Set column names
colnames(table_step6) <- c(
  "Trial",
  "ROX strategy",
  "Time point",
  "Usual care, % (95% CI)",
  "ROX strategy, % (95% CI)",
  "Risk difference, percentage points (95% CI)",
  "Risk ratio (95% CI)"
)

# Create flextable
ft_step6 <- flextable(table_step6) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1:3, align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  # Merge cells for Trial column
  merge_at(i = 1:3, j = 1) %>%
  merge_at(i = 4:6, j = 1) %>%
  merge_at(i = 7:9, j = 1) %>%
  # Merge cells for ROX strategy column
  merge_at(i = 1:3, j = 2) %>%
  merge_at(i = 4:6, j = 2) %>%
  merge_at(i = 7:9, j = 2) %>%
  # Vertical alignment for merged cells

  valign(j = 1:2, valign = "center", part = "body") %>%
  # Add horizontal lines between trials
  hline(i = c(3, 6), border = fp_border(color = "gray70", width = 1)) %>%
  add_footer_lines(
    values = c(
      "Sensitivity Analysis: Excluding respiratory variables (SpO2, FiO2, respiratory rate) from outcome model",
      "*Primary outcome time point.",
      "Time-varying ROX strategy: thresholds were 2.85 (hours 1-5), 3.47 (hours 6-11), and 3.85 (>=12 hours).",
      "Estimates were obtained from models adjusted for prespecified baseline covariates and hourly time-varying covariates (excluding respiratory variables); 95% confidence intervals were obtained using 500 bootstrap replicates.")
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")


ft_step6

# Save Step 6 table to Word document
path1_step6 <- paste0("outputs/ROX/sensitivity/table_multiple_timepoints_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.docx")

ft_step6 %>%
  autofit() %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = path1_step6)

cat("\nMultiple timepoints table saved:", path1_step6, "\n")


path2_step6 <- paste0("C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Table 8_ROX_strategy_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.docx")

ft_step6 |> autofit() |>
  set_table_properties(layout = "autofit") |>
  save_as_docx(path = path2_step6)

cat("\nSensitivity analysis table saved:", path1, "\n")



## Step 7: Forest Plot -----
# Create forest plot for manuscript

# Extract numeric values for forest plot
# Function to extract RR estimates and CIs from bootstrap
extract_rr_for_plot <- function(time_h, graph_data, risk0_boot, risk1_boot) {
  idx <- which(graph_data$time_0 == time_h)
  rr_boot <- risk1_boot[, idx] / risk0_boot[, idx]

  list(
    estimate = mean(rr_boot, na.rm = TRUE),
    lower = quantile(rr_boot, 0.025, na.rm = TRUE),
    upper = quantile(rr_boot, 0.975, na.rm = TRUE)
  )
}

# Extract RR for all strategies and time points
forest_data <- data.frame(
  Strategy = rep(c("ROX < 3.85 vs Usual Care", "ROX < 4.88 vs Usual Care", "Time-varying vs Usual Care"), each = 3),
  Time_point = rep(c("Day 7", "Day 14", "Day 30"), 3),
  RR = NA,
  RR_lower = NA,
  RR_upper = NA
)

# Fill in ROX 3.85
for (i in 1:3) {
  res <- extract_rr_for_plot(time_points[i], graph_3.85, risk0_boot_3.85, risk1_boot_3.85)
  forest_data$RR[i] <- res$estimate
  forest_data$RR_lower[i] <- res$lower
  forest_data$RR_upper[i] <- res$upper
}

# Fill in ROX 4.88
for (i in 1:3) {
  res <- extract_rr_for_plot(time_points[i], graph_4.88, risk0_boot_4.88, risk1_boot_4.88)
  forest_data$RR[i + 3] <- res$estimate
  forest_data$RR_lower[i + 3] <- res$lower
  forest_data$RR_upper[i + 3] <- res$upper
}

# Fill in Complex
for (i in 1:3) {
  res <- extract_rr_for_plot(time_points[i], graph_complex, risk0_boot_complex, risk1_boot_complex)
  forest_data$RR[i + 6] <- res$estimate
  forest_data$RR_lower[i + 6] <- res$lower
  forest_data$RR_upper[i + 6] <- res$upper
}

# Order: Strategy groups from top to bottom (ROX 3.85 at top), within each: Day 7, 14, 30
forest_data$Strategy <- factor(forest_data$Strategy,
                                levels = c("ROX < 3.85 vs Usual Care", "ROX < 4.88 vs Usual Care", "Time-varying vs Usual Care"))
forest_data$Time_point <- factor(forest_data$Time_point,
                                  levels = c("Day 7", "Day 14", "Day 30"))

# Create y position (1 = bottom, 9 = top)
forest_data <- forest_data %>%
  arrange(desc(Strategy), desc(Time_point)) %>%
  mutate(y_pos = row_number())

# Format RR text for display
forest_data$RR_text <- sprintf("%.2f (%.2f, %.2f)",
                                forest_data$RR,
                                forest_data$RR_lower,
                                forest_data$RR_upper)

# For Strategy label: show only for middle row of each group
# y_pos 7-9: ROX < 3.85 (top), y_pos 4-6: ROX < 4.88 (middle), y_pos 1-3: Time-varying (bottom)
forest_data$Strategy_label <- ""
forest_data$Strategy_label[forest_data$y_pos == 8] <- "ROX < 3.85 vs Usual Care"
forest_data$Strategy_label[forest_data$y_pos == 5] <- "ROX < 4.88 vs Usual Care"
forest_data$Strategy_label[forest_data$y_pos == 2] <- "Time-varying vs Usual Care"

# Define x-axis range for forest plot area
x_min <- 0.2
x_max <- 1.5
text_strategy_x <- -0.65
text_day_x <- -0.25
text_rr_x <- -0.05

# Create forest plot with left-side labels
forest_plot <- ggplot(forest_data, aes(y = y_pos)) +

  # Reference line at RR = 1
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.5) +

  # Confidence intervals (color by strategy)
  geom_errorbarh(aes(xmin = RR_lower, xmax = RR_upper, x = RR, color = Strategy),
                 height = 0.25, linewidth = 0.8) +

  # Point estimates (diamond shape, color by strategy)
  geom_point(aes(x = RR, color = Strategy), shape = 18, size = 4) +

  # Color scale
  scale_color_manual(
    values = c("ROX < 3.85 vs Usual Care" = col_navy,
               "ROX < 4.88 vs Usual Care" = col_light_blue,
               "Time-varying vs Usual Care" = col_dark_mauve),
    guide = "none"
  ) +

  # Strategy labels (left)
  geom_text(aes(x = text_strategy_x, label = Strategy_label),
            hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +

  # Time point labels
  geom_text(aes(x = text_day_x, label = as.character(Time_point)),
            hjust = 0, size = 3.2, family = "Helvetica") +

  # RR estimate text
  geom_text(aes(x = text_rr_x, label = RR_text),
            hjust = 0, size = 3, family = "Helvetica") +

  # Horizontal separators between strategies
  geom_hline(yintercept = c(3.5, 6.5), linetype = "solid", color = "gray70", linewidth = 0.4) +

  # Scale

  scale_x_continuous(
    limits = c(text_strategy_x - 0.1, x_max),
    breaks = seq(0.4, 1.4, by = 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.3, 9.7),
    expand = c(0, 0)
  ) +

  # Labels
  labs(
    x = "Risk Ratio",
    y = NULL
  ) +

  # Add column headers
  annotate("text", x = text_strategy_x, y = 9.5, label = "Strategy",
           hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +
  annotate("text", x = text_day_x, y = 9.5, label = "Time point",
           hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +
  annotate("text", x = text_rr_x, y = 9.5, label = "RR (95% CI)",
           hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +

  # Theme with Helvetica font
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 10)),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(15, 15, 10, 15),
    axis.line.x = element_line(color = "black", linewidth = 0.5)
  ) +

  # Add vertical line to separate text from plot
  geom_vline(xintercept = x_min, linetype = "solid", color = "gray50", linewidth = 0.3) +

  # Clip off to show only forest plot area x-axis
  coord_cartesian(clip = "off")

print(forest_plot)

# Save forest plot
forest_path1 <- paste0("outputs/ROX/sensitivity/Fig_forest_plot_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png")

ggsave(forest_path1, forest_plot, width = 10, height = 6, dpi = 300)
cat("\nForest plot saved:", forest_path1, "\n")

## Step 8: Forest Plot for Risk Difference -----
# Create forest plot for Risk Difference

# Function to extract RD estimates and CIs from bootstrap
extract_rd_for_plot <- function(time_h, graph_data, risk0_boot, risk1_boot) {
  idx <- which(graph_data$time_0 == time_h)
  rd_boot <- (risk1_boot[, idx] - risk0_boot[, idx]) * 100  # Convert to percentage points

  list(
    estimate = mean(rd_boot, na.rm = TRUE),
    lower = quantile(rd_boot, 0.025, na.rm = TRUE),
    upper = quantile(rd_boot, 0.975, na.rm = TRUE)
  )
}

# Extract RD for all strategies and time points
forest_data_rd <- data.frame(
  Strategy = rep(c("ROX < 3.85 vs Usual Care", "ROX < 4.88 vs Usual Care", "Time-varying vs Usual Care"), each = 3),
  Time_point = rep(c("Day 7", "Day 14", "Day 30"), 3),
  RD = NA,
  RD_lower = NA,
  RD_upper = NA
)

# Fill in ROX 3.85
for (i in 1:3) {
  res <- extract_rd_for_plot(time_points[i], graph_3.85, risk0_boot_3.85, risk1_boot_3.85)
  forest_data_rd$RD[i] <- res$estimate
  forest_data_rd$RD_lower[i] <- res$lower
  forest_data_rd$RD_upper[i] <- res$upper
}

# Fill in ROX 4.88
for (i in 1:3) {
  res <- extract_rd_for_plot(time_points[i], graph_4.88, risk0_boot_4.88, risk1_boot_4.88)
  forest_data_rd$RD[i + 3] <- res$estimate
  forest_data_rd$RD_lower[i + 3] <- res$lower
  forest_data_rd$RD_upper[i + 3] <- res$upper
}

# Fill in Complex
for (i in 1:3) {
  res <- extract_rd_for_plot(time_points[i], graph_complex, risk0_boot_complex, risk1_boot_complex)
  forest_data_rd$RD[i + 6] <- res$estimate
  forest_data_rd$RD_lower[i + 6] <- res$lower
  forest_data_rd$RD_upper[i + 6] <- res$upper
}

# Order: Strategy groups from top to bottom (ROX 3.85 at top), within each: Day 7, 14, 30
forest_data_rd$Strategy <- factor(forest_data_rd$Strategy,
                                   levels = c("ROX < 3.85 vs Usual Care", "ROX < 4.88 vs Usual Care", "Time-varying vs Usual Care"))
forest_data_rd$Time_point <- factor(forest_data_rd$Time_point,
                                     levels = c("Day 7", "Day 14", "Day 30"))

# Create y position (1 = bottom, 9 = top)
forest_data_rd <- forest_data_rd %>%
  arrange(desc(Strategy), desc(Time_point)) %>%
  mutate(y_pos = row_number())

# Format RD text for display
forest_data_rd$RD_text <- sprintf("%.1f (%.1f, %.1f)",
                                   forest_data_rd$RD,
                                   forest_data_rd$RD_lower,
                                   forest_data_rd$RD_upper)

# For Strategy label: show only for middle row of each group
forest_data_rd$Strategy_label <- ""
forest_data_rd$Strategy_label[forest_data_rd$y_pos == 8] <- "ROX < 3.85 vs Usual Care"
forest_data_rd$Strategy_label[forest_data_rd$y_pos == 5] <- "ROX < 4.88 vs Usual Care"
forest_data_rd$Strategy_label[forest_data_rd$y_pos == 2] <- "Time-varying vs Usual Care"

# Define x-axis range for forest plot area
x_min_rd <- -20
x_max_rd <- 10
text_strategy_x_rd <- -55
text_day_x_rd <- -38
text_rd_x_rd <- -28

# Create forest plot for RD
forest_plot_rd <- ggplot(forest_data_rd, aes(y = y_pos)) +

  # Reference line at RD = 0
geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +

  # Confidence intervals (color by strategy)
  geom_errorbarh(aes(xmin = RD_lower, xmax = RD_upper, x = RD, color = Strategy),
                 height = 0.25, linewidth = 0.8) +

  # Point estimates (diamond shape, color by strategy)
  geom_point(aes(x = RD, color = Strategy), shape = 18, size = 4) +

  # Color scale
  scale_color_manual(
    values = c("ROX < 3.85 vs Usual Care" = col_navy,
               "ROX < 4.88 vs Usual Care" = col_light_blue,
               "Time-varying vs Usual Care" = col_dark_mauve),
    guide = "none"
  ) +

  # Strategy labels (left)
  geom_text(aes(x = text_strategy_x_rd, label = Strategy_label),
            hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +

  # Time point labels
  geom_text(aes(x = text_day_x_rd, label = as.character(Time_point)),
            hjust = 0, size = 3.2, family = "Helvetica") +

  # RD estimate text
  geom_text(aes(x = text_rd_x_rd, label = RD_text),
            hjust = 0, size = 3, family = "Helvetica") +

  # Horizontal separators between strategies
  geom_hline(yintercept = c(3.5, 6.5), linetype = "solid", color = "gray70", linewidth = 0.4) +

  # Scale
  scale_x_continuous(
    limits = c(text_strategy_x_rd - 3, x_max_rd),
    breaks = seq(-20, 10, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0.3, 9.7),
    expand = c(0, 0)
  ) +

  # Labels
  labs(
    x = "Risk Difference (percentage points)",
    y = NULL
  ) +

  # Add column headers
  annotate("text", x = text_strategy_x_rd, y = 9.5, label = "Strategy",
           hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +
  annotate("text", x = text_day_x_rd, y = 9.5, label = "Time point",
           hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +
  annotate("text", x = text_rd_x_rd, y = 9.5, label = "RD (95% CI)",
           hjust = 0, size = 3.5, fontface = "bold", family = "Helvetica") +

  # Theme with Helvetica font
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 11, margin = margin(t = 10)),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(15, 15, 10, 15),
    axis.line.x = element_line(color = "black", linewidth = 0.5)
  ) +

  # Add vertical line to separate text from plot
  geom_vline(xintercept = x_min_rd, linetype = "solid", color = "gray50", linewidth = 0.3) +

  # Clip off to show only forest plot area x-axis
  coord_cartesian(clip = "off")

print(forest_plot_rd)

# Save forest plot for RD
forest_rd_path1 <- paste0("outputs/ROX/sensitivity/Fig_forest_plot_RD_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.png")

ggsave(forest_rd_path1, forest_plot_rd, width = 10, height = 6, dpi = 300)
cat("\nForest plot (RD) saved:", forest_rd_path1, "\n")

## Step 9: E-value Sensitivity Analysis for Primary Outcome (30-Day Mortality) -----
# E-value quantifies the minimum strength of association that an unmeasured confounder
# would need to have with both the treatment and outcome to fully explain away the observed effect

# Load EValue package
library(EValue)

# Calculate E-values for 30-day mortality (primary outcome)
cat("\n=== E-value Sensitivity Analysis for Primary Outcome (30-Day Mortality) ===\n\n")

# ROX < 3.85 vs Usual Care
rr_3.85_point <- rr_3.85_est
rr_3.85_lower_ci <- rr_3.85_ci[1]
rr_3.85_upper_ci <- rr_3.85_ci[2]

evalue_3.85 <- evalues.RR(est = rr_3.85_point,
                           lo = rr_3.85_lower_ci,
                           hi = rr_3.85_upper_ci,
                           true = 1)
evalue_3.85_point <- evalue_3.85["E-values", "point"]
evalue_3.85_ci_lower <- evalue_3.85["E-values", "lower"]
evalue_3.85_ci_upper <- evalue_3.85["E-values", "upper"]

cat("ROX < 3.85 vs Usual Care:\n")
cat(sprintf("  Risk Ratio: %.2f (95%% CI: %.2f to %.2f)\n", rr_3.85_point, rr_3.85_lower_ci, rr_3.85_upper_ci))
cat(sprintf("  E-value for point estimate: %.2f\n", evalue_3.85_point))
cat(sprintf("  E-value for lower CI: %s\n", ifelse(is.na(evalue_3.85_ci_lower), "NA", sprintf("%.2f", evalue_3.85_ci_lower))))
cat(sprintf("  E-value for upper CI: %s\n\n", ifelse(is.na(evalue_3.85_ci_upper), "NA", sprintf("%.2f", evalue_3.85_ci_upper))))
print(evalue_3.85)

# ROX < 4.88 vs Usual Care
rr_4.88_point <- rr_4.88_est
rr_4.88_lower_ci <- rr_4.88_ci[1]
rr_4.88_upper_ci <- rr_4.88_ci[2]

evalue_4.88 <- evalues.RR(est = rr_4.88_point,
                           lo = rr_4.88_lower_ci,
                           hi = rr_4.88_upper_ci,
                           true = 1)
evalue_4.88_point <- evalue_4.88["E-values", "point"]
evalue_4.88_ci_lower <- evalue_4.88["E-values", "lower"]
evalue_4.88_ci_upper <- evalue_4.88["E-values", "upper"]

cat("\nROX < 4.88 vs Usual Care:\n")
cat(sprintf("  Risk Ratio: %.2f (95%% CI: %.2f to %.2f)\n", rr_4.88_point, rr_4.88_lower_ci, rr_4.88_upper_ci))
cat(sprintf("  E-value for point estimate: %.2f\n", evalue_4.88_point))
cat(sprintf("  E-value for lower CI: %s\n", ifelse(is.na(evalue_4.88_ci_lower), "NA", sprintf("%.2f", evalue_4.88_ci_lower))))
cat(sprintf("  E-value for upper CI: %s\n\n", ifelse(is.na(evalue_4.88_ci_upper), "NA", sprintf("%.2f", evalue_4.88_ci_upper))))
print(evalue_4.88)

# Time-varying ROX vs Usual Care
rr_complex_point <- rr_complex_est
rr_complex_lower_ci <- rr_complex_ci[1]
rr_complex_upper_ci <- rr_complex_ci[2]

evalue_complex <- evalues.RR(est = rr_complex_point,
                              lo = rr_complex_lower_ci,
                              hi = rr_complex_upper_ci,
                              true = 1)
evalue_complex_point <- evalue_complex["E-values", "point"]
evalue_complex_ci_lower <- evalue_complex["E-values", "lower"]
evalue_complex_ci_upper <- evalue_complex["E-values", "upper"]

cat("\nTime-varying ROX vs Usual Care:\n")
cat(sprintf("  Risk Ratio: %.2f (95%% CI: %.2f to %.2f)\n", rr_complex_point, rr_complex_lower_ci, rr_complex_upper_ci))
cat(sprintf("  E-value for point estimate: %.2f\n", evalue_complex_point))
cat(sprintf("  E-value for lower CI: %s\n", ifelse(is.na(evalue_complex_ci_lower), "NA", sprintf("%.2f", evalue_complex_ci_lower))))
cat(sprintf("  E-value for upper CI: %s\n\n", ifelse(is.na(evalue_complex_ci_upper), "NA", sprintf("%.2f", evalue_complex_ci_upper))))
print(evalue_complex)

# Create E-value summary table
# Helper function to format E-value (handle NA)
format_evalue <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.2f", x))
}

evalue_table <- data.frame(
  Strategy = c("ROX < 3.85", "ROX < 4.88", "Time-varying ROX"),

  RR = c(
    sprintf("%.2f (%.2f to %.2f)", rr_3.85_point, rr_3.85_lower_ci, rr_3.85_upper_ci),
    sprintf("%.2f (%.2f to %.2f)", rr_4.88_point, rr_4.88_lower_ci, rr_4.88_upper_ci),
    sprintf("%.2f (%.2f to %.2f)", rr_complex_point, rr_complex_lower_ci, rr_complex_upper_ci)
  ),

  Evalue_point = c(
    format_evalue(evalue_3.85_point),
    format_evalue(evalue_4.88_point),
    format_evalue(evalue_complex_point)
  ),

  Evalue_CI_lower = c(
    format_evalue(evalue_3.85_ci_lower),
    format_evalue(evalue_4.88_ci_lower),
    format_evalue(evalue_complex_ci_lower)
  ),

  Evalue_CI_upper = c(
    format_evalue(evalue_3.85_ci_upper),
    format_evalue(evalue_4.88_ci_upper),
    format_evalue(evalue_complex_ci_upper)
  )
)

colnames(evalue_table) <- c(
  "ROX Strategy",
  "Risk Ratio (95% CI)",
  "E-value (point estimate)",
  "E-value (lower CI)",
  "E-value (upper CI)"
)

# Create flextable for E-values
ft_evalue <- flextable(evalue_table) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  add_footer_lines(
    values = c(
      "Sensitivity Analysis: Excluding respiratory variables (SpO2, FiO2, respiratory rate) from outcome model",
      "E-value: The minimum strength of association on the risk ratio scale that an unmeasured confounder would need to have with both the treatment and the outcome to fully explain away the observed treatment-outcome association.",
      "For protective effects (RR < 1), the E-value is calculated using the inverse of the RR.",
      "NA indicates that the CI bound is on the opposite side of the null (RR = 1) from the point estimate, so no E-value is computed.",
      "For RR < 1, the upper CI bound is closer to the null value of 1.",
      paste0("Grace period: ", GRACE_HOURS, " hour(s).")
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

ft_evalue

# Save E-value table
evalue_path <- paste0("outputs/ROX/sensitivity/S.Tab_Evalue_sensitivity_no_resp_spo2_fio2_grace", GRACE_HOURS, "hr.docx")

ft_evalue %>%
  autofit() %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = evalue_path)

cat("\nE-value sensitivity analysis table saved:", evalue_path, "\n")


## Step 10: Summary -----
cat("\n=== Analysis Complete ===\n")
cat("Sensitivity Analysis: Excluding respiratory variables (SpO2, FiO2, respiratory rate)\n")
cat("Grace period:", GRACE_HOURS, "hour(s)\n")
cat("Bootstrap resamples:", n_boot_3.85, "\n\n")

cat("Output files saved to: outputs/ROX/sensitivity/\n")
cat("- Cumulative incidence plots\n")
cat("- Individual trial plots\n")
cat("- Comparison tables\n")
cat("- Forest plots (RR and RD)\n")
cat("- E-value sensitivity analysis\n")
