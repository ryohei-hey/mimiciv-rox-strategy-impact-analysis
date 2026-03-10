# ============================================================================
# Summary Analysis with Rubin's Rule Pooling for Multiple Imputation (m=5)
# ============================================================================
# Purpose: Pool bootstrap results from m=5 multiply imputed datasets
# Method: Rubin's Rule for combining estimates and variances
# Grace period: 1 hour
# Strategies: ROX < 3.85, ROX < 4.88, Time-varying ROX
# ============================================================================

# Setup memory
rm(list = ls())
gc()
gc()

# Load libraries -----
library(pacman)
p_load(tidyverse, here, data.table, boot, lubridate, skimr, naniar, gtsummary, EValue,
       patchwork, ggsci, parallel, speedglm, doParallel, foreach,
       officer, flextable)

# Disable printing results in scientific notation
options(scipen = 999)

# Parameters -----
GRACE_HOURS <- 1
N_IMPUTATIONS <- 5
N_BOOTSTRAP <- 200

cat("\n============================================================================\n")
cat("  Multiple Imputation Summary Analysis with Rubin's Rule (m=5)\n")
cat("============================================================================\n")
cat(sprintf("Grace Period: %d hour(s)\n", GRACE_HOURS))
cat(sprintf("Number of Imputations: %d\n", N_IMPUTATIONS))
cat(sprintf("Bootstrap per Imputation: %d\n", N_BOOTSTRAP))
cat(sprintf("Total Bootstrap Samples: %d\n", N_IMPUTATIONS * N_BOOTSTRAP))
cat("============================================================================\n\n")

# ============================================================================
# Function: Apply Rubin's Rule for pooling
# ============================================================================
# Pools estimates and variances from multiple imputations
# theta_m: vector of point estimates from each imputation
# var_m: vector of within-imputation variances from each imputation
# Returns: list with pooled estimate, total variance, and 95% CI

rubins_rule <- function(theta_m, var_m, m = 5, alpha = 0.05) {
  # Pooled estimate (average across imputations)
  theta_pooled <- mean(theta_m, na.rm = TRUE)

  # Within-imputation variance (average of variances)
  W <- mean(var_m, na.rm = TRUE)

  # Between-imputation variance
  B <- var(theta_m, na.rm = TRUE)

  # Total variance (Rubin's Rule formula)
  T_var <- W + (1 + 1/m) * B

  # Standard error
  SE <- sqrt(T_var)

  # Degrees of freedom (Barnard-Rubin adjustment)
  df <- (m - 1) * (1 + W / ((1 + 1/m) * B))^2

  # Critical value from t-distribution
  t_crit <- qt(1 - alpha/2, df)

  # 95% Confidence interval
  ci_lower <- theta_pooled - t_crit * SE
  ci_upper <- theta_pooled + t_crit * SE

  return(list(
    estimate = theta_pooled,
    variance = T_var,
    se = SE,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    df = df
  ))
}

# ============================================================================
# Load Bootstrap Data from All Imputations
# ============================================================================

cat("=== Loading Bootstrap Data ===\n")

# Initialize lists to store bootstrap data from all imputations
boot_data_3.85_all <- list()
boot_data_4.88_all <- list()
boot_data_complex_all <- list()

# Load data from each imputation
for (imp in 1:N_IMPUTATIONS) {
  cat(sprintf("Loading imputation %d/%d...\n", imp, N_IMPUTATIONS))

  boot_data_3.85_all[[imp]] <- readRDS(sprintf("outputs/ROX/bootdata/bootstrap_results_3.85_imp%d.rds", imp))
  boot_data_4.88_all[[imp]] <- readRDS(sprintf("outputs/ROX/bootdata/bootstrap_results_4.88_imp%d.rds", imp))
  boot_data_complex_all[[imp]] <- readRDS(sprintf("outputs/ROX/bootdata/bootstrap_results_complex_imp%d.rds", imp))
}

cat("\nAll data loaded successfully\n")
cat(sprintf("Total datasets loaded: %d per strategy\n\n", N_IMPUTATIONS))

# ============================================================================
# Function: Extract Risk Matrices from Bootstrap Results
# ============================================================================

extract_risk_matrices <- function(boot_data) {
  n_boot <- length(boot_data)
  n_time <- nrow(boot_data[[1]])

  risk0_boot <- matrix(NA, nrow = n_boot, ncol = n_time)
  risk1_boot <- matrix(NA, nrow = n_boot, ncol = n_time)

  for (i in 1:n_boot) {
    if (is.data.frame(boot_data[[i]])) {
      risk0_boot[i, ] <- boot_data[[i]]$risk0
      risk1_boot[i, ] <- boot_data[[i]]$risk1
    }
  }

  list(risk0 = risk0_boot, risk1 = risk1_boot, time_0 = boot_data[[1]]$time_0)
}

# ============================================================================
# Apply Rubin's Rule to Pool Results Across Imputations
# ============================================================================

cat("=== Applying Rubin's Rule to Pool Results ===\n")

# Function to pool results for a specific strategy
pool_strategy_results <- function(boot_data_all, strategy_name) {
  cat(sprintf("\nPooling results for %s strategy...\n", strategy_name))

  # Extract risk matrices for each imputation
  risk_matrices_list <- lapply(boot_data_all, extract_risk_matrices)

  # Get dimensions
  n_time <- length(risk_matrices_list[[1]]$time_0)
  time_0 <- risk_matrices_list[[1]]$time_0

  # Initialize arrays to store estimates and variances from each imputation
  risk0_est_imp <- matrix(NA, nrow = N_IMPUTATIONS, ncol = n_time)
  risk1_est_imp <- matrix(NA, nrow = N_IMPUTATIONS, ncol = n_time)
  risk0_var_imp <- matrix(NA, nrow = N_IMPUTATIONS, ncol = n_time)
  risk1_var_imp <- matrix(NA, nrow = N_IMPUTATIONS, ncol = n_time)

  # Calculate estimates and variances for each imputation
  for (imp in 1:N_IMPUTATIONS) {
    # Point estimates (mean across bootstrap samples within imputation)
    risk0_est_imp[imp, ] <- colMeans(risk_matrices_list[[imp]]$risk0, na.rm = TRUE)
    risk1_est_imp[imp, ] <- colMeans(risk_matrices_list[[imp]]$risk1, na.rm = TRUE)

    # Within-imputation variances (variance across bootstrap samples)
    risk0_var_imp[imp, ] <- apply(risk_matrices_list[[imp]]$risk0, 2, var, na.rm = TRUE)
    risk1_var_imp[imp, ] <- apply(risk_matrices_list[[imp]]$risk1, 2, var, na.rm = TRUE)
  }

  # Apply Rubin's Rule to pool across imputations for each time point
  risk0_pooled <- rep(NA, n_time)
  risk1_pooled <- rep(NA, n_time)
  risk0_lower <- rep(NA, n_time)
  risk0_upper <- rep(NA, n_time)
  risk1_lower <- rep(NA, n_time)
  risk1_upper <- rep(NA, n_time)

  for (t in 1:n_time) {
    # Pool risk0
    pooled_r0 <- rubins_rule(risk0_est_imp[, t], risk0_var_imp[, t], m = N_IMPUTATIONS)
    risk0_pooled[t] <- pooled_r0$estimate
    risk0_lower[t] <- pooled_r0$ci_lower
    risk0_upper[t] <- pooled_r0$ci_upper

    # Pool risk1
    pooled_r1 <- rubins_rule(risk1_est_imp[, t], risk1_var_imp[, t], m = N_IMPUTATIONS)
    risk1_pooled[t] <- pooled_r1$estimate
    risk1_lower[t] <- pooled_r1$ci_lower
    risk1_upper[t] <- pooled_r1$ci_upper
  }

  # Create pooled graph data frame
  graph_pooled <- data.frame(
    time_0 = time_0,
    time_days = time_0 / 24,
    risk0 = risk0_pooled,
    risk1 = risk1_pooled,
    risk0_lower = risk0_lower,
    risk0_upper = risk0_upper,
    risk1_lower = risk1_lower,
    risk1_upper = risk1_upper
  )

  # Also keep bootstrap samples from all imputations for risk ratio/difference calculations
  # Combine all bootstrap samples across imputations
  all_risk0_boot <- do.call(rbind, lapply(risk_matrices_list, function(x) x$risk0))
  all_risk1_boot <- do.call(rbind, lapply(risk_matrices_list, function(x) x$risk1))

  list(
    graph = graph_pooled,
    risk0_boot = all_risk0_boot,
    risk1_boot = all_risk1_boot,
    risk0_est_imp = risk0_est_imp,
    risk1_est_imp = risk1_est_imp,
    risk0_var_imp = risk0_var_imp,
    risk1_var_imp = risk1_var_imp
  )
}

# Pool results for each strategy
results_3.85 <- pool_strategy_results(boot_data_3.85_all, "ROX < 3.85")
results_4.88 <- pool_strategy_results(boot_data_4.88_all, "ROX < 4.88")
results_complex <- pool_strategy_results(boot_data_complex_all, "Time-varying ROX")

# Extract pooled graphs
graph_3.85 <- results_3.85$graph
graph_4.88 <- results_4.88$graph
graph_complex <- results_complex$graph

# Extract combined bootstrap matrices (for RD/RR calculations)
risk0_boot_3.85 <- results_3.85$risk0_boot
risk1_boot_3.85 <- results_3.85$risk1_boot
risk0_boot_4.88 <- results_4.88$risk0_boot
risk1_boot_4.88 <- results_4.88$risk1_boot
risk0_boot_complex <- results_complex$risk0_boot
risk1_boot_complex <- results_complex$risk1_boot

n_boot_combined <- N_IMPUTATIONS * N_BOOTSTRAP

cat("\n✓ Pooling complete\n")
cat(sprintf("Pooled estimates calculated for %d time points\n", nrow(graph_3.85)))
cat(sprintf("Combined bootstrap samples: %d\n\n", n_boot_combined))

# ============================================================================
# Create Combined Visualization
# ============================================================================

cat("=== Creating Visualizations ===\n")

# Grace period label for legend
grace_label <- ifelse(GRACE_HOURS == 0, "", paste0(", ", GRACE_HOURS, "hr grace"))

# Color palette
col_black <- "#000000"
col_navy <- "#0c385c"
col_light_blue <- "#1770b8"
col_dark_mauve <- "#6b8ac5"

# Combined plot with confidence intervals
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
plot_filename <- paste0("outputs/ROX/Fig.cumulative_incidence_all_strategies_m5_grace", GRACE_HOURS, "hr.png")
ggsave(plot_filename, plot_combined, width = 12, height = 7, dpi = 300)

cat("\n✓ Combined plot saved:", plot_filename, "\n")

# ============================================================================
# Create Comparison Table
# ============================================================================

cat("\n=== Creating Comparison Tables ===\n")

# Helper function to calculate RD with Rubin's Rule pooling
calc_rd_rubins <- function(time_h, risk0_est_imp, risk1_est_imp, risk0_var_imp, risk1_var_imp, graph_data) {
  idx <- which(graph_data$time_0 == time_h)

  # Calculate RD for each imputation
  rd_est_imp <- (risk1_est_imp[, idx] - risk0_est_imp[, idx]) * 100

  # Variance of RD for each imputation (using variance addition rule)
  rd_var_imp <- (risk0_var_imp[, idx] + risk1_var_imp[, idx]) * 100^2

  # Apply Rubin's Rule
  pooled <- rubins_rule(rd_est_imp, rd_var_imp, m = N_IMPUTATIONS)

  list(
    est = pooled$estimate,
    lower = pooled$ci_lower,
    upper = pooled$ci_upper
  )
}

# Calculate RD for each strategy and time point using Rubin's Rule
rd_3.85_7d <- calc_rd_rubins(168, results_3.85$risk0_est_imp, results_3.85$risk1_est_imp,
                              results_3.85$risk0_var_imp, results_3.85$risk1_var_imp, graph_3.85)
rd_3.85_14d <- calc_rd_rubins(336, results_3.85$risk0_est_imp, results_3.85$risk1_est_imp,
                               results_3.85$risk0_var_imp, results_3.85$risk1_var_imp, graph_3.85)
rd_3.85_30d <- calc_rd_rubins(720, results_3.85$risk0_est_imp, results_3.85$risk1_est_imp,
                               results_3.85$risk0_var_imp, results_3.85$risk1_var_imp, graph_3.85)

rd_4.88_7d <- calc_rd_rubins(168, results_4.88$risk0_est_imp, results_4.88$risk1_est_imp,
                              results_4.88$risk0_var_imp, results_4.88$risk1_var_imp, graph_4.88)
rd_4.88_14d <- calc_rd_rubins(336, results_4.88$risk0_est_imp, results_4.88$risk1_est_imp,
                               results_4.88$risk0_var_imp, results_4.88$risk1_var_imp, graph_4.88)
rd_4.88_30d <- calc_rd_rubins(720, results_4.88$risk0_est_imp, results_4.88$risk1_est_imp,
                               results_4.88$risk0_var_imp, results_4.88$risk1_var_imp, graph_4.88)

rd_complex_7d <- calc_rd_rubins(168, results_complex$risk0_est_imp, results_complex$risk1_est_imp,
                                 results_complex$risk0_var_imp, results_complex$risk1_var_imp, graph_complex)
rd_complex_14d <- calc_rd_rubins(336, results_complex$risk0_est_imp, results_complex$risk1_est_imp,
                                  results_complex$risk0_var_imp, results_complex$risk1_var_imp, graph_complex)
rd_complex_30d <- calc_rd_rubins(720, results_complex$risk0_est_imp, results_complex$risk1_est_imp,
                                  results_complex$risk0_var_imp, results_complex$risk1_var_imp, graph_complex)

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

  # Risk Difference (3.85 vs Usual Care) - Rubin's Rule
  RD_3.85 = sprintf("%.1f (%.1f to %.1f)",
                    c(rd_3.85_7d$est, rd_3.85_14d$est, rd_3.85_30d$est),
                    c(rd_3.85_7d$lower, rd_3.85_14d$lower, rd_3.85_30d$lower),
                    c(rd_3.85_7d$upper, rd_3.85_14d$upper, rd_3.85_30d$upper)),

  # Risk Difference (4.88 vs Usual Care) - Rubin's Rule
  RD_4.88 = sprintf("%.1f (%.1f to %.1f)",
                    c(rd_4.88_7d$est, rd_4.88_14d$est, rd_4.88_30d$est),
                    c(rd_4.88_7d$lower, rd_4.88_14d$lower, rd_4.88_30d$lower),
                    c(rd_4.88_7d$upper, rd_4.88_14d$upper, rd_4.88_30d$upper)),

  # Risk Difference (Complex vs Usual Care) - Rubin's Rule
  RD_Complex = sprintf("%.1f (%.1f to %.1f)",
                       c(rd_complex_7d$est, rd_complex_14d$est, rd_complex_30d$est),
                       c(rd_complex_7d$lower, rd_complex_14d$lower, rd_complex_30d$lower),
                       c(rd_complex_7d$upper, rd_complex_14d$upper, rd_complex_30d$upper)),

  # Risk Ratio (3.85 vs Usual Care) - using combined bootstrap
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

  # Risk Ratio (4.88 vs Usual Care) - using combined bootstrap
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

  # Risk Ratio (Complex vs Usual Care) - using combined bootstrap
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
      "CI indicates confidence interval; RD, risk difference; RR, risk ratio.",
      "Time-varying strategy: ROX < 2.85 (hours 1-5), ROX < 3.47 (hours 6-11), ROX < 3.85 (hours 12+).",
      "Values are percentages with 95% confidence intervals.",
      "Risk differences are in percentage points.",
      paste0("Grace period: ", GRACE_HOURS, " hour(s)."),
      paste0("Bootstrap resamples: ", n_boot_combined)
    )
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

ft_comparison

# Save to Word document
doc_comparison <- read_docx()

doc_comparison <- doc_comparison %>%
  body_add_par(paste0("Table. Comparison of Cumulative Mortality Outcomes Across ROX-Guided Intubation Strategies (Grace Period: ", GRACE_HOURS, " hour)"),
               style = "heading 1") %>%
  body_add_flextable(ft_comparison)

doc_filename <- paste0("outputs/ROX/mortality_comparison_all_strategies_m5_grace", GRACE_HOURS, "hr.docx")
print(doc_comparison, target = doc_filename)

cat("\n✓ Comparison table saved:", doc_filename, "\n")

# ============================================================================
# Step 6: Supplementary Table - Multiple time points
# ============================================================================
# Table with 7-day, 14-day, and 30-day mortality for each strategy

cat("\n=== Creating Step 6 Table (Multiple Time Points) ===\n")

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
results_step6_3.85 <- lapply(time_points, function(t) {
  extract_results(t, graph_3.85, risk0_boot_3.85, risk1_boot_3.85)
})

# ROX 4.88
results_step6_4.88 <- lapply(time_points, function(t) {
  extract_results(t, graph_4.88, risk0_boot_4.88, risk1_boot_4.88)
})

# Complex
results_step6_complex <- lapply(time_points, function(t) {
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
    "Time-varying ROX strategy§", "", ""
  ),

  Outcome = rep(time_labels, 3),

  Usual_care = c(
    results_step6_3.85[[1]]$risk0, results_step6_3.85[[2]]$risk0, results_step6_3.85[[3]]$risk0,
    results_step6_4.88[[1]]$risk0, results_step6_4.88[[2]]$risk0, results_step6_4.88[[3]]$risk0,
    results_step6_complex[[1]]$risk0, results_step6_complex[[2]]$risk0, results_step6_complex[[3]]$risk0
  ),

  ROX_risk = c(
    results_step6_3.85[[1]]$risk1, results_step6_3.85[[2]]$risk1, results_step6_3.85[[3]]$risk1,
    results_step6_4.88[[1]]$risk1, results_step6_4.88[[2]]$risk1, results_step6_4.88[[3]]$risk1,
    results_step6_complex[[1]]$risk1, results_step6_complex[[2]]$risk1, results_step6_complex[[3]]$risk1
  ),

  Risk_difference = c(
    results_step6_3.85[[1]]$rd, results_step6_3.85[[2]]$rd, results_step6_3.85[[3]]$rd,
    results_step6_4.88[[1]]$rd, results_step6_4.88[[2]]$rd, results_step6_4.88[[3]]$rd,
    results_step6_complex[[1]]$rd, results_step6_complex[[2]]$rd, results_step6_complex[[3]]$rd
  ),

  Risk_ratio = c(
    results_step6_3.85[[1]]$rr, results_step6_3.85[[2]]$rr, results_step6_3.85[[3]]$rr,
    results_step6_4.88[[1]]$rr, results_step6_4.88[[2]]$rr, results_step6_4.88[[3]]$rr,
    results_step6_complex[[1]]$rr, results_step6_complex[[2]]$rr, results_step6_complex[[3]]$rr
  )
)

# Set column names
colnames(table_step6) <- c(
  "Trial",
  "ROX strategy",
  "Time point",
  "Usual care†, % (95% CI)‡",
  "ROX strategy, % (95% CI)‡",
  "Risk difference, percentage points (95% CI)‡",
  "Risk ratio (95% CI)‡"
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
      "*Primary outcome time point.",
      "†The usual care strategy was identical across emulated target trials and estimated from the same underlying cohort.",
      paste0("‡Estimates were obtained from models adjusted for prespecified baseline covariates and hourly time-varying covariates; 95% confidence intervals were obtained using ", N_BOOTSTRAP, " bootstrap replicates pooled from ", N_IMPUTATIONS, " multiply imputed datasets."),
      "§Time-varying ROX strategy: thresholds were 2.85 (hours 1–5), 3.47 (hours 6–11), and 3.85 (≥12 hours).")
  ) %>%
  fontsize(part = "footer", size = 8) %>%
  italic(part = "footer")

ft_step6

# Save Step 6 table to Word document
path1_step6 <- paste0("outputs/ROX/table3_multiple_timepoints_m5_grace", GRACE_HOURS, "hr.docx")

ft_step6 %>%
  autofit() %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = path1_step6)


path2_step6 <- paste0("C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Table 9_ROX_strategy_MI_grace", GRACE_HOURS, "hr.docx")

ft_step6 %>%
  autofit() %>%
  set_table_properties(layout = "autofit") %>%
  save_as_docx(path = path1_step6)

cat("\n✓ Step 6 Table saved:", path1_step6, "\n")

# Also save to manuscript folder if it exists
if (dir.exists(dirname(path2_step6))) {
  ft_step6 %>%
    autofit() %>%
    set_table_properties(layout = "autofit") %>%
    save_as_docx(path = path2_step6)
  cat("✓ Step 6 Table also saved:", path2_step6, "\n")
}

cat("\n✓ Step 6 Table saved:", path1_step6, "\n")

# ============================================================================
# Analysis Complete
# ============================================================================

cat("\n============================================================================\n")
cat("  Analysis Complete\n")
cat("============================================================================\n")
cat("Summary:\n")
cat(sprintf("- Pooled %d multiply imputed datasets using Rubin's Rule\n", N_IMPUTATIONS))
cat(sprintf("- Combined %d bootstrap samples (100 per imputation)\n", n_boot_combined))
cat(sprintf("- Grace period: %d hour(s)\n", GRACE_HOURS))
cat("\nOutput files:\n")
cat("- Combined cumulative incidence plot\n")
cat("- Comparison table with all strategies\n")
cat("\nNote: Results incorporate both within-imputation and between-imputation variability\n")
cat("      using Rubin's Rule for proper uncertainty quantification.\n")
cat("============================================================================\n\n")

