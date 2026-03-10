#setup memory
rm(list = ls())
gc()
gc()


# load library -----
library(pacman)

p_load(tidyverse, here, data.table,boot,lubridate, skimr, naniar,gtsummary,
       patchwork, ggsci, parallel,speedglm, doParallel, foreach,
       officer,flextable)

# load data
#dt <- readRDS("./data/ROX/imputed_data_for_analysis.rds")
dt <-readRDS("./data/ROX/imputed_data_for_analysis_m5_imp1.rds")
dt<-as.data.table(dt)
dt<-dt |> filter(elig_1==1) 

# Disable printing results in scientific notation
options(scipen=999)

# Set time point for estimation of risks (time 0:720 hours; 30 days)
K <- 721


### Censoring function for "intubate when ROX < cutoff" strategy WITH GRACE PERIOD
index<-3.85

glimpse(dt)

arm1.censor.grace <- function(rox, treat_cum, treat, highflow_lag, cutoff = index,
                              grace_hours = 1){
  n <- length(rox)
  my.censor <- rep(0, n)
  first.cens <- NA

  treat_cum_lag <- c(0, treat_cum[1:(n-1)])
  first_below_cutoff <- NA

  for(i in 1:n){
    # 挿管前 かつ 前時点でHFNC使用中のみチェック（roxがlag1なので時点を揃える）
    if(!is.na(treat_cum_lag[i]) && treat_cum_lag[i] == 0 &&
       !is.na(highflow_lag[i]) && highflow_lag[i] == 1){

      # === 非遵守1: ROX < cutoffへの対応（GRACE PERIOD付き） ===
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

      # === 非遵守2: ROX >= cutoffなのに挿管した ===
      if(!is.na(rox[i]) && !is.na(treat[i]) &&
         rox[i] >= cutoff && treat[i] == 1){
        first.cens <- i
        break
      }
    }
    # highflow_lag=0 の時点ではループは継続するが、遵守判定はスキップ
    # → 検閲は発生しない、フォローアップ継続
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

# boot strap ----

# 患者IDリスト
dt_ids <- data.frame(id = unique(dt$id))

clone.boot.graph <- function(data, indices, grace_hours = 1){
  
  # Select individuals into each bootstrapped sample
  boot.ids <- data.frame(id = data$id[indices])
  boot.ids$bid <- 1:nrow(boot.ids)
  
  # Subset person-time data to individuals selected into the bootstrapped sample
  d <- merge(boot.ids, dt, by = "id")
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
      race_grouped + antibiotic_use + 
      elixhauser_vanwalraven + hr_at_hfnc_start +
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


## Step 1: Do boot strap ----
# 利用可能なCPUコア数を確認
n_cores <- detectCores()
cat("利用可能なコア数:", n_cores, "\n")

# 並列処理に使うコア数（全体の80%程度推奨）
use_cores <- max(1, floor(n_cores * 0.8))
cat("使用するコア数:", use_cores, "\n")

# Grace period設定
GRACE_HOURS <- 0  # ここを変更することで感度分析が可能（0, 1, 2, 3など）

# 並列処理準備
cl <- makeCluster(use_cores)
registerDoParallel(cl)

clusterExport(cl, c("dt", "dt_ids", "arm1.censor.grace", "arm1.grace.flag", "clone.boot.graph", "GRACE_HOURS", "index"))
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
})

cat("\n=== 本番ブートストラップ（R=500, GRACE PERIOD=", GRACE_HOURS, "時間） ===\n")
cat("累積発生率曲線用\n")
cat("開始時刻:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

start_time <- Sys.time()
set.seed(12345)

boot_graphs_list <- foreach(
  i = 1:500,　# 500
  .packages = c("data.table", "dplyr"),
  .errorhandling = "pass"
) %dopar% {
  set.seed(12345 + i)
  indices <- sample(1:nrow(dt_ids), replace = TRUE)
  clone.boot.graph(dt_ids, indices, grace_hours = GRACE_HOURS)
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time

cat("\n完了時刻:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("実行時間:", round(as.numeric(elapsed_time, units = "mins"), 1), "分\n\n")

stopCluster(cl)

# 保存
output_filename <- paste0("outputs/ROX/bootdata/bootstrap_graphs_",index,"_r500_grace", GRACE_HOURS, "hr.rds")
saveRDS(boot_graphs_list, output_filename)
cat("✓ 結果を保存:", output_filename, "\n")

# READ DATA
# パラメータを設定
#index <- 3.85
#GRACE_HOURS <- 1

# ファイル名を生成
#input_filename <- paste0("outputs/ROX/bootdata/bootstrap_graphs_", index, "_r500_grace", GRACE_HOURS, "hr.rds")

# 読み込み
#boot_graphs_list <- readRDS(input_filename)

## Step 2: Calc 95% CI -----

n_boot <- length(boot_graphs_list)
n_time <- nrow(boot_graphs_list[[1]])

risk0_boot <- matrix(NA, nrow = n_boot, ncol = n_time)
risk1_boot <- matrix(NA, nrow = n_boot, ncol = n_time)

for (i in 1:n_boot) {
  if (is.data.frame(boot_graphs_list[[i]])) {
    risk0_boot[i, ] <- boot_graphs_list[[i]]$risk0
    risk1_boot[i, ] <- boot_graphs_list[[i]]$risk1
  }
}

time_0 <- boot_graphs_list[[1]]$time_0

max(time_0)

graph <- data.frame(
  time_0 = time_0,
  time_days = time_0 / 24,  # 時間から日へ変換
  
  # 点推定値
  risk0 = colMeans(risk0_boot, na.rm = TRUE),
  risk1 = colMeans(risk1_boot, na.rm = TRUE),
  
  # 95%信頼区間
  risk0_lower = apply(risk0_boot, 2, quantile, 0.025, na.rm = TRUE),
  risk0_upper = apply(risk0_boot, 2, quantile, 0.975, na.rm = TRUE),
  
  risk1_lower = apply(risk1_boot, 2, quantile, 0.025, na.rm = TRUE),
  risk1_upper = apply(risk1_boot, 2, quantile, 0.975, na.rm = TRUE)
)

## Step 3: Graph ------

plot.ipw <- ggplot(graph, aes(x = time_days)) +
  
  geom_ribbon(aes(ymin = risk0_lower, ymax = risk0_upper),
              fill = "#E7B800", alpha = 0.3) +
  geom_ribbon(aes(ymin = risk1_lower, ymax = risk1_upper),
              fill = "#2E9FDF", alpha = 0.3) +
  
  geom_line(aes(y = risk0, color = "Usual Care"), linewidth = 1.5) + 
  geom_line(aes(y = risk1, color = "ROX Strategy"), linewidth = 1.5) +
  
  xlab("Days") + 
  ylab("Cumulative Incidence of 30-Day Mortality") + 
  
  scale_x_continuous(limits = c(0, 30.5), 
                     breaks = seq(0, 30, by = 5)) + 
  scale_y_continuous(
    limits = c(0, 0.50), 
    breaks = seq(0, 0.50, by = 0.10),
    labels = scales::percent_format(accuracy = 1)
  ) + 
  
  theme_minimal() + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.3, 0.85),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(fill = "white", color = "gray80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  
  scale_color_manual(
    values = c("Usual Care" = "#E7B800", 
               "ROX Strategy" = "#2E9FDF"),
    labels = c("Usual Care" ="Usual Care", 
               "ROX Strategy" = paste0("ROX-guided intubation (ROX < ",index, ", ", GRACE_HOURS, " hr grace)"))
  )

print(plot.ipw)

plot_filename <- paste0("outputs/ROX/Fig.cumulative_incidence_r500_",index,"_30d_with_CI_grace", GRACE_HOURS, "hr.png")
ggsave(plot_filename, plot.ipw, width = 10, height = 7, dpi = 300)


# データの範囲を確認
nrow(graph)

# 最初の数行を詳しく見る
head(graph, 10)

# グラフが正しく表示されているか確認
print(plot.ipw)

# 実際にグラフが表示されていれば問題ありません
## Step 4: 30 day mortality -----

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

# Calculate RD for each time point using bootstrap percentile method
rd_7d <- calc_rd_bootstrap(168, risk0_boot, risk1_boot, graph)
rd_14d <- calc_rd_bootstrap(336, risk0_boot, risk1_boot, graph)
rd_30d <- calc_rd_bootstrap(720, risk0_boot, risk1_boot, graph)

# Complete results table with Risk Ratio
complete_results <- data.frame(
  Time_Point = c("7-Day", "14-Day", "30-Day"),

  Usual_Care = sprintf("%.1f (%.1f-%.1f)",
                       c(graph$risk0[graph$time_0 == 168] * 100,
                         graph$risk0[graph$time_0 == 336] * 100,
                         graph$risk0[graph$time_0 == 720] * 100),
                       c(graph$risk0_lower[graph$time_0 == 168] * 100,
                         graph$risk0_lower[graph$time_0 == 336] * 100,
                         graph$risk0_lower[graph$time_0 == 720] * 100),
                       c(graph$risk0_upper[graph$time_0 == 168] * 100,
                         graph$risk0_upper[graph$time_0 == 336] * 100,
                         graph$risk0_upper[graph$time_0 == 720] * 100)),

  ROX_Strategy = sprintf("%.1f (%.1f-%.1f)",
                         c(graph$risk1[graph$time_0 == 168] * 100,
                           graph$risk1[graph$time_0 == 336] * 100,
                           graph$risk1[graph$time_0 == 720] * 100),
                         c(graph$risk1_lower[graph$time_0 == 168] * 100,
                           graph$risk1_lower[graph$time_0 == 336] * 100,
                           graph$risk1_lower[graph$time_0 == 720] * 100),
                         c(graph$risk1_upper[graph$time_0 == 168] * 100,
                           graph$risk1_upper[graph$time_0 == 336] * 100,
                           graph$risk1_upper[graph$time_0 == 720] * 100)),

  # Risk Difference - Bootstrap percentile method
  Risk_Difference = sprintf("%.1f (%.1f to %.1f)",
                            c(rd_7d$est, rd_14d$est, rd_30d$est),
                            c(rd_7d$lower, rd_14d$lower, rd_30d$lower),
                            c(rd_7d$upper, rd_14d$upper, rd_30d$upper)),

  Risk_Ratio = sprintf("%.2f (%.2f-%.2f)",
                       sapply(c(168, 336, 720), function(h) {
                         idx <- which(graph$time_0 == h)
                         mean(risk1_boot[, idx] / risk0_boot[, idx], na.rm = TRUE)
                       }),
                       sapply(c(168, 336, 720), function(h) {
                         idx <- which(graph$time_0 == h)
                         quantile(risk1_boot[, idx] / risk0_boot[, idx], 0.025, na.rm = TRUE)
                       }),
                       sapply(c(168, 336, 720), function(h) {
                         idx <- which(graph$time_0 == h)
                         quantile(risk1_boot[, idx] / risk0_boot[, idx], 0.975, na.rm = TRUE)
                       }))
)

colnames(complete_results) <- c(
  "Time Point",
  "Usual Care",
  "ROX-Guided Strategy",
  "Risk Difference",
  "Risk Ratio"
)

# 論文スタイルのflextable
ft_paper <- flextable(complete_results) %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "center", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 10, part = "all") %>%
  add_header_row(
    values = c("", "Cumulative Mortality, % (95% CI)", "", ""),
    colwidths = c(1, 2, 1, 1)
  ) %>%
  align(part = "header", i = 1, align = "center") %>%
  bold(part = "header") %>%
  add_footer_lines(
    values = c(
      "CI indicates confidence interval.",
      "Values are percentages with 95% confidence intervals in parentheses.",
      "Risk differences are in percentage points.",
      paste0("Grace period: ", GRACE_HOURS, " hour(s)"),
      paste0("Bootstrap resamples: ", n_boot)
    )
  ) %>%
  fontsize(part = "footer", size = 9) %>%
  italic(part = "footer")

ft_paper

# 論文スタイルWord文書
doc_paper <- read_docx()

doc_paper <- doc_paper %>%
  body_add_par(paste0("Table 1. Cumulative Mortality Outcomes at 7, 14, and 30 Days (Grace Period: ", GRACE_HOURS, " hour)"), 
               style = "heading 1") %>%
  body_add_flextable(ft_paper)

doc_filename <- paste0("outputs/ROX/mortality_results_7_14_30day_",index,"_grace", GRACE_HOURS, "hr.docx")
print(doc_paper, target = doc_filename)

cat("\n✓ 論文スタイル版を保存:", doc_filename, "\n")

