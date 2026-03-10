#setup memory
rm(list = ls())
gc()
gc()

# load library-----
library(pacman)
p_load(tidyverse, here, lubridate, skimr, naniar, gtsummary,
       patchwork, ggsci, ggstream, ggalluvial, networkD3, htmlwidgets)

# load data-----
cat("=== 多重補完データの読み込み ===\n")

# 方法1: リストから読み込み（推奨）
dt_list <- readRDS("./data/ROX/imputed_cohort_mice_m5_list.rds")

cat(sprintf("読み込んだ補完データセット数: %d\n", length(dt_list)))
cat(sprintf("各データセットのサイズ: %s行 × %d列\n",
            format(nrow(dt_list[[1]]), big.mark = ","),
            ncol(dt_list[[1]])))

# 各補完データセットに共通の処理を適用する関数-----
wrangle_imputed_data <- function(dt, imputation_number = NULL) {

  if(!is.null(imputation_number)) {
    cat(sprintf("\n=== 補完データセット%d の処理開始 ===\n", imputation_number))
  }

  # time filtering
  dt <- dt %>% filter(time >= 0)

  # Intubation ----
  ## intubation_start, intubation_start_cum, and intubation ----

  dt <- dt %>%
    group_by(subject_id) %>%
    arrange(subject_id, hr) %>%
    mutate(
      # time > 0の範囲内で初めてinvasive == 1になった時点
      invasive_start = if_else(
        time > 0 & invasive == 1 &
          lag(cummax(if_else(time > 0, invasive, 0)), default = 0) == 0,
        1, 0
      ),
      # 元のinvasiveをバックアップしてから置き換え
      invasive_original = invasive,
      invasive = if_else(
        time == 0 & invasive == 1,
        lag(invasive, n = 1, default = 0),
        invasive),

      # 元のnoninvasiveをバックアップしてから置き換え
      noninvasive_original = noninvasive,
      noninvasive = if_else(
        time == 0 & noninvasive == 1,
        lag(noninvasive, n = 1, default = 0),
        noninvasive),

      # _flag: time >= 0以降で一度でもxxx == 1になったか
      invasive_flag = if_else(any(time >= 0 & invasive == 1, na.rm = TRUE), 1, 0),
      noninvasive_flag = if_else(any(time >= 0 & noninvasive == 1, na.rm = TRUE), 1, 0),
      highflow_flag = if_else(any(time >= 0 & highflow == 1, na.rm = TRUE), 1, 0),

      # cum
      invasive_start_cum = if_else(cumsum(invasive_start) >= 1, 1, 0)
    ) %>%
    ungroup()

  # time > 0で invasive == 1 なら highflow を強制的に 0 にする
  dt <- dt %>%
    mutate(
      # invasive使用中は他の呼吸サポートを0にする
      highflow = if_else(time > 0 & invasive == 1, 0, highflow),
      noninvasive = if_else(time > 0 & invasive == 1, 0, noninvasive),
      noninvasive = if_else(time > 0 & highflow == 1, 0, noninvasive)
    )

  # oxygen device trend ----

  # 1. 死亡患者を特定し、死亡時点の情報を取得
  death_info <- dt %>%
    filter(death_event == 1) %>%
    dplyr::select(subject_id, death_time = time,
           age, female, race_grouped, admission_location_grouped,
           intime, outtime, admittime, dischtime,
           baseline_sofa, elixhauser_vanwalraven, sepsis3,
           invasive_flag, vasopressor, los,
           death_datetime, death_30day, death_location)

  # 2. 死亡後の時間点を作成（death_time+1 から 720まで）
  death_extended <- death_info %>%
    group_by(subject_id) %>%
    reframe(
      time = (death_time + 1):720,
      death_time = death_time,
      across(c(age, female, race_grouped, admission_location_grouped,
               intime, outtime, admittime, dischtime,
               baseline_sofa, elixhauser_vanwalraven, sepsis3,
               invasive_flag, vasopressor, los,
               death_datetime, death_30day, death_location), first)
    ) %>%
    mutate(
      # 死亡後の状態
      death_event = 0,
      invasive = 0,
      noninvasive = 0,
      highflow = 0,
      # その他のバイタルはNA
      heart_rate = NA_real_,
      sbp = NA_real_,
      dbp = NA_real_,
      mbp = NA_real_,
      resp_rate = NA_real_,
      temperature = NA_real_,
      spo2 = NA_real_,
      fio2 = NA_real_,
      gcs = NA_real_,
      rox = NA_real_
    )

  # 3. 元のdtと結合
  dt_extended <- bind_rows(dt, death_extended) %>%
    arrange(subject_id, time)

  # 4. respiratory_supportを再計算
  dt_extended <- dt_extended %>%
    group_by(subject_id) %>%
    mutate(
      # 死亡時点を特定
      death_occurred = any(death_event == 1),
      death_time_point = if_else(death_occurred,
                                 time[which(death_event == 1)[1]],
                                 NA_real_),
      # 死亡後はずっとDead
      after_death = death_occurred & time > death_time_point,

      respiratory_support = case_when(
        after_death ~ "Dead",
        death_event == 1 ~ "Dead",
        invasive == 1 ~ "Invasive ventilation",
        noninvasive == 1 ~ "Non-invasive ventilation",
        highflow == 1 ~ "High-flow oxygen",
        TRUE ~ "Normal/No oxygen"
      ),
      respiratory_support = factor(respiratory_support,
                                   levels = c("Normal/No oxygen",
                                              "High-flow oxygen",
                                              "Non-invasive ventilation",
                                              "Invasive ventilation",
                                              "Dead"))
    ) %>%
    ungroup()

  dt <- dt_extended

  # Treatment Limitation ----
  dt <- dt %>%
    mutate(
      treatment_limitation = if_else(fullcode == 1, 0, 1)
    )

  # Lag all time-varying covariates ----
  lag_vars <- c("heart_rate", "resp_rate", "sbp", "dbp", "mbp",
                "spo2", "temperature", "fio2", "glucose", "gcs",
                "ph", "po2", "pco2",
                "treatment_limitation",
                "vasopressor", "crrt", "sofa_24hours")

  dt <- dt %>%
    group_by(subject_id) %>%
    mutate(across(all_of(lag_vars),
                  ~{
                    lagged <- lag(.)
                    coalesce(lagged, .)
                  })) %>%
    ungroup()

  # GCS category ----
  dt <- dt %>%
    mutate(gcs_category = case_when(
      gcs >= 13 ~ "13 to 15",
      gcs >= 9 ~ "9 to 12",
      gcs < 9 ~ "< 9",
      .default = NA_character_
    )) %>%
    mutate(gcs_category = factor(gcs_category,
                                 levels = c("13 to 15",
                                            "9 to 12", "< 9")))

  # ROX ----
  dt <- dt %>%
    mutate(
      rox = (spo2 / (fio2 / 100)) / resp_rate
    )

  # vasopressor_use ----
  dt <- dt %>%
    group_by(subject_id) %>%
    arrange(subject_id, hr) %>%
    mutate(
      vasopressor_use = if_else(
        any(time >= 0 & vasopressor == 1, na.rm = TRUE), 1, 0
      )
    ) %>%
    ungroup()

  # Censor ----
  dt <- dt %>%
    group_by(subject_id) %>%
    arrange(subject_id, hr) %>%
    mutate(
      censor = if_else(
        invasive_start_cum == 0 &
          highflow == 1 &
          lead(highflow, default = 0) == 0 &
          lead(noninvasive, default = 0) == 1,
        1, 0
      )
    ) %>%
    ungroup()

  # censor後のデータを削除
  dt <- dt %>%
    group_by(subject_id) %>%
    arrange(subject_id, hr) %>%
    mutate(
      censor_occurred = cumsum(censor) > 0,
      keep = !lag(censor_occurred, default = FALSE)
    ) %>%
    filter(keep) %>%
    dplyr::select(-censor_occurred, -keep) %>%
    ungroup()

  # Data preparation ----
  df <- dt
  df$id <- as.numeric(as.factor(df$subject_id))

  # Change Variable Name ----
  df$cal_time <- df$time
  df$cal_timesqr <- df$time^2
  df$death <- df$death_event
  df$treat <- df$invasive_start
  df$treat_cum <- df$invasive_start_cum

  df <- df %>%
    group_by(id) %>%
    arrange(id, hr) %>%
    mutate(
      treat_lag1 = lag(invasive_start, n = 1, default = 0),
      treat_cum_lag1 = lag(invasive_start_cum, n = 1, default = 0)
    ) %>%
    ungroup()

  df$time <- df$cal_time
  df$timesqr <- df$cal_timesqr

  df <- df %>%
    group_by(id) %>%
    mutate(highflow_lag = lag(highflow, n = 1,
                              default = first(highflow))) %>%
    ungroup()

  df <- df %>% dplyr::select(
    id, elig_1, elig_2, cal_time, cal_timesqr,
    hr_at_hfnc_start,anchor_year_group, 
    age, female, medicare_medicaid, language_grouped, Married, race_grouped,
    highflow, highflow_lag, invasive, noninvasive,
    heart_rate, mbp, temperature,
    resp_rate, spo2, fio2, rox,
    gcs_category, ph, po2, pco2, sofa_24hours, elixhauser_vanwalraven,
    sepsis3, antibiotic_use,vasopressor, crrt,
    treatment_limitation,
    icu_discharge_event,
    hosp_discharge_event,
    death_30day,
    death,
    censor,
    treat,
    treat_lag1,
    treat_cum,
    treat_cum_lag1,
    time,
    timesqr
  )

  if(!is.null(imputation_number)) {
    cat(sprintf("✓ 補完データセット%d の処理完了: %s行\n",
                imputation_number, format(nrow(df), big.mark = ",")))
  }

  return(df)
}

# 各補完データセットに処理を適用-----
cat("\n=== 全補完データセットの処理開始 ===\n")

processed_list <- list()

for(i in 1:length(dt_list)) {
  processed_list[[i]] <- wrangle_imputed_data(dt_list[[i]], imputation_number = i)
}

# 基本統計の確認-----
cat("\n=== 処理後の基本統計 ===\n")

for(i in 1:length(processed_list)) {
  cat(sprintf("\n補完データセット%d:\n", i))
  cat(sprintf("  サンプルサイズ: %s行\n", format(nrow(processed_list[[i]]), big.mark = ",")))
  cat(sprintf("  患者数: %d\n", n_distinct(processed_list[[i]]$id)))
  cat(sprintf("  挿管イベント数: %d\n", sum(processed_list[[i]]$treat)))
  cat(sprintf("  死亡イベント数: %d\n", sum(processed_list[[i]]$death)))
  cat(sprintf("  Censorイベント数: %d\n", sum(processed_list[[i]]$censor)))
}

# データの一致性確認（各補完データセットの構造が同じか）-----
cat("\n=== データ構造の確認 ===\n")

structure_check <- sapply(2:length(processed_list), function(i) {
  all(names(processed_list[[1]]) == names(processed_list[[i]])) &&
    nrow(processed_list[[1]]) == nrow(processed_list[[i]])
})

if(all(structure_check)) {
  cat("✓ 全補完データセットの構造が一致しています\n")
} else {
  cat("⚠️ 警告: 一部のデータセットで構造が異なります\n")
}

# 可視化（補完データセット1のみ使用）-----
cat("\n=== 可視化の作成（補完データセット1のみ）===\n")

# 時系列拡張したデータを使用（dt_extendedを再作成）
dt_plot <- dt_list[[1]] %>%
  filter(time >= 0) %>%
  # invasive_start等を再計算（簡略版）
  group_by(subject_id) %>%
  arrange(subject_id, hr) %>%
  mutate(
    invasive_start = if_else(
      time > 0 & invasive == 1 &
        lag(cummax(if_else(time > 0, invasive, 0)), default = 0) == 0,
      1, 0
    ),
    invasive_original = invasive,
    invasive = if_else(
      time == 0 & invasive == 1,
      lag(invasive, n = 1, default = 0),
      invasive),
    noninvasive_original = noninvasive,
    noninvasive = if_else(
      time == 0 & noninvasive == 1,
      lag(noninvasive, n = 1, default = 0),
      noninvasive),
    invasive_start_cum = if_else(cumsum(invasive_start) >= 1, 1, 0)
  ) %>%
  ungroup() %>%
  mutate(
    highflow = if_else(time > 0 & invasive == 1, 0, highflow),
    noninvasive = if_else(time > 0 & invasive == 1, 0, noninvasive),
    noninvasive = if_else(time > 0 & highflow == 1, 0, noninvasive)
  )

# 死亡情報
death_info_plot <- dt_plot %>%
  filter(death_event == 1) %>%
  dplyr::select(subject_id, death_time = time,
         age, female, race_grouped, admission_location_grouped,
         intime, outtime, admittime, dischtime,
         baseline_sofa, elixhauser_vanwalraven, sepsis3,
         invasive_flag = invasive_start_cum, vasopressor, los,
         death_datetime, death_30day, death_location)

death_extended_plot <- death_info_plot %>%
  group_by(subject_id) %>%
  reframe(
    time = (death_time + 1):720,
    death_time = death_time,
    across(c(age, female, race_grouped, admission_location_grouped,
             intime, outtime, admittime, dischtime,
             baseline_sofa, elixhauser_vanwalraven, sepsis3,
             invasive_flag, vasopressor, los,
             death_datetime, death_30day, death_location), first)
  ) %>%
  mutate(
    death_event = 0,
    invasive = 0,
    noninvasive = 0,
    highflow = 0
  )

dt_extended_plot <- bind_rows(dt_plot, death_extended_plot) %>%
  arrange(subject_id, time)

dt_extended_plot <- dt_extended_plot %>%
  group_by(subject_id) %>%
  mutate(
    death_occurred = any(death_event == 1),
    death_time_point = if_else(death_occurred,
                               time[which(death_event == 1)[1]],
                               NA_real_),
    after_death = death_occurred & time > death_time_point,

    respiratory_support = case_when(
      after_death ~ "Dead",
      death_event == 1 ~ "Dead",
      invasive == 1 ~ "Invasive ventilation",
      noninvasive == 1 ~ "Non-invasive ventilation",
      highflow == 1 ~ "High-flow oxygen",
      TRUE ~ "Normal/No oxygen"
    ),
    respiratory_support = factor(respiratory_support,
                                 levels = c("Normal/No oxygen",
                                            "High-flow oxygen",
                                            "Non-invasive ventilation",
                                            "Invasive ventilation",
                                            "Dead"))
  ) %>%
  ungroup()

# 時間ごとの分布
support_by_time <- dt_extended_plot %>%
  filter(time >= 0, time <= 720) %>%
  group_by(time, respiratory_support) %>%
  summarise(n = n(), .groups = "drop")


# Save data ----
cat("\n=== データ保存中 ===\n")

# 1. 各補完データセットを個別に保存
for(i in 1:length(processed_list)) {
  output_file <- sprintf("./data/ROX/imputed_data_for_analysis_m5_imp%d.rds", i)
  saveRDS(processed_list[[i]], output_file)
  cat(sprintf("✓ 保存: %s (%.1f MB)\n",
              output_file,
              file.info(output_file)$size / 1024^2))
}

# 2. リストとして保存
output_list <- "./data/ROX/imputed_data_for_analysis_m5_list.rds"
saveRDS(processed_list, output_list)
cat(sprintf("✓ リスト保存: %s (%.1f MB)\n",
            output_list,
            file.info(output_list)$size / 1024^2))

# サマリー-----
cat("\n========================================\n")
cat("✓✓✓ 多重補完データの処理完了！ ✓✓✓\n")
cat("========================================\n")
cat(sprintf("処理したデータセット数: %d\n", length(processed_list)))
cat(sprintf("各データセットのサンプルサイズ: %s行\n",
            format(nrow(processed_list[[1]]), big.mark = ",")))
cat(sprintf("変数数: %d\n", ncol(processed_list[[1]])))

cat("\n保存ファイル:\n")
cat("  個別ファイル:\n")
for(i in 1:length(processed_list)) {
  cat(sprintf("    - imputed_data_for_analysis_m5_imp%d.rds\n", i))
}
cat(sprintf("  リストファイル: %s\n", output_list))

cat("\n次のステップ:\n")
cat("  # 個別ファイルから読み込み:\n")
cat("  df1 <- readRDS('./data/ROX/imputed_data_for_analysis_m5_imp1.rds')\n\n")
cat("  # リストから読み込み:\n")
cat("  df_list <- readRDS('./data/ROX/imputed_data_for_analysis_m5_list.rds')\n")
cat("  df1 <- df_list[[1]]\n\n")
cat("  # 多重補完データでの分析例:\n")
cat("  # 各データセットで分析を実行し、結果をプール\n")
cat("  results <- lapply(df_list, function(df) {\n")
cat("    # 分析コード\n")
cat("  })\n")
cat("========================================\n")
