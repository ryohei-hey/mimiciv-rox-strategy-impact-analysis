library(duckdb)
library(DBI)

# パス設定
data_dir <- "C:/Users/ryohe/Dropbox (個人)/Research/mimic_projects/physionet.org/files/mimiciv/3.1"
db_path <- "data/mimic.duckdb"

# 既存のデータベースファイルを削除（クリーンスタート）
if (file.exists(db_path)) {
  cat("既存のデータベースファイルを削除中...\n")
  file.remove(db_path)
}

# DuckDB接続を作成
con <- dbConnect(duckdb::duckdb(), dbdir = db_path)

# hospフォルダのCSVファイルをインポート
hosp_dir <- file.path(data_dir, "hosp")
hosp_files <- list.files(hosp_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

cat("\nhospフォルダのテーブルをインポート中...\n")
cat(sprintf("見つかったファイル数: %d\n\n", length(hosp_files)))

for (file in hosp_files) {
  table_name <- gsub("\\.csv\\.gz$", "", basename(file))
  cat(sprintf("  - %s をインポート中...\n", table_name))
  
  tryCatch({
    # サンプルサイズを増やして型推定の精度を向上
    query <- sprintf(
      "CREATE TABLE hosp_%s AS SELECT * FROM read_csv_auto('%s', 
       sample_size=-1, 
       ignore_errors=false,
       nullstr=['', 'NA', 'NULL', '___', 'NaN'])",
      table_name,
      normalizePath(file, winslash = "/")
    )
    dbExecute(con, query)
    cat(sprintf("    ✓ 完了\n"))
  }, error = function(e) {
    cat(sprintf("    ✗ エラー: %s\n", substr(conditionMessage(e), 1, 100)))
    cat("    → all_varcharモードで再試行中...\n")
    
    # エラーが発生した場合は、全てVARCHARとして読み込む
    query_varchar <- sprintf(
      "CREATE TABLE hosp_%s AS SELECT * FROM read_csv_auto('%s', 
       all_varchar=1)",
      table_name,
      normalizePath(file, winslash = "/")
    )
    dbExecute(con, query_varchar)
    cat(sprintf("    ✓ VARCHAR型で完了\n"))
  })
}

# icuフォルダのCSVファイルをインポート
icu_dir <- file.path(data_dir, "icu")
icu_files <- list.files(icu_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)

cat("\nicuフォルダのテーブルをインポート中...\n")
cat(sprintf("見つかったファイル数: %d\n\n", length(icu_files)))

for (file in icu_files) {
  table_name <- gsub("\\.csv\\.gz$", "", basename(file))
  cat(sprintf("  - %s をインポート中...\n", table_name))
  
  tryCatch({
    query <- sprintf(
      "CREATE TABLE icu_%s AS SELECT * FROM read_csv_auto('%s', 
       sample_size=-1, 
       ignore_errors=false,
       nullstr=['', 'NA', 'NULL', '___', 'NaN'])",
      table_name,
      normalizePath(file, winslash = "/")
    )
    dbExecute(con, query)
    cat(sprintf("    ✓ 完了\n"))
  }, error = function(e) {
    cat(sprintf("    ✗ エラー: %s\n", substr(conditionMessage(e), 1, 100)))
    cat("    → all_varcharモードで再試行中...\n")
    
    query_varchar <- sprintf(
      "CREATE TABLE icu_%s AS SELECT * FROM read_csv_auto('%s', 
       all_varchar=1)",
      table_name,
      normalizePath(file, winslash = "/")
    )
    dbExecute(con, query_varchar)
    cat(sprintf("    ✓ VARCHAR型で完了\n"))
  })
}

# テーブル一覧を確認
cat("\n" , rep("=", 50), "\n", sep="")
cat("作成されたテーブル:\n")
cat(rep("=", 50), "\n", sep="")
tables <- dbListTables(con)
print(tables)

# 各テーブルの行数を表示
cat("\n", rep("=", 50), "\n", sep="")
cat("各テーブルの行数:\n")
cat(rep("=", 50), "\n", sep="")
for (table in tables) {
  count <- dbGetQuery(con, sprintf("SELECT COUNT(*) as n FROM %s", table))
  cat(sprintf("  %-30s: %15s 行\n", table, format(count$n, big.mark = ",")))
}

# 接続を閉じる
dbDisconnect(con, shutdown = TRUE)

cat("\n", rep("=", 50), "\n", sep="")
cat("✅ DuckDBの作成が完了しました！\n")
cat(rep("=", 50), "\n", sep="")
cat(sprintf("データベースファイル: %s\n", normalizePath(db_path)))
cat(sprintf("ファイルサイズ: %.2f MB\n", file.size(db_path) / 1024^2))
