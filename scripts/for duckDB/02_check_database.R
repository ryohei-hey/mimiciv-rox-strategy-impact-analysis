library(duckdb)
library(DBI)

# データベースに接続（読み取り専用）
con <- dbConnect(duckdb::duckdb(), dbdir = "data/mimic.duckdb", read_only = TRUE)

# テーブル一覧
cat("=== テーブル一覧 ===\n")
tables <- dbListTables(con)
print(tables)

# 各テーブルの行数を確認
cat("\n=== 各テーブルの行数 ===\n")
for (table in tables) {
  count <- dbGetQuery(con, sprintf("SELECT COUNT(*) as n FROM %s", table))
  cat(sprintf("%-30s: %15s 行\n", table, format(count$n, big.mark = ",")))
}

# サンプルデータの表示（例：patients）
if ("hosp_patients" %in% tables) {
  cat("\n=== hosp_patients のサンプルデータ（5行）===\n")
  sample_data <- dbGetQuery(con, "SELECT * FROM hosp_patients LIMIT 5")
  print(sample_data)
}

# テーブルの構造を確認（例：admissions）
if ("hosp_admissions" %in% tables) {
  cat("\n=== hosp_admissions のカラム情報 ===\n")
  columns <- dbGetQuery(con, "DESCRIBE hosp_admissions")
  print(columns)
}

dbDisconnect(con, shutdown = TRUE)
