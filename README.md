# Effect of a ROX Index–Guided Intubation Strategy on Outcomes in Patients Receiving High-Flow Nasal Cannula: A Target Trial Emulation

> Analysis code for [Journal Name]

## Overview

This repository contains the analysis code for a **target trial emulation** evaluating whether **ROX Index–guided intubation strategies** improve outcomes in ICU patients receiving **high-flow nasal cannula (HFNC)**. The study was designed and reported in accordance with the **TARGET** (Transparent Reporting of Observational Studies Emulating a Target Trial) guideline.

Using a **clone–censor–weight** approach, we emulate three separate two-arm target trials, each comparing a ROX-guided intubation strategy with usual care at HFNC initiation.

**Data source**: [MIMIC-IV Clinical Database v3.1](https://physionet.org/content/mimiciv/3.1/)

## Methods Summary

We emulated three two-arm target trials using MIMIC-IV v3.1, each comparing usual care with a ROX-guided intubation strategy among ICU patients on HFNC. ROX thresholds were < 3.85 (primary), < 4.88, and a time-varying cutoff. Per-protocol effects on 30-day mortality (primary outcome) and 7/14-day mortality (secondary) were estimated via the clone–censor–weight approach with time-varying inverse probability weighting and 500 bootstrap resamples. For full details, see the accompanying manuscript.

## Data Access

The analysis uses the **MIMIC-IV v3.1** database, which requires credentialed access through PhysioNet.

1. Complete the required CITI training course
2. Sign the data use agreement on [PhysioNet](https://physionet.org/content/mimiciv/3.1/)
3. Download the MIMIC-IV v3.1 CSV files

**This repository contains analysis scripts only.** Data files and outputs are not included.

## Requirements

- **R** (>= 4.5.0)
- **DuckDB** (for database storage and querying)

### Key R Packages

| Category                 | Packages                                                    |
| ------------------------ | ----------------------------------------------------------- |
| Data manipulation        | `tidyverse`, `data.table`, `lubridate`, `here`              |
| Database                 | `duckdb`, `DBI`                                             |
| Imputation               | `mice`                                                      |
| Statistical analysis     | `boot`, `speedglm`, `EValue`                                |
| Parallel computing       | `parallel`, `doParallel`, `foreach`                         |
| Visualization            | `patchwork`, `ggsci`, `ggalluvial`, `ggstream`, `networkD3` |
| Tables & reporting       | `gtsummary`, `flextable`, `officer`                         |
| Missing data diagnostics | `naniar`, `skimr`                                           |
| Utilities                | `pacman`, `tictoc`, `benchmarkme`                           |

Install all packages at once:

```r
install.packages(c(
  "tidyverse", "data.table", "lubridate", "here",
  "duckdb", "DBI",
  "mice", "naniar", "skimr",
  "boot", "speedglm", "EValue",
  "parallel", "doParallel", "foreach",
  "patchwork", "ggsci", "ggalluvial", "ggstream", "networkD3", "htmlwidgets",
  "gtsummary", "flextable", "officer",
  "pacman", "tictoc", "benchmarkme"
))
```

## Repository Structure

```
scripts/
├── for duckDB/                  # Database creation from MIMIC CSV files
├── for person_hourly_data_ROX/  # Quarto docs for hourly data extraction
├── Main Analysis/               # Primary analysis pipeline
└── Sensitivity Analyses/        # Robustness checks
```

## Analysis Pipeline

Run the scripts in the following order to reproduce the analysis.

### Step 1: Database Setup

```
scripts/for duckDB/01_create_duckdb.R
```

Creates a DuckDB database from the raw MIMIC-IV v3.1 CSV files for efficient SQL querying. See the [MIMIC-IV DuckDB setup guide](https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iv/concepts_duckdb) for instructions.

### Step 2: Data Wrangling

```
scripts/Main Analysis/01_wangling.R
```

Extracts and structures hourly ICU data for the HFNC cohort, including vital signs, laboratory values, and treatment variables.

### Step 3: Multiple Imputation

```
scripts/Main Analysis/02_Imputation_m5.R
```

Performs multiple imputation using MICE (m = 5) to handle missing time-varying covariates.

### Step 4: Post-Imputation Processing

```
scripts/Main Analysis/03_Wrangling_after_imputation_m5.R
```

Processes imputed datasets, derives analysis variables (e.g., ROX index), and prepares the final analytic cohort.

### Step 5: Censoring & Cloning

```
scripts/Main Analysis/04_01_censor_table_grace0to2_with_niv.R
```

Implements the clone–censor–weight framework: clones each patient into treatment strategy arms, applies censoring rules with grace periods (0–2 hours), and generates censoring indicators.

Visualization scripts:
- `04_02_censor_plot_grace0to2.R` — Censoring pattern plots
- `04_03_censor_plot_cumulative_status.R` — Cumulative treatment status over time

### Step 6: Bootstrap Clone–Censor–Weighting Analysis

Three ROX cutoff strategies are evaluated, each via 500 bootstrap iterations with stabilized inverse probability weighting:

| Script                               | ROX Cutoff           |
| ------------------------------------ | -------------------- |
| `05_1_Bootstrap_CCW_3.85.R`          | Fixed cutoff at 3.85 |
| `05_2_Bootstrap_CCW_4.88.R`          | Fixed cutoff at 4.88 |
| `05_3_Bootstrap_CCW_complexcutoff.R` | Time-varying cutoff  |

### Step 7: Results Summary

```
scripts/Main Analysis/05_4_summary 3.85_4.88_complex_grace1hr.R
```

Pools bootstrap results across cutoff strategies, generates cumulative incidence curves, risk difference estimates, and publication-ready figures.

## Sensitivity Analysis Scripts

All sensitivity analysis scripts are in `scripts/Sensitivity Analyses/`:

| Analysis                             | Description                                                                          |
| ------------------------------------ | ------------------------------------------------------------------------------------ |
| **Grace period = 1 hour**            | `05_1–3_*_grace0.R` — Requiring intubation within 1 hour of threshold breach         |
| **Alternative eligibility criteria** | `05_1–3_*_elig2.R` — Excluding patients intubated within 1 hour of HFNC initiation   |
| **Excluding respiratory variables**  | `05_1–3_*_no_resp_*.R` — Removes respiratory rate, SpO2, and FiO2 from weight models |
| **Multiple imputation (m = 5)**      | `05_01–03_*_m5.R` — 200 bootstrap resamples per dataset, combined with Rubin's rules |

Summary scripts for each sensitivity analysis:
- `05_06_summary*_grace0hr.R`
- `05_07_summary*_elig2_grace1hr.R`
- `05_08_summary*_no_resp_spo2_fio2_grace1hr.R`
- `05_09_summary*_grace1hr_m5.R`

## Citation

If you use this code, please cite the accompanying manuscript (citation details will be added upon publication).

## License

This project is licensed under the [MIT License](LICENSE).
