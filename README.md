This repository provides reproducible code for an impact analysis framework that evaluates prediction scores as clinical strategies (model-guided policies) using target trial emulation with time-updated ICU data.
We demonstrate the approach by comparing ROX-guided intubation strategies versus usual care after HFNC initiation in MIMIC-IV.

# Aim
The aim is to estimate the causal effect of a ROX-guided (simple prediction model-guided) strategy compared with usual care on 30-day mortality among adults who have been initiated on high flow nasal cannula (HFNC) for respiratory failure.

# Introduction  
The ROX index is a widely validated tool for predicting High flow nasal canula (HFNC) failure and has been proposed to support decisions regarding early intubation for patients with respiratory failure. The ROX index is defined as (SpO₂/FiO₂)/respiratory rate, and its predictive performance has been externaly validated across diverse diseases and care settings. However, whether using the ROX index to guide intubation decisions improves patient outcomes remains unknown. Prediction alone does not imply causation, and prognostic accuracy does not guarantee clinical utility. For use in routine care, prediction models require impact analysis, which evaluates whether model-guided decisions improve outcomes compared with usual practice. A target trial emulation (TTE) is well suited to this question because the treatment strategies—ROX-guided intubation versus usual care—can be clearly defined at HFNC initiation and aligned with routinely collected ICU data, thereby avoiding design-related biases such as immortal time bias. As no randomized trial has evaluated a ROX-guided intubation strategy among HFNC recipients, robust causal evidence derived from observational data is needed.

# PICO
P: ICU patients on HFNC  
I: ROX guided intubation strategy  
C: Usual care (physician decsision making)   
O: 30 day mortality  
Follow up start/end: HFNC initiated (not ICU admission) to 30 days (720 hr)  
    
# Analysis Steps  
- Step0: Construct DuckDB (MIMIC-IV v3.1)  
- Step1: Extracting person hourly data
- Step2: Data cleaning  
- Step3: Imputation  
- Step4: Clone Censor Weighting analysis   

# Step0  
## Construct DuckDB (MIMIC-IV v3.1)
Please see instruction in below URL to build local DuckDB database from MIMIC-IV CSV/Parquet files.

https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iv/concepts_duckdb

# Step1  
## Goal - Person-Hourly Data Structure  
  
This project uses a **person-hourly data structure**, where each row represents one patient at one specific hour. This format enables time-to-event analysis with time-varying covariates.

### Example Data

| stay_id | hr | time | heart_rate | spo2 | fio2 | highflow |
|---------|-----|------|------------|------|------|----------|
| 101 | 0 | -6 | 88 | 94 | 40 | 0 |
| 101 | 1 | -5 | 90 | 93 | 40 | 0 |
| 101 | ... | ... | ... | ... | ... | ... |
| 101 | 6 | **0** | 92 | 91 | 60 | **1** |
| 101 | 7 | 1 | 88 | 93 | 55 | 1 |

*Note: `time = 0` indicates HFNC initiation*

### Two Time Axes

| Variable | Reference Point | Description |
|----------|-----------------|-------------|
| `hr` | ICU admission | Hours since ICU admission (always starts at 0) |
| `time` | HFNC initiation | Hours relative to HFNC start (negative = before, positive = after) |

### Advantages

- Suitable for **time-to-event analysis** (survival analysis)
- Captures **time-varying covariates** (vitals, labs, treatments) at each hour
- Enables **target trial emulation** with time-dependent exposures

### Final Dataset Size

~1,771 patients × ~700 hours ≈ **1+ million rows**
## HTML files  
If you see files through html, please click below URL  
  
- 01_Basic_hourly_data  
https://www.dropbox.com/scl/fi/d6b122q29q1ix5m8ozl9u/01_Basic_hourly_data.html?rlkey=pi460cscx7bap6xlbxp4ckqm0&st=803cxi45&dl=0  
  
- 02_time_varying_variables  
https://www.dropbox.com/scl/fi/3fwcnznqbz2dc6eyymjt1/02_time_varying_variables.html?rlkey=v6u4wemluzk9zzs1amsm5iuql&st=seee0p2r&dl=0  
  
03_final_integration  
https://www.dropbox.com/scl/fi/9u9wqh4jy1m1n199g5qif/03_final_integration.html?rlkey=g5owfpbw34fieio6w1s1mcuap&st=3jnu3wa8&dl=0  
  
## File Contets
## Overview

A three-part data processing pipeline for Target Trial Emulation of High Flow Nasal Cannula (HFNC) treatment outcomes in ICU patients. Processes MIMIC-IV database using DuckDB with 30-day mortality as the primary outcome.

---
🗄️ Data Infrastructure
**MIMIC-IV Database
**
- Downloaded from PhysioNet (credentialed access required)
- Version: MIMIC-IV v3.1
- Modules used: hosp, icu

**DuckDB Setup
**
- Local DuckDB database built from MIMIC-IV CSV/Parquet files
- Database path: data/mimic.duckdb
- Enables fast SQL queries without BigQuery costs
- Tables imported include:

  - icu_icustays, icu_chartevents, icu_outputevents
  - hosp_admissions, hosp_patients
  - Derived tables from MIMIC-IV official concepts
---

## 📁 01_Basic_hourly_data.qmd

### Purpose
Cohort selection and construction of basic person-hourly data structure

### Main Processing Steps

| Step | Description |
|------|-------------|
| 1-3 | Base table preparation (ICU stays, ventilator data) |
| 4 | **Cohort selection** (inclusion/exclusion criteria) |
| 5 | Time array generation (ICU admission to 720h post-HFNC) |
| 6 | Add `time` variable (relative to HFNC initiation) |
| 7 | Calculate death/discharge outcomes |
| 8 | Censor data after death |

### Inclusion Criteria
- HFNC initiated  
- First ICU admission only  
- Age 18-85 years  
- Non-neurological ICU  

### Exclusion Criteria
1. HFNC started >1 hours before ICU admission (allow patients who start HFNC therapy in ED or are initiated on HFNC to be transferred to ICU) 
2. Invasive ventilation before HFNC (Post extubation use were excluded)
3. DNR/DNI/CMO status before HFNC initiation (Any patients with treatment limitation was excluded)
4. Death before HFNC initiation

### Output Tables
- `hfnc_cohort_final` - Final cohort (~1,771 patients)
- `person_hourly_censored` - Person-hourly data (censored)

---

## 📁 02_time_varying_variables.qmd

### Purpose
Extract time-varying clinical variables and integrate into person-hourly data

### Main Processing Steps

| Step | Variable Category | Forward Fill |
|------|-------------------|--------------|
| 9 | Ventilator data (highflow, invasive, etc.) | - |
| 10 | FiO2 (derived_bg + chartevents) | 6 hours |
| 11 | Vasopressors & inotropes | - |
| 12 | Vital signs (HR, BP, SpO2, etc.) | 24 hours |
| 13 | O2 flow | 6 hours |
| 14 | ABG, Chemistry, CBC | 6-48 hours |
| 15 | CRRT | - |
| 16 | GCS | 24 hours |
| 17 | Suspicion of Infection | - |
| 18 | SOFA score calculation | - |
| 19 | Sepsis-3 diagnosis | - |

### Data Source Priority
- **FiO2**: Prioritize derived_bg (ABG measurements), supplement with chartevents
- **Vitals**: Use derived_vitalsign
- **Labs**: Use derived_chemistry, derived_complete_blood_count

### Forward Fill Strategy
| Variable Type | Time Limit |
|---------------|------------|
| Vital signs | 24 hours |
| Blood gas | 6 hours |
| Labs (Chem, CBC) | 48 hours |
| GCS | 24 hours |

### Output Tables
- `person_hourly_with_gcs` - Fully integrated data
- `derive_sofa` - SOFA scores (hourly)
- `derive_sepsis3` - Sepsis-3 diagnosis

---

## 📁 03_final_integration.qmd (Planned)

### Purpose
Final dataset integration, data quality checks, and export

### Planned Processing Steps

1. **Data Integration**
   - Merge File 01 + File 02
   - Define baseline variables (at time=0)
   - Calculate derived variables (ROX index, etc.)

2. **Data Quality Checks**
   - Missing data patterns
   - Outlier detection
   - Time series consistency

3. **Export**
   - Final dataset (CSV/Parquet)
   - Data dictionary
   - Analysis-ready subsets

---

## 🔗 File Dependencies

```
01_Basic_hourly_data.qmd
    ↓
    ├── hfnc_cohort_final
    ├── person_hourly_censored
    └── icu_stays_base
          ↓
02_time_varying_variables.qmd
    ↓
    ├── person_hourly_with_gcs (all variables)
    ├── derive_sofa
    └── derive_sepsis3
          ↓
03_final_integration.qmd
    ↓
    └── Final analysis dataset
```

---

## 📊 Cohort Flow Summary

```
Total ICU stays: ~94,000
    ↓
HFNC used: ~3,500
    ↓
After inclusion criteria: ~2,500
    ↓
After exclusion criteria: ~1,771
    ↓
Final cohort: 1,771 patients
(30-day mortality: ~32%)
```
