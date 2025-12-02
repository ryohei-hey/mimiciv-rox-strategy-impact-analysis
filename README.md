# Impact-Analaysis-of-simple-prediction-model
# Title: Impact study of ROX index (simple prediction model) using MIMIC-IV. 

# Aim
The aim is to estimate the causal effect of a ROX-guided (simple prediction model-guided) strategy compared with usual care on 30-day mortality among adults who have been initiated on high flow nasal cannula (HFNC) for respiratory failure.

# Introduction  
The ROX index is a widely validated tool for predicting HFNC failure and has been proposed to support decisions regarding early intubation for patients with respiratory failure. The ROX index is defined as (SpO₂/FiO₂)/respiratory rate, and its predictive performance has been externaly validated across diverse diseases and care settings. However, whether using the ROX index to guide intubation decisions improves patient outcomes remains unknown. Prediction alone does not imply causation, and prognostic accuracy does not guarantee clinical utility. For use in routine care, prediction models require impact analysis, which evaluates whether model-guided decisions improve outcomes compared with usual practice. A target trial emulation (TTE) is well suited to this question because the treatment strategies—ROX-guided intubation versus usual care—can be clearly defined at HFNC initiation and aligned with routinely collected ICU data, thereby avoiding design-related biases such as immortal time bias. As no randomized trial has evaluated a ROX-guided intubation strategy among HFNC recipients, robust causal evidence derived from observational data is needed.

# PICO
P:
I:
C:
O: 
  
# Analysis Steps  
- Step1: Extracting person hourly data from MIMIC-IV
- Step2: DataWrangling  
- Step3: Imputation  
- Step4: Clone Censor Weighting analysis   
  

# File Contets
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
- First ICU admission only
- Age 18-85 years
- Exclude Neuro ICU patients

### Exclusion Criteria
1. HFNC started >1 hours before ICU admission
2. Invasive ventilation before HFNC
3. DNR/DNI/CMO status before HFNC initiation
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
