# Variable Extraction Methods for HFNC Study Dataset

## Overview

The final dataset contains **124 variables** organized into the following categories:

| Category | Number of Variables |
|----------|---------------------|
| Identifiers & Time Axes | 3 |
| Static: Basic IDs | 2 |
| Static: Patient Demographics | 6 |
| Static: Social Information | 2 |
| Static: Admission Information | 4 |
| Static: ICU Information | 2 |
| Static: Elixhauser Comorbidities | 31 |
| Static: Elixhauser Scores | 3 |
| Time-varying: HFNC-related | 4 |
| Time-varying: Mortality-related | 8 |
| Time-varying: Discharge Events | 5 |
| Time-varying: Outcome | 1 |
| Time-varying: Ventilation | 3 |
| Time-varying: Vital Signs | 11 |
| Time-varying: Respiratory Settings | 4 |
| Time-varying: Vasopressors | 6 |
| Time-varying: Blood Gas Analysis | 7 |
| Time-varying: Chemistry | 6 |
| Time-varying: Complete Blood Count | 4 |
| Time-varying: Therapeutic Intervention | 1 |
| Time-varying: SOFA Score | 8 |
| **Total** | **124** |

---

## Variable Extraction Methods (For Methods Section)

### Table 1. Identifiers and Time Axes

| Variable | Description | Source | Definition |
|----------|-------------|--------|------------|
| stay_id | ICU stay identifier | icu_icustays | Primary key for each ICU admission |
| hr | Hours since ICU admission | Calculated | FLOOR((charttime - intime) / 3600); hr=0 represents ICU admission |
| time | Hours since HFNC initiation | Calculated | hr - hr_at_hfnc_start; negative values indicate pre-HFNC period |

---

### Table 2. Static Variables: Patient Demographics

| Variable | Description | Source Table | Definition/Notes |
|----------|-------------|--------------|------------------|
| subject_id | Patient identifier | hosp_patients | Unique patient ID across admissions |
| hadm_id | Hospital admission ID | hosp_admissions | Unique admission ID |
| gender | Sex | hosp_patients | M=Male, F=Female |
| anchor_age | Age at anchor year | hosp_patients | Ages ≥89 years recorded as 91 for privacy protection |
| anchor_year | Anchor year | hosp_patients | De-identified reference year |
| anchor_year_group | Anchor year range | hosp_patients | Categorical grouping of anchor years |
| race | Race/ethnicity | hosp_admissions | Self-reported race/ethnicity |
| marital_status | Marital status | hosp_admissions | MARRIED, SINGLE, WIDOWED, etc. |
| insurance | Insurance type | hosp_admissions | Private, Medicare, Medicaid, etc. |
| language | Primary language | hosp_admissions | Primary language spoken |

---

### Table 3. Static Variables: Admission and ICU Information

| Variable | Description | Source Table | Definition/Notes |
|----------|-------------|--------------|------------------|
| intime | ICU admission time | icu_icustays | Timestamp of ICU entry (hr=0 reference point) |
| outtime | ICU discharge time | icu_icustays | Timestamp of ICU exit |
| admittime | Hospital admission time | hosp_admissions | Timestamp of hospital entry |
| dischtime | Hospital discharge time | hosp_admissions | Timestamp of hospital exit |
| admission_type | Admission type | hosp_admissions | EMERGENCY, URGENT, ELECTIVE, etc. |
| admission_location | Admission source | hosp_admissions | EMERGENCY ROOM, TRANSFER, etc. |
| discharge_location | Discharge destination | hosp_admissions | HOME, REHAB, SNF, etc. |
| los | ICU length of stay | Calculated | (outtime - intime) in days |
| first_careunit | First ICU type | icu_icustays | CCU, MICU, SICU, etc. |
| last_careunit | Last ICU type | icu_icustays | CCU, MICU, SICU, etc. |

---

### Table 4. Static Variables: Elixhauser Comorbidities (31 variables)

| Variable | Description | ICD Codes | Definition |
|----------|-------------|-----------|------------|
| CONGESTIVE_HEART_FAILURE | Congestive heart failure | ICD-9: 398.91, 402.01, 428.x; ICD-10: I09.9, I11.0, I50.x | Binary (0/1) |
| CARDIAC_ARRHYTHMIAS | Cardiac arrhythmias | ICD-9: 426.x, 427.x; ICD-10: I44.x, I47.x, I49.x | Binary (0/1) |
| VALVULAR_DISEASE | Valvular disease | ICD-9: 093.2, 394-397, 424.x; ICD-10: I05-I08, I34-I39 | Binary (0/1) |
| PULMONARY_CIRCULATION | Pulmonary circulation disorders | ICD-9: 415.0, 416.x, 417.x; ICD-10: I26.x, I27.x, I28.x | Binary (0/1) |
| PERIPHERAL_VASCULAR | Peripheral vascular disease | ICD-9: 440.x, 441.x, 443.x; ICD-10: I70.x, I71.x, I73.x | Binary (0/1) |
| HYPERTENSION | Hypertension (uncomplicated and complicated) | ICD-9: 401.x-405.x; ICD-10: I10-I16 | Binary (0/1) |
| PARALYSIS | Paralysis | ICD-9: 342.x, 343.x, 344.x; ICD-10: G04.1, G80-G83 | Binary (0/1) |
| OTHER_NEUROLOGICAL | Other neurological disorders | ICD-9: 330-337, 340-345; ICD-10: G10-G13, G20-G26, G35-G37, G40-G47 | Binary (0/1) |
| CHRONIC_PULMONARY | Chronic pulmonary disease | ICD-9: 490-496, 500-505, 506.4; ICD-10: J40-J47, J60-J67, J68.4 | Binary (0/1) |
| DIABETES_UNCOMPLICATED | Diabetes without chronic complications | ICD-9: 250.0-250.3; ICD-10: E10.0-E10.1, E11.0-E11.1 | Binary (0/1) |
| DIABETES_COMPLICATED | Diabetes with chronic complications | ICD-9: 250.4-250.9; ICD-10: E10.2-E10.9, E11.2-E11.9 | Binary (0/1) |
| HYPOTHYROIDISM | Hypothyroidism | ICD-9: 243-244; ICD-10: E00-E03 | Binary (0/1) |
| RENAL_FAILURE | Renal failure | ICD-9: 585-586, 588; ICD-10: N17-N19, N25 | Binary (0/1) |
| LIVER_DISEASE | Liver disease | ICD-9: 070.x, 570-573; ICD-10: B18.x, K70-K77 | Binary (0/1) |
| PEPTIC_ULCER | Peptic ulcer disease excluding bleeding | ICD-9: 531-534; ICD-10: K25-K28 | Binary (0/1) |
| AIDS | AIDS/HIV | ICD-9: 042-044; ICD-10: B20-B24 | Binary (0/1) |
| LYMPHOMA | Lymphoma | ICD-9: 200-202, 203.0; ICD-10: C81-C85, C88, C90.0, C96 | Binary (0/1) |
| METASTATIC_CANCER | Metastatic cancer | ICD-9: 196-199; ICD-10: C77-C80 | Binary (0/1) |
| SOLID_TUMOR | Solid tumor without metastasis | ICD-9: 140-172, 174-195; ICD-10: C00-C75 | Binary (0/1) |
| RHEUMATOID_ARTHRITIS | Rheumatoid arthritis/collagen vascular diseases | ICD-9: 446, 701.0, 710-714, 719.3, 720, 725; ICD-10: L94.x, M05-M06, M08, M12, M30-M36, M45 | Binary (0/1) |
| COAGULOPATHY | Coagulopathy | ICD-9: 286, 287.1, 287.3-287.5; ICD-10: D65-D68, D69.1, D69.3-D69.6 | Binary (0/1) |
| OBESITY | Obesity | ICD-9: 278.0; ICD-10: E66 | Binary (0/1) |
| WEIGHT_LOSS | Weight loss | ICD-9: 260-263, 783.2; ICD-10: E40-E46, R63.4 | Binary (0/1) |
| FLUID_ELECTROLYTE | Fluid and electrolyte disorders | ICD-9: 276; ICD-10: E86-E87 | Binary (0/1) |
| BLOOD_LOSS_ANEMIA | Blood loss anemia | ICD-9: 280.0; ICD-10: D50.0 | Binary (0/1) |
| DEFICIENCY_ANEMIAS | Deficiency anemias | ICD-9: 280.1-281.9; ICD-10: D50.8-D53 | Binary (0/1) |
| ALCOHOL_ABUSE | Alcohol abuse | ICD-9: 291, 303, 305.0; ICD-10: F10, E52, G62.1, etc. | Binary (0/1) |
| DRUG_ABUSE | Drug abuse | ICD-9: 292, 304, 305.2-305.9; ICD-10: F11-F16, F18-F19, etc. | Binary (0/1) |
| PSYCHOSES | Psychoses | ICD-9: 293-298; ICD-10: F20, F22-F25, F28-F29, etc. | Binary (0/1) |
| DEPRESSION | Depression | ICD-9: 296.2-296.3, 296.5, 300.4, 309, 311; ICD-10: F32-F33, etc. | Binary (0/1) |

**Elixhauser Composite Scores (3 variables):**
- elixhauser_vanwalraven: van Walraven weighted score for mortality prediction
- elixhauser_SID29: SID29 weighted score for length of stay prediction
- elixhauser_SID30: SID30 weighted score

**Note:** Primary diagnosis (seq_num = 1) was excluded from comorbidity assessment per standard methodology.

---

### Table 5. HFNC-Related Variables

| Variable | Description | Source | Definition |
|----------|-------------|--------|------------|
| first_hfnc_start | First HFNC start time | icu_chartevents | Timestamp when itemid=226732 AND value='High flow nasal cannula' first recorded |
| hr_at_hfnc_start | Hour at HFNC initiation | Calculated | FLOOR((first_hfnc_start - intime) / 3600) |
| original_hfnc_start | Original HFNC start (unadjusted) | icu_chartevents | Pre-adjustment HFNC start time |
| original_hr_at_hfnc_start | Original hour at HFNC start | Calculated | Pre-adjustment hr value |

**HFNC Definition:**
- Source: icu_chartevents, itemid = 226732 (Oxygen Delivery Device)
- Value: 'High flow nasal cannula' (exact string match; excludes "High flow neb")
- Continuity: Records within 24 hours were considered a continuous HFNC session

---

### Table 6. Time-Varying Variables: Ventilation Status

| Variable | Description | Source | itemid | Definition | Forward Fill |
|----------|-------------|--------|--------|------------|--------------|
| invasive | Invasive mechanical ventilation | icu_procedureevents | 225792 | Binary (0/1) | No |
| noninvasive | Non-invasive ventilation (BiPAP/CPAP) | icu_procedureevents | 225794 | Binary (0/1) | No |
| highflow | High flow nasal cannula | icu_chartevents | 226732 | Binary (0/1); value='High flow nasal cannula' | No |

---

### Table 7. Time-Varying Variables: Vital Signs

| Variable | Description | Source | itemid | Unit | Valid Range | Aggregation | Forward Fill |
|----------|-------------|--------|--------|------|-------------|-------------|--------------|
| heart_rate | Heart rate | icu_chartevents | 220045 | bpm | 0 < x < 300 | Mean | 24 hours |
| sbp | Systolic BP (arterial) | icu_chartevents | 220050, 220179, 225309 | mmHg | 0 < x < 400 | Mean | 24 hours |
| dbp | Diastolic BP (arterial) | icu_chartevents | 220051, 220180, 225310 | mmHg | 0 < x < 300 | Mean | 24 hours |
| mbp | Mean arterial pressure (arterial) | icu_chartevents | 220052, 220181, 225312 | mmHg | 0 < x < 300 | Mean | 24 hours |
| sbp_ni | Systolic BP (non-invasive) | icu_chartevents | 220179 | mmHg | 0 < x < 400 | Mean | 24 hours |
| dbp_ni | Diastolic BP (non-invasive) | icu_chartevents | 220180 | mmHg | 0 < x < 300 | Mean | 24 hours |
| mbp_ni | Mean arterial pressure (non-invasive) | icu_chartevents | 220181 | mmHg | 0 < x < 300 | Mean | 24 hours |
| resp_rate | Respiratory rate | icu_chartevents | 220210, 224690 | breaths/min | 0 < x < 70 | Mean | 24 hours |
| temperature | Body temperature | icu_chartevents | 223761 (°F), 223762 (°C) | °C | 10 < x < 50 (°C); 70 < x < 120 (°F) | Mean | 24 hours |
| spo2 | Oxygen saturation | icu_chartevents | 220277 | % | 0 < x ≤ 100 | Mean | 24 hours |
| glucose | Blood glucose (fingerstick) | icu_chartevents | 220621, 225664, 226537 | mg/dL | > 0 | Mean | 24 hours |

**Temperature Conversion:** Fahrenheit values (70-120°F) were converted to Celsius using (°F - 32) / 1.8

---

### Table 8. Time-Varying Variables: Respiratory Settings

| Variable | Description | Source | itemid | Unit | Valid Range | Aggregation | Forward Fill |
|----------|-------------|--------|--------|------|-------------|-------------|--------------|
| fio2 | Fraction of inspired oxygen | icu_chartevents / derived_bg | 223835 | fraction (0.21-1.0) | 21 ≤ x ≤ 100 (%) | Mean | 6 hours |
| o2_flow | Oxygen flow rate (primary) | icu_chartevents | 223834 | L/min | > 0 | Mean | 6 hours |
| o2_flow_additional | Additional O2 flow | icu_chartevents | 227287 | L/min | > 0 | Mean | 6 hours |
| o2_flow_total | Total O2 flow | Calculated | - | L/min | - | Sum | 6 hours |

**FiO2 Processing:**
- Values 0 < x ≤ 1: Multiplied by 100 to convert to percentage
- Values 1 < x < 21: Excluded as invalid
- Values 21 ≤ x ≤ 100: Used as-is
- FiO2 values populated only when ventilator (invasive/noninvasive/highflow) is active

---

### Table 9. Time-Varying Variables: Vasopressors and Inotropes

| Variable | Description | Source | itemid | Unit | Aggregation | Forward Fill |
|----------|-------------|--------|--------|------|-------------|--------------|
| norepinephrine | Norepinephrine infusion | icu_inputevents | 221906 | mcg/kg/min | Mean | 4 hours |
| epinephrine | Epinephrine infusion | icu_inputevents | 221289 | mcg/kg/min | Mean | 4 hours |
| dopamine | Dopamine infusion | icu_inputevents | 221662 | mcg/kg/min | Mean | 4 hours |
| phenylephrine | Phenylephrine infusion | icu_inputevents | 221749, 229617 | mcg/kg/min | Mean | 4 hours |
| dobutamine | Dobutamine infusion | icu_inputevents | 221653 | mcg/kg/min | Mean | 4 hours |
| vasopressor | Any vasopressor use | Calculated | - | Binary (0/1) | MAX | No |

**Note:** vasopressor = MAX(norepinephrine, epinephrine, dopamine, phenylephrine); dobutamine excluded from vasopressor composite as it is an inotrope.

---

### Table 10. Time-Varying Variables: Arterial Blood Gas Analysis

| Variable | Description | Source | itemid | Unit | Valid Range | Aggregation | Forward Fill |
|----------|-------------|--------|--------|------|-------------|-------------|--------------|
| ph | Arterial pH | hosp_labevents (derived_bg) | 50820 | - | - | Last | 6 hours |
| po2 | Partial pressure of oxygen | hosp_labevents (derived_bg) | 50821 | mmHg | - | Last | 6 hours |
| pco2 | Partial pressure of CO2 | hosp_labevents (derived_bg) | 50818 | mmHg | - | Last | 6 hours |
| lactate | Arterial lactate | hosp_labevents (derived_bg) | 50813 | mmol/L | - | Last | 6 hours |
| pao2fio2ratio | P/F ratio | Calculated | - | mmHg | - | Last | 6 hours |
| baseexcess | Base excess | hosp_labevents (derived_bg) | 50802 | mEq/L | - | Last | 6 hours |
| bicarbonate | Bicarbonate | hosp_labevents (derived_bg) | 50803 | mEq/L | - | Last | 6 hours |

**Specimen Filter:** Only arterial blood samples (specimen = 'ART.') were included.

---

### Table 11. Time-Varying Variables: Chemistry Panel

| Variable | Description | Source | itemid | Unit | Aggregation | Forward Fill |
|----------|-------------|--------|--------|------|-------------|--------------|
| creatinine | Serum creatinine | hosp_labevents (derived_chemistry) | 50912 | mg/dL | Last | 24 hours |
| sodium | Serum sodium | hosp_labevents (derived_chemistry) | 50983 | mEq/L | Last | 24 hours |
| potassium | Serum potassium | hosp_labevents (derived_chemistry) | 50971 | mEq/L | Last | 24 hours |
| chloride | Serum chloride | hosp_labevents (derived_chemistry) | 50902 | mEq/L | Last | 24 hours |
| bun | Blood urea nitrogen | hosp_labevents (derived_chemistry) | 51006 | mg/dL | Last | 24 hours |
| glucose_chem | Serum glucose (chemistry) | hosp_labevents (derived_chemistry) | 50931 | mg/dL | Last | 24 hours |

---

### Table 12. Time-Varying Variables: Complete Blood Count

| Variable | Description | Source | itemid | Unit | Aggregation | Forward Fill |
|----------|-------------|--------|--------|------|-------------|--------------|
| wbc | White blood cell count | hosp_labevents (derived_cbc) | 51301 | K/μL | Last | 48 hours |
| hemoglobin | Hemoglobin | hosp_labevents (derived_cbc) | 51222 | g/dL | Last | 48 hours |
| hematocrit | Hematocrit | hosp_labevents (derived_cbc) | 51221 | % | Last | 48 hours |
| platelet | Platelet count | hosp_labevents (derived_cbc) | 51265 | K/μL | Last | 48 hours |

---

### Table 13. Time-Varying Variables: Therapeutic Interventions

| Variable | Description | Source | itemid | Definition |
|----------|-------------|--------|--------|------------|
| crrt | Continuous renal replacement therapy | icu_chartevents | 224144 (Blood Flow), 227519 (Fluid Removal) | Binary (0/1); 1 if any CRRT parameter recorded |

---

### Table 14. Time-Varying Variables: SOFA Score and Sepsis-3

| Variable | Description | Definition | Valid Range |
|----------|-------------|------------|-------------|
| sofa_24hours | Total SOFA score | Sum of 6 organ component scores based on worst values in rolling 24-hour window | 0-24 |
| respiration_24hours | SOFA respiratory component | Based on P/F ratio: 0 (≥400), 1 (<400), 2 (<300), 3 (<200 with MV), 4 (<100 with MV) | 0-4 |
| coagulation_24hours | SOFA coagulation component | Based on platelets (K/μL): 0 (≥150), 1 (<150), 2 (<100), 3 (<50), 4 (<20) | 0-4 |
| liver_24hours | SOFA liver component | Based on bilirubin (mg/dL): 0 (<1.2), 1 (1.2-1.9), 2 (2.0-5.9), 3 (6.0-11.9), 4 (≥12) | 0-4 |
| cardiovascular_24hours | SOFA cardiovascular component | Based on MAP and vasopressor requirements | 0-4 |
| cns_24hours | SOFA neurological component | Based on GCS: 0 (15), 1 (13-14), 2 (10-12), 3 (6-9), 4 (<6) | 0-4 |
| renal_24hours | SOFA renal component | Based on creatinine (mg/dL): 0 (<1.2), 1 (1.2-1.9), 2 (2.0-3.4), 3 (3.5-4.9), 4 (≥5.0) | 0-4 |
| sepsis3 | Sepsis-3 diagnosis | Binary (0/1); 1 if suspected infection + SOFA ≥2 point increase from baseline | 0/1 |

**SOFA Calculation:** Based on MIMIC-IV official implementation (sofa.sql). Each component uses the worst value within a rolling 24-hour window.

**Sepsis-3 Criteria:** Based on MIMIC-IV official implementation (sepsis3.sql). Requires:
1. Suspected infection: Antibiotic administration and culture order within 24-hour window
2. SOFA score ≥2 points increase from baseline (assumed 0)

---

### Table 15. Mortality and Discharge Variables

| Variable | Description | Source | Definition |
|----------|-------------|--------|------------|
| deathtime_inhosp | In-hospital death time | hosp_admissions (deathtime) | Timestamp of in-hospital death |
| dod | Date of death | hosp_patients | Death date (includes out-of-hospital deaths) |
| death_datetime | Combined death datetime | Calculated | COALESCE(deathtime_inhosp, dod) |
| death_location | Location of death | Calculated | 'inhosp' or 'out_of_hosp' or NULL |
| time_to_death | Hours from ICU admission to death | Calculated | (death_datetime - intime) in hours |
| death_30day | 30-day mortality | Calculated | Binary (0/1); 1 if death within 720 hours of ICU admission |
| hr_at_death | Hour of death | Calculated | FLOOR((death_datetime - intime) / 3600) |
| death_event | Death event indicator | Calculated | Binary (0/1); 1 at the hour when death occurred |
| hr_at_icu_discharge | Hour of ICU discharge | Calculated | FLOOR((outtime - intime) / 3600) |
| hr_at_hosp_discharge | Hour of hospital discharge | Calculated | FLOOR((dischtime - intime) / 3600) |
| icu_discharge_event | ICU discharge indicator | Calculated | Binary (0/1); 1 when hr ≥ hr_at_icu_discharge |
| hosp_discharge_event | Hospital discharge indicator | Calculated | Binary (0/1); 1 when hr ≥ hr_at_hosp_discharge |
| discharge_outcome | Combined discharge outcome | Calculated | Binary (0/1); MAX(icu_discharge_event, hosp_discharge_event) |

---

## Outlier Handling Summary

### Table 16. Plausible Value Ranges for Outlier Exclusion

| Variable | Lower Bound | Upper Bound | Unit | Notes |
|----------|-------------|-------------|------|-------|
| heart_rate | >0 | <300 | bpm | Values outside range excluded |
| sbp | >0 | <400 | mmHg | Values outside range excluded |
| dbp | >0 | <300 | mmHg | Values outside range excluded |
| mbp | >0 | <300 | mmHg | Values outside range excluded |
| resp_rate | >0 | <70 | breaths/min | Values outside range excluded |
| temperature (°F) | >70 | <120 | °F | Converted to Celsius |
| temperature (°C) | >10 | <50 | °C | Values outside range excluded |
| spo2 | >0 | ≤100 | % | Values outside range excluded |
| fio2 | ≥21 | ≤100 | % | Values 1-20 excluded as invalid |

**FiO2 Special Processing:**
- Values 0 < x ≤ 1: Converted to percentage (×100)
- Values 1 < x < 21: Excluded as measurement errors
- Values 21 ≤ x ≤ 100: Used as recorded

---

## Forward Fill Strategy

### Table 17. Forward Fill Duration by Variable Category

| Variable Category | Forward Fill Duration | Rationale |
|-------------------|----------------------|-----------|
| Vital signs | 24 hours | Measurements typically performed every 1-4 hours in ICU |
| Arterial blood gas | 6 hours | Respiratory status can change rapidly; frequent reassessment expected |
| Chemistry panel | 24 hours | Routine daily laboratory measurements |
| Complete blood count | 48 hours | Less frequently measured; values stable over longer periods |
| Vasopressors | 4 hours | Infusion rates adjusted frequently |
| FiO2 / O2 flow | 6 hours | Respiratory settings changed based on patient status |

**Implementation:** Forward filling was implemented using a group counter (g-counter) approach with window functions, where each new measurement resets the group and values are carried forward within each group up to the maximum duration.

---

## Data Processing Notes

1. **Observation Window:** ICU admission (hr=0) to 720 hours (30 days) post-ICU admission or death, whichever occurred first

2. **Post-Death Data Censoring:** All observations after the hour of death were excluded

3. **HFNC Timing Adjustment:** For patients who initiated HFNC before ICU admission, hr_at_hfnc_start was adjusted to 0

4. **Missing Value Handling:** 
   - Continuous variables: NULL values retained
   - Binary indicators (invasive, noninvasive, highflow, vasopressor, crrt): NULL → 0 via COALESCE

5. **Aggregation Methods:**
   - Vital signs: Mean value within each hour
   - Laboratory values: Last value within each hour
   - Binary indicators: MAX within each hour
   - SOFA scores: MAX over rolling 24-hour window

6. **Derived Tables:** Blood gas (derived_bg), chemistry (derived_chemistry), and complete blood count (derived_cbc) tables were created following MIMIC-IV official SQL implementations

---

## Data Source Summary

| Source Table | Variables Extracted |
|--------------|---------------------|
| hosp_patients | subject_id, gender, anchor_age, anchor_year, anchor_year_group, dod |
| hosp_admissions | hadm_id, admittime, dischtime, deathtime, admission_type, admission_location, discharge_location, insurance, language, marital_status, race |
| icu_icustays | stay_id, first_careunit, last_careunit, intime, outtime, los |
| icu_chartevents | Vital signs, FiO2, O2 flow, HFNC status, CRRT, GCS components |
| icu_procedureevents | Invasive/noninvasive ventilation episodes |
| icu_inputevents | Vasopressors and inotropes |
| hosp_labevents | ABG, chemistry, CBC (via derived tables) |
| hosp_diagnoses_icd | Elixhauser comorbidities (ICD-9/10 codes) |

---

## Software and Analysis Environment

- **Database:** DuckDB (local) / BigQuery (MIMIC-IV 3.1)
- **Programming:** R with Quarto
- **Key Packages:** duckdb, tidyverse
- **MIMIC-IV Version:** 3.1
