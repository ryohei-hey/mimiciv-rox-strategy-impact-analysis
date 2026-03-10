#setup memory
rm(list = ls())
gc()
gc()

setwd("C:/Users/ryohe/Dropbox (個人)/Research/mimic_projects/mimic_analysis")

# load library-----
library(pacman)
p_load(tidyverse,here, lubridate,skimr,naniar,gtsummary,flextable,
       slider)

#load data
dt<-read_csv("output/final_dataset.csv") 

# Data check-----
## Unique person and first admission----
unique(dt$stay_id) |> length()
unique(dt$subject_id) |> length()
unique(dt$hadm_id) |> length()

## death check----

## Check basic information of death-related variables
dt %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  summarise(
    total_patients = n(),
    in_hospital_deaths = sum(!is.na(deathtime_inhosp)),
    post_discharge_deaths = sum(death_location == "post_discharge", na.rm = TRUE),
    in_hospital_deaths_location = sum(death_location == "in_hospital", na.rm = TRUE),
    no_death = sum(is.na(death_datetime)),
    death_within_30days = sum(death_30day == 1, na.rm = TRUE)
  )

## Compare hr_at_death and max_hr for deceased patients
death_patients_check <- dt %>%
  filter(!is.na(death_datetime)) %>%  # Deceased patients only
  group_by(stay_id) %>%
  summarise(
    hr_at_death = first(hr_at_death),
    max_time = max(time, na.rm = TRUE),
    min_time = min(time, na.rm = TRUE),
    n_records = n(),
    death_location = first(death_location),
    .groups = "drop"
  ) %>%
  mutate(
    data_after_death = max_time > hr_at_death,  # Whether there is data after death
    hours_after_death = max_time - hr_at_death  # How many hours of data after death
  )

# Number of patients with data after death (OK if FALSE is 100%)
death_patients_check %>%
  count(data_after_death) %>%
  mutate(percentage = round(n / sum(n) * 100, 2))


# Variable wrangling----

# Age at ICU admission
dt <- dt %>%
  mutate(
    # Extract year from intime
    intime_year = year(intime),
    # Calculate difference from anchor_year
    year_diff = intime_year - anchor_year,
    # Calculate age at ICU admission by adding year difference to anchor_age
    age = anchor_age + year_diff
  )

## gender-> INTEGER（M=1, F=0）----
dt <- dt %>%
  mutate(
    female = case_when(
      gender == "M" ~ 0L,
      gender == "F" ~ 1L,
      TRUE ~ NA_integer_
    )
  )


## insurance----
dt <- dt %>%
  mutate(medicare_medicaid = case_when(
    insurance == 'Medicaid' | insurance == 'Medicare' ~ 1,
    is.na(insurance) ~ NA_real_,
    TRUE ~ 0)) 


dt <- dt %>%
  mutate(
    insurance_grouped = case_when(
      insurance == "Medicare" ~ "Medicare",
      insurance == "Medicaid" ~ "Medicaid",
      insurance %in% c("Private", "Other", "No charge") ~ "Private/Other",
      is.na(insurance) ~ NA_character_,
    ) %>% factor(levels = c("Medicare", "Medicaid", "Private/Other"))
  )

## Language----
dt <- dt %>%
  mutate(
    language_grouped = case_when(
      language == "English" ~ "English",
      language == "Spanish" ~ "Spanish",
      is.na(language) ~ NA_character_,
      !is.na(language) & language != "English" & language != "Spanish" ~ "Other",
    ) %>% factor(levels = c("English", "Spanish", "Other"))
  )


## marital_status----
dt<-dt |> mutate(
  Married = case_when(
    marital_status=="MARRIED" ~1,
    is.na(marital_status) ~ NA_real_,
    TRUE~0
  )
)

dt <- dt %>%
  mutate(
    marital_status = case_when(
      marital_status == "MARRIED" ~ "Married",
      marital_status == "SINGLE" ~ "Single",
      marital_status == "DIVORCED" ~ "Divorced",
      marital_status == "WIDOWED" ~ "Widowed",
      is.na(marital_status) ~ NA_character_,
    ) %>% factor(levels = c("Married", "Single", "Divorced", "Widowed"))
  )

## admission type ----
dt <- dt %>%
  mutate(
    admission_type_grouped = case_when(
      admission_type %in% c("EW EMER.", "DIRECT EMER.", "URGENT") ~ "Emergency/Urgent",
      admission_type %in% c("SURGICAL SAME DAY ADMISSION", "ELECTIVE") ~ "Planned Surgery",
      admission_type == "OBSERVATION ADMIT" ~ "Observation",
      is.na(admission_type) ~ NA_character_,
    ) %>% factor(levels = c("Emergency/Urgent", "Planned Surgery", "Observation"))
  )

## race ----
dt <- dt %>%
  mutate(
    race_grouped = case_when(
      # White (combine all White categories)
      str_detect(race, "WHITE") ~ "White",

      # Black/African American (combine all Black categories)
      str_detect(race, "BLACK") ~ "Black",

      # Hispanic/Latino (combine all Hispanic categories)
      str_detect(race, "HISPANIC|LATINO") ~ "Hispanic/Latino",

      # Asian (combine all Asian categories and Pacific Islander)
      str_detect(race, "ASIAN|HAWAIIAN|PACIFIC") ~ "Asian/Pacific Islander",

      # Unknown (UNKNOWN, UNABLE TO OBTAIN, DECLINED)
      race %in% c("UNKNOWN", "UNABLE TO OBTAIN", "PATIENT DECLINED TO ANSWER") ~ NA_character_,

      # Other (American Indian, all others)
      TRUE ~ "Other"
    ) %>% factor(levels = c("White", "Black", "Hispanic/Latino", 
                            "Asian/Pacific Islander", "Other"))
  )

## admission_location----
dt <- dt %>%
  mutate(
    admission_location_grouped = case_when(
      # Transfer from other hospital
      admission_location == "TRANSFER FROM HOSPITAL" ~ "Transfer from Hospital",

      # Referral from healthcare provider (physician/clinic)
      admission_location %in% c("PHYSICIAN REFERRAL", "CLINIC REFERRAL") ~ "Physician/Clinic Referral",

      # Emergency room / Walk-in
      admission_location %in% c("EMERGENCY ROOM", "WALK-IN/SELF REFERRAL") ~ "Emergency/Walk-in",

      # Transfer from facility (SNF)
      admission_location == "TRANSFER FROM SKILLED NURSING FACILITY" ~ "From Nursing Facility",

      # Internal transfer (procedure room, OR, PACU, psychiatry, etc.)
      admission_location %in% c("PROCEDURE SITE", "AMBULATORY SURGERY TRANSFER",
                                "PACU", "INTERNAL TRANSFER TO OR FROM PSYCH") ~ "Internal Transfer",

      # Unknown
      admission_location == "INFORMATION NOT AVAILABLE" ~ NA_character_
    ) %>% factor(levels = c("Transfer from Hospital", "Physician/Clinic Referral", 
                            "Emergency/Walk-in", "From Nursing Facility", 
                            "Internal Transfer"))
  )


## Blood Pressure -----
dt %>%
  summarise(
    # Invasive blood pressure missing count
    sbp_missing = sum(is.na(sbp)),
    dbp_missing = sum(is.na(dbp)),
    mbp_missing = sum(is.na(mbp)),
    # Non-invasive blood pressure missing count
    sbp_ni_missing = sum(is.na(sbp_ni)),
    dbp_ni_missing = sum(is.na(dbp_ni)),
    mbp_ni_missing = sum(is.na(mbp_ni)),
    # Total rows
    total_rows = n()
  ) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "count") %>%
  mutate(percentage = round(count / max(count) * 100, 1))

# 2. Calculate and impute mbp
dt <- dt %>%
  mutate(
    # Calculate mbp from sbp and dbp if missing (mbp = dbp + (sbp - dbp)/3)
    mbp = case_when(
      !is.na(mbp) ~ mbp,  # Keep mbp if available
      !is.na(sbp) & !is.na(dbp) ~ dbp + (sbp - dbp) / 3,  # Calculate
      TRUE ~ NA_real_
    ),

    # Calculate mbp_ni from sbp_ni and dbp_ni if missing
    mbp_ni = case_when(
      !is.na(mbp_ni) ~ mbp_ni,  # Keep mbp_ni if available
      !is.na(sbp_ni) & !is.na(dbp_ni) ~ dbp_ni + (sbp_ni - dbp_ni) / 3,  # Calculate
      TRUE ~ NA_real_
    )
  )

# 3. Use non-invasive values if invasive measurements are missing
dt <- dt %>%
  mutate(
    sbp = if_else(is.na(sbp) & !is.na(sbp_ni), sbp_ni, sbp),
    dbp = if_else(is.na(dbp) & !is.na(dbp_ni), dbp_ni, dbp),
    mbp = if_else(is.na(mbp) & !is.na(mbp_ni), mbp_ni, mbp)
  )

# mbp again
dt <- dt %>%
  mutate(
    # Calculate mbp from sbp and dbp if missing (mbp = dbp + (sbp - dbp)/3)
    mbp = case_when(
      !is.na(mbp) ~ mbp,  # Keep mbp if available
      !is.na(sbp) & !is.na(dbp) ~ dbp + (sbp - dbp) / 3,  # Calculate
      TRUE ~ NA_real_
    ),

    # Calculate mbp_ni from sbp_ni and dbp_ni if missing
    mbp_ni = case_when(
      !is.na(mbp_ni) ~ mbp_ni,  # Keep mbp_ni if available
      !is.na(sbp_ni) & !is.na(dbp_ni) ~ dbp_ni + (sbp_ni - dbp_ni) / 3,  # Calculate
      TRUE ~ NA_real_
    )
  )

#CHeck
dt %>%
  summarise(
    sbp_missing_after = sum(is.na(sbp)),
    dbp_missing_after = sum(is.na(dbp)),
    mbp_missing_after = sum(is.na(mbp)),
    total_rows = n()
  ) %>%
  pivot_longer(cols = ends_with("after"), names_to = "variable", values_to = "count") %>%
  mutate(percentage = round(count / max(count) * 100, 1))


## fio2-----
#fix some missing values and fio2 var

dt |> dplyr::select(fio2,o2_flow,o2_flow_total,glucose) |> skim()

dt <- dt %>%
  mutate(fio2 = case_when(
 
      #fio2 miscoded as o2_flow
    is.na(fio2) & highflow==0 & noninvasive==0 & invasive==0 & o2_flow > 25 & o2_flow <= 100 ~ o2_flow,
    
    #Reliability of methods to estimate the fraction of inspired oxygen in patients with acute respiratory failure breathing through non-rebreather reservoir bag oxygen mask.Thorax 2020;75:805–807
    is.na(fio2) & highflow==0 & noninvasive==0 & invasive==0 & o2_flow >= 0 & o2_flow <= 25 ~ 21+(o2_flow*3),
    is.na(fio2) & (highflow==1 | noninvasive==1 | invasive==1) & o2_flow >= 0 ~ pmin(21+(o2_flow*3), 100),
    TRUE ~ fio2)) %>%
  mutate(glucose = case_when(glucose >= 500000 ~ NA_real_, TRUE ~ glucose))


# SOFA ----

# 1. Initialize with ICU admission values + Forward fill
dt <- dt %>%
  arrange(stay_id, hr) %>%
  group_by(stay_id) %>%
  mutate(
    # Initialize each component with ICU admission values
    across(c(respiration, coagulation, liver, cardiovascular, cns, renal),
           ~{
             # Get value at hr=0, if not available use first non-NA value, otherwise 0
             baseline_val <- if (any(hr == 0 & !is.na(.))) {
               .[hr == 0 & !is.na(.)][1]
             } else {
               first_non_na <- first(na.omit(.))
               # If vector length is 0 (all NA), initialize with 0
               if (length(first_non_na) == 0) {
                 0
               } else {
                 first_non_na
               }
             }
             # Replace NA with baseline value
             if_else(is.na(.), baseline_val, .)
           })
  ) %>%
  # Then forward fill (just in case)
  fill(respiration, coagulation, liver, cardiovascular, cns, renal,
       .direction = "down") %>%
  # 2. Rolling 24-hour maximum SOFA
  mutate(
    respiration_24h_max = slide_index_dbl(
      respiration, hr,
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
      .before = 23, .after = 0,
      .complete = FALSE
    ),
    coagulation_24h_max = slide_index_dbl(
      coagulation, hr,
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
      .before = 23, .after = 0,
      .complete = FALSE
    ),
    liver_24h_max = slide_index_dbl(
      liver, hr,
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
      .before = 23, .after = 0,
      .complete = FALSE
    ),
    cardiovascular_24h_max = slide_index_dbl(
      cardiovascular, hr,
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
      .before = 23, .after = 0,
      .complete = FALSE
    ),
    cns_24h_max = slide_index_dbl(
      cns, hr,
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
      .before = 23, .after = 0,
      .complete = FALSE
    ),
    renal_24h_max = slide_index_dbl(
      renal, hr,
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE),
      .before = 23, .after = 0,
      .complete = FALSE
    ),
    sofa_24h_max = if_else(
      is.na(respiration_24h_max) | is.na(coagulation_24h_max) |
        is.na(liver_24h_max) | is.na(cardiovascular_24h_max) |
        is.na(cns_24h_max) | is.na(renal_24h_max),
      NA_real_,
      respiration_24h_max + coagulation_24h_max + liver_24h_max +
        cardiovascular_24h_max + cns_24h_max + renal_24h_max
    )
  ) %>%
  ungroup()

# 3. Calculate Baseline SOFA (maximum value during hours 0-23)
baseline_sofa_df <- dt %>%
  filter(time >= 0 & time <= 23) %>%
  group_by(stay_id) %>%
  summarise(
    respiration_baseline = max(respiration, na.rm = TRUE),
    coagulation_baseline = max(coagulation, na.rm = TRUE),
    liver_baseline = max(liver, na.rm = TRUE),
    cardiovascular_baseline = max(cardiovascular, na.rm = TRUE),
    cns_baseline = max(cns, na.rm = TRUE),
    renal_baseline = max(renal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # Convert -Inf to NA (when all values are NA)
    dplyr::across(ends_with("_baseline"), ~if_else(is.infinite(.), NA_real_, .)),

    # Baseline SOFA total
    baseline_sofa = if_else(
      is.na(respiration_baseline) | is.na(coagulation_baseline) |
        is.na(liver_baseline) | is.na(cardiovascular_baseline) |
        is.na(cns_baseline) | is.na(renal_baseline),
      NA_real_,
      respiration_baseline + coagulation_baseline + liver_baseline +
        cardiovascular_baseline + cns_baseline + renal_baseline
    )
  )

# 4. Join to original data
dt <- dt %>%
  left_join(baseline_sofa_df, by = "stay_id")

# Check
dt %>%
  dplyr::select(stay_id, hr, 
         respiration, respiration_24h_max, respiration_baseline,
         sofa_24h_max, baseline_sofa) %>%
  filter(stay_id == first(stay_id)) %>%
  head(30)


dt |> filter(time==0) |> dplyr::select(baseline_sofa,
                                respiration,
                                cardiovascular,
                                renal,
                                liver,
                                coagulation,
                                cns) |> skim()


# Carry forward variable values-----

dt <- dt %>%
  group_by(subject_id) %>%
  fill(heart_rate,resp_rate,sbp,dbp,mbp,spo2,temperature,fio2,
       glucose,gcs,ph,po2,pco2) |> 
  ungroup()



# Missing check

dt |> dplyr::select(stay_id,subject_id,hadm_id,hr,time,
             intime,outtime,admittime,dischtime,deathtime,
             anchor_age) |> skim()

dt |> dplyr::select(first_hfnc_start,hr_at_hfnc_start,
             deathtime_inhosp,dod,death_datetime,death_location,
             time_to_death,death_30day,
             hr_at_death,hr_at_icu_discharge,hr_at_hosp_discharge) |> skim()

dt |> dplyr::select(death_event,icu_discharge_event,hosp_discharge_event,
             gender,anchor_year,anchor_year_group) |> skim()

dt |> dplyr::select(gender,female,insurance_grouped,race_grouped,marital_status,
             language_grouped,admission_type_grouped,admission_location_grouped,
             discharge_location) |> skim()

dt |> dplyr::select(first_careunit,last_careunit,los,elixhauser_vanwalraven) |> skim()

dt |> dplyr::select(discharge_outcome,invasive,noninvasive,highflow) |> skim()

dt |> filter(time>=-1)|> dplyr::select(sbp,dbp,mbp,resp_rate,temperature,spo2,glucose,gcs,fio2,o2_flow) |> skim()

dt |> filter(time>=-1)|> dplyr::select(
  ph,po2,pco2,
  vasopressor,crrt,sofa_24hours,sepsis3) |> skim()

# Code status -----

#All patients full
dt |> filter(time==0 & fullcode == 0) 

dt |> filter(time==0 & fullcode == 1 & dnr ==1) 

dt |> filter(fullcode == 1 & (dni==1 | dnr ==1 | dnr_dni==1 | cmo ==1)) |> dplyr::select(stay_id,time,invasive,fullcode,dni,dnr,dnr_dni,cmo,death_event) 

#dt |> filter(fullcode == 1 & dnr ==1) |> dplyr::select(stay_id,time,invasive,fullcode,dni,dnr,dnr_dni,cmo,death_event) |> view()


dt |> filter(dni == 1 & dnr ==1 & dnr_dni==1) |> dplyr::select(stay_id,time,invasive,fullcode,dni,dnr,dnr_dni,cmo,death_event) 

# Fix: Change fullcode to 0 when fullcode==1 & dnr==1
dt$fullcode[dt$fullcode == 1 & dt$dnr == 1] <- 0

# Fix 2: Change dni and dnr to 0 when dni==1 & dnr==1 & dnr_dni==1
# First get the indices of rows that match the condition
idx_to_fix <- which(dt$dni == 1 & dt$dnr == 1 & dt$dnr_dni == 1)
# Change dni and dnr to 0 simultaneously for the matching rows
dt$dni[idx_to_fix] <- 0
dt$dnr[idx_to_fix] <- 0

dt<-dt |> mutate(code_sum=dni+dnr+dnr_dni+cmo+fullcode)

dt |> filter(cmo == 1 ) |> dplyr::select(stay_id) |> unique()

# Zero expected
dt |> filter(code_sum > 2 ) |> dplyr::select(stay_id) |> unique()

# code sum >1
dt |> filter(code_sum > 1 ) |> dplyr::select(stay_id) |> unique()

# All sum>1 are cmo + alpha
#dt |> filter(stay_id %in% c(33376820, 34250110, 34478728, 35234112, 35373167, 35393635, 36375397, 37881237) & code_sum>1 ) |> dplyr::select(stay_id,time,invasive,fullcode,dni,dnr,dnr_dni,cmo,death_event) |> view()

#dt |> filter(cmo == 1 ) |> dplyr::select(stay_id,time,invasive,fullcode,dni,dnr,dnr_dni,cmo,death_event) |> view()

dt |> filter(fullcode == 1 & (dni==1 | dnr ==1 | dnr_dni==1 | cmo ==1)) |> dplyr::select(stay_id,time,invasive,fullcode,dni,dnr,dnr_dni,cmo,death_event) 


# Elig flag --------
# Step 1: Create elig_1 at patient level
patient_elig_1 <- dt %>%
  group_by(subject_id) %>%
  slice(1) %>%  # Get first row for each patient
  ungroup() %>%
  mutate(
    elig_1 = if_else(!is.na(hr_at_hfnc_start) & hr_at_hfnc_start <= 168, 1, 0)
  ) %>%
  dplyr::select(subject_id, elig_1)

# Step 2: Join to original dt
dt <- dt %>%
  left_join(patient_elig_1, by = "subject_id")

# Check
cat("=== elig_1 distribution ===\n")
dt %>%
  distinct(subject_id, elig_1) %>%
  count(elig_1) %>%
  mutate(
    label = case_when(
      elig_1 == 1 ~ "Eligible (HFNC ≤ 168hr)",
      elig_1 == 0 ~ "Not eligible (HFNC > 168hr or NA)",
      TRUE ~ "Other"
    ),
    percentage = sprintf("%.1f%%", n / sum(n) * 100)
  )

# Check relationship between hr_at_hfnc_start and elig_1
cat("\n=== Relationship between hr_at_hfnc_start and elig_1 ===\n")
dt %>%
  distinct(subject_id, hr_at_hfnc_start, elig_1) %>%
  group_by(elig_1) %>%
  summarise(
    n = n(),
    mean_hr = mean(hr_at_hfnc_start, na.rm = TRUE),
    median_hr = median(hr_at_hfnc_start, na.rm = TRUE),
    min_hr = min(hr_at_hfnc_start, na.rm = TRUE),
    max_hr = max(hr_at_hfnc_start, na.rm = TRUE),
    .groups = "drop"
  )

# Check data structure (verify all rows have same elig_1 for each patient)
cat("\n=== Check elig_1 in long format (first 20 rows) ===\n")
dt %>%
  dplyr::select(subject_id, time, hr_at_hfnc_start, elig_1) %>%
  head(20)

# Check if elig_1 is consistent for each patient
consistency_check <- dt %>%
  group_by(subject_id) %>%
  summarise(
    n_rows = n(),
    n_unique_elig = n_distinct(elig_1),
    .groups = "drop"
  ) %>%
  filter(n_unique_elig > 1)

if(nrow(consistency_check) == 0) {
  cat("\n OK: All rows have the same elig_1 value for each patient\n")
} else {
  cat("\n Warning: Some patients have different elig_1 values across rows\n")
  print(consistency_check)
}

# Example of filtering by elig_1
cat("\n=== Number of patients and observations with elig_1 = 1 ===\n")
dt %>%
  filter(elig_1 == 1) %>%
  summarise(
    n_patients = n_distinct(subject_id),
    n_observations = n()
  )

## elig_2
## Patients intubated within 1 hour

# Step 1: Create elig_2 at patient level
patient_elig_2 <- dt %>%
  filter(time == 0) |> 
  mutate(
    elig_2 = case_when(elig_1 == 1 & invasive == 0 ~ 1,
                       TRUE ~ 0)
  ) %>%
  dplyr::select(subject_id, elig_2)

# Step 2: Join to original dt
dt <- dt %>%
  left_join(patient_elig_2, by = "subject_id")

# Check
dt |> filter(time == 0 & invasive==1) |> 
  dplyr::select(subject_id,elig_1,elig_2)

dt |> filter(time == 0 & invasive==1) |> nrow()
dt |>   dplyr::select(subject_id)|> unique()|> nrow()
dt |> filter(elig_2==1) |>   dplyr::select(subject_id)|> unique() |> nrow()
dt %>%
  filter(elig_1 == 1) |>   dplyr::select(subject_id)|> unique() |> nrow()

# select valiables----

#remove unnecessary variables

dt$anchor_year_group<-factor(dt$anchor_year_group,
                                levels = c(
                                  "2008 - 2010",
                                  "2011 - 2013",
                                  "2014 - 2016",
                                  "2017 - 2019",
                                  "2020 - 2022")
                                )



df<-dt

dt<-df |> dplyr::select(
  subject_id,elig_1,elig_2,hr,time,intime,outtime,admittime,dischtime,
  hr_at_hfnc_start,
  age,anchor_year_group,female,
  medicare_medicaid,
  language_grouped,
  Married,
  race_grouped,admission_location_grouped,
  
  # ventilation
  invasive,noninvasive,highflow,
  
  #vital
  heart_rate,sbp,dbp,mbp,resp_rate,temperature,spo2,
  fio2,gcs,
  
  #labo
  ph,po2,pco2,glucose,lactate,
  
  # score
  baseline_sofa,
  sofa_24hours,
  elixhauser_vanwalraven,
  ## a comorbidity summary score based on 30 comorbidities PMID: 25769057

    # condition 
  sepsis3,
  
  # drug
  antibiotic_use,
  norepinephrine,epinephrine,dopamine,phenylephrine,
  dobutamine,vasopressor,
  
  # treat
  crrt,
  
  # code
  fullcode,dnr,dni,dnr_dni,cmo,
  
  #outcome
  icu_discharge_event,# When discharge, then 1
  hosp_discharge_event,# When discharge, then 1
  los,
  death_datetime, #including dod (00:00)
  death_30day, # If death, all row is 1
  death_event, # When death, then 1)
  death_location)


saveRDS(dt, file = "./data/ROX/preIM_cohort.rds")


# Create table 1
# temporary lag
#lag all time-varying covariates since modelling will depend on past covariate values
lag_vars <- c("heart_rate", "resp_rate", "sbp", "dbp", "mbp", 
              "spo2", "temperature", "fio2", "glucose", "gcs", 
              "ph", "po2", "pco2","antibiotic_use",
              "vasopressor","crrt","sofa_24hours")

cohort <- dt %>%  
  group_by(subject_id) %>%
  mutate(dplyr::across(all_of(lag_vars), 
                ~{
                  lagged <- lag(.)
                  coalesce(lagged, .)
                })) %>%
  ungroup()

# Check: first row should have original value, not NA
cohort %>%
  group_by(subject_id) %>%
  slice(1:3) %>%
  dplyr::select(subject_id, hr,time,heart_rate, sbp, gcs)

cohort <- cohort %>%
  mutate(
    # ROX Index = (SpO2 / FiO2) / Respiratory Rate
    # Divide FiO2 by 100 since it's in percentage
    rox = (spo2 / (fio2 / 100)) / resp_rate
  )

# ROX index distribution
summary(cohort$rox)



# GCS category
cohort <-cohort |> 
  mutate(gcs_category = case_when(
    gcs >= 13 ~ "13 to 15",
    gcs >= 9 ~ "9 to 12",
    gcs < 9 ~ "< 9",
    .default = NA_character_
  )) |> 
  mutate(gcs_category = factor(gcs_category, 
                               levels = c("13 to 15", 
                                          "9 to 12", "< 9")))

  
# Table 1----------

baseline <-cohort |> filter(time==0)
baseline |> nrow() # 1755

baseline<-baseline |> filter(elig_1==1)
baseline |> nrow() # 1651


table1 <- baseline %>%
  dplyr::select(
    # Demographics
    age, female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    
    # Vital signs (baseline)
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    
    # Laboratory values (baseline)
    ph, po2, pco2, #glucose, lactate,
    
    # Severity scores
    baseline_sofa, elixhauser_vanwalraven,
    
    # Condition
    antibiotic_use,
    #sepsis3,
    
    # Vasopressor use
    vasopressor,
    
    # Treatment
    crrt,
    
    # hr
    hr_at_hfnc_start
    
  ) %>%
  tbl_summary(missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      # Demographics
      age ~ "Age, years",
      female ~ "Female",
      medicare_medicaid ~ "Medicare/Medicaid",
      language_grouped ~ "Language",
      Married ~ "Married",
      race_grouped ~ "Race/Ethnicity",
      
      # Vital signs
      heart_rate ~ "Heart rate, bpm",
      sbp ~ "Systolic blood pressure, mmHg",
      dbp ~ "Diastolic blood pressure, mmHg",
      mbp ~ "Mean blood pressure, mmHg",
      resp_rate ~ "Respiratory rate, /min",
      temperature ~ "Temperature, °C",
      spo2 ~ "SpO2, %",
      fio2 ~ "FiO2, %",
      rox  ~ "ROX index",
      gcs_category ~ "Glasgow Coma Scale",
      
      # Laboratory
      ph ~ "pH",
      po2 ~ "PaO2, mmHg",
      pco2 ~ "PaCO2, mmHg",
      #glucose ~ "Glucose, mg/dL",
      #lactate ~ "Lactate, mmol/L",
      
      # Scores
      baseline_sofa ~ "SOFA score, first 24 hours",
      elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
      
      # Condition
      antibiotic_use ~ "Antibiotic use",
      
      # Treatment
      vasopressor ~ "Vasopressor use",
      crrt ~ "CRRT",
      
      # hour
      hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
    ),
    missing_text = "Missing"
  ) %>%
  add_overall() %>% 
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**Missing**") %>%
  bold_labels() %>%
  modify_header(
    stat_0 ~ "**Overall**\nN = {N}",
    n ~ "**Missing, N (%)**"
  ) %>% 
  modify_caption("**Table 1. Baseline Characteristics of Study Population**") %>%
  as_flex_table() %>% 
  flextable::fontsize(size = 9, part = "body") %>%
  flextable::fontsize(size = 10, part = "header") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::align(align = "left", j = 1, part = "all") %>%
  flextable::align(align = "center", j = 2:3, part = "all") %>%
  autofit() %>% 
  set_table_properties(
    width = 1,
    layout = "autofit"
  ) %>% 
  flextable::add_footer_lines(
    values = c("SOFA: Sequential Organ Failure Assessment.Values represent maximum scores observed during the first 24 hours of HFNC initiation.; SpO2: Peripheral oxygen saturation; FiO2: Fraction of inspired oxygen; CRRT: Continuous renal replacement therapy.", "* Owing to the cloning step in the clone-censor-weight approach, patients' characteristics are identical at baseline for usual care and ROX-guided intubation strategy groups. For detailed explanation of the clone-censor-weight method, see methods section and supplementary methods."
    )
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::italic(part = "footer")

  
# Display table
table1

# Export to Word
table1 |>
  flextable::save_as_docx(path = "./outputs/ROX/table1.docx")

# Export to Word
table1 |>
  flextable::save_as_docx(path = "C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/Table2_patient_characteristics_singlecolumn.docx")


# Export to PowerPoint
#flextable::save_as_pptx(table1, path = "./outputs/ROX/table1.pptx")

# Table 1 clone------
# Clone table ------

# Step 1: Create Overall table (original baseline, N=1,651)
tbl_overall <- baseline %>%
  dplyr::select(
    age, female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    ph, po2, pco2,
    baseline_sofa, elixhauser_vanwalraven,
    antibiotic_use,
    #sepsis3,
    vasopressor, crrt,
    hr_at_hfnc_start
  ) %>%
  tbl_summary(
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age, years",
      female ~ "Female",
      medicare_medicaid ~ "Medicare/Medicaid",
      language_grouped ~ "Language",
      Married ~ "Married",
      race_grouped ~ "Race/Ethnicity",
      heart_rate ~ "Heart rate, bpm",
      sbp ~ "Systolic blood pressure, mmHg",
      dbp ~ "Diastolic blood pressure, mmHg",
      mbp ~ "Mean blood pressure, mmHg",
      resp_rate ~ "Respiratory rate, /min",
      temperature ~ "Temperature, °C",
      spo2 ~ "SpO2, %",
      fio2 ~ "FiO2, %",
      rox  ~ "ROX index",
      gcs_category ~ "Glasgow Coma Scale",
      ph ~ "pH",
      po2 ~ "PaO2, mmHg",
      pco2 ~ "PaCO2, mmHg",
      baseline_sofa ~ "SOFA score, first 24 hours",
      elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
      antibiotic_use ~ "Antibiotic use",
      #sepsis3 ~ "Sepsis-3",
      vasopressor ~ "Vasopressor use",
      crrt ~ "CRRT",
      hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
    )
  ) %>%
  modify_header(stat_0 ~ "**Overall**\n(N = {N})")

# Step 2: Create cloned table for treatment strategies (each arm N=1,651)
baseline_clone <- baseline %>%
  mutate(arm = "Usual Care") %>%
  bind_rows(
    baseline %>% mutate(arm = "ROX-guided Intubation (ROX <3.85)")
  ) %>%
  mutate(
    arm = factor(arm, 
                 levels = c("Usual Care", "ROX-guided Intubation (ROX <3.85)"))
  )

tbl_arms <- baseline_clone %>%
  dplyr::select(
    arm,
    age, female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    ph, po2, pco2,
    baseline_sofa, elixhauser_vanwalraven,
    antibiotic_use,
    #sepsis3,
    vasopressor, crrt,
    hr_at_hfnc_start
  ) %>%
  tbl_summary(
    by = arm,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age, years",
      female ~ "Female",
      medicare_medicaid ~ "Medicare/Medicaid",
      language_grouped ~ "Language",
      Married ~ "Married",
      race_grouped ~ "Race/Ethnicity",
      heart_rate ~ "Heart rate, bpm",
      sbp ~ "Systolic blood pressure, mmHg",
      dbp ~ "Diastolic blood pressure, mmHg",
      mbp ~ "Mean blood pressure, mmHg",
      resp_rate ~ "Respiratory rate, /min",
      temperature ~ "Temperature, °C",
      spo2 ~ "SpO2, %",
      fio2 ~ "FiO2, %",
      rox  ~ "ROX index",
      gcs_category ~ "Glasgow Coma Scale",
      ph ~ "pH",
      po2 ~ "PaO2, mmHg",
      pco2 ~ "PaCO2, mmHg",
      baseline_sofa ~ "SOFA score, first 24 hours",
      elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
      antibiotic_use ~ "Antibiotic use",
      #      sepsis3 ~ "Sepsis-3",
      vasopressor ~ "Vasopressor use",
      crrt ~ "CRRT",
      hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
    )
  ) %>%
  modify_header(
    stat_1 ~ "**Usual Care**\n(N = {n})*",
    stat_2 ~ "**ROX-guided Intubation**\n(ROX <3.85) (N = {n})*"
  )

# Step 3: Merge tables
table1_clone <- tbl_merge(
  tbls = list(tbl_overall, tbl_arms),
  tab_spanner = c(NA, "**Treatment Strategy**")
) %>%
  bold_labels() %>%
  as_flex_table() %>% 
  flextable::fontsize(size = 9, part = "body") %>%
  flextable::fontsize(size = 10, part = "header") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::align(align = "left", j = 1, part = "all") %>%
  flextable::align(align = "center", j = 2:4, part = "all") %>%
  autofit() %>% 
  set_table_properties(
    width = 1,
    layout = "autofit"
  ) %>% 
  flextable::add_footer_lines(
    values = c(
      "SOFA: Sequential Organ Failure Assessment (maximum scores during first 24 hours of HFNC); SpO2: Peripheral oxygen saturation; FiO2: Fraction of inspired oxygen; CRRT: Continuous renal replacement therapy; ROX index: (SpO2/FiO2)/respiratory rate.",
      "* Owing to the cloning step in the clone-censor-weight approach, patients' characteristics are identical at baseline for usual care and ROX-guided intubation strategy groups. For detailed explanation of the clone-censor-weight method, see Methods section and Supplementary Methods."
    )
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::italic(part = "footer")

# Display table
table1_clone

# Export to Word
table1_clone |>
  flextable::save_as_docx(path = "./outputs/ROX/table1_clone.docx")

# Export to Word
table1_clone |>
  flextable::save_as_docx(path = "C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/Table2_patient_characteristics_clone.docx")

# All clone table 1 -------
# Step 1: Create Overall table (original baseline, N=1,651)
tbl_overall <- baseline %>%
  dplyr::select(
    # Demographics
    age, female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    # Vital signs
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    # Laboratory
    ph, po2, pco2,
    # Severity scores
    baseline_sofa, elixhauser_vanwalraven,
    # Condition
    #sepsis3,
    antibiotic_use,
    # Treatment
    vasopressor, crrt,
    # hour
    hr_at_hfnc_start
  ) %>%
  tbl_summary(
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age, years",
      female ~ "Female",
      medicare_medicaid ~ "Medicare/Medicaid",
      language_grouped ~ "Language",
      Married ~ "Married",
      race_grouped ~ "Race/Ethnicity",
      heart_rate ~ "Heart rate, bpm",
      sbp ~ "Systolic blood pressure, mmHg",
      dbp ~ "Diastolic blood pressure, mmHg",
      mbp ~ "Mean blood pressure, mmHg",
      resp_rate ~ "Respiratory rate, /min",
      temperature ~ "Temperature, °C",
      spo2 ~ "SpO2, %",
      fio2 ~ "FiO2, %",
      rox  ~ "ROX index",
      gcs_category ~ "Glasgow Coma Scale",
      ph ~ "pH",
      po2 ~ "PaO2, mmHg",
      pco2 ~ "PaCO2, mmHg",
      baseline_sofa ~ "SOFA score, first 24 hours",
      elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
      antibiotic_use ~ "Antibiotic use",
      #sepsis3 ~ "Sepsis-3",
      vasopressor ~ "Vasopressor use",
      crrt ~ "CRRT",
      hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
    )
  ) %>%
  modify_header(stat_0 ~ "**Overall**\n(N = {N})")

# Step 2: Create cloned table for treatment strategies (each arm N=1,651)
baseline_clone_full <- baseline %>%
  mutate(arm = "Usual Care") %>%
  bind_rows(
    baseline %>% mutate(arm = "ROX <3.85"),
    baseline %>% mutate(arm = "ROX <4.88"),
    baseline %>% mutate(arm = "Time-varying")
  ) %>%
  mutate(
    arm = factor(arm, 
                 levels = c("Usual Care", "ROX <3.85", "ROX <4.88", "Time-varying"))
  )

tbl_arms <- baseline_clone_full %>%
  dplyr::select(
    arm,
    age, female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    ph, po2, pco2,
    baseline_sofa, elixhauser_vanwalraven,
    antibiotic_use,
    #    sepsis3,
    vasopressor, crrt,
    hr_at_hfnc_start
  ) %>%
  tbl_summary(
    by = arm,
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age, years",
      female ~ "Female",
      medicare_medicaid ~ "Medicare/Medicaid",
      language_grouped ~ "Language",
      Married ~ "Married",
      race_grouped ~ "Race/Ethnicity",
      heart_rate ~ "Heart rate, bpm",
      sbp ~ "Systolic blood pressure, mmHg",
      dbp ~ "Diastolic blood pressure, mmHg",
      mbp ~ "Mean blood pressure, mmHg",
      resp_rate ~ "Respiratory rate, /min",
      temperature ~ "Temperature, °C",
      spo2 ~ "SpO2, %",
      fio2 ~ "FiO2, %",
      rox  ~ "ROX index",
      gcs_category ~ "Glasgow Coma Scale",
      ph ~ "pH",
      po2 ~ "PaO2, mmHg",
      pco2 ~ "PaCO2, mmHg",
      baseline_sofa ~ "SOFA score, first 24 hours",
      elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
      antibiotic_use ~ "Antibiotic use",
      #      sepsis3 ~ "Sepsis-3",
      vasopressor ~ "Vasopressor use",
      crrt ~ "CRRT",
      hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
    )
  ) %>%
  modify_header(
    stat_1 ~ "**Usual Care**\n(N = {n})*",
    stat_2 ~ "**ROX <3.85**\n(N = {n})*",
    stat_3 ~ "**ROX <4.88**\n(N = {n})*",
    stat_4 ~ "**Time-varying**\n(N = {n})*"
  )

# Step 3: Create Missing table (from original baseline)
tbl_missing <- baseline %>%
  dplyr::select(
    age, female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    ph, po2, pco2,
    baseline_sofa, elixhauser_vanwalraven,
    antibiotic_use,
    #sepsis3,
    vasopressor, crrt,
    hr_at_hfnc_start
  ) %>%
  tbl_summary(
    missing = "no",
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age ~ "Age, years",
      female ~ "Female",
      medicare_medicaid ~ "Medicare/Medicaid",
      language_grouped ~ "Language",
      Married ~ "Married",
      race_grouped ~ "Race/Ethnicity",
      heart_rate ~ "Heart rate, bpm",
      sbp ~ "Systolic blood pressure, mmHg",
      dbp ~ "Diastolic blood pressure, mmHg",
      mbp ~ "Mean blood pressure, mmHg",
      resp_rate ~ "Respiratory rate, /min",
      temperature ~ "Temperature, °C",
      spo2 ~ "SpO2, %",
      fio2 ~ "FiO2, %",
      rox  ~ "ROX index",
      gcs_category ~ "Glasgow Coma Scale",
      ph ~ "pH",
      po2 ~ "PaO2, mmHg",
      pco2 ~ "PaCO2, mmHg",
      baseline_sofa ~ "SOFA score, first 24 hours",
      elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
      antibiotic_use ~ "Antibiotic use",
      #sepsis3 ~ "Sepsis-3",
      vasopressor ~ "Vasopressor use",
      crrt ~ "CRRT",
      hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
    )
  ) %>%
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**Missing**\nN (%)") %>%
  modify_column_hide(stat_0)  # Hide the stat column, keep only N

# Step 4: Merge tables
stable1_clone <- tbl_merge(
  tbls = list(tbl_overall, tbl_arms, tbl_missing),
  tab_spanner = c(NA, "**Treatment Strategy**", NA)
) %>%
  bold_labels() %>%
  as_flex_table() %>% 
  flextable::fontsize(size = 8, part = "body") %>%
  flextable::fontsize(size = 9, part = "header") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::align(align = "left", j = 1, part = "all") %>%
  flextable::align(align = "center", j = 2:7, part = "all") %>%
  autofit() %>% 
  set_table_properties(
    width = 1,
    layout = "autofit"
  ) %>% 
  flextable::add_footer_lines(
    values = c(
      "SOFA: Sequential Organ Failure Assessment (maximum scores during first 24 hours of HFNC); SpO2: Peripheral oxygen saturation; FiO2: Fraction of inspired oxygen; CRRT: Continuous renal replacement therapy; ROX index: (SpO2/FiO2)/respiratory rate.",
      "* Owing to the cloning step in the clone-censor-weight approach, patients' characteristics are identical at baseline for all treatment strategy groups. For detailed explanation of the clone-censor-weight method, see Methods section and Supplementary Methods.",
      "Time-varying strategy: ROX <2.85 (hours 1 to 5), ROX <3.47 (hours 6 to 11), ROX <3.85 (≥12 hours)."
    )
  ) %>%
  flextable::fontsize(size = 7, part = "footer") %>%
  flextable::italic(part = "footer")

# Display table
stable1_clone

# Save to Word document
# Export to Word
stable1_clone|>
  flextable::save_as_docx(path = "./outputs/ROX/S.Table 3_patient_characteristic_all_clone.docx")

stable1_clone |> 
  flextable::save_as_docx(path = "C:/Users/ryohe/Dropbox (個人)/Research/2025ROX_ImpactAnalysis/Manuscript/Table and figures/S.Table 3_patient_characteristic_all_clone.docx")


## table 1 by age group  ----

table1_age <- baseline %>%
  dplyr::select(
    # Demographics
    age, anchor_year_group,female, medicare_medicaid, language_grouped, 
    Married, race_grouped,
    
    # Vital signs (baseline)
    heart_rate, sbp, dbp, mbp, resp_rate, 
    temperature, spo2, fio2, rox, gcs_category,
    
    # Laboratory values (baseline)
    ph, po2, pco2, #glucose, lactate,
    
    # Severity scores
    baseline_sofa, elixhauser_vanwalraven,
    
    # Condition
    antibiotic_use,
#    sepsis3,
    
    # Vasopressor use
    vasopressor,
    
    # Treatment
    crrt,
    
    # hr
    hr_at_hfnc_start
    
  ) %>%
  tbl_summary(by=anchor_year_group,
              missing = "no",
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75})",
                all_categorical() ~ "{n} ({p}%)"
              ),
              digits = all_continuous() ~ 1,
              label = list(
                # Demographics
                age ~ "Age, years",
                female ~ "Female",
                medicare_medicaid ~ "Medicare/Medicaid",
                language_grouped ~ "Language",
                Married ~ "Married",
                race_grouped ~ "Race/Ethnicity",
                
                # Vital signs
                heart_rate ~ "Heart rate, bpm",
                sbp ~ "Systolic blood pressure, mmHg",
                dbp ~ "Diastolic blood pressure, mmHg",
                mbp ~ "Mean blood pressure, mmHg",
                resp_rate ~ "Respiratory rate, /min",
                temperature ~ "Temperature, °C",
                spo2 ~ "SpO2, %",
                fio2 ~ "FiO2, %",
                rox  ~ "ROX index",
                gcs_category ~ "Glasgow Coma Scale",
                
                # Laboratory
                ph ~ "pH",
                po2 ~ "PaO2, mmHg",
                pco2 ~ "PaCO2, mmHg",
                #glucose ~ "Glucose, mg/dL",
                #lactate ~ "Lactate, mmol/L",
                
                # Scores
                baseline_sofa ~ "SOFA score, first 24 hours",
                elixhauser_vanwalraven ~ "Elixhauser comorbidity score",
                
                # Condition
                antibiotic_use ~ "Antibiotic use",
                
                #sepsis3 ~ "Sepsis-3",
                
                # Treatment
                vasopressor ~ "Vasopressor use",
                crrt ~ "CRRT",
                
                # hour
                hr_at_hfnc_start ~ "Time from ICU admission to HFNC initiation, hours"
              ),
              missing_text = "Missing"
  ) %>%
  add_overall() %>% 
  add_n(statistic = "{N_miss} ({p_miss}%)", 
        col_label = "**Missing**") %>%
  bold_labels() %>%
  modify_header(
    stat_0 ~ "**Overall**\nN = {N}",
    n ~ "**Missing, N (%)**"
  ) %>% 
  #modify_caption("**Table 1. Baseline Characteristics of Study Population**") %>%
  as_flex_table() %>% 
  flextable::fontsize(size = 9, part = "body") %>%
  flextable::fontsize(size = 10, part = "header") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  flextable::align(align = "left", j = 1, part = "all") %>%
  flextable::align(align = "center", j = 2:3, part = "all") %>%
  autofit() %>% 
  set_table_properties(
    width = 1,
    layout = "autofit"
  ) %>% 
  flextable::add_footer_lines(
    values = c("SOFA: Sequential Organ Failure Assessment.Values represent maximum scores observed during the first 24 hours of HFNC initiation.; SpO2: Peripheral oxygen saturation; FiO2: Fraction of inspired oxygen; CRRT: Continuous renal replacement therapy.", "* Owing to the cloning step in the clone-censor-weight approach, patients' characteristics are identical at baseline for usual care and ROX-guided intubation strategy groups. For detailed explanation of the clone-censor-weight method, see methods section and supplementary methods."
    )
  ) %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  flextable::italic(part = "footer")


# Display table
table1_age

# Export to Word
table1 |>
  flextable::save_as_docx(path = "./outputs/ROX/table1 by age.docx")



# Other outcomes -----
# Calculate Hospital LOS and mortality
baseline <- baseline %>%
  mutate(
    # Hospital LOS
    hospital_los = as.numeric(difftime(dischtime, admittime, units = "days")),

    # Time from HFNC start to death (hours)
    hfnc_to_death_hours = as.numeric(difftime(death_datetime, intime + lubridate::dhours(hr_at_hfnc_start), units = "hours")),

    # 7-day mortality (from HFNC start)
    death_7day = if_else(
      !is.na(hfnc_to_death_hours) & hfnc_to_death_hours <= 168, 1, 0
    ),

    # 14-day mortality (from HFNC start)
    death_14day = if_else(
      !is.na(hfnc_to_death_hours) & hfnc_to_death_hours <= 336, 1, 0
    ),

    # ICU mortality: Death during ICU stay
    icu_mortality = if_else(
      !is.na(death_datetime) & death_datetime >= intime & death_datetime <= outtime,
      1, 0
    ),

    # Hospital mortality: Death during hospitalization
    hospital_mortality = if_else(
      !is.na(death_datetime) & death_datetime >= admittime & death_datetime <= dischtime,
      1, 0
    )
  )


# Select outcome variables and create summary
outcome_table <- baseline %>%
  dplyr::select(los, hospital_los, icu_mortality, hospital_mortality, death_7day, death_14day, death_30day) %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{median} [{p25} - {p75}]",
      all_dichotomous() ~ "{n} ({p}%)"
    ),
    label = list(
      los ~ "ICU LOS, days",
      hospital_los ~ "Hospital LOS, days",
      icu_mortality ~ "ICU mortality",
      hospital_mortality ~ "Hospital mortality",
      death_7day ~ "7-day mortality",
      death_14day ~ "14-day mortality",
      death_30day ~ "30-day mortality"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_dichotomous() ~ c(0, 1)
    )
  ) %>%
  modify_header(label ~ "**Outcome**") %>%
  bold_labels() |> as_flex_table() %>% autofit()

outcome_table

# Export to Word
outcome_table |>
  flextable::save_as_docx(path = "./outputs/ROX/outcometable.docx")


## Outcome table by age -----
outcome_table_age <- baseline %>%
  dplyr::select(anchor_year_group,los, hospital_los, icu_mortality, hospital_mortality, death_30day) %>%
  tbl_summary(by=anchor_year_group,
    statistic = list(
      all_continuous() ~ "{median} [{p25} - {p75}]",
      all_dichotomous() ~ "{n} ({p}%)"
    ),
    label = list(
      los ~ "ICU LOS, days",
      hospital_los ~ "Hospital LOS, days",
      icu_mortality ~ "ICU mortality",
      hospital_mortality ~ "Hospital mortality",
      death_30day ~ "30-day mortality"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_dichotomous() ~ c(0, 1)
    )
  ) %>%
  modify_header(label ~ "**Outcome**") %>%
  bold_labels() |> as_flex_table() %>% autofit()

outcome_table_age

# Export to Word
outcome_table |>
  flextable::save_as_docx(path = "./outputs/ROX/outcometable by age.docx")


print("all code run")
