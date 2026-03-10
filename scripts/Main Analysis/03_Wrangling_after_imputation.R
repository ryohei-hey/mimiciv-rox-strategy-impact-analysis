#setup memory
rm(list = ls())
gc()
gc()

# load library-----
library(pacman)
p_load(tidyverse, here, lubridate, skimr, naniar,gtsummary,
       patchwork, ggsci, ggstream,ggalluvial,networkD3,htmlwidgets)



glimpse(dt)


dt <- readRDS("./data/ROX/imputed_cohort_mice.rds")

skim(dt)
glimpse(dt)


# time triming
dt<-dt |> filter(time>=0)

# Intubation ----
## intubation_start, intubation_start_cum, and intubation ----

dt <- dt %>%
  group_by(subject_id) %>%
  arrange(subject_id, hr) %>%
  mutate(
    # First time invasive == 1 within time > 0 range
    invasive_start = if_else(
      time > 0 & invasive == 1 &
        lag(cummax(if_else(time > 0, invasive, 0)), default = 0) == 0,
      1, 0
    ),
    # Backup original invasive before replacement
    invasive_original = invasive,  # Save original value
    invasive = if_else(
      time == 0 & invasive == 1,
      lag(invasive, n = 1, default = 0),
      invasive),

    # Backup original noninvasive before replacement
    noninvasive_original = noninvasive,  # Save original value
    noninvasive = if_else(
      time == 0 & noninvasive == 1,
      lag(noninvasive, n = 1, default = 0),
      noninvasive),

    # _flag: Whether xxx == 1 at least once after time >= 0
    invasive_flag = if_else(any(time >= 0 & invasive == 1, na.rm = TRUE), 1, 0),
    noninvasive_flag = if_else(any(time >= 0 & noninvasive == 1, na.rm = TRUE), 1, 0),
    highflow_flag = if_else(any(time >= 0 & highflow == 1, na.rm = TRUE), 1, 0),
    
    # cum
    invasive_start_cum = if_else(cumsum(invasive_start) >= 1, 1, 0)
  ) %>%
  ungroup()


## Check ----
# Pattern 1: Patients with invasive == 1 at time == 0 should become 0
dt %>%
  filter(subject_id %in% (dt %>% filter(time == 0, invasive == 1) %>% 
                            pull(subject_id))) %>%
  dplyr::select(subject_id, hr, time, highflow, invasive, invasive_start, invasive_start_cum,
         invasive_flag, 
         death_event) %>%
  filter(time %in% c(-1, 0, 1, 2, 3))

# Pattern 2: Patients who first become invasive == 1 at time > 0 should become 0
dt %>%
  filter(subject_id %in% (dt %>% 
                            group_by(subject_id) %>% 
                            filter(all(time <= 0 | invasive_original == 0 | time <= 0), 
                                   any(time > 0 & invasive_original == 1)) %>% 
                            pull(subject_id) %>% head(3))) %>%
  dplyr::select(subject_id, hr, time, invasive_original, invasive_start) %>%
  filter(time %in% c(-1:5))

# Check if all invasive_start values are 0 at time <= 0
dt %>%
  filter(time <= 0) %>%
  summarise(any_invasive_start = sum(invasive_start))

# Check patients who had invasive==1 at time==0
dt %>%
  filter(subject_id %in% (
    dt %>% filter(invasive_original == 1, time == 0) %>% pull(subject_id) %>% head(5)
  )) %>%
  dplyr::select(subject_id, hr, time, highflow,invasive,invasive_original) %>%
  filter(time %in% c(-2:5))



# Force highflow to 0 if invasive == 1 at time > 0
dt <- dt %>%
  mutate(
    # Set other respiratory support to 0 during invasive ventilation
    highflow = if_else(time > 0 & invasive == 1, 0, highflow),
    noninvasive = if_else(time > 0 & invasive == 1, 0, noninvasive),
    noninvasive = if_else(time > 0 & highflow == 1, 0, noninvasive)
  )


# No more overlapping devices now
dt |> filter(highflow==1 & invasive==1) %>%
  dplyr::select(subject_id, hr, time, highflow,invasive,invasive_original) 

dt |> filter(noninvasive==1 & invasive==1) %>%
  dplyr::select(subject_id, hr, time, highflow,invasive,invasive_original) 

dt |> filter(noninvasive==1 & highflow==1) %>%
  dplyr::select(subject_id, hr, time, highflow,noninvasive,invasive) 


#dt |> filter(noninvasive_flag==1 & noninvasive==1) |> dplyr::select(subject_id,time,highflow,noninvasive,invasive, invasive_original,invasive_flag,death_event) |> view()


# oxygen device trend ----

# 1. Identify deceased patients and get information at time of death
death_info <- dt %>%
  filter(death_event == 1) %>%
  dplyr::select(subject_id, death_time = time, 
         # Retain fixed values (patient attributes)
         age, female, race_grouped, admission_location_grouped,
         intime, outtime, admittime, dischtime,
         baseline_sofa, elixhauser_vanwalraven, sepsis3,
         invasive_flag, vasopressor, los,  # intubation → invasive_flag
         death_datetime, death_30day, death_location)


# 2. Create time points after death (from death_time+1 to 720)
death_extended <- death_info %>%
  group_by(subject_id) %>%
  reframe(
    time = (death_time + 1):720,
    death_time = death_time,
    # Repeat patient attributes
    across(c(age, female, race_grouped, admission_location_grouped,
             intime, outtime, admittime, dischtime,
             baseline_sofa, elixhauser_vanwalraven, sepsis3,
             invasive_flag, vasopressor, los,  # intubation → invasive_flag
             death_datetime, death_30day, death_location), first)
  ) %>%
  mutate(
    # Post-death status
    death_event = 0,  # death_event is 1 only at time of death
    invasive = 0,
    noninvasive = 0,
    highflow = 0,
    # Other vitals are NA
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

# 3. Combine with original dt
dt_extended <- bind_rows(dt, death_extended) %>%
  arrange(subject_id, time)

# 4. Recalculate respiratory_support
dt_extended <- dt_extended %>%
  group_by(subject_id) %>%
  mutate(
    # Identify time of death
    death_occurred = any(death_event == 1),
    death_time_point = if_else(death_occurred,
                               time[which(death_event == 1)[1]],
                               NA_real_),
    # Dead status continues after death
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

# Distribution by time
support_by_time <- dt_extended %>%
  filter(time >= 0, time <= 720) %>%
  group_by(time, respiratory_support) %>%
  summarise(n = n(), .groups = "drop")

## plot ----
# Display in actual numbers
support_by_time %>%
  ggplot(aes(x = time, y = n, fill = respiratory_support)) +
  geom_area() +
  scale_fill_manual(values = c(
    "Normal/No oxygen" = "lightgreen",
    "High-flow oxygen" = "yellow",
    "Non-invasive ventilation" = "orange",
    "Invasive ventilation" = "red",
    "Dead" = "black"
  )) +
  labs(title = "Number of Patients by Respiratory Support Over Time",
       x = "Time from HFNC initiation (hours)",
       y = "Number of patients",
       fill = "Status") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Distribution by time
support_by_time <- dt_extended %>%
  filter(time >= 0, time <= 720) %>%
  group_by(time, respiratory_support) %>%
  summarise(n = n(), .groups = "drop")


# Colorblind-friendly d3 palette
p1<-support_by_time %>%
  ggplot(aes(x = time, y = n, fill = respiratory_support)) +
  geom_area(alpha = 0.8) +
  scale_fill_d3(palette = "category10") +
  labs(x = "Time from HFNC initiation (hours)",
       y = "Number of patients",
       fill = "") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10),
        axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(color = "black"),
        panel.grid.major.y = element_line(color = "gray90"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

p1

# Save plot
ggsave(
  filename = here("outputs/ROX", "S.Fig.StackedArea_respiratory_support.png"),
  plot = p1,
  width = 10,
  height = 6,
  dpi = 300
)


# Bin time and track state transitions (0-72 hours, every 12 hours)
dt_extended %>%
  filter(time >= 0, time <= 72) %>%
  mutate(time_bin = cut(time, 
                        breaks = c(0, 12, 24, 36, 48, 60, 72),
                        labels = c("0-12h", "12-24h", "24-36h", "36-48h", 
                                   "48-60h", "60-72h"),
                        include.lowest = TRUE)) %>%
  group_by(subject_id, time_bin) %>%
  slice(1) %>%
  ungroup() %>%
  count(time_bin, respiratory_support) %>%
  ggplot(aes(x = time_bin, stratum = respiratory_support, 
             alluvium = respiratory_support, y = n, fill = respiratory_support)) +
  geom_flow(alpha = 0.6) +
  geom_stratum() +
  scale_fill_npg() +
  labs(x = "Time period",
       y = "Number of patients",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Display both proportion and actual numbers
plot_data <- dt_extended %>%
  filter(time >= 0, time <= 72) %>%
  mutate(time_bin = cut(time, 
                        breaks = c(0, 12, 24, 36, 48, 60, 72),
                        labels = c("0-12h", "12-24h", "24-36h", "36-48h", 
                                   "48-60h", "60-72h"),
                        include.lowest = TRUE)) %>%
  group_by(subject_id, time_bin) %>%
  slice(1) %>%
  ungroup() %>%
  count(time_bin, respiratory_support) %>%
  group_by(time_bin) %>%
  mutate(
    proportion = n / sum(n) * 100,
    label = sprintf("%.1f%%\n(n=%d)", proportion, n)
  ) %>%
  ungroup()

p1 <- plot_data %>%
  ggplot(aes(x = time_bin, stratum = respiratory_support, 
             alluvium = respiratory_support, y = proportion, 
             fill = respiratory_support)) +
  geom_flow(alpha = 0.6) +
  geom_stratum(alpha = 0.8) +
  scale_fill_npg() +
  labs(x = "Time period",
       y = "Proportion (%)",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

p1

# Save plot
ggsave(
  filename = here("outputs/ROX", "S.Fig.alluvial_respiratory_support.png"),
  plot = p1,
  width = 10,
  height = 6,
  dpi = 300
)

# Code Status change -----

library(networkD3)
library(dplyr)
library(tidyr)

# Step 1: Define code status
dt_code <- dt %>%
  mutate(
    code_status = case_when(
      death_event == 1 ~ "Death",
      dnr_dni == 1 ~ "DNR/DNI",
      dnr == 1 ~ "DNR",
      dni == 1 ~ "DNI",
      fullcode == 1 ~ "Full Code",
      cmo == 1 ~ "CMO",
      TRUE ~ "Unknown"
    )
  )

# Step 2: Identify first death time for each patient
first_death <- dt_code %>%
  filter(code_status == "Death") %>%
  group_by(subject_id) %>%
  summarise(death_time = min(time), .groups = "drop")

# Step 3: Extend data so "Death" status continues until 720 after death
dt_extended <- dt_code %>%
  left_join(first_death, by = "subject_id") %>%
  filter(is.na(death_time) | time <= death_time) %>%
  dplyr::select(subject_id, time, code_status, death_time)

# Add "Death" rows from death time to 720 for deceased patients
death_extensions <- first_death %>%
  filter(death_time < 720) %>%
  rowwise() %>%
  reframe(
    subject_id = subject_id,
    time = seq(death_time + 1, 720, by = 1),
    code_status = "Death",
    death_time = death_time
  )

# Combine original data with extended data
dt_complete <- bind_rows(dt_extended, death_extensions) %>%
  arrange(subject_id, time)

# Step 4: Get status at specific time points (0, 24, 48, 72, 96, 120, 144 hours)
time_points <- c(0, 24, 48, 72, 96, 120, 144)

# Get status at each time point (most recent status before that point)
time_states <- dt_complete %>%
  filter(time <= 144) %>%
  crossing(target_time = time_points) %>%
  filter(time <= target_time) %>%
  group_by(subject_id, target_time) %>%
  arrange(desc(time)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(subject_id, time = target_time, code_status) %>%
  mutate(time_label = paste0(time, "h"))

# Data verification
cat("=== Code status distribution at each time point ===\n")
time_states %>%
  count(time_label, code_status) %>%
  group_by(time_label) %>%
  mutate(
    total = sum(n),
    percentage = sprintf("%.1f%%", n / sum(n) * 100)
  ) %>%
  arrange(time_label, desc(n)) %>%
  print()

# Step 5: Create transitions between consecutive time points
transitions <- time_states %>%
  arrange(subject_id, time) %>%
  group_by(subject_id) %>%
  mutate(
    next_status = lead(code_status),
    next_time = lead(time),
    next_time_label = lead(time_label)
  ) %>%
  filter(!is.na(next_status)) %>%
  ungroup()

# Aggregate link data
links <- transitions %>%
  count(time_label, code_status, next_time_label, next_status, name = "value") %>%
  filter(value > 0) %>%
  mutate(
    source = paste(time_label, code_status, sep = "|||"),
    target = paste(next_time_label, next_status, sep = "|||")
  )

# Create node list
time_order <- c("0h", "24h", "48h", "72h", "96h", "120h", "144h")
status_order <- c("Full Code", "DNI", "DNR", "DNR/DNI", "CMO", "Death", "Unknown")

all_nodes <- unique(c(links$source, links$target))

nodes <- data.frame(
  name = all_nodes
) %>%
  mutate(
    time_label = sub("\\|\\|\\|.*", "", name),
    status = sub(".*\\|\\|\\|", "", name),
    time_num = match(time_label, time_order),
    status_num = match(status, status_order)
  ) %>%
  arrange(time_num, status_num) %>%
  mutate(node_id = row_number() - 1)

# Set colors for nodes
node_colors <- data.frame(
  status = c("Full Code", "DNI", "DNR", "DNR/DNI", "CMO", "Death", "Unknown"),
  color = c("#4DBBD5FF", "#00A087FF", "#E64B35FF", "#F39B7FFF", 
            "#8491B4FF", "#3C5488FF", "#B09C85FF")
)

nodes <- nodes %>%
  left_join(node_colors, by = "status") %>%
  mutate(display_name = paste0(time_label, ": ", status))

# Convert link source and target to node IDs
links_indexed <- links %>%
  left_join(nodes %>% dplyr::select(name, node_id), by = c("source" = "name")) %>%
  rename(source_id = node_id) %>%
  left_join(nodes %>% dplyr::select(name, node_id), by = c("target" = "name")) %>%
  rename(target_id = node_id) %>%
  dplyr::select(source_id, target_id, value)

# Create Sankey diagram
sn <- sankeyNetwork(
  Links = links_indexed,
  Nodes = nodes,
  Source = "source_id",
  Target = "target_id",
  Value = "value",
  NodeID = "display_name",
  NodeGroup = "status",
  colourScale = JS(sprintf(
    'd3.scaleOrdinal()
      .domain([%s])
      .range([%s])',
    paste0('"', node_colors$status, '"', collapse = ", "),
    paste0('"', node_colors$color, '"', collapse = ", ")
  )),
  fontSize = 12,
  nodeWidth = 15,
  nodePadding = 20,
  height = 900,
  width = 1600,
  iterations = 0,
  sinksRight = TRUE
)

sn



# Save as HTML file (preserve interactive figure)
saveWidget(sn, "outputs/ROX/S.Fig.code_status_sankey.html", selfcontained = TRUE)


# Code status and intubation -----

# Step 1: Define code status
dt_code <- dt %>%
  mutate(
    code_status = case_when(
      death_event == 1 ~ "Death",
      dnr_dni == 1 ~ "DNR/DNI",
      dnr == 1 ~ "DNR",
      dni == 1 ~ "DNI",
      fullcode == 1 ~ "Full Code",
      cmo == 1 ~ "CMO",
      TRUE ~ "Unknown"
    )
  )

# Step 2: Identify intubation and death times for each patient (avoid warnings)
patient_events <- dt %>%
  group_by(subject_id) %>%
  summarise(
    intubation_time = {
      intub_times <- time[invasive_start == 1 & !is.na(invasive_start)]
      if(length(intub_times) > 0) min(intub_times) else NA_real_
    },
    death_time = {
      death_times <- time[death_event == 1 & !is.na(death_event)]
      if(length(death_times) > 0) min(death_times) else NA_real_
    },
    .groups = "drop"
  )

# Step 3: Extend data so "Death" status continues until 720 after death
dt_extended <- dt_code %>%
  left_join(patient_events, by = "subject_id") %>%
  filter(is.na(death_time) | time <= death_time) %>%
  dplyr::select(subject_id, time, code_status, intubation_time, death_time)

# Add "Death" rows from death time to 720 for deceased patients
death_extensions <- patient_events %>%
  filter(!is.na(death_time), death_time < 720) %>%
  rowwise() %>%
  reframe(
    subject_id = subject_id,
    time = seq(death_time + 1, 720, by = 1),
    code_status = "Death",
    intubation_time = intubation_time,
    death_time = death_time
  )

# Combine original data with extended data
dt_complete <- bind_rows(dt_extended, death_extensions) %>%
  arrange(subject_id, time)

# Step 4: Get status at specific time points (0, 24, 48, 72, 96, 120, 144 hours)
time_points <- c(0, 24, 48, 72, 96, 120, 144)

# Get status at each time point (code status + intubation status)
time_states_combined <- dt_complete %>%
  filter(time <= 144) %>%
  crossing(target_time = time_points) %>%
  filter(time <= target_time) %>%
  group_by(subject_id, target_time) %>%
  arrange(desc(time)) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(subject_id, time = target_time, code_status, intubation_time, death_time) %>%
  mutate(
    # Determine intubation status
    is_intubated = !is.na(intubation_time) & time >= intubation_time,
    # Combination of code status + intubation status
    combined_status = case_when(
      code_status == "Death" ~ "Death",
      is_intubated ~ paste0(code_status, " (Intubated)"),
      !is_intubated ~ paste0(code_status, " (Non-intubated)")
    ),
    time_label = paste0(time, "h"),
    # Base status for sorting
    base_status = code_status,
    intubation_status = case_when(
      code_status == "Death" ~ "Death",
      is_intubated ~ "Intubated",
      TRUE ~ "Non-intubated"
    )
  )

# Data verification (count including time column)
cat("=== Status distribution at each time point (code status + intubation) ===\n")
time_states_combined %>%
  count(time, time_label, combined_status) %>%
  group_by(time_label) %>%
  mutate(
    total = sum(n),
    percentage = sprintf("%.1f%%", n / sum(n) * 100)
  ) %>%
  arrange(time, desc(n)) %>%
  dplyr::select(-time) %>%
  print(n = 100)

# Step 5: Create transitions between consecutive time points
transitions_combined <- time_states_combined %>%
  arrange(subject_id, time) %>%
  group_by(subject_id) %>%
  mutate(
    next_status = lead(combined_status),
    next_time = lead(time),
    next_time_label = lead(time_label)
  ) %>%
  filter(!is.na(next_status)) %>%
  ungroup()

# Aggregate link data
links_combined <- transitions_combined %>%
  count(time_label, combined_status, next_time_label, next_status, name = "value") %>%
  filter(value > 0) %>%
  mutate(
    source = paste(time_label, combined_status, sep = "|||"),
    target = paste(next_time_label, next_status, sep = "|||")
  )

# Create node list (structured order)
time_order <- c("0h", "24h", "48h", "72h", "96h", "120h", "144h")
# Base status order
base_status_order <- c("Full Code", "DNI", "DNR", "DNR/DNI", "CMO", "Unknown", "Death")
# Intubation status order (within each base status)
intubation_order <- c("Non-intubated", "Intubated", "Death")

all_nodes <- unique(c(links_combined$source, links_combined$target))

nodes_combined <- data.frame(
  name = all_nodes
) %>%
  mutate(
    time_label = sub("\\|\\|\\|.*", "", name),
    combined_status = sub(".*\\|\\|\\|", "", name),
    # Extract base status and intubation status
    base_status = case_when(
      combined_status == "Death" ~ "Death",
      TRUE ~ sub(" \\(.*\\)", "", combined_status)
    ),
    intubation_status = case_when(
      combined_status == "Death" ~ "Death",
      grepl("Intubated\\)", combined_status) ~ "Intubated",
      TRUE ~ "Non-intubated"
    ),
    time_num = match(time_label, time_order),
    base_status_num = match(base_status, base_status_order),
    intubation_num = match(intubation_status, intubation_order)
  ) %>%
  arrange(time_num, base_status_num, intubation_num) %>%
  mutate(node_id = row_number() - 1)

# Set colors for nodes (same base status has similar colors, Intubated is darker)
node_colors <- data.frame(
  combined_status = c(
    "Full Code (Non-intubated)", "Full Code (Intubated)",
    "DNI (Non-intubated)", "DNI (Intubated)",
    "DNR (Non-intubated)", "DNR (Intubated)",
    "DNR/DNI (Non-intubated)", "DNR/DNI (Intubated)",
    "CMO (Non-intubated)", "CMO (Intubated)",
    "Unknown (Non-intubated)", "Unknown (Intubated)",
    "Death"
  ),
  color = c(
    "#80D0E8", "#4DBBD5",  # Full Code: light blue -> dark blue
    "#66C9A8", "#00A087",  # DNI: light green -> dark green
    "#F49999", "#E64B35",  # DNR: light red -> dark red
    "#F9BDAB", "#F39B7F",  # DNR/DNI: light orange -> dark orange
    "#A5A5C5", "#8491B4",  # CMO: light purple -> dark purple
    "#C8B9A0", "#B09C85",  # Unknown: light gray -> dark gray
    "#3C5488"              # Death: dark blue-purple
  )
)

nodes_combined <- nodes_combined %>%
  left_join(node_colors, by = "combined_status") %>%
  mutate(display_name = paste0(time_label, ": ", combined_status))

# Convert link source and target to node IDs
links_combined_indexed <- links_combined %>%
  left_join(nodes_combined %>% dplyr::select(name, node_id), by = c("source" = "name")) %>%
  rename(source_id = node_id) %>%
  left_join(nodes_combined %>% dplyr::select(name, node_id), by = c("target" = "name")) %>%
  rename(target_id = node_id) %>%
  dplyr::select(source_id, target_id, value)

# Create Sankey diagram
sn_combined <- sankeyNetwork(
  Links = links_combined_indexed,
  Nodes = nodes_combined,
  Source = "source_id",
  Target = "target_id",
  Value = "value",
  NodeID = "display_name",
  NodeGroup = "combined_status",
  colourScale = JS(sprintf(
    'd3.scaleOrdinal()
      .domain([%s])
      .range([%s])',
    paste0('"', node_colors$combined_status, '"', collapse = ", "),
    paste0('"', node_colors$color, '"', collapse = ", ")
  )),
  fontSize = 11,
  nodeWidth = 15,
  nodePadding = 15,
  height = 800,
  width = 1600,
  iterations = 0,
  sinksRight = TRUE
)

sn_combined

# Save
saveWidget(sn_combined, "outputs/ROX/code_status_intubation_combined_sankey.html", selfcontained = TRUE)

# Statistics
cat("\n=== Intubation rate by code status (by time point) ===\n")
time_states_combined %>%
  filter(code_status != "Death") %>%
  group_by(time, time_label, base_status) %>%
  summarise(
    total = n(),
    intubated = sum(is_intubated),
    intubation_rate = sprintf("%.1f%%", intubated / total * 100),
    .groups = "drop"
  ) %>%
  arrange(time, base_status) %>%
  dplyr::select(-time) %>%
  print(n = 100)

# Treatment Limitation ------------
# Full Code = 0 (no treatment limitation), others (DNI/DNR/DNR-DNI/CMO) = 1 (treatment limitation)
dt <- dt %>%
  mutate(
    treatment_limitation = if_else(fullcode == 1, 0, 1)
  )



# Lag all time-varying covariates ----
## since modelling will depend on past covariate values

lag_vars <- c("heart_rate", "resp_rate", "sbp", "dbp", "mbp", 
              "spo2", "temperature", "fio2", "glucose", "gcs", 
              "ph", "po2", "pco2",
              "treatment_limitation",
              "vasopressor","crrt","sofa_24hours")

dt <- dt %>% 
  group_by(subject_id) %>%
  mutate(across(all_of(lag_vars), 
                ~{
                  lagged <- lag(.)
                  coalesce(lagged, .)
                })) %>%
  ungroup()

# Check: first row should have original value, not NA
dt %>%
  group_by(subject_id) %>%
  slice(1:3) %>%
  dplyr::select(subject_id, hr,time,heart_rate, sbp, gcs,treatment_limitation)


# GCS category

dt <-dt |> 
  mutate(gcs_category = case_when(
    gcs >= 13 ~ "13 to 15",
    gcs >= 9 ~ "9 to 12",
    gcs < 9 ~ "< 9",
    .default = NA_character_
  )) |> 
  mutate(gcs_category = factor(gcs_category, 
                               levels = c("13 to 15", 
                                          "9 to 12", "< 9")))


barplot(table(dt$gcs_category))


## ROX ----

dt |> dplyr::select(spo2,fio2,resp_rate) |> skim()

dt <- dt %>%
  mutate(
    # ROX Index = (SpO2 / FiO2) / Respiratory Rate
    # Divide FiO2 by 100 since it's in percentage
    rox = (spo2 / (fio2 / 100)) / resp_rate
  )

# ROX index distribution
summary(dt$rox)

# Manual calculation check with sample data
dt %>%
  filter(!is.na(rox)) %>%
  dplyr::select(subject_id, time, spo2, fio2, resp_rate, rox) %>%
  head(10)

# Check missing values
dt %>%
  summarise(
    n_total = n(),
    n_rox_missing = sum(is.na(rox)),
    pct_missing = round(100 * n_rox_missing / n_total, 1)
  )

# ROX index histogram
dt %>%
  filter(!is.na(rox), rox < 50) %>%  # 外れ値除外
  ggplot(aes(x = rox)) +
  geom_histogram(bins = 50) +
  labs(title = "ROX Index Distribution",
       x = "ROX Index",
       y = "Count")

dt |> filter(rox>50) |> dplyr::select(spo2,fio2,resp_rate,rox)


# Extract data at time points where invasive_start == 1
invasive_start_data <- dt %>%
  filter(invasive_start == 1)

summary(invasive_start_data$rox)
summary(invasive_start_data$time)

# 1. ROX distribution
p_rox <- invasive_start_data %>%
  ggplot(aes(x = rox)) +
  geom_histogram(binwidth = 1 , 
                 fill = "white",
                 color = "black"
                 ) +
  labs(x = "ROX index at intubation",
       y = "Number of patients") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  coord_cartesian(xlim = c(0, 25))  # データを捨てずにズーム

p_rox

# 1. Time distribution
p1 <- invasive_start_data %>%
  ggplot(aes(x = time)) +
  geom_histogram(binwidth = 6, fill = "steelblue", color = "white") +
  labs(title = "Distribution of Time at Intubation Start",
       x = "Time (hours from ICU admission)",
       y = "Count") +
  theme_minimal()
p1


# Display excluding ROX outliers
# Calculate mean and median of ROX
rox_stats <- invasive_start_data %>%
  filter(!is.na(rox), rox < 30) %>%
  summarise(
    mean_rox = mean(rox),
    median_rox = median(rox)
  )
max_count <- invasive_start_data %>%
  filter(!is.na(rox), rox < 30) %>%
  pull(rox) %>%
  hist(breaks = 30, plot = FALSE) %>%
  .$counts %>%
  max()

p2 <- invasive_start_data %>%
  filter(!is.na(rox), rox < 30) %>%
  ggplot(aes(x = rox)) +
  geom_histogram(bins = 30, fill = "coral", color = "white") +
  # Existing threshold values
  geom_vline(xintercept = 3.85, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 4.88, linetype = "dashed", color = "red", linewidth = 1) +
  # Mean
  geom_vline(xintercept = rox_stats$mean_rox, linetype = "solid", color = "darkgreen", linewidth = 1) +
  # Median
  geom_vline(xintercept = rox_stats$median_rox, linetype = "dotted", color = "purple", linewidth = 1) +
  # Labels (y-coordinate placed as proportion of maximum)
  annotate("text", x = 3.85, y = max_count * 0.95, 
           label = "ROX = 3.85", 
           vjust = 0, hjust = 1.1, color = "blue", size = 3.5) +
  annotate("text", x = 4.88, y = max_count * 0.90, 
           label = "ROX = 4.88", 
           vjust = 0, hjust = 1.1, color = "red", size = 3.5) +
  annotate("text", x = rox_stats$mean_rox, y = max_count * 0.85, 
           label = sprintf("Mean = %.2f", rox_stats$mean_rox), 
           vjust = 0, hjust = -0.1, color = "darkgreen", size = 3.5) +
  annotate("text", x = rox_stats$median_rox, y = max_count * 0.80, 
           label = sprintf("Median = %.2f", rox_stats$median_rox), 
           vjust = 0, hjust = -0.1, color = "purple", size = 3.5) +
  labs(title = "Distribution of ROX Index at Intubation Start",
       x = "ROX Index",
       y = "Count") +
  theme_minimal()

p2


dt %>%
  filter(invasive_start == 1, !is.na(rox)) %>%
  mutate(
    rox_category = factor(
      case_when(
        rox < 3.85 ~ "ROX < 3.85",
        rox >= 3.85 & rox < 4.88 ~ "3.85 ≤ ROX < 4.88",
        rox >= 4.88 ~ "ROX ≥ 4.88"
      ),
      levels = c("ROX < 3.85", "3.85 ≤ ROX < 4.88", "ROX ≥ 4.88")
    )
  ) %>%
  dplyr::select(rox_category) %>%
  tbl_summary(
    label = list(rox_category ~ "ROX Category at Intubation")
  ) %>%
  bold_labels()


# ROX when not intubated
p3<-dt %>%
  filter(!is.na(rox), rox < 30) %>%
  filter(time>=0 & invasive_start_cum == 0) |> 
  ggplot(aes(x = rox)) +
  geom_histogram(bins = 30, fill = "coral", color = "white") +
  geom_vline(xintercept = 3.85, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 4.88, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 3.85, y = Inf, label = "3.85", 
           vjust = 1.5, hjust = -0.2, color = "blue", size = 4) +
  annotate("text", x = 4.88, y = Inf, label = "4.88", 
           vjust = 1.5, hjust = -0.2, color = "red", size = 4) +
  labs(title = "Distribution of ROX Index at non intubation",
       x = "ROX Index",
       y = "Count") +
  theme_minimal()

p3

dt %>%
  filter(time>=0 & invasive_start_cum == 0) |> 
  mutate(
    rox_category = factor(
      case_when(
        rox < 3.85 ~ "ROX < 3.85",
        rox >= 3.85 & rox < 4.88 ~ "3.85 ≤ ROX < 4.88",
        rox >= 4.88 ~ "ROX ≥ 4.88"
      ),
      levels = c("ROX < 3.85", "3.85 ≤ ROX < 4.88", "ROX ≥ 4.88")
    )
  ) %>%
  dplyr::select(rox_category) %>%
  tbl_summary(
    label = list(rox_category ~ "ROX Category at Intubation")
  ) %>%
  bold_labels()

# Back to Back plot ----

# Create histogram data
breaks_rox <- seq(0, 30, by = 1)

# Histogram at intubation (proportion)
intubated_hist <- dt %>%
  filter(!is.na(rox), rox < 30, invasive_start == 1) %>%
  mutate(rox_bin = cut(rox, breaks = breaks_rox, include.lowest = TRUE)) %>%
  count(rox_bin) %>%
  mutate(
    group = "Intubated",
    total = sum(n),
    proportion = n / total * 100,
    count = -proportion,  # Make negative (display on lower side)
    # Calculate midpoint from rox_bin
    rox_midpoint = as.numeric(gsub("\\(|\\]|\\[", "", sub(",.*", "", as.character(rox_bin)))) + 0.5
  )

# Histogram when not intubated (proportion)
non_intubated_hist <- dt %>%
  filter(!is.na(rox), rox < 30, time >= 0, invasive_start_cum == 0) %>%
  mutate(rox_bin = cut(rox, breaks = breaks_rox, include.lowest = TRUE)) %>%
  count(rox_bin) %>%
  mutate(
    group = "Non-intubated",
    total = sum(n),
    proportion = n / total * 100,
    count = proportion,  # Keep positive (display on upper side)
    # Calculate midpoint from rox_bin
    rox_midpoint = as.numeric(gsub("\\(|\\]|\\[", "", sub(",.*", "", as.character(rox_bin)))) + 0.5
  )

# Check sample sizes
cat("Intubated n =", unique(intubated_hist$total), "\n")
cat("Non-intubated n =", unique(non_intubated_hist$total), "\n")

# Combine
combined_hist <- bind_rows(intubated_hist, non_intubated_hist)

# Vertical plot (proportion)
p_pyramid_proportion <- ggplot(combined_hist, aes(x = rox_midpoint, y = count, fill = group)) +
  geom_col(width = 0.9, alpha = 0.8) +
  geom_vline(xintercept = 3.85, linetype = "dashed", color = "blue", linewidth = 1) +
  geom_vline(xintercept = 4.88, linetype = "dashed", color = "red", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  annotate("text", x = 3.85, y = max(combined_hist$count) * 0.9, 
           label = "ROX = 3.85", color = "blue", size = 3.5, hjust = -0.1) +
  annotate("text", x = 4.88, y = max(combined_hist$count) * 0.9, 
           label = "ROX = 4.88", color = "red", size = 3.5, hjust = -0.1) +
  scale_fill_manual(
    values = c("Intubated" = "coral", "Non-intubated" = "lightblue"),
    labels = c(
      "Intubated" = sprintf("Intubated (n=%d)", unique(intubated_hist$total)),
      "Non-intubated" = sprintf("Non-intubated (n=%d)", unique(non_intubated_hist$total))
    )
  ) +
  scale_y_continuous(
    labels = function(x) paste0(abs(x), "%"),
    name = "Proportion (%)"
  ) +
  labs(x = "ROX Index",
       fill = "",
       title = "Distribution of ROX Index: Intubated vs Non-intubated") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black"),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

p_pyramid_proportion



# Why they were intubated 10 > rox
invasive_start_data |> filter(rox>10) |> dplyr::select(
  time,rox,fio2,spo2,resp_rate,mbp,sbp,heart_rate,gcs
)


# HFNC initiated after ICU day 7
dt |> filter(hr>7*24 & time==0) |> dplyr::select(hr) |> summary()


# data sparse
skim(dt)

dt <- dt %>%
  group_by(subject_id) %>%
  arrange(subject_id, hr) %>%
  mutate(
    # Whether vasopressor == 1 at least once after time >= 0
    vasopressor_use = if_else(
      any(time >= 0 & vasopressor == 1, na.rm = TRUE), 1, 0
    )
  ) %>%
  ungroup()


# Censor
dt <- dt %>%
  group_by(subject_id) %>%
  arrange(subject_id, hr) %>%
  mutate(
    # Censor determination:
    # Currently on HFNC, HFNC ends at next time point and NIV starts
    censor = if_else(
      invasive_start_cum == 0 &           # Not yet intubated
        highflow == 1 &                     # Currently on HFNC
        lead(highflow, default = 0) == 0 &  # HFNC ends at next time point
        lead(noninvasive, default = 0) == 1,  # NIV starts at next time point
      1, 0
    )
  ) %>%
  ungroup()

# Check cases where censor==1
censor_cases <- dt %>%
  filter(censor == 1) %>%
  dplyr::select(subject_id, time, highflow, noninvasive, 
         invasive_start_cum, censor)

print(censor_cases)
cat("\nTotal censored observations:", nrow(censor_cases), "\n")
cat("Unique censored patients:", n_distinct(censor_cases$subject_id), "\n")

# Check time series for sample patient
if(nrow(censor_cases) > 0) {
  sample_id <- censor_cases$subject_id[1]
  dt %>%
    filter(subject_id == sample_id) %>%
    filter(time >= min(censor_cases$time[censor_cases$subject_id == sample_id]) - 5 &
             time <= min(censor_cases$time[censor_cases$subject_id == sample_id]) + 5) %>%
    dplyr::select(subject_id, time, highflow, noninvasive, invasive, invasive_start_cum, censor)
}

# Backup original data
dt_original <- dt
#dt<-dt_original

# Remove data after censor
dt <- dt %>%
  group_by(subject_id) %>%
  arrange(subject_id, hr) %>%
  mutate(
    censor_occurred = cumsum(censor) > 0,
    keep = !lag(censor_occurred, default = FALSE)  # Keep the censor row itself
  ) %>%
  filter(keep) %>%
  dplyr::select(-censor_occurred, -keep) %>%
  ungroup()

# Verification
# Check censored patients
censor_patients <- dt %>%
  filter(censor == 1) %>%
  distinct(subject_id) %>%
  pull(subject_id)

cat("\nNumber of censored patients:", length(censor_patients), "\n")


# Check if last observation has censor==1
dt %>%
  filter(subject_id %in% censor_patients) %>%
  group_by(subject_id) %>%
  slice_tail(n = 1) %>%
  count(censor)

cat("Original data rows:", nrow(dt_original), "\n")
cat("Filtered data rows:", nrow(dt), "\n")
cat("Removed rows:", nrow(dt_original) - nrow(dt), "\n")


# Check if data after censor is deleted for sample patient
if(length(censor_patients) > 0) {
  sample_id <- censor_patients[1]

  cat("\nSample patient", sample_id, "data:\n")
  dt %>%
    filter(subject_id == sample_id) %>%
    dplyr::select(subject_id, time, highflow, noninvasive, invasive, censor) %>%
    tail(10) %>%
    print()
}

# Check if last observation has censor==1
dt %>%
  filter(subject_id %in% censor_patients) %>%
  group_by(subject_id) %>%
  slice_tail(n = 1) %>%
  count(censor)

# Data preparation
df<-dt
df$id<-as.numeric(as.factor(df$subject_id))

# Change Variable Name -----

df$cal_time<-df$time
df$cal_timesqr<-df$time^2
df$death<-df$death_event
df$treat<-df$invasive_start
df$treat_cum<-df$invasive_start_cum

df <- df %>%
  group_by(id) %>%
  arrange(id, hr) %>%
  mutate(
    treat_lag1 = lag(invasive_start, n = 1, default = 0),
    treat_cum_lag1 = lag(invasive_start_cum, n = 1, default = 0)
  ) %>%
  ungroup()

glimpse(df)


df$time<-df$cal_time
df$timesqr<-df$cal_timesqr

# Using dplyr
df <- df %>%
  group_by(id) %>%
  mutate(highflow_lag = lag(highflow, n = 1, 
                            default = first(highflow))) %>%
  ungroup()

df<-df |> dplyr::select(id,elig_1,elig_2,cal_time,cal_timesqr,
             hr_at_hfnc_start,
             age,anchor_year_group,
             female,medicare_medicaid,language_grouped,Married,race_grouped,
             highflow,highflow_lag,invasive,noninvasive,
             heart_rate, mbp, temperature, 
             resp_rate, spo2,fio2,rox,
             gcs_category,ph,po2,pco2, sofa_24hours,elixhauser_vanwalraven,
             sepsis3,vasopressor,crrt,
             treatment_limitation,
             icu_discharge_event,
             hosp_discharge_event,
             death_30day,
             death,
             censor,
             treat, # invasive_start 
             treat_lag1, #invasive_start,
             treat_cum, # invasive_start_cum,
             treat_cum_lag1,
             time,
             timesqr)

# Save data ----
saveRDS(df, "./data/ROX/imputed_data_for_analysis.rds")

