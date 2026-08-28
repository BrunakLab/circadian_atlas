# 01_data_normalization.R
# Performs the steps indicated in the Methods section "Statistical and population-level
# normalization". The processing of the raw LAB table described in the
# "Data Preprocessing" section are not included. The output files are the
# 24-hour means for every biomarker.
# Author: Álvaro G. León

# Input: demo_lab_data_clean.tsv 
# Output: 1.normalized_hourly_data_population.tsv
# Output: 1.normalized_hourly_data_age.tsv


library(tidyverse)

# 1. Loading the cleaned LAB data  -------
# The original data at this point has been cleaned to only retain 
# numeric values, and NPU codes have been harmonized into standard names
# indicated in the column "component_simple_lookup"

lab <- read_tsv("simulated_input_data/demo_lab_data_clean.tsv")


# 2. Tukey's fences ------
# Removing outliers
quantiles_tbl <- lab %>%
  group_by(component_simple_lookup, unit_clean, lab_id) %>%
  summarise(
    q1 = quantile(shown_clean, 0.25),
    q3 = quantile(shown_clean, 0.75),
    .groups = "drop"
  ) %>%
  mutate(
    iqr = q3 - q1,
    lower_fence = q1 - 3 * iqr,
    upper_fence = q3 + 3 * iqr
  )

lab <- lab %>%
  left_join(quantiles_tbl,
            by = c("component_simple_lookup", "unit_clean", "lab_id")) %>%
  filter(
    shown_clean >= lower_fence,
    shown_clean <= upper_fence
  ) %>%
  ungroup() %>%
  select(-q1, -q3, -iqr, -lower_fence, -upper_fence) 


# Creating age groups 
lab <- lab %>%
  mutate(
    age = as.numeric((as.Date(date) - as.Date(DOB)) / 365),
    age_group = case_when(
      age < 10 ~ 1,
      age >= 10 & age < 20 ~ 2,
      age >= 20 & age < 30 ~ 3,
      age >= 30 & age < 40 ~ 4,
      age >= 40 & age < 50 ~ 5,
      age >= 50 & age < 60 ~ 6,
      age >= 60 & age < 70 ~ 7,
      age >= 70 & age < 80 ~ 8,
      age >= 80 & age < 90 ~ 9,
      age >= 90 & age < 100 ~ 10,
      age >= 100 ~ 11
    )
  )


# Changing the hours "24" to convert them to "00"
lab <- lab %>%
  mutate(
    hour_rounded = case_when(
      hour_rounded == 24 ~ 0,       
      TRUE             ~ hour_rounded  
    )
  )


# 3. Intra-individual averaging -----
# if pid has more than one test per HOUR, compute mean value

lab <- lab %>% 
  mutate(group = paste(pid, component_simple_lookup, 
                       hour_rounded, unit_clean, lab_id,
                       age_group))

# the slightly cryptic code here is due to technical limitations when
# working with database tables, which was the original format of the data
lab_ref <- lab

lab <- lab_ref %>%
  group_by(group) %>%
  summarise(shown_clean = mean(shown_clean, na.rm = TRUE), .groups = "drop") %>%
  inner_join(
    lab_ref %>% select(pid, hour_rounded, component_simple_lookup,  group, lab_id, 
                       unit_clean, sex, age_group),
    by = "group"
  ) %>%
  distinct(pid, hour_rounded, component_simple_lookup, shown_clean, unit_clean, 
           lab_id, sex, age_group)



# Info: number of unique biomarkers: 1880



# 4. Centering of values -----------------
lab <- lab %>% 
  mutate(group = paste(component_simple_lookup, 
                       unit_clean, lab_id, sex,
                       age_group))



lab_ref <- lab

lab <- lab_ref %>%
  group_by(group) %>%
  summarise(
    group_n = n(),
    group_mean = mean(shown_clean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(lab_ref, by = "group") %>%
  
  # filter valid groups with at least 500 observations
  # filter(group_n >= 500) %>%
  
  # NOTE: In the demo, the group n has been reduced to 50 
  filter(group_n >= 50) %>%
  
  # compute normalized values
  mutate(normal = (shown_clean / group_mean) - 1)






# 5. Computing the 24-hour means for every test  ------

df_testing <- lab %>%
  mutate(
    normal = case_when(
      abs(shown_clean - group_mean) < 1e-8 ~ 0,
      is.infinite((normal)) & (normal) > 0 ~ 1,
      is.infinite((normal)) & (normal) < 0 ~ -1,
      TRUE ~ (normal)
    )
  ) 

# Population 24-hour means (all strata included)
df_mean_population <- df_testing %>%
  group_by(component_simple_lookup, hour_rounded) %>%
  summarise(
    mean_value = mean(normal, na.rm = TRUE),
    number_entries = n(),
    number_unique_patients = n_distinct(pid),
    .groups = "drop"
  ) %>% collect()



write_tsv(df_mean_population, "simulated_output_data/1.normalized_hourly_data_population.tsv")


# Age-stratified 24-hour means 
df_mean_age <- df_testing %>%
  group_by(component_simple_lookup, hour_rounded,age_group) %>%
  summarise(
    mean_value = mean(normal, na.rm = TRUE),
    number_entries = n(),
    number_unique_patients = n_distinct(pid),
    .groups = "drop"
  )


write_tsv(df_mean_age, "simulated_output_data/1.normalized_hourly_data_age.tsv")




