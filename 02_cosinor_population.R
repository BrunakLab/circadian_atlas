# 02_cosinor_population.R
# This script applies the robustness filters on the 24-hour population data
# and fits the cosinor model. The output is a table with the circadian parameters
# of each biomarker.
# Author: Álvaro G. León

# Note: some filters have been changed in the demo compared to the original code
# due to the lack of coverage of the simulated data (indicated by "NOTE")

# Input: 1.normalized_hourly_data_population.tsv
# Output: 2.cosinor_results_population.tsv

library(tidyverse)
source("functions.R")


# Data loading ----
df <- read_tsv("simulated_output_data/1.normalized_hourly_data_population.tsv")


# Step 1: Robustness filters  ----
####### Removing hour data produced by fewer than 100 unique patients ######
# NOTE: in the demo, number of unique patients per hour has been reduced 
# from 100 to 15
# patient_variable = 100
patient_variable = 15

df <- df %>%
  filter(number_unique_patients >= patient_variable) 


# Info: Number of unique tests remaining: 504 



####### Filtering by max hour gap of 5 hours  ######

# It keeps only components whose largest gap between observed
# hours (0–23) does not exceed 5.

max_gap_threshold <-  5

valid_coverage <- df %>%
  group_by(component_simple_lookup) %>%
  summarize(
    hours_with_data = n_distinct(hour_rounded),
    max_gap = max_hour_gap(hour_rounded),
    .groups = "drop"
  ) %>%
  filter(max_gap <= max_gap_threshold)


df <- df %>%
  semi_join(valid_coverage, by = "component_simple_lookup")

# Info: Number of unique tests remaining: 186 



# Step 2: Apply weights based on number of number of observations (entries) ----

hour_weights <- df %>%
  group_by(component_simple_lookup) %>%
  mutate(weight = number_entries / sum(number_entries)) %>%
  ungroup()  %>%
  select(-number_entries, -number_unique_patients)


df_weighted <- df %>%
  inner_join(hour_weights)


# Step 3: Fitting the cosinor model -------

# The goal here is to:
# Produce the predicted 24-hour rhythms for each biomarker (curves)
# Obtain the circadian parameters of each biomarker (params)
# Merge them into a combined dataset (curves)

fits <- df_weighted %>%
  group_by(component_simple_lookup) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$component_simple_lookup))) %>%
  map(purrr::safely(fit_linear_cosinor)) 


curves <- fits %>%
  keep(~ is.null(.x$error)) %>%
  map("result") %>%
  map_dfr("fitted_full", .id = "component_simple_lookup") %>% 
  rename(hour_rounded = time)


params <- fits %>%
  keep(~ is.null(.x$error)) %>%
  map("result") %>%
  map_dfr("params", .id = "component_simple_lookup")


# FDR correction (BH)
params <- params %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    significant = ifelse(p_adj < 0.05, TRUE, FALSE)
  )


curves <- full_join(curves, params, by = "component_simple_lookup")

### Adding number of patient information back into each component-hour #####
curves <- curves %>%
  left_join(df, by = c("component_simple_lookup", "hour_rounded")) %>% 
  mutate(
    rel_amp = amplitude / (1 + M),
    peak = M+amplitude,
    trough = M-amplitude,
  )


write_tsv(curves, "simulated_output_data/2.cosinor_results_population.tsv")



