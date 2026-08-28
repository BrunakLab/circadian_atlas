# 03_cosinor_age.R
# This script filters the 24-hour age-stratified data to retain the biomarkers that
# were part of the population cosinor model, and fits an age-stratified model to them.
# The output is a table with the circadian parameters of each biomarker per age group.
# Author: Álvaro G. León

# Input: 1.normalized_hourly_data_age.tsv
# Input: 2.cosinor_results_population.tsv
# Output: 3.cosinor_results_age.tsv

library(tidyverse)
source("functions.R")


## Data loading -------------

### Loading cosinor results from population cohort #####
# This is used to filter the tests to be analysed in the age cohorts
curves_population <- read_tsv( "simulated_output_data/2.cosinor_results_population.tsv")


### Loading normalized data from age cohorts #####
df <- read_tsv("simulated_output_data/1.normalized_hourly_data_age.tsv")


# Step 1: Restrict age lab tests to those of population cosinor #####

# In the age sub-analysis, we do not apply the patient or hour gap filters
# We restrict the tests to those that passed the filters in the population
# analysis

df <- df %>%
  filter(component_simple_lookup %in% curves_population$component_simple_lookup)


# Step 2: Apply weights based on number of number of observations (entries) ----

hour_weights <- df %>%
  group_by(component_simple_lookup,age_group) %>%
  mutate(weight = number_entries / sum(number_entries)) %>%
  ungroup()  %>%
  select(-number_entries, -number_unique_patients)

df_weighted <- df %>%
  inner_join(hour_weights, by = c("hour_rounded", "component_simple_lookup", "mean_value", "age_group"))




# Step 3: Fitting the cosinor model -------

# The goal here is to:
# Produce the predicted 24-hour rhythms for each biomarker (curves)
# Obtain the circadian parameters of each biomarker (params)
# Merge them into a combined dataset (curves)

# Fit to the cosinor
fits <- df_weighted %>%
  group_by(component_simple_lookup, age_group) %>%
  group_split() %>%
  set_names(map_chr(., ~ paste(unique(.x$component_simple_lookup),
                               unique(.x$age_group), sep = "_"))) %>%
  map(purrr::safely(fit_linear_cosinor))

params <- fits %>%
  keep(~ is.null(.x$error)) %>%
  map("result") %>%
  map_dfr("params", .id = "id")


# FDR correction (BH)
params <- params %>%
  separate(id, into = c("component_simple_lookup", "age_group"), sep = "_", remove = FALSE) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    significant = ifelse(p_adj < 0.05, TRUE, FALSE)
  )

curves <- fits %>%
  keep(~ is.null(.x$error)) %>%
  map("result") %>%
  map_dfr("fitted_full", .id = "id") %>%
  rename(hour_rounded = time) %>% 
  separate(id, into = c("component_simple_lookup", "age_group"), sep = "_", remove = FALSE)


curves <- full_join(curves, params, by = c("component_simple_lookup", "age_group", "id"))

df <- df %>%
  mutate(age_group = as.character(age_group))

curves <- curves %>%
  left_join(df, by = c("component_simple_lookup","age_group", "hour_rounded"))


write_tsv(curves, "simulated_output_data/3.cosinor_results_age.tsv")

