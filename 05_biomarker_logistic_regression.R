# 05_biomarker_logistic_regression.R
# This script fits a logistic model of 28-day mortality for each biomarker
# based on the classification categories (TN, TP, FN, FP) derived from the 
# concordance of static vs dynamic ranges, with age and sex as confounders.
# Author: Álvaro G. León

# Input: demo_mortality_data.tsv
# Output: 5.logistic_regression_results.tsv


library(tidyverse)
library(fixest)
library(parallel)
source("functions.R")


## Data loading and processing -------------

### Loading #####

df_mortality <- read_tsv("simulated_input_data/demo_mortality_data.tsv")



### Obtaining relative standard error for each classification #####

quality_map <- df_mortality %>%
  group_by(component_simple_lookup, classification, mortality_28) %>%
  summarize(n = n(), .groups = "drop") %>%
  collect() %>%
  # Creates columns like FN_0, FN_1, TN_0, TN_1, etc.
  pivot_wider(
    names_from = c(classification, mortality_28),
    values_from = n,
    values_fill = 0
  ) %>%
  # 3. Calculate RSEs for every category 
  mutate(
    rse_FN = ifelse(FN_1 > 0, 1 / sqrt(FN_1), Inf),
    rse_FP = ifelse(FP_1 > 0, 1 / sqrt(FP_1), Inf),
    rse_TP = ifelse(TP_1 > 0, 1 / sqrt(TP_1), Inf),
    rse_TN = ifelse(TN_1 > 0, 1 / sqrt(TN_1), Inf)
  )

# Filtering to Relative standard error below 0.5.
components_to_run <- quality_map %>%
  filter(rse_FN < 0.5) %>%
  pull(component_simple_lookup)



# Running the logistic regression in parallel ----
no_cores <- 6
print(paste("Launching parallel DB modeling on", no_cores, "cores..."))

final_list <- mclapply(components_to_run, run_parallel_db_model,
                       data = df_mortality, mc.cores = no_cores)

print("Combining results...")
combined_results <- bind_rows(final_list)

write_tsv(combined_results, "simulated_output_data/5.logistic_regression_results.tsv")


