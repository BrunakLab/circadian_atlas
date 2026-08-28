# 04_acrophase_age_correlations.R
# This script computes the Mardia linear-circular correlation between age group
# and acrophase for every biomarker, as well as the median shift in acrophase across
# age groups. 
# Author: Álvaro G. León

# Note: some filters have been changed in the demo compared to the original code
# due to the lack of coverage of the simulated data (indicated by "NOTE")

# Input: 3.cosinor_results_age.tsv
# Output: 4.acrophase_age_correlation.tsv


library(tidyverse)
library(broom)
source("functions.R")


# Data loading
curves <- read_tsv("simulated_output_data/3.cosinor_results_age.tsv")

# Choosing only significant biomarker-age circadian parameters ----
# NOTE: this step is skipped in the demo data due to its limitations
curves <- curves %>% filter(!is.na(significant))
# curves <- curves %>% filter(significant == TRUE)

curves <- curves %>%
  mutate(age_group = as.numeric(age_group))

# Getting only unique biomarker-age circadian parameters -----
curves_unique <- curves %>%
  select(-id,-hour_rounded, -fitted, -mean_value, -number_entries, -number_unique_patients) %>%
  distinct() 



# Step 1: Producing the correlation and acrophase velocity per biomarker - age group ------
component_names <- unique(curves_unique$component_simple_lookup)
n_comp <- length(component_names)

analysis_results_speed <- data.frame(
  component_simple_lookup = character(n_comp),
  Acro_R2 = numeric(n_comp),
  Acro_P = numeric(n_comp),
  Acro_Is_Sig = logical(n_comp),
  Acro_Velocity = numeric(n_comp),
  Acro_Flow = character(n_comp),      
  N_Unique_Age_Groups = integer(n_comp),
  stringsAsFactors = FALSE
)


set.seed(123)

for (i in 1:n_comp) {
  comp <- component_names[i]
  sub_data <- curves_unique[curves_unique$component_simple_lookup == comp, ]
  sub_data <- sub_data[order(as.numeric(as.character(sub_data$age_group))), ]
  
  current_ages <- as.numeric(as.character(sub_data$age_group))
  current_hrs  <- as.numeric(sub_data$acrophase)
  valid_idx <- complete.cases(current_ages, current_hrs)
  
  current_ages <- current_ages[valid_idx]; current_hrs <- current_hrs[valid_idx]
  n_obs <- length(current_ages)
  
  if (n_obs < 3 || sd(current_hrs) == 0) {
    analysis_results_speed[i, ] <- list(comp, NA, NA, FALSE, NA, "Stable", n_obs)
    
    next
  }
  
  # 1. ACROPHASE VELOCITY & FLOW 
  diff_hrs  <- diff(current_hrs)
  diff_ages <- diff(current_ages)
  
  adj_diffs_acro <- ifelse(diff_hrs > 12, diff_hrs - 24, 
                           ifelse(diff_hrs < -12, diff_hrs + 24, diff_hrs))
  
  rates_per_decade_acro <- (adj_diffs_acro / diff_ages) #normalization to manage non-contiguous gaps
  avg_shift_acro <- median(rates_per_decade_acro)
  acro_flow <- if(avg_shift_acro > 0) "Later" else if(avg_shift_acro < 0) "Earlier" else "Stable"
  
  # 2. OBSERVED STATS
  obs_acro_R2 <- calculate_mardia_R2(current_ages, current_hrs)

  # 3. PERMUTATIONS
  n_perm <- 5000
  acro_exceed <- 0 
  for(j in 1:n_perm) {
    sh_hrs <- sample(current_hrs) 
    if(calculate_mardia_R2(current_ages, sh_hrs) >= obs_acro_R2) acro_exceed <- acro_exceed + 1
  }
  
  # 4. SIGNIFICANCE
  ap <- (acro_exceed + 1) / (n_perm + 1)

  # 6. SAVE
  analysis_results_speed[i, ] <- list(
    comp, obs_acro_R2, ap, ap < 0.05, avg_shift_acro, acro_flow, n_obs
  )
  
  # Progress tracker
  if(i %% 20 == 0) cat("Processed", i, "of", n_comp, "components...\n")
}



# FDR for the p values
analysis_results_speed$Acro_P_adj <- p.adjust(analysis_results_speed$Acro_P, method = "fdr")
analysis_results_speed$Acro_Is_Sig <- analysis_results_speed$Acro_P_adj < 0.05
 
analysis_results_speed <- analysis_results_speed %>% 
  relocate(Acro_P_adj, .after="Acro_P")


write_tsv(analysis_results_speed, "simulated_output_data/4.acrophase_age_correlation.tsv")


