# 06_mortality_metaregression.R
# This script performs a weighted linear regression across biomarkers to assess
# the influence of circadian amplitude in 28-day mortality
# Author: Álvaro G. León

# Input: 2.cosinor_results_population.tsv
# Input: 5.logistic_regression_results.tsv
# Output: 6.summary_data_statistics_meta_regression.tsv

library(tidyverse)
library(broom)

## Data loading and processing -------------
### Loading population circadian data #####
curves_population <- read_tsv( "simulated_output_data/2.cosinor_results_population.tsv")%>% 
  select(component_simple_lookup,amplitude, acrophase) %>% 
  distinct()

### Loading logistic regression results #####
df <- read_tsv("simulated_output_data/5.logistic_regression_results.tsv") %>%
  distinct() %>% 
  filter(term %in% c("classificationFN")) 


### Processing and merging #####
df <- df %>%
  select(component_simple_lookup, term, estimate, std.error, conf.low, conf.high,p.value) %>%
  pivot_wider(
    names_from = term, 
    values_from = c(p.value,estimate, std.error, conf.low, conf.high),
    names_glue = "{term}_{.value}"
  )

merged <- df %>%
  left_join(curves_population, by = "component_simple_lookup") %>% 
  mutate(amplitude = scale(amplitude)) %>% distinct()


# Meta regression --------
model_fn <- lm(log(classificationFN_estimate) ~ amplitude, 
               data = merged, 
               weights = 1 / (classificationFN_std.error^2))

summary <- broom::tidy(model_fn, conf.int = TRUE) %>%
  mutate(df_residual = model_fn$df.residual)

write_tsv(summary, paste0("simulated_output_data/6.summary_data_statistics_meta_regression.tsv"))

