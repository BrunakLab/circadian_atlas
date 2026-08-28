# 00_generate_demo_data_lab_clean.R
# Generates a simulated dataset with the same columns as the cleaned LAB table. 
# Author: Álvaro G. León
# Input: none
# Output: demo_lab_data_clean.tsv

set.seed(1)

N_ROWS     <- 12000
N_PATIENTS <- 1500

components <- c("Cortisol - Plasma", "Glucose - Plasma", "Leukocytes - Blood", "Creatinine - Plasma")
units      <- c("nmol/L", "mmol/L", "10^9/L", "umol/L")
labs       <- c("LAB01", "LAB02")

patients <- data.frame(
  pid = sprintf("PT%05d", seq_len(N_PATIENTS)),
  DOB = as.Date("1970-01-01") - sample(0:14600, N_PATIENTS, replace = TRUE),
  sex = sample(c("F", "M"), N_PATIENTS, replace = TRUE),
  stringsAsFactors = FALSE
)

i <- sample(N_PATIENTS, N_ROWS, replace = TRUE)
k <- sample(seq_along(components), N_ROWS, replace = TRUE)

# Non-uniform sampling
hour_weights <- c(rep(1, 7),     # 00:00-06:00  night
                  rep(12, 9),    # 07:00-15:00  daytime peak
                  rep(4, 4),     # 16:00-19:00  evening
                  rep(1, 4))     # 20:00-23:00  night

demo <- data.frame(
  pid                     = patients$pid[i],
  DOB                     = patients$DOB[i],
  sex                     = patients$sex[i],
  date                    = as.Date("2020-01-01") + sample(0:365, N_ROWS, replace = TRUE),
  hour_rounded            = as.numeric(sample(0:23, N_ROWS, replace = TRUE,
                                              prob = hour_weights)),
  component_simple_lookup = components[k],
  shown_clean             = round(runif(N_ROWS, 1, 100), 1),
  unit_clean              = units[k],
  lab_id                  = sample(labs, N_ROWS, replace = TRUE),
  stringsAsFactors = FALSE
)


write.table(demo, "simulated_input_data/demo_lab_data_clean.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)

