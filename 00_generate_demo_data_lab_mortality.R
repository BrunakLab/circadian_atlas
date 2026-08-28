# 00_generate_demo_data_lab_mortality.R
# Generates a simulated dataset of lab results with a classification for 28-day mortality.
# Author: Álvaro G. León
# Input: none
# Output: demo_mortality_data.tsv

set.seed(1)

N_ROWS     <- 12000
N_PATIENTS <- 1500
N_VISITS   <- 3000    

components <- c("Cortisol - Plasma", "Glucose - Plasma", "Leukocytes - Blood", "Creatinine - Plasma")

# Relative frequencies of the four classification outcomes.
W_TN <- 70; W_TP <- 25; W_FN <- 10; W_FP <- 5


# --- Patients ------------------------------------------------------------
patients <- data.frame(
  pid = sprintf("PT%05d", seq_len(N_PATIENTS)),
  age = round(runif(N_PATIENTS, 18, 95), 1),
  sex = sample(c("F", "M"), N_PATIENTS, replace = TRUE),
  stringsAsFactors = FALSE
)

# --- Visits (one row per patient-day) ------------------------------------
visits <- data.frame(
  pid  = sprintf("PT%05d", sample(N_PATIENTS, N_VISITS, replace = TRUE)),
  date = as.Date("2020-01-01") + sample(0:365, N_VISITS, replace = TRUE),
  stringsAsFactors = FALSE
)
visits <- visits[!duplicated(paste(visits$pid, visits$date)), ]

# Mortality is a property of the patient-day, so every test taken that day
# carries the same value.
visits$mortality_28 <- rbinom(nrow(visits), 1, (W_TP + W_FN) / (W_TN + W_TP + W_FN + W_FP))

# --- Rows ----------------------------------------------------------------
i <- sample(nrow(visits), N_ROWS, replace = TRUE)

demo <- data.frame(
  pid                     = visits$pid[i],
  date                    = visits$date[i],
  component_simple_lookup = sample(components, N_ROWS, replace = TRUE),
  shown_clean             = round(runif(N_ROWS, 1, 100), 1),
  classification          = NA_character_,
  mortality_28            = as.numeric(visits$mortality_28[i]),
  age                     = patients$age[match(visits$pid[i], patients$pid)],
  sex                     = patients$sex[match(visits$pid[i], patients$pid)],
  stringsAsFactors = FALSE
)


demo$classification <- sample(c("TN", "TP", "FN", "FP"), nrow(demo),
                              replace = TRUE,
                              prob = c(W_TN, W_TP, W_FN, W_FP))

demo <- demo[sample(nrow(demo)), ]
rownames(demo) <- NULL


write.table(demo, "simulated_input_data/demo_mortality_data.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)



