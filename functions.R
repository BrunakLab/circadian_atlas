# functions.R
# Shared functions for the processing and circadian analysis of laboratory data 
# Author: Álvaro G. León

# ------------------------------------------------------------------------------
# 1. Data coverage
# ------------------------------------------------------------------------------
# Largest gap in hourly coverage across the day
#
# Returns the size of the longest run of consecutive missing hours, treating
# the day as circular so that a gap spanning midnight is counted. Used to
# exclude laboratory tests with insufficient coverage before cosinor fitting.
#
# @param hours  numeric vector of hours of day (0-23) with observations
# @return       number of missing hours in the largest gap; 24 if fewer than
#               two distinct hours are present
max_hour_gap <- function(hours) {
  
  hours <- sort(unique(hours))
  
  if (length(hours) <= 1) return(24)   # if only one hour, gap = full day
  
  diffs <- diff(hours)
  
  # include circular gap from last hour to first (wrap around 24)
  diffs <- c(diffs, 24 - (max(hours) - min(hours)))
  
  gap <- max(diffs)
  
  return(gap - 1)
}


# ------------------------------------------------------------------------------
# 2. Cosinor model
# ------------------------------------------------------------------------------

# Fit a weighted linear cosinor model to hourly means
#
# Fits value ~ cos(2*pi*t/period) + sin(2*pi*t/period) by weighted least
# squares and returns the cosinor parameters with 95% confidence intervals
# and the Zero-Amplitude Test statistics.
#
# Confidence intervals use t for the mesor and amplitude, and F for the acrophase (Bingham's method). 
# The amplitude standard error (SE) is obtained by the delta method.
#
# @param df         data frame for a single laboratory test, containing the
#                   time and value columns plus a `weight` column
# @param time_col   name of the hour-of-day column (here, hour_rounded)
# @param value_col  name of the value column (here, mean_value)
# @param period     rhythm period in hours (here, assumed 24)
# @return list with elements:
#   \describe{
#     \item{model}{the fitted lm object}
#     \item{params}{one-row tibble of parameters, CIs and test statistics}
#     \item{fitted}{observed-hour fitted values and residuals}
#     \item{fitted_full}{fitted curve over all 24 hours}
#   }

fit_linear_cosinor <- function(df, time_col = "hour_rounded", value_col = "mean_value", period = 24) {
  time  <- df[[time_col]]
  value <- df[[value_col]]
  w     <- df$weight
  
  # Linear terms
  x_coord <- cos(2 * pi * time / period)
  z_coord <- sin(2 * pi * time / period)
  
  # Weighted Linear Model
  model <- lm(value ~ x_coord + z_coord, weights = w)
  
  # Extract model info
  s      <- summary(model)
  coeffs <- s$coefficients
  M_hat  <- coeffs[1, 1]
  beta   <- coeffs[2, 1]
  gamma  <- coeffs[3, 1]
  
  # Derived Params
  A_hat   <- sqrt(beta^2 + gamma^2)
  phi_hat <- (atan2(gamma, beta) * period / (2 * pi)) %% period
  
  # Significance & R2
  f_stats <- s$fstatistic
  p_val   <- pf(f_stats[1], f_stats[2], f_stats[3], lower.tail = FALSE)
  R2      <- s$r.squared
  
  # Variances for CI 
  V      <- vcov(model)
  v_b    <- V[2, 2]
  v_g    <- V[3, 3]
  cov_bg <- V[2, 3]
  
  t_crit <- qt(0.975, model$df.residual)
  
  # Mesor SE 
  se_mesor <- sqrt(V[1, 1])
  
  # Amplitude SE (delta method) 
  se_amp  <- sqrt((beta^2 * v_b + gamma^2 * v_g + 2 * beta * gamma * cov_bg) / A_hat^2)
  
  # Bingham Acrophase CI 
  f_crit   <- qf(0.95, 2, model$df.residual)
  D <- (beta^2 * v_g + gamma^2 * v_b - 2 * beta * gamma * cov_bg) - 
    (f_crit * 2 * (v_b * v_g - cov_bg^2))
  
  acro_low <- NA; acro_high <- NA
  if (D >= 0) {
    delta_phi <- (period / (2 * pi)) * asin(sqrt(2 * f_crit * (v_b * v_g - cov_bg^2) / (A_hat^2 * (v_b + v_g))))
    acro_low  <- (phi_hat - delta_phi) %% period
    acro_high <- (phi_hat + delta_phi) %% period
  }
  
  # Generate fitted hourly values
  hours_full <- tibble(time = 0:23) %>%
    mutate(fitted = M_hat + A_hat * cos(2 * pi * (time - phi_hat) / period))
  
  list(
    model = model,
    params = tibble(
      M            = M_hat,
      se_mesor     = se_mesor,
      mesor_low    = M_hat - t_crit * se_mesor,
      mesor_high   = M_hat + t_crit * se_mesor,
      amplitude    = A_hat,
      se_amplitude = se_amp,
      amp_low      = A_hat - t_crit * se_amp,
      amp_high     = A_hat + t_crit * se_amp,
      acrophase    = phi_hat,
      acro_low     = acro_low,
      acro_high    = acro_high,
      f_statistic  = unname(f_stats[1]),
      df_num       = unname(f_stats[2]),
      df_den       = unname(f_stats[3]),
      p_value      = unname(p_val),
      r_squared    = R2
      ),
    fitted = tibble(
      !!time_col := time,
      fitted     = fitted(model),
      residual   = resid(model)
    ),
    fitted_full = hours_full
  )
}



# ------------------------------------------------------------------------------
# 3. Acrophase and age
# ------------------------------------------------------------------------------

# Mardia's circular-linear correlation coefficient
#
# Computes the squared circular-linear correlation between a linear variable
# (age group) and a circular one (acrophase), from the Pearson correlations
# between the linear variable and the sine and cosine components of the angle.
# We assume a 24-hour period.
#
# @param lin_vec   linear variable
# @param circ_vec  circular variable, in hours (0-24), not radians
# @return          R-squared, constrained to [0, 1]
# 

calculate_mardia_R2 <- function(lin_vec, circ_vec) {
  
  # Convert 24h hours to radians
  theta <- (circ_vec / 24) * 2 * pi
  y_cos <- cos(theta)
  y_sin <- sin(theta)
  
  # Pearson correlations
  r_xc <- cor(lin_vec, y_cos)
  r_xs <- cor(lin_vec, y_sin)
  r_cs <- cor(y_cos, y_sin)
  
  # Check: (prevents Inf/NaN for small N)
  denom <- 1 - r_cs^2
  if (is.na(denom) || denom < 1e-10) { denom <- 1e-10 }
  
  R2 <- (r_xc^2 + r_xs^2 - 2 * r_xc * r_xs * r_cs) / denom
  
  # Check: R2 must stay between 0 and 1
  return(min(max(R2, 0), 1))
}


# ------------------------------------------------------------------------------
# 4. Mortality models
# ------------------------------------------------------------------------------

# Fit the 28-day mortality logistic model for a single biomarker
#
# Opens a read-only DuckDB connection to the database containing the dynamic ranges
# pulls the rows for one biomarker, fits a logistic model of 28-day mortality 
# on the four-way classification with age and sex as adjustments and standard 
# errors clustered by patient, and returns tidied coefficients exponentiated to odds ratios.
#
#
# NOTE: in the original code, the function reads a table for a given biomarker
# and collects it. Therefore, the code related to DuckDB functionality
# is commented out in its respectives lines.
#
#
# @param x     biomarker identifier (`component_simple_lookup`)
# @param data  data frame of mortality data to draw the biomarker rows from
#              (in the demo, the contents of demo_mortality_data.tsv)
# @return      tibble of tidied coefficients, or NULL if the fit failed

run_parallel_db_model <- function(x, data) {
  message(paste0("[", Sys.time(), "] Starting: ", x))

  # con <- dbConnect(duckdb(), dbdir = db_string, read_only = TRUE)

  result <- try({

    # df <- tbl(con, "lab_dynamic_ranges_dod") %>%
    #  filter(component_simple_lookup == x) %>%
    #  collect()

    # In the demo version, it draws from the demo_mortality_data tsv file
    df <- data %>%
      filter(component_simple_lookup == x)
    
    message(paste0("[", Sys.time(), "] ", x, ": pulled ", nrow(df), " rows"))
    
    # Data prep
    df$classification <- fct_relevel(factor(df$classification), "TN")
    
    # covariates
    covs <- c("classification", "age", "sex")
    form <- as.formula(paste("mortality_28 ~", paste(covs, collapse = " + ")))
    
    # E. Modeling
    model <- feglm(
      form,
      data    = df,
      family  = "logit",
      cluster = ~pid
    )
    
    # convergence check
    if (isFALSE(model$convStatus)) {
      warning(paste0(x, ": model ran but did not converge (convStatus = FALSE)"))
    }
    
    tidy_result <- broom::tidy(model, conf.int = TRUE) %>%
      mutate(
        estimate  = exp(estimate),
        conf.low  = exp(conf.low),
        conf.high = exp(conf.high),
        component_simple_lookup = x
      )
    
    tidy_result
  }, silent = TRUE)   
    
  # dbDisconnect(con, shutdown = TRUE)
  
  # success messages
  if (inherits(result, "try-error")) {
    message(paste0("[", Sys.time(), "] FAILED: ", x, " — ", conditionMessage(attr(result, "condition"))))
    return(NULL)
  }
  
  message(paste0("[", Sys.time(), "] Finished: ", x))
  return(result)
}














