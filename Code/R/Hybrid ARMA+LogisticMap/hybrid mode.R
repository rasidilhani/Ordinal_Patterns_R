library(StatOrdPattHxC)
library(dplyr)
library(writexl)

# ==============================================================================
# PARAMETERS
# ==============================================================================
set.seed(1234567890, kind = "Mersenne-Twister")
BETA <- 1.5               # Parameter for Rényi / Tsallis
N    <- c(1000, 5000)       # Sample sizes
D    <- 4                 # Embedding dimension
R    <- 50                # Number of repetitions
r_logistic <- 3.8         # Logistic map parameter (chaotic regime)
r_mix <- 0.5              # Mixing parameter for hybrid models

# ==============================================================================
# DEFINE ALL ARMA MODEL CONFIGURATIONS
# ==============================================================================
arma_models <- list(
  ARMA22 = list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), name = "ARMA(2,2)"),
  ARMA11 = list(ar = c(0.8), ma = c(0.8), name = "ARMA(1,1)"),
  AR2    = list(ar = c(0.1, 0.8), ma = NULL, name = "AR(2)"),
  AR1    = list(ar = c(-0.1), ma = NULL, name = "AR(1)"),
  MA2    = list(ar = NULL, ma = c(0.1, -0.8), name = "MA(2)"),
  MA1    = list(ar = NULL, ma = c(0.8), name = "MA(1)")
)

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Jensen–Shannon divergence between two discrete pmfs p, q
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12)/(m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12)/(m + 1e-12))))
  return(js)
}

# Generalized statistical complexity based on a given entropy value
GeneralizedComplexity <- function(prob, entropy_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * entropy_value)
}

# Fisher-based statistical complexity
FisherBasedComplexity <- function(prob, fisher_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * fisher_value)
}

# ==============================================================================
# PURE LOGISTIC MAP
# ==============================================================================
generate_logistic_map <- function(n, r_log = 3.6, burnin = 500) {
  x <- numeric(n + burnin)
  x[1] <- 0.5  # Initial condition
  
  for (i in 2:(n + burnin)) {
    x[i] <- r_log * x[i-1] * (1 - x[i-1])
  }
  
  # Remove burnin period
  x <- x[(burnin + 1):(n + burnin)]
  return(x)
}

# ==============================================================================
# PURE ARMA MODEL
# ==============================================================================
generate_arma <- function(n, ar, ma) {
  arma_model <- list(ar = ar, ma = ma)
  arma_ts <- arima.sim(model = arma_model, n = n)
  return(as.numeric(arma_ts))
}

# ==============================================================================
# HYBRID MODEL: r*Logistic + (1-r)*ARMA
# ==============================================================================
generate_hybrid_logistic_arma <- function(n, r, ar, ma, burnin = 500) {
  # Generate ARMA component
  arma_model <- list(ar = ar, ma = ma)
  arma_ts <- arima.sim(model = arma_model, n = n + burnin)
  print(arma_ts)
  
  # Normalize ARMA to [0,1] for mixing with logistic map
  arma_norm <- (arma_ts - min(arma_ts)) / (max(arma_ts) - min(arma_ts))
  
  # Initialize logistic map
  x <- numeric(n + burnin)
  x[1] <- 0.5  # Initial condition
  
  # Generate hybrid time series
  for (i in 2:(n + burnin)) {
    logistic_component <- r_logistic * x[i-1] * (1 - x[i-1])
    x[i] <- r * logistic_component + (1 - r) * arma_norm[i]
    # Keep in [0,1] range
    x[i] <- pmax(0, pmin(1, x[i]))
  }
  
  # Remove burnin period
  x <- x[(burnin + 1):(n + burnin)]
  
  return(x)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================
generate_model_data <- function(n, D, R, model_type, ar = NULL, ma = NULL, 
                                r = NULL, beta = BETA) {
  results   <- list()
  ts_store  <- list()
  z_alpha   <- qnorm(1 - 0.05 / 2)   # z for 95% CI
  
  for (rep in 1:R) {
    cat(sprintf("Processing %s, n=%d, repetition %d/%d\n", model_type, n, rep, R))
    
    #-----------------------------
    # Generate time series based on model type
    #-----------------------------
    if (model_type == "Logistic") {
      ts_data <- generate_logistic_map(n = n, r_log = r_logistic)
    } else if (grepl("^ARMA|^AR|^MA", model_type) && !grepl("Hybrid", model_type)) {
      ts_data <- generate_arma(n = n, ar = ar, ma = ma)
    } else {  # Hybrid models
      ts_data <- generate_hybrid_logistic_arma(n = n, r = r, ar = ar, ma = ma)
    }
    
    ProbTS  <- OPprob(ts_data, emb = D)
    
    #-----------------------------
    # Entropies
    #-----------------------------
    Hs <- HShannon(ProbTS)
    Hr <- HRenyi(ProbTS, beta = beta)
    Ht <- HTsallis(ProbTS, beta = beta)
    Hf <- HFisher(ProbTS)
    
    #-----------------------------
    # Complexities
    #-----------------------------
    C_Shannon <- StatComplexity(ProbTS)
    C_Renyi   <- GeneralizedComplexity(ProbTS, Hr)
    C_Tsallis <- GeneralizedComplexity(ProbTS, Ht)
    C_Fisher  <- FisherBasedComplexity(ProbTS, Hf)
    
    #-----------------------------
    # Variances of entropies via sigma2q
    #-----------------------------
    Var_HS <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HR <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "R", beta = beta))
    Var_HT <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "T", beta = beta))
    Var_HF <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "F"))
    
    #-----------------------------
    # Shannon complexity variance (distance-based)
    #-----------------------------
    Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(ProbTS, n - (D-1)))
    a_ratio <- Var_HS / Var_HI
    Var_CI  <- suppressWarnings(varC(ProbTS, n - (D-1)))
    Var_CS  <- a_ratio * Var_CI
    
    #-----------------------------
    # Semi-lengths for entropies & Shannon complexity
    #-----------------------------
    Semi_HS <- ifelse(!is.finite(Var_HS) | Var_HS <= 0, NA,
                      sqrt(Var_HS) / sqrt(n - (D-1)) * z_alpha)
    Semi_HR <- ifelse(!is.finite(Var_HR) | Var_HR <= 0, NA,
                      sqrt(Var_HR) / sqrt(n - (D-1)) * z_alpha)
    Semi_HT <- ifelse(!is.finite(Var_HT) | Var_HT <= 0, NA,
                      sqrt(Var_HT) / sqrt(n - (D-1)) * z_alpha)
    Semi_HF <- ifelse(!is.finite(Var_HF) | Var_HF <= 0, NA,
                      sqrt(Var_HF) / sqrt(n - (D-1)) * z_alpha)
    Semi_CS <- ifelse(!is.finite(Var_CS) | Var_CS <= 0, NA,
                      sqrt(Var_CS) / sqrt(n - (D-1)) * z_alpha)
    
    #-----------------------------
    # Store results
    #-----------------------------
    results[[rep]] <- data.frame(
      Model = model_type,
      n = n,
      r_mix = ifelse(is.null(r), NA, r),
      rep = rep,
      H_Shannon = Hs,
      H_Renyi = Hr,
      H_Tsallis = Ht,
      H_Fisher = Hf,
      C_Shannon = C_Shannon,
      C_Renyi = C_Renyi,
      C_Tsallis = C_Tsallis,
      C_Fisher = C_Fisher,
      Var_HS = Var_HS,
      Var_HR = Var_HR,
      Var_HT = Var_HT,
      Var_HF = Var_HF,
      Var_CS = Var_CS,
      Semi_HS = Semi_HS,
      Semi_HR = Semi_HR,
      Semi_HT = Semi_HT,
      Semi_HF = Semi_HF,
      Semi_CS = Semi_CS
    )
    
    # Store time series
    ts_store[[rep]] <- ts_data
  }
  
  return(list(results = do.call(rbind, results), 
              timeseries = ts_store))
}

# ==============================================================================
# RUN ANALYSIS FOR ALL 13 MODELS
# ==============================================================================
all_results <- list()
all_timeseries <- list()

for (n_size in N) {
  cat(sprintf("\n========== SAMPLE SIZE n = %d ==========\n", n_size))
  
  # 1. Pure ARMA models (6 models)
  for (model_name in names(arma_models)) {
    model_info <- arma_models[[model_name]]
    
    cat(sprintf("\n--- Processing %s ---\n", model_info$name))
    output <- generate_model_data(
      n = n_size,
      D = D,
      R = R,
      model_type = model_info$name,
      ar = model_info$ar,
      ma = model_info$ma,
      r = NULL,
      beta = BETA
    )
    
    key <- paste0(model_info$name, "_n", n_size)
    all_results[[key]] <- output$results
    all_timeseries[[key]] <- output$timeseries
  }
  
  # 2. Pure Logistic Map (1 model)
  cat(sprintf("\n--- Processing Logistic Map ---\n"))
  output <- generate_model_data(
    n = n_size,
    D = D,
    R = R,
    model_type = "Logistic",
    ar = NULL,
    ma = NULL,
    r = NULL,
    beta = BETA
  )
  
  key <- paste0("Logistic_n", n_size)
  all_results[[key]] <- output$results
  all_timeseries[[key]] <- output$timeseries
  
  # 3. Hybrid models: Logistic + ARMA (6 models)
  for (model_name in names(arma_models)) {
    model_info <- arma_models[[model_name]]
    hybrid_name <- paste0("Hybrid_", model_info$name)
    
    cat(sprintf("\n--- Processing %s ---\n", hybrid_name))
    output <- generate_model_data(
      n = n_size,
      D = D,
      R = R,
      model_type = hybrid_name,
      ar = model_info$ar,
      ma = model_info$ma,
      r = r_mix,
      beta = BETA
    )
    
    key <- paste0(hybrid_name, "_n", n_size)
    all_results[[key]] <- output$results
    all_timeseries[[key]] <- output$timeseries
  }
}

# ==============================================================================
# SAVE RESULTS TO SPECIFIED LOCATION
# ==============================================================================
# Define base path
base_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Hybrid ARMA+LogisticMap_data/D4_Data_r_3.6"

# Create directory if it doesn't exist
if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
  cat(sprintf("Created directory: %s\n", base_path))
}

# ==============================================================================
# SAVE HC RESULTS TO EXCEL (All models in one file, different sheets)
# ==============================================================================
summary_path <- file.path(base_path, "All_Models_HC_Results_D4_r_3_6.xlsx")
write_xlsx(all_results, path = summary_path)
cat(sprintf("\nHC Results saved to: %s\n", summary_path))
cat(sprintf("Total sheets: %d\n", length(all_results)))

# ==============================================================================
# COMBINE ALL SHEETS INTO ONE SINGLE SHEET
# ==============================================================================
cat("\n========== COMBINING ALL SHEETS INTO ONE ==========\n")

# Combine all results into one dataframe
combined_all_results <- bind_rows(all_results, .id = "Sheet_Name")

# Reorder columns to put Sheet_Name and Model first
combined_all_results <- combined_all_results %>%
  select(Sheet_Name, Model, n, r_mix, rep, everything())

# Save combined results to Excel
combined_path <- file.path(base_path, "Combined_All_Models_HC_Results_D4_r_3_6.xlsx")
write_xlsx(list(All_Models_Combined = combined_all_results), path = combined_path)

cat(sprintf("Combined HC Results saved to: %s\n", combined_path))
cat(sprintf("Total rows: %d\n", nrow(combined_all_results)))
cat(sprintf("Models included: %d\n", length(unique(combined_all_results$Model))))

# Optional: Save as CSV as well for easier analysis
combined_csv_path <- file.path(base_path, "Combined_All_Models_HC_Results_D4_r_3_6.csv")
write.csv(combined_all_results, file = combined_csv_path, row.names = FALSE)
cat(sprintf("Combined HC Results also saved as CSV to: %s\n", combined_csv_path))


# ==============================================================================
# SAVE HC RESULTS TO EXCEL (All models in one file, different sheets)
# ==============================================================================
summary_path <- file.path(base_path, "All_Models_HC_Results_D4_r_3_6.xlsx")
write_xlsx(all_results, path = summary_path)
cat(sprintf("\nHC Results saved to: %s\n", summary_path))
cat(sprintf("Total sheets: %d\n", length(all_results)))


# ==============================================================================
# SAVE TIME SERIES AS CSV FILES
# ==============================================================================
# Create subdirectory for time series
ts_dir <- file.path(base_path, "All_Models_TimeSeries_D4_r_3_6")
if (!dir.exists(ts_dir)) {
  dir.create(ts_dir)
}

# Save each time series as a separate CSV file
for (model_key in names(all_timeseries)) {
  ts_list <- all_timeseries[[model_key]]
  
  for (rep in 1:R) {
    # Create filename with model name
    filename <- file.path(ts_dir, sprintf("%s_rep%03d.csv", model_key, rep))
    
    # Save time series as CSV with time index
    ts_df <- data.frame(
      time = 1:length(ts_list[[rep]]),
      value = ts_list[[rep]]
    )
    write.csv(ts_df, file = filename, row.names = FALSE)
  }
  
  cat(sprintf("Saved %d time series for %s\n", R, model_key))
}

cat(sprintf("\nAll files saved successfully!\n"))
cat(sprintf("Location: %s\n", base_path))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
cat("\n========== SUMMARY STATISTICS ==========\n")

# Combine all results into one dataframe for summary
combined_results <- bind_rows(all_results)

# Summary by Model and Sample Size
summary_stats <- combined_results %>%
  group_by(Model, n) %>%
  summarise(
    Mean_H_Shannon = mean(H_Shannon, na.rm = TRUE),
    SD_H_Shannon = sd(H_Shannon, na.rm = TRUE),
    Mean_C_Shannon = mean(C_Shannon, na.rm = TRUE),
    SD_C_Shannon = sd(C_Shannon, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# Save summary statistics
summary_stats_path <- file.path(base_path, "Summary_Statistics_D4_r_3_6.xlsx")
write_xlsx(list(Summary = summary_stats), path = summary_stats_path)
cat(sprintf("\nSummary statistics saved to: %s\n", summary_stats_path))

#----------------------------------------------------------------------------
# The new code below is for generating the ordinal pattern distributions 
# and plotting them. It reads the time series data for each model, 
# calculates the ordinal pattern probabilities, and creates faceted bar 
# plots of the distributions for each model. It also creates a combined plot 
# with all models overlaid for comparison. Finally, it generates a summary 
# table of the ordinal pattern statistics for each model.
#----------------------------------------------------------------------------
library(StatOrdPattHxC)
library(dplyr)
library(writexl)

# ==============================================================================
# PARAMETERS
# ==============================================================================
set.seed(1234567890, kind = "Mersenne-Twister")
BETA <- 1.5               # Parameter for Rényi / Tsallis
N    <- c(1000, 5000)     # Sample sizes
D    <- 3                 # Embedding dimension
R    <- 50                # Number of repetitions
r1 <- 3.8                 # Logistic map parameter (chaotic regime) - for original models
r2 <- 3.6                 # Logistic map parameter (different regime) - for new model
alpha1 <- 0.5             # Mixing parameter for hybrid ARMA models
alpha2 <- 0.7             # Mixing parameter for logistic+sine model
f <- 0.05                 # frequency for sine wave (f = 1/T)

# ==============================================================================
# DEFINE ALL ARMA MODEL CONFIGURATIONS
# ==============================================================================
arma_models <- list(
  ARMA22 = list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), name = "ARMA(2,2)"),
  ARMA11 = list(ar = c(0.8), ma = c(0.8), name = "ARMA(1,1)"),
  AR2    = list(ar = c(0.1, 0.8), ma = NULL, name = "AR(2)"),
  AR1    = list(ar = c(-0.1), ma = NULL, name = "AR(1)"),
  MA2    = list(ar = NULL, ma = c(0.1, -0.8), name = "MA(2)"),
  MA1    = list(ar = NULL, ma = c(0.8), name = "MA(1)")
)

# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Jensen–Shannon divergence between two discrete pmfs p, q
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12)/(m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12)/(m + 1e-12))))
  return(js)
}

# Generalized statistical complexity based on a given entropy value
GeneralizedComplexity <- function(prob, entropy_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * entropy_value)
}

# Fisher-based statistical complexity
FisherBasedComplexity <- function(prob, fisher_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * fisher_value)
}

# ==============================================================================
# PURE LOGISTIC MAP
# ==============================================================================
generate_logistic_map <- function(n, r_log = 3.8, burnin = 500) {
  x <- numeric(n + burnin)
  x[1] <- 0.5  # Initial condition
  
  for (i in 2:(n + burnin)) {
    x[i] <- r_log * x[i-1] * (1 - x[i-1])
  }
  
  # Remove burnin period
  x <- x[(burnin + 1):(n + burnin)]
  return(x)
}

# ==============================================================================
# SINE WAVE FUNCTION
# ==============================================================================
generate_sine_wave <- function(n, f = 0.05) {
  t <- 1:n
  x <- sin(2 * pi * f * t)
  # Normalize to [0, 1] range
  x <- (x + 1) / 2
  return(x)
}

# ==============================================================================
# COMBINED LOGISTIC + SINE MODEL
# ==============================================================================
generate_logistic_sine_combined <- function(n, alpha = 0.7, r_log = 3.6, f = 0.05, burnin = 500) {
  # Generate logistic map component
  logistic <- generate_logistic_map(n = n, r_log = r_log, burnin = burnin)
  
  # Generate sine wave component
  sine <- generate_sine_wave(n = n, f = f)
  
  # Combine: alpha*logistic + (1-alpha)*sine
  combined <- alpha * logistic + (1 - alpha) * sine
  
  return(combined)
}

# ==============================================================================
# PURE ARMA MODEL
# ==============================================================================
generate_arma <- function(n, ar, ma) {
  arma_model <- list(ar = ar, ma = ma)
  arma_ts <- arima.sim(model = arma_model, n = n)
  return(as.numeric(arma_ts))
}

# ==============================================================================
# HYBRID MODEL: alpha*Logistic + (1-alpha)*ARMA
# ==============================================================================
generate_hybrid_logistic_arma <- function(n, alpha, ar, ma, r_log = 3.8, burnin = 500) {
  # Generate ARMA component
  arma_model <- list(ar = ar, ma = ma)
  arma_ts <- arima.sim(model = arma_model, n = n + burnin)
  
  # Normalize ARMA to [0,1] for mixing with logistic map
  arma_norm <- (arma_ts - min(arma_ts)) / (max(arma_ts) - min(arma_ts))
  
  # Initialize logistic map
  x <- numeric(n + burnin)
  x[1] <- 0.5  # Initial condition
  
  # Generate hybrid time series
  for (i in 2:(n + burnin)) {
    logistic_component <- r_log * x[i-1] * (1 - x[i-1])
    x[i] <- alpha * logistic_component + (1 - alpha) * arma_norm[i]
    # Keep in [0,1] range
    x[i] <- pmax(0, pmin(1, x[i]))
  }
  
  # Remove burnin period
  x <- x[(burnin + 1):(n + burnin)]
  
  return(x)
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================
generate_model_data <- function(n, D, R, model_type, ar = NULL, ma = NULL, 
                                alpha = NULL, beta = BETA, r_log = NULL, 
                                f = NULL) {
  results   <- list()
  ts_store  <- list()
  z_alpha   <- qnorm(1 - 0.05 / 2)   # z for 95% CI
  
  for (rep in 1:R) {
    cat(sprintf("Processing %s, n=%d, repetition %d/%d\n", model_type, n, rep, R))
    
    #-----------------------------
    # Generate time series based on model type
    #-----------------------------
    if (model_type == "Logistic") {
      # Original logistic with r=3.8
      ts_data <- generate_logistic_map(n = n, r_log = r1)
      
    } else if (model_type == "Logistic_r3_6") {
      # NEW: Logistic with r=3.6
      ts_data <- generate_logistic_map(n = n, r_log = r2)
      
    } else if (model_type == "Sine_Wave") {
      # NEW: Sine wave
      ts_data <- generate_sine_wave(n = n, f = f)
      
    } else if (model_type == "Logistic_Sine_Combined") {
      # NEW: Combined logistic (r=3.6) + sine
      ts_data <- generate_logistic_sine_combined(n = n, alpha = alpha, 
                                                 r_log = r2, f = f)
      
    } else if (grepl("^ARMA|^AR|^MA", model_type) && !grepl("Hybrid", model_type)) {
      # Pure ARMA models
      ts_data <- generate_arma(n = n, ar = ar, ma = ma)
      
    } else if (grepl("Hybrid", model_type)) {
      # Hybrid ARMA + Logistic models (use r1 = 3.8)
      ts_data <- generate_hybrid_logistic_arma(n = n, alpha = alpha, ar = ar, ma = ma, r_log = r1)
      
    } else {
      stop(sprintf("Unknown model type: %s", model_type))
    }
    
    ProbTS  <- OPprob(ts_data, emb = D)
    
    #-----------------------------
    # Entropies
    #-----------------------------
    Hs <- HShannon(ProbTS)
    Hr <- HRenyi(ProbTS, beta = beta)
    Ht <- HTsallis(ProbTS, beta = beta)
    Hf <- HFisher(ProbTS)
    
    #-----------------------------
    # Complexities
    #-----------------------------
    C_Shannon <- StatComplexity(ProbTS)
    C_Renyi   <- GeneralizedComplexity(ProbTS, Hr)
    C_Tsallis <- GeneralizedComplexity(ProbTS, Ht)
    C_Fisher  <- FisherBasedComplexity(ProbTS, Hf)
    
    #-----------------------------
    # Variances of entropies via sigma2q
    #-----------------------------
    Var_HS <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HR <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "R", beta = beta))
    Var_HT <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "T", beta = beta))
    Var_HF <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "F"))
    
    #-----------------------------
    # Shannon complexity variance (distance-based)
    #-----------------------------
    Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(ProbTS, n - (D-1)))
    a_ratio <- Var_HS / Var_HI
    Var_CI  <- suppressWarnings(varC(ProbTS, n - (D-1)))
    Var_CS  <- a_ratio * Var_CI
    
    #-----------------------------
    # Semi-lengths for entropies & Shannon complexity
    #-----------------------------
    Semi_HS <- ifelse(!is.finite(Var_HS) | Var_HS <= 0, NA,
                      sqrt(Var_HS) / sqrt(n - (D-1)) * z_alpha)
    Semi_HR <- ifelse(!is.finite(Var_HR) | Var_HR <= 0, NA,
                      sqrt(Var_HR) / sqrt(n - (D-1)) * z_alpha)
    Semi_HT <- ifelse(!is.finite(Var_HT) | Var_HT <= 0, NA,
                      sqrt(Var_HT) / sqrt(n - (D-1)) * z_alpha)
    Semi_HF <- ifelse(!is.finite(Var_HF) | Var_HF <= 0, NA,
                      sqrt(Var_HF) / sqrt(n - (D-1)) * z_alpha)
    Semi_CS <- ifelse(!is.finite(Var_CS) | Var_CS <= 0, NA,
                      sqrt(Var_CS) / sqrt(n - (D-1)) * z_alpha)
    
    #-----------------------------
    # Store results
    #-----------------------------
    results[[rep]] <- data.frame(
      Model = model_type,
      n = n,
      r_mix = ifelse(is.null(alpha), NA, alpha),
      rep = rep,
      H_Shannon = Hs,
      H_Renyi = Hr,
      H_Tsallis = Ht,
      H_Fisher = Hf,
      C_Shannon = C_Shannon,
      C_Renyi = C_Renyi,
      C_Tsallis = C_Tsallis,
      C_Fisher = C_Fisher,
      Var_HS = Var_HS,
      Var_HR = Var_HR,
      Var_HT = Var_HT,
      Var_HF = Var_HF,
      Var_CS = Var_CS,
      Semi_HS = Semi_HS,
      Semi_HR = Semi_HR,
      Semi_HT = Semi_HT,
      Semi_HF = Semi_HF,
      Semi_CS = Semi_CS
    )
    
    # Store time series
    ts_store[[rep]] <- ts_data
  }
  
  return(list(results = do.call(rbind, results), 
              timeseries = ts_store))
}

# ==============================================================================
# RUN ANALYSIS FOR ALL 16 MODELS (13 original + 3 new)
# ==============================================================================
all_results <- list()
all_timeseries <- list()

for (n_size in N) {
  cat(sprintf("\n========== SAMPLE SIZE n = %d ==========\n", n_size))
  
  # 1. Pure ARMA models (6 models)
  for (model_name in names(arma_models)) {
    model_info <- arma_models[[model_name]]
    
    cat(sprintf("\n--- Processing %s ---\n", model_info$name))
    output <- generate_model_data(
      n = n_size,
      D = D,
      R = R,
      model_type = model_info$name,
      ar = model_info$ar,
      ma = model_info$ma,
      beta = BETA
    )
    
    key <- paste0(model_info$name, "_n", n_size)
    all_results[[key]] <- output$results
    all_timeseries[[key]] <- output$timeseries
  }
  
  # 2. Pure Logistic Map r1=3.8 (1 model)
  cat(sprintf("\n--- Processing Logistic Map (r=%.1f) ---\n", r1))
  output <- generate_model_data(
    n = n_size,
    D = D,
    R = R,
    model_type = "Logistic",
    beta = BETA
  )
  
  key <- paste0("Logistic_n", n_size)
  all_results[[key]] <- output$results
  all_timeseries[[key]] <- output$timeseries
  
  # 3. Hybrid models: alpha1*Logistic(r1) + (1-alpha1)*ARMA (6 models)
  for (model_name in names(arma_models)) {
    model_info <- arma_models[[model_name]]
    hybrid_name <- paste0("Hybrid_", model_info$name)
    
    cat(sprintf("\n--- Processing %s (alpha=%.1f, r=%.1f) ---\n", hybrid_name, alpha1, r1))
    output <- generate_model_data(
      n = n_size,
      D = D,
      R = R,
      model_type = hybrid_name,
      ar = model_info$ar,
      ma = model_info$ma,
      alpha = alpha1,  # Use alpha1 for hybrid ARMA models
      beta = BETA
    )
    
    key <- paste0(hybrid_name, "_n", n_size)
    all_results[[key]] <- output$results
    all_timeseries[[key]] <- output$timeseries
  }
  
  # ===========================================================================
  # NEW MODELS (3 additional models)
  # ===========================================================================
  
  # 4. NEW: Logistic Map with r2=3.6
  cat(sprintf("\n--- Processing Logistic Map (r=%.1f) ---\n", r2))
  output <- generate_model_data(
    n = n_size,
    D = D,
    R = R,
    model_type = "Logistic_r3_6",
    r_log = r2,
    beta = BETA
  )
  
  key <- paste0("Logistic_r3_6_n", n_size)
  all_results[[key]] <- output$results
  all_timeseries[[key]] <- output$timeseries
  
  # 5. NEW: Sine Wave 
  cat(sprintf("\n--- Processing Sine Wave (f=%.2f) ---\n", f))
  output <- generate_model_data(
    n = n_size,
    D = D,
    R = R,
    model_type = "Sine_Wave",
    f = f,
    beta = BETA
  )
  
  key <- paste0("Sine_Wave_n", n_size)
  all_results[[key]] <- output$results
  all_timeseries[[key]] <- output$timeseries
  
  # 6. NEW: Combined alpha2*Logistic(r2) + (1-alpha2)*Sine
  cat(sprintf("\n--- Processing Logistic+Sine Combined (alpha=%.1f, r=%.1f, f=%.2f) ---\n", 
              alpha2, r2, f))
  output <- generate_model_data(
    n = n_size,
    D = D,
    R = R,
    model_type = "Logistic_Sine_Combined",
    alpha = alpha2,   # Use alpha2 for logistic+sine model
    f = f,
    beta = BETA
  )
  
  key <- paste0("Logistic_Sine_Combined_n", n_size)
  all_results[[key]] <- output$results
  all_timeseries[[key]] <- output$timeseries
}

# ==============================================================================
# SAVE RESULTS TO SPECIFIED LOCATION
# ==============================================================================
base_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Hybrid model_data/D3_Data"

if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
  cat(sprintf("Created directory: %s\n", base_path))
}

# ==============================================================================
# SAVE HC RESULTS TO EXCEL
# ==============================================================================
summary_path <- file.path(base_path, paste0("All_Models_HC_Results_D", D, ".xlsx"))
write_xlsx(all_results, path = summary_path)
cat(sprintf("\nHC Results saved to: %s\n", summary_path))
cat(sprintf("Total sheets: %d\n", length(all_results)))

# ==============================================================================
# COMBINE ALL SHEETS INTO ONE
# ==============================================================================
cat("\n========== COMBINING ALL SHEETS INTO ONE ==========\n")

combined_all_results <- bind_rows(all_results, .id = "Sheet_Name")

combined_all_results <- combined_all_results %>%
  select(Sheet_Name, Model, n, r_mix, rep, everything())

combined_path <- file.path(base_path, paste0("Combined_All_Models_HC_Results_D", D, ".xlsx"))
write_xlsx(list(All_Models_Combined = combined_all_results), path = combined_path)

cat(sprintf("Combined HC Results saved to: %s\n", combined_path))
cat(sprintf("Total rows: %d\n", nrow(combined_all_results)))
cat(sprintf("Models included: %d\n", length(unique(combined_all_results$Model))))

combined_csv_path <- file.path(base_path, paste0("Combined_All_Models_HC_Results_D", D, ".csv"))
write.csv(combined_all_results, file = combined_csv_path, row.names = FALSE)
cat(sprintf("Combined HC Results also saved as CSV to: %s\n", combined_csv_path))

# ==============================================================================
# SAVE TIME SERIES AS CSV FILES
# ==============================================================================
ts_dir <- file.path(base_path, paste0("All_Models_TimeSeries_D", D))
if (!dir.exists(ts_dir)) {
  dir.create(ts_dir)
}

for (model_key in names(all_timeseries)) {
  ts_list <- all_timeseries[[model_key]]
  
  for (rep in 1:R) {
    filename <- file.path(ts_dir, sprintf("%s_rep%03d.csv", model_key, rep))
    
    ts_df <- data.frame(
      time = 1:length(ts_list[[rep]]),
      value = ts_list[[rep]]
    )
    write.csv(ts_df, file = filename, row.names = FALSE)
  }
  
  cat(sprintf("Saved %d time series for %s\n", R, model_key))
}

cat(sprintf("\nAll files saved successfully!\n"))
cat(sprintf("Location: %s\n", base_path))

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
cat("\n========== SUMMARY STATISTICS ==========\n")

combined_results <- bind_rows(all_results)

summary_stats <- combined_results %>%
  group_by(Model, n) %>%
  summarise(
    Mean_H_Shannon = mean(H_Shannon, na.rm = TRUE),
    SD_H_Shannon = sd(H_Shannon, na.rm = TRUE),
    Mean_C_Shannon = mean(C_Shannon, na.rm = TRUE),
    SD_C_Shannon = sd(C_Shannon, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

summary_stats_path <- file.path(base_path, paste0("Summary_Statistics_16_Models_D", D, ".xlsx"))
write_xlsx(list(Summary = summary_stats), path = summary_stats_path)
cat(sprintf("\nSummary statistics saved to: %s\n", summary_stats_path))


#-----------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------
# This script generate a ordinal pattern for selected models and sample sizes, and saves 
#them as CSV files for later use in plotting.
#install.packages("here")
library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

# ==============================================================================
# PARAMETERS
# ==============================================================================
D <- 3  # Embedding dimension

ts_dir <- here("Data", "Hybrid model_data", "D3_Data", "All_Models_TimeSeries_D3")
output_dir <- here("Plots", "Hybrid Model plots", "D3 plots")

# ==============================================================================
# DEFINE MODEL ORDER, COLORS, AND SHAPES (CONSISTENT FOR ALL FUTURE WORK)
# ==============================================================================

# Define the 16 models in the order you want them displayed
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", 
  "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

# FIXED COLOR PALETTE (use these colors consistently in all future plots)
model_colors <- c(
  "ARMA(2,2)" = "#D55E00",   # Vermillion
  "ARMA(1,1)" = "#0072B2",   # Blue
  "AR(2)" = "#009E73",       # Bluish green
  "AR(1)" = "#CC79A7",       # Reddish purple
  "MA(2)" = "#E69F00",       # Orange
  "MA(1)" = "#56B4E9",       # Sky blue
  "Logistic" = "#000000",    # Black
  
  "Hybrid_ARMA(2,2)" = "#8B4513", # Dark brown
  "Hybrid_ARMA(1,1)" = "#4B0082", # Indigo
  "Hybrid_AR(2)" = "#696969",     # Dark gray
  "Hybrid_AR(1)" = "#20B2AA",     # Light sea green
  "Hybrid_MA(2)" = "#B22222",     # Firebrick
  "Hybrid_MA(1)" = "#4169E1",     # Royal blue
  
  "Logistic_r3_6" = "#2F4F4F",    # Slate gray
  "Sine_Wave" = "#228B22",        # Forest green
  "Logistic_Sine_Combined" = "#8A2BE2" # Blue violet
)


# FIXED SHAPE PALETTE (use these shapes consistently in all future plots)
model_shapes <- c(
  "ARMA(2,2)" = 16,
  "ARMA(1,1)" = 17,
  "AR(2)" = 15,
  "AR(1)" = 18,
  "MA(2)" = 3,
  "MA(1)" = 7,
  "Logistic" = 8,
  
  "Hybrid_ARMA(2,2)" = 1,
  "Hybrid_ARMA(1,1)" = 2,
  "Hybrid_AR(2)" = 0,
  "Hybrid_AR(1)" = 5,
  "Hybrid_MA(2)" = 6,
  "Hybrid_MA(1)" = 4,
  
  "Logistic_r3_6" = 9,
  "Sine_Wave" = 10,
  "Logistic_Sine_Combined" = 11
)


# ==============================================================================
# READ TIME SERIES AND CALCULATE ORDINAL PATTERN PROBABILITIES
# ==============================================================================
all_op_probs <- list()

for (model in model_names) {
  # Construct filename
  filename <- file.path(ts_dir, paste0(model, "_n5000_rep001.csv"))
  
  # Check if file exists
  if (!file.exists(filename)) {
    cat(sprintf("Warning: File not found - %s\n", filename))
    next
  }
  
  # Read time series
  ts_data <- read.csv(filename)
  ts_values <- ts_data$value
  
  # Calculate ordinal pattern probabilities
  op_prob <- OPprob(ts_values, emb = D)
  
  # Store results
  all_op_probs[[model]] <- data.frame(
    Model = model,
    Pattern = 1:length(op_prob),
    Probability = op_prob
  )
  
  cat(sprintf("Processed: %s (n = %d)\n", model, length(ts_values)))
}

# Combine all results
combined_op <- bind_rows(all_op_probs)

# Factor the Model variable to control plotting order
combined_op$Model <- factor(combined_op$Model, levels = model_names)

# ==============================================================================
# CREATE FACETED ORDINAL PATTERNS PLOT
# ==============================================================================

# Calculate number of possible patterns
n_patterns <- factorial(D)

# Create the faceted plot
p_faceted <- ggplot(combined_op, aes(x = Pattern, y = Probability, fill = Model)) +
  geom_bar(stat = "identity", position = "identity", width = 0.8) +
  facet_wrap(~ Model, ncol = 3, scales = "free_y") +
  theme_bw(base_family = "serif", base_size = 13) +
  labs(
    title = paste0("Ordinal Pattern Distributions (D = ", D, ", n = 500)"),
    #subtitle = "First repetition of each model",
    x = "Ordinal Pattern",
    y = "Probability"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, n_patterns, by = 4)) +
  scale_fill_manual(values = model_colors)

# Save the faceted plot
ggsave(
  filename = file.path(output_dir, "Ordinal_Patterns_Faceted_D3.png"),
  plot = p_faceted,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "Ordinal_Patterns_Faceted_D3.pdf"),
  plot = p_faceted,
  width = 14,
  height = 10
)

cat(sprintf("\nFaceted plot saved to: %s\n", 
            file.path(output_dir, "Ordinal_Patterns_Faceted_D3.png")))

# Display the plot
print(p_faceted)

# ==============================================================================
# CREATE COMBINED PLOT (All models overlaid with lines and points)
# ==============================================================================

p_combined <- ggplot(combined_op, aes(x = Pattern, y = Probability, 
                                      color = Model, shape = Model, group = Model)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_family = "serif", base_size = 13) +
  labs(
    title = paste0("Ordinal Pattern Distributions - All Models (D = ", D, ", n = 500)"),
    #subtitle = "First repetition of each model",
    x = "Ordinal Pattern",
    y = "Probability",
    color = "Model",
    shape = "Model"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.8, "cm"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, n_patterns, by = 2)) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  guides(
    color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1))
  )

ggsave(
  filename = file.path(output_dir, "Ordinal_Patterns_Hybrid_Models_D3.pdf"),
  plot = p_combined,
  width = 14,
  height = 8
)

print(p_combined)

# ==============================================================================
# CREATE SUMMARY TABLE OF ORDINAL PATTERNS
# ==============================================================================
summary_table <- combined_op %>%
  group_by(Model) %>%
  summarise(
    N_Nonzero_Patterns = sum(Probability > 0),
    Total_Patterns = n_patterns,
    Coverage_Percent = round(100 * sum(Probability > 0) / n_patterns, 2),
    Max_Probability = round(max(Probability), 4),
    Min_Nonzero_Prob = round(min(Probability[Probability > 0]), 6),
    Shannon_Entropy = round(-sum(ifelse(Probability > 0, 
                                        Probability * log(Probability), 0)), 4),
    .groups = 'drop'
  )

print(summary_table)

# Save summary table
write.csv(summary_table, 
          file = file.path(output_dir, "Ordinal_Patterns_Summary_D3.csv"),
          row.names = FALSE)

cat(sprintf("\nSummary table saved to: %s\n", 
            file.path(output_dir, "Ordinal_Patterns_Summary_D3.csv")))

# ==============================================================================
# SAVE COLOR AND SHAPE REFERENCE FOR FUTURE USE
# ==============================================================================

# Create a reference table
#color_shape_reference <- data.frame(
#  Model = model_names,
#  Color = model_colors[model_names],
#  Color_Name = c("Red", "Blue", "Green", "Purple", "Orange", "Yellow",
#                 "Black", "Pink", "Brown", "Gray", "Teal", "Coral", "Light Blue"),
#  Shape = model_shapes[model_names],
#  Shape_Name = c("Filled Circle", "Filled Triangle", "Filled Square", "Filled Diamond",
#                 "Filled Circle", "Filled Triangle", "Star",
#                 "Open Circle", "Open Triangle", "Open Square", "Open Diamond",
#                 "Inverted Triangle", "Plus")
#)

#----------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# This script generates time series plots versus HC plane. It reads the time series data for each model, 
# read the corresponding HC values, and creates scatter plots of entropy vs. complexity, colored by 
# model type. 

library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(gridExtra)
library(grid)
library(gtable)
library(here)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 3  # Embedding dimension

# Paths
ts_dir <- here("Data", "Hybrid model_data", "D4_Data", "All_Models_TimeSeries_D4")
hc_data_path <- here("Data", "Hybrid model_data", "D4_Data", "Combined_All_Models_HC_Results_D3.xlsx")
output_dir <- here("Plots", "Hybrid Model plots", "D4 plots")
# ==============================================================================
# DEFINE MODEL ORDER, COLORS, AND SHAPES (SAME AS BEFORE)
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

model_colors <- c(
  "ARMA(2,2)" = "#D55E00",   # Vermillion
  "ARMA(1,1)" = "#0072B2",   # Blue
  "AR(2)" = "#009E73",       # Bluish green
  "AR(1)" = "#CC79A7",       # Reddish purple
  "MA(2)" = "#E69F00",       # Orange
  "MA(1)" = "#56B4E9",       # Sky blue
  "Logistic" = "#000000",    # Black
  
  "Hybrid_ARMA(2,2)" = "#8B4513", # Dark brown
  "Hybrid_ARMA(1,1)" = "#4B0082", # Indigo
  "Hybrid_AR(2)" = "#696969",     # Dark gray
  "Hybrid_AR(1)" = "#20B2AA",     # Light sea green
  "Hybrid_MA(2)" = "#B22222",     # Firebrick
  "Hybrid_MA(1)" = "#4169E1",     # Royal blue
  
  "Logistic_r3_6" = "#2F4F4F",    # Slate gray
  "Sine_Wave" = "#228B22",        # Forest green
  "Logistic_Sine_Combined" = "#8A2BE2" # Blue violet
)


model_shapes <- c(
  "ARMA(2,2)" = 16,
  "ARMA(1,1)" = 17,
  "AR(2)" = 15,
  "AR(1)" = 18,
  "MA(2)" = 3,
  "MA(1)" = 7,
  "Logistic" = 8,
  
  "Hybrid_ARMA(2,2)" = 1,
  "Hybrid_ARMA(1,1)" = 2,
  "Hybrid_AR(2)" = 0,
  "Hybrid_AR(1)" = 5,
  "Hybrid_MA(2)" = 6,
  "Hybrid_MA(1)" = 4,
  
  "Logistic_r3_6" = 9,
  "Sine_Wave" = 10,
  "Logistic_Sine_Combined" = 11
)


# ==============================================================================
# LOAD HC DATA
# ==============================================================================
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Filter for n=5000 and first repetition
hc_data <- hc_data_all %>%
  filter(n == 5000, rep == 1) %>%
  select(Model, H_Shannon, C_Shannon, Semi_HS, Semi_CS)

# Rename columns for convenience
names(hc_data) <- c("Model", "H_val", "C_val", "Semi_H", "Semi_C")

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

# Sort by model order
hc_data <- hc_data %>% arrange(Model)

print("HC Data loaded:")
print(hc_data)

# ==============================================================================
# LOAD LinfLsup BOUNDARIES
# ==============================================================================
data("LinfLsup")

# OPTION 1: FULL BOUNDARY CURVE (recommended - shows complete theoretical boundaries)
# Filter by Dimension (the LinfLsup uses "Dimension" column as character)
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# OPTION 2: SUBSET BOUNDARY (use only if data points are in limited H range)
# Uncomment these lines if you want to use subset instead:
# H_range <- range(hc_data$H_val, na.rm = TRUE)
# H_buffer <- 0.1  # Add 10% buffer
# H_min <- max(0, H_range[1] - H_buffer)
# H_max <- min(1, H_range[2] + H_buffer)
# 
# LinfLsup_subset <- LinfLsup_full %>%
#   filter(H >= H_min & H <= H_max)
# 
# cat("\nUsing SUBSET boundary curve\n")
# cat("Data H range:", H_range, "\n")
# cat("Subset H range:", c(H_min, H_max), "\n")

# Choose which boundary to use (comment/uncomment as needed)
boundary_data <- LinfLsup_full      # Use this for FULL boundary
# boundary_data <- LinfLsup_subset  # Use this for SUBSET boundary

# Separate upper and lower boundaries
boundary_lower <- boundary_data %>% filter(Side == "Lower")
boundary_upper <- boundary_data %>% filter(Side == "Upper")

# ==============================================================================
# CREATE HC PLANE PLOT WITH CONFIDENCE INTERVALS (NO LEGEND)
# ==============================================================================

# Prepare data for error bars (only where Semi values are valid)
hc_data_ci <- hc_data %>%
  mutate(
    has_H_error = !is.na(Semi_H) & Semi_H > 0 & is.finite(Semi_H),
    has_C_error = !is.na(Semi_C) & Semi_C > 0 & is.finite(Semi_C)
  )

# Create the plot with SEPARATE upper and lower boundaries
p_hc <- ggplot() +
  # Add LOWER boundary (dashed gray)
  geom_line(data = boundary_lower, aes(x = H, y = C), 
            color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
  
  # Add UPPER boundary (dashed gray)
  geom_line(data = boundary_upper, aes(x = H, y = C), 
            color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
  
  # Add horizontal error bars (H_Shannon confidence intervals)
  geom_errorbarh(data = hc_data_ci %>% filter(has_H_error),
                 aes(y = C_val, xmin = H_val - Semi_H, xmax = H_val + Semi_H, 
                     color = Model),
                 height = 0.005, size = 0.8, alpha = 0.8) +
  
  # Add vertical error bars (C_Shannon confidence intervals)
  geom_errorbar(data = hc_data_ci %>% filter(has_C_error),
                aes(x = H_val, ymin = C_val - Semi_C, ymax = C_val + Semi_C, 
                    color = Model),
                width = 0.01, size = 0.8, alpha = 0.8) +
  
  # Add points
  geom_point(data = hc_data, aes(x = H_val, y = C_val, color = Model, shape = Model),
             size = 2, stroke = 1) +
  
  # Styling
  theme_bw(base_family = "serif", base_size = 13) +
  labs(
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(p_hc)

ggsave(
  filename = file.path(output_dir, "HC_Plane_with_CI_D3.pdf"),
  plot = p_hc,
  width = 8,
  height = 8,
  bg = "white"
)

# ==============================================================================
# READ TIME SERIES DATA
# ==============================================================================
all_ts_data <- list()

for (model in model_names) {
  filename <- file.path(ts_dir, paste0(model, "_n5000_rep001.csv"))
  
  if (!file.exists(filename)) {
    cat(sprintf("Warning: File not found - %s\n", filename))
    next
  }
  
  ts_data <- read.csv(filename)
  all_ts_data[[model]] <- ts_data
  
  cat(sprintf("Loaded time series: %s\n", model))
}

# ==============================================================================
# CREATE INDIVIDUAL TIME SERIES PLOTS WITH MODEL NAMES
# ==============================================================================

create_small_ts_plot <- function(model_name, ts_data, color) {
  # Limit to first 500 points for visibility
  plot_data <- ts_data[1:min(500, nrow(ts_data)), ]
  
  p <- ggplot(plot_data, aes(x = time, y = value)) +
    geom_line(color = color, size = 0.6) +
    theme_bw(base_family = "serif", base_size = 9) +
    labs(title = model_name, x = NULL, y = NULL) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, 
                                color = color),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2, color = "gray90"),
      plot.margin = unit(c(4, 4, 4, 4), "pt"),
      panel.border = element_rect(color = color, size = 1.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  
  return(p)
}

# Create all time series plots
ts_plots <- list()
for (model in model_names) {
  if (model %in% names(all_ts_data)) {
    ts_plots[[model]] <- create_small_ts_plot(
      model, 
      all_ts_data[[model]], 
      model_colors[model]
    )
  }
}

# ==============================================================================
# CREATE COMBINED LAYOUT: HC PLANE SURROUNDED BY TIME SERIES
# ==============================================================================

# Define layout matrix - CORRECTED
# You have 17 plots total: plots 1-16 are time series, plot 17 is HC plane
layout_matrix <- rbind(
  c(NA,  1,  2,  3,  4),    # Row 1: Top 5 time series
  c(5,  17, 17, 17, 6),    # Row 2: TS left, HC center (3 cols), TS right
  c(7,  17, 17, 17, 8),    # Row 3: TS left, HC center (3 cols), TS right
  c(9, 17, 17, 17, 10),   # Row 4: TS left, HC center (3 cols), TS right
  c(11, 17, 17, 17, 12),   # Row 5: TS left, HC center (3 cols), TS right
  c(NA, 13, 14, 15, 16)    # Row 6: Bottom 3 time series
)


# Create list of all plots in order
all_plots <- c(
  ts_plots[model_names[1:16]],
  list(p_hc)
)

# Create combined plot with MATCHING dimensions
combined_plot <- grid.arrange(
  grobs = all_plots,
  layout_matrix = layout_matrix,
  widths = c(1.2, 1.2, 1.2, 1.2, 1.2),      # 5 columns: narrow sides, wide center
  heights = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2)      # 6 rows: narrow top/bottom, tall center
)

ggsave(
  filename = file.path(output_dir, "HC_Plane_with_TimeSeries_Surrounded_D3.pdf"),
  plot = combined_plot,
  width = 20,
  height = 16,
  bg = "white"
)

# ==============================================================================
# CREATE SUMMARY INFORMATION
# ==============================================================================

summary_info <- hc_data %>%
  mutate(
    Color = model_colors[as.character(Model)],
    Has_H_CI = !is.na(Semi_H) & Semi_H > 0,
    Has_C_CI = !is.na(Semi_C) & Semi_C > 0,
    H_CI_lower = ifelse(Has_H_CI, H_val - Semi_H, NA),
    H_CI_upper = ifelse(Has_H_CI, H_val + Semi_H, NA),
    C_CI_lower = ifelse(Has_C_CI, C_val - Semi_C, NA),
    C_CI_upper = ifelse(Has_C_CI, C_val + Semi_C, NA)
  ) %>%
  select(Model, Color, H_val, C_val, Has_H_CI, H_CI_lower, H_CI_upper, 
         Has_C_CI, C_CI_lower, C_CI_upper)

write.csv(summary_info,
          file = file.path(output_dir, "HC_Plane_Summary_with_CI_D3.csv"),
          row.names = FALSE)

print(summary_info)

#-----------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#This script generates scatter plots of entropy vs. complexity for various ARMA models,
#using extended entropy measures (Shannon, Rényi, Tsallis, Fisher).
#======================================================================
#----------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(here)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension
#BETA <- 1.5  # Beta parameter used

# Paths
hc_data_path <- here("Data", "Hybrid model_data", "D3_Data", "Combined_All_Models_HC_Results_D3.xlsx")
output_dir <- here("Plots", "Hybrid Model plots", "D3 plots")
# Sample sizes to process
sample_sizes <- c(1000, 5000)

# ==============================================================================
# DEFINE MODEL ORDER, COLORS, AND SHAPES
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)


model_colors <- c(
  "ARMA(2,2)" = "#D55E00",   # Vermillion
  "ARMA(1,1)" = "#0072B2",   # Blue
  "AR(2)" = "#009E73",       # Bluish green
  "AR(1)" = "#CC79A7",       # Reddish purple
  "MA(2)" = "#E69F00",       # Orange
  "MA(1)" = "#56B4E9",       # Sky blue
  "Logistic" = "#000000",    # Black
  
  "Hybrid_ARMA(2,2)" = "#8B4513", # Dark brown
  "Hybrid_ARMA(1,1)" = "#4B0082", # Indigo
  "Hybrid_AR(2)" = "#696969",     # Dark gray
  "Hybrid_AR(1)" = "#20B2AA",     # Light sea green
  "Hybrid_MA(2)" = "#B22222",     # Firebrick
  "Hybrid_MA(1)" = "#4169E1",     # Royal blue
  
  "Logistic_r3_6" = "#2F4F4F",    # Slate gray
  "Sine_Wave" = "#228B22",        # Forest green
  "Logistic_Sine_Combined" = "#8A2BE2" # Blue violet
)


model_shapes <- c(
  "ARMA(2,2)" = 16,
  "ARMA(1,1)" = 17,
  "AR(2)" = 15,
  "AR(1)" = 18,
  "MA(2)" = 3,
  "MA(1)" = 7,
  "Logistic" = 8,
  
  "Hybrid_ARMA(2,2)" = 1,
  "Hybrid_ARMA(1,1)" = 2,
  "Hybrid_AR(2)" = 0,
  "Hybrid_AR(1)" = 5,
  "Hybrid_MA(2)" = 6,
  "Hybrid_MA(1)" = 4,
  
  "Logistic_r3_6" = 9,
  "Sine_Wave" = 10,
  "Logistic_Sine_Combined" = 11
)


# ==============================================================================
# LOAD HC DATA (ALL REPETITIONS)
# ==============================================================================
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Check available sample sizes
available_n <- unique(hc_data_all$n)
cat("Available sample sizes in data:", available_n, "\n")

# Verify we have the expected sample sizes
if (!all(sample_sizes %in% available_n)) {
  cat("WARNING: Not all expected sample sizes found in data!\n")
  cat("Expected:", sample_sizes, "\n")
  cat("Found:", available_n, "\n")
  sample_sizes <- intersect(sample_sizes, available_n)
  cat("Will process:", sample_sizes, "\n")
}

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon, Semi_HS, Semi_CS,
         H_Renyi, C_Renyi, Semi_HR, 
         H_Tsallis, C_Tsallis, Semi_HT,
         H_Fisher, C_Fisher, Semi_HF)

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

cat("\nTotal data points loaded:", nrow(hc_data), "\n")
cat("Sample sizes found:", unique(hc_data$n), "\n")
cat("Repetitions per sample size:\n")
print(table(hc_data$n))

# ==============================================================================
# LOAD LinfLsup BOUNDARIES (for Shannon only)
# ==============================================================================
data("LinfLsup")

# Full boundary curve
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# Separate upper and lower boundaries
boundary_lower <- LinfLsup_full %>% filter(Side == "Lower")
boundary_upper <- LinfLsup_full %>% filter(Side == "Upper")

cat("\nBoundary curves prepared for Shannon entropy\n")
cat("Boundary points - Lower:", nrow(boundary_lower), "Upper:", nrow(boundary_upper), "\n")

# ==============================================================================
# ENTROPY CONFIGURATIONS
# ==============================================================================
entropy_configs <- list(
  Shannon = list(
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_HS",
    Semi_C_col = "Semi_CS",
    add_bounds = TRUE,
    title = "Shannon Entropy-Complexity Plane",
    xlab = expression(italic(H)[Shannon]),
    ylab = expression(italic(C)[Shannon]),
    xlim = c(0, 1),
    ylim = c(0, 0.5)
  ),
  Renyi = list(
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_HR",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Rényi Entropy-Complexity Plane",
    xlab = expression(italic(H)[Rényi]),
    ylab = expression(italic(C)[Rényi]),
    xlim = NULL,
    ylim = NULL
  ),
  Tsallis = list(
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_HT",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Tsallis Entropy-Complexity Plane",
    xlab = expression(italic(H)[Tsallis]),
    ylab = expression(italic(C)[Tsallis]),
    xlim = NULL,
    ylim = NULL
  ),
  Fisher = list(
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_HF",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Fisher Entropy-Complexity Plane",
    xlab = expression(italic(H)[Fisher]),
    ylab = expression(italic(C)[Fisher]),
    xlim = NULL,
    ylim = NULL
  )
)

# ==============================================================================
# CREATE HC PLANE PLOTS FOR EACH SAMPLE SIZE AND ENTROPY TYPE
# ==============================================================================

for (n_val in sample_sizes) {
  
  # Create subdirectory for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n========== Processing Sample Size n = %d ==========\n", n_val))
  
  # Filter data for this sample size
  hc_data_n <- hc_data %>%
    filter(n == n_val)
  
  cat(sprintf("Data points for n=%d: %d\n", n_val, nrow(hc_data_n)))
  cat(sprintf("Models: %d, Max repetitions: %d\n", 
              length(unique(hc_data_n$Model)),
              max(hc_data_n$rep, na.rm = TRUE)))
  
  for (ent_name in names(entropy_configs)) {
    
    cfg <- entropy_configs[[ent_name]]
    
    cat(sprintf("\n=== Creating %s HC plane plot for n=%d ===\n", ent_name, n_val))
    
    # Extract relevant columns
    H_col <- cfg$H_col
    C_col <- cfg$C_col
    Semi_H_col <- cfg$Semi_H_col
    Semi_C_col <- cfg$Semi_C_col
    
    # Check if columns exist
    if (!all(c(H_col, C_col) %in% names(hc_data_n))) {
      cat(sprintf("⚠️ Skipping %s - columns not found\n", ent_name))
      next
    }
    
    # Prepare data with renamed columns for easier plotting
    plot_data <- hc_data_n %>%
      select(Model, rep, all_of(c(H_col, C_col))) %>%
      rename(H_val = !!H_col, C_val = !!C_col)
    
    # Add semi-length columns if they exist
    if (!is.na(Semi_H_col) && Semi_H_col %in% names(hc_data_n)) {
      plot_data$Semi_H <- hc_data_n[[Semi_H_col]]
    } else {
      plot_data$Semi_H <- NA
    }
    
    if (!is.na(Semi_C_col) && Semi_C_col %in% names(hc_data_n)) {
      plot_data$Semi_C <- hc_data_n[[Semi_C_col]]
    } else {
      plot_data$Semi_C <- NA
    }
    
    # Remove rows with NA in H or C
    plot_data <- plot_data %>%
      filter(!is.na(H_val), !is.na(C_val), is.finite(H_val), is.finite(C_val))
    
    if (nrow(plot_data) == 0) {
      cat(sprintf("⚠️ Skipping %s - no valid data\n", ent_name))
      next
    }
    
    cat(sprintf("Valid data points: %d\n", nrow(plot_data)))
    
    # Prepare error bar data (only for first repetition to avoid clutter)
    plot_data_ci <- plot_data %>%
      filter(rep == 1) %>%
      mutate(
        has_H_error = !is.na(Semi_H) & Semi_H > 0 & is.finite(Semi_H),
        has_C_error = !is.na(Semi_C) & Semi_C > 0 & is.finite(Semi_C)
      )
    
    # Create the plot
    p <- ggplot()
    
    # Add boundaries (Shannon only)
    if (cfg$add_bounds) {
      p <- p +
        geom_line(data = boundary_lower, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
        geom_line(data = boundary_upper, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8)
    }
    
    # Add error bars (only for rep 1 to avoid visual clutter)
    # Horizontal error bars (Entropy confidence intervals)
    if (any(plot_data_ci$has_H_error, na.rm = TRUE)) {
      p <- p +
        geom_errorbarh(data = plot_data_ci %>% filter(has_H_error),
                       aes(y = C_val, xmin = H_val - Semi_H, xmax = H_val + Semi_H, 
                           color = Model),
                       height = 0.002, size = 0.6, alpha = 0.5, show.legend = FALSE)
      cat(sprintf("  Added H error bars for %d models\n", 
                  sum(plot_data_ci$has_H_error, na.rm = TRUE)))
    } else {
      cat("  No H error bars (Semi_H not available)\n")
    }
    
    # Vertical error bars (Complexity confidence intervals)
    if (any(plot_data_ci$has_C_error, na.rm = TRUE)) {
      p <- p +
        geom_errorbar(data = plot_data_ci %>% filter(has_C_error),
                      aes(x = H_val, ymin = C_val - Semi_C, ymax = C_val + Semi_C, 
                          color = Model),
                      width = 0.01, size = 0.6, alpha = 0.5, show.legend = FALSE)
      cat(sprintf("  Added C error bars for %d models\n", 
                  sum(plot_data_ci$has_C_error, na.rm = TRUE)))
    } else {
      cat("  No C error bars (Semi_C not available)\n")
    }
    
    # Add points (all repetitions) - COMBINED color and shape in ONE aesthetic
    p <- p +
      geom_point(data = plot_data, 
                 aes(x = H_val, y = C_val, 
                     color = Model, 
                     shape = Model,
                     fill = Model),  # Add fill for shapes 21-25
                 size = 3, stroke = 1.2, alpha = 0.7) +
      
      # Styling - SINGLE COMBINED LEGEND
      theme_bw(base_family = "serif", base_size = 13) +
      labs(
        title = cfg$title,
        subtitle = paste0("D = ", D, ", n = ", n_val, ", β = ", BETA, 
                          " (", max(plot_data$rep, na.rm = TRUE), " repetitions)"),
        x = cfg$xlab,
        y = cfg$ylab
      ) +
      
      # Only use guide for color, suppress for shape and fill to get single legend
      scale_color_manual(
        values = model_colors, 
        name = "Model",
        guide = guide_legend(
          override.aes = list(
            size = 4,        # Larger points in legend
            alpha = 1,       # Full opacity in legend
            stroke = 1.5     # Thicker border in legend
          )
        )
      ) +
      scale_shape_manual(values = model_shapes, guide = "none") +
      scale_fill_manual(values = model_colors, guide = "none") +
      
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.7, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.3),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      )
    
    # Set axis limits if specified
    if (!is.null(cfg$xlim)) {
      p <- p + scale_x_continuous(limits = cfg$xlim, 
                                  breaks = seq(cfg$xlim[1], cfg$xlim[2], length.out = 5))
    }
    if (!is.null(cfg$ylim)) {
      p <- p + scale_y_continuous(limits = cfg$ylim, 
                                  breaks = seq(cfg$ylim[1], cfg$ylim[2], length.out = 5))
    }
    
    # Save plot
    filename_pdf <- file.path(n_dir, paste0("HC_Plane_", ent_name, "_D", D, "_n", n_val, ".pdf"))
    
      ggsave(
      filename = filename_pdf,
      plot = p,
      width = 12,
      height = 8,
      bg = "white"
    )
    
    cat(sprintf("✅ %s HC plane plot saved for n=%d\n", ent_name, n_val))
    
    # Print the plot
    print(p)
  }
}  
#----------------------------------------------------------------------------------------------
# This script has no CI in all scatter plots. It generates scatter plots of entropy vs. complexity 
#for various ARMA models, using extended entropy measures (Shannon, Rényi, Tsallis, Fisher). 
library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(here)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 3  # Embedding dimension
#BETA <- 1.5  # Beta parameter used

# Paths
hc_data_path <- here("Data", "Hybrid model_data", "D3_Data", "Combined_All_Models_HC_Results_D3.xlsx")
output_dir <- here("Plots", "Hybrid Model plots", "D3 plots")

# Sample sizes to process
sample_sizes <- c(1000, 5000)

# ==============================================================================
# DEFINE MODEL ORDER, COLORS, AND SHAPES
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

model_colors <- c(
  "ARMA(2,2)" = "#D55E00",   # Vermillion
  "ARMA(1,1)" = "#0072B2",   # Blue
  "AR(2)" = "#009E73",       # Bluish green
  "AR(1)" = "#CC79A7",       # Reddish purple
  "MA(2)" = "#E69F00",       # Orange
  "MA(1)" = "#56B4E9",       # Sky blue
  "Logistic" = "#000000",    # Black
  
  "Hybrid_ARMA(2,2)" = "#8B4513", # Dark brown
  "Hybrid_ARMA(1,1)" = "#4B0082", # Indigo
  "Hybrid_AR(2)" = "#696969",     # Dark gray
  "Hybrid_AR(1)" = "#20B2AA",     # Light sea green
  "Hybrid_MA(2)" = "#B22222",     # Firebrick
  "Hybrid_MA(1)" = "#4169E1",     # Royal blue
  
  "Logistic_r3_6" = "#2F4F4F",    # Slate gray
  "Sine_Wave" = "#228B22",        # Forest green
  "Logistic_Sine_Combined" = "#8A2BE2" # Blue violet
)


model_shapes <- c(
  "ARMA(2,2)" = 16,
  "ARMA(1,1)" = 17,
  "AR(2)" = 15,
  "AR(1)" = 18,
  "MA(2)" = 3,
  "MA(1)" = 7,
  "Logistic" = 8,
  
  "Hybrid_ARMA(2,2)" = 1,
  "Hybrid_ARMA(1,1)" = 2,
  "Hybrid_AR(2)" = 0,
  "Hybrid_AR(1)" = 5,
  "Hybrid_MA(2)" = 6,
  "Hybrid_MA(1)" = 4,
  
  "Logistic_r3_6" = 9,
  "Sine_Wave" = 10,
  "Logistic_Sine_Combined" = 11
)


# ==============================================================================
# LOAD HC DATA (ALL REPETITIONS)
# ==============================================================================
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Check available sample sizes
available_n <- unique(hc_data_all$n)
cat("Available sample sizes in data:", available_n, "\n")

# Verify we have the expected sample sizes
if (!all(sample_sizes %in% available_n)) {
  cat("WARNING: Not all expected sample sizes found in data!\n")
  cat("Expected:", sample_sizes, "\n")
  cat("Found:", available_n, "\n")
  sample_sizes <- intersect(sample_sizes, available_n)
  cat("Will process:", sample_sizes, "\n")
}

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon, Semi_HS, Semi_CS,
         H_Renyi, C_Renyi, Semi_HR, 
         H_Tsallis, C_Tsallis, Semi_HT,
         H_Fisher, C_Fisher, Semi_HF)

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

cat("\nTotal data points loaded:", nrow(hc_data), "\n")
cat("Sample sizes found:", unique(hc_data$n), "\n")
cat("Repetitions per sample size:\n")
print(table(hc_data$n))

# ==============================================================================
# LOAD LinfLsup BOUNDARIES (for Shannon only)
# ==============================================================================
data("LinfLsup")

# Full boundary curve
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# Separate upper and lower boundaries
boundary_lower <- LinfLsup_full %>% filter(Side == "Lower")
boundary_upper <- LinfLsup_full %>% filter(Side == "Upper")

cat("\nBoundary curves prepared for Shannon entropy\n")
cat("Boundary points - Lower:", nrow(boundary_lower), "Upper:", nrow(boundary_upper), "\n")

# ==============================================================================
# ENTROPY CONFIGURATIONS
# ==============================================================================
entropy_configs <- list(
  Shannon = list(
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_HS",
    Semi_C_col = "Semi_CS",
    add_bounds = TRUE,
    title = "Shannon Entropy-Complexity Plane",
    xlab = expression(italic(H)[Shannon]),
    ylab = expression(italic(C)[Shannon]),
    xlim = c(0, 1),
    ylim = c(0, 0.5)
  ),
  Renyi = list(
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_HR",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Rényi Entropy-Complexity Plane",
    xlab = expression(italic(H)[Rényi]),
    ylab = expression(italic(C)[Rényi]),
    xlim = NULL,
    ylim = NULL
  ),
  Tsallis = list(
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_HT",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Tsallis Entropy-Complexity Plane",
    xlab = expression(italic(H)[Tsallis]),
    ylab = expression(italic(C)[Tsallis]),
    xlim = NULL,
    ylim = NULL
  ),
  Fisher = list(
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_HF",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Fisher Entropy-Complexity Plane",
    xlab = expression(italic(H)[Fisher]),
    ylab = expression(italic(C)[Fisher]),
    xlim = NULL,
    ylim = NULL
  )
)

# ==============================================================================
# CREATE HC PLANE PLOTS FOR EACH SAMPLE SIZE AND ENTROPY TYPE
# ==============================================================================

for (n_val in sample_sizes) {
  
  # Create subdirectory for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n========== Processing Sample Size n = %d ==========\n", n_val))
  
  # Filter data for this sample size
  hc_data_n <- hc_data %>%
    filter(n == n_val)
  
  cat(sprintf("Data points for n=%d: %d\n", n_val, nrow(hc_data_n)))
  cat(sprintf("Models: %d, Max repetitions: %d\n", 
              length(unique(hc_data_n$Model)),
              max(hc_data_n$rep, na.rm = TRUE)))
  
  for (ent_name in names(entropy_configs)) {
    
    cfg <- entropy_configs[[ent_name]]
    
    cat(sprintf("\n=== Creating %s HC plane plot for n=%d ===\n", ent_name, n_val))
    
    # Extract relevant columns
    H_col <- cfg$H_col
    C_col <- cfg$C_col
    Semi_H_col <- cfg$Semi_H_col
    Semi_C_col <- cfg$Semi_C_col
    
    # Check if columns exist
    if (!all(c(H_col, C_col) %in% names(hc_data_n))) {
      cat(sprintf("⚠️ Skipping %s - columns not found\n", ent_name))
      next
    }
    
    # Prepare data with renamed columns for easier plotting
    plot_data <- hc_data_n %>%
      select(Model, rep, all_of(c(H_col, C_col))) %>%
      rename(H_val = !!H_col, C_val = !!C_col)
    
    # Remove rows with NA in H or C
    plot_data <- plot_data %>%
      filter(!is.na(H_val), !is.na(C_val), is.finite(H_val), is.finite(C_val))
    
    if (nrow(plot_data) == 0) {
      cat(sprintf("⚠️ Skipping %s - no valid data\n", ent_name))
      next
    }
    
    cat(sprintf("Valid data points: %d\n", nrow(plot_data)))
    
    # Create the plot
    p <- ggplot()
    
    # Add boundaries (Shannon only)
    if (cfg$add_bounds) {
      p <- p +
        geom_line(data = boundary_lower, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
        geom_line(data = boundary_upper, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8)
    }
    
    # Add points (all repetitions) - NO CONFIDENCE INTERVALS
    p <- p +
      geom_point(data = plot_data, 
                 aes(x = H_val, y = C_val, 
                     color = Model, 
                     shape = Model,
                     fill = Model),  # Add fill for shapes 21-25
                 size = 3, stroke = 1.2, alpha = 0.7) +
      
      # Styling - SINGLE COMBINED LEGEND
      theme_bw(base_family = "serif", base_size = 13) +
      labs(
        subtitle = paste0("D = ", D, ", n = ", n_val),
        x = cfg$xlab,
        y = cfg$ylab
      ) +
      
      # Only use guide for color, suppress for shape and fill to get single legend
      scale_color_manual(
        values = model_colors, 
        name = "Model",
        guide = guide_legend(
          override.aes = list(
            size = 4,        # Larger points in legend
            alpha = 1,       # Full opacity in legend
            stroke = 1.5     # Thicker border in legend
          )
        )
      ) +
      scale_shape_manual(values = model_shapes, guide = "none") +
      scale_fill_manual(values = model_colors, guide = "none") +
      
      theme(
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.7, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.3),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      )
    
    # Set axis limits if specified
    if (!is.null(cfg$xlim)) {
      p <- p + scale_x_continuous(limits = cfg$xlim, 
                                  breaks = seq(cfg$xlim[1], cfg$xlim[2], length.out = 5))
    }
    if (!is.null(cfg$ylim)) {
      p <- p + scale_y_continuous(limits = cfg$ylim, 
                                  breaks = seq(cfg$ylim[1], cfg$ylim[2], length.out = 5))
    }
    
    # Save plot
   filename_pdf <- file.path(n_dir, paste0("HC_Plane_", ent_name, "_D", D, "_n", n_val, ".pdf"))
    
    ggsave(
      filename = filename_pdf,
      plot = p,
      width = 12,
      height = 8,
      bg = "white"
    )
    
    cat(sprintf("✅ %s HC plane plot saved for n=%d\n", ent_name, n_val))
    
    # Print the plot
    print(p)
  }
}

#---------------------------------------------------------------------------------------------
# This script generates a feature analysis of the models. It includes correlation analysis and 
#prepares data for multinomial logistic regression. It uses the same features as the previous 
#script (H and C for Shannon, Rényi, Tsallis, Fisher, and Disequilibrium).
#----------------------------------------------------------------------------------------------
library(nnet)
library(caret)
library(corrplot)
library(readxl)
library(dplyr)
library(tidyr)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension

# Paths
hc_data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Hybrid ARMA+LogisticMap_data/D4_Data/Combined_All_Models_HC_Results_D4.xlsx"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Hybrid_Analysis/Multinomial_Logistic_Regression"

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sample sizes to process
sample_sizes <- c(1000, 5000)

# ==============================================================================
# DEFINE MODEL NAMES
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data...\n")
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon,
         H_Renyi, C_Renyi,
         H_Tsallis, C_Tsallis,
         H_Fisher, C_Fisher,
         all_of(disequil_col))

# Rename for consistency
names(hc_data)[names(hc_data) == disequil_col] <- "Disequilibrium"

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)
# ==============================================================================
# DEFINE FEATURES
# ==============================================================================
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "Disequilibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

print(features)

# ==============================================================================
# ANALYSIS FOR EACH SAMPLE SIZE
# ==============================================================================

for (n_val in sample_sizes) {
  
  cat(sprintf("\n\n================================================================================\n"))
  cat(sprintf("ANALYSIS FOR SAMPLE SIZE n = %d\n", n_val))
  cat(sprintf("================================================================================\n\n"))
  
  # Create subdirectory for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # PREPARE DATA
  # --------------------------------------------------------------------------
  
  analysis_data <- hc_data %>%
    filter(n == n_val) %>%
    select(Model, all_of(features)) %>%
    na.omit()  # Remove any rows with missing values
  
  cat(sprintf("Data after filtering and removing NA: %d rows\n", nrow(analysis_data)))
  cat(sprintf("Models in dataset: %d\n", length(unique(analysis_data$Model))))
  
  # Check class distribution
  cat("\nClass distribution:\n")
  print(table(analysis_data$Model))
  
  # --------------------------------------------------------------------------
  # CORRELATION ANALYSIS
  # --------------------------------------------------------------------------
  
  cat("\n### CORRELATION ANALYSIS ###\n")
  
  # Calculate correlation matrix
  cor_matrix <- cor(analysis_data[, features], use = "complete.obs")
  
  cat("\nCorrelation Matrix:\n")
  print(round(cor_matrix, 3))
  
  # Save correlation matrix as CSV
  write.csv(round(cor_matrix, 3),
            file = file.path(n_dir, paste0("Correlation_Matrix_n", n_val, ".csv")),
            row.names = TRUE)
  
  # Save correlation plot as PDF
  pdf(file.path(n_dir, paste0("Correlation_Matrix_n", n_val, ".pdf")),
      width = 10, height = 10)
  corrplot(cor_matrix, 
           method = "color", 
           addCoef.col = "black",
           tl.col = "black",
           tl.srt = 45,
           number.cex = 0.7,
           title = paste0("Feature Correlation Matrix (n = ", n_val, ")"),
           mar = c(0, 0, 2, 0))
  dev.off()
  
  cat(sprintf("\n✅ Correlation matrix saved as PDF\n"))
  
  # Identify highly correlated pairs (|r| > 0.8)
  high_cor <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("\nHighly correlated feature pairs (|r| > 0.8):\n")
    for (i in 1:nrow(high_cor)) {
      row_idx <- high_cor[i, 1]
      col_idx <- high_cor[i, 2]
      if (row_idx < col_idx) {  # Avoid duplicates
        cat(sprintf("  %s <-> %s: %.3f\n", 
                    rownames(cor_matrix)[row_idx],
                    colnames(cor_matrix)[col_idx],
                    cor_matrix[row_idx, col_idx]))
      }
    }
  } else {
    cat("\nNo highly correlated feature pairs (|r| > 0.8) found.\n")
  }
  
  # --------------------------------------------------------------------------
  # FIT MULTINOMIAL LOGISTIC REGRESSION MODEL
  # --------------------------------------------------------------------------
  
  cat("\n### MULTINOMIAL LOGISTIC REGRESSION MODEL ###\n")
  
  # Fit the model
  cat("\nFitting multinomial logistic regression model...\n")
  
  set.seed(42)  # For reproducibility
  
  multinom_model <- multinom(
    Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
      C_Shannon + C_Renyi + C_Tsallis + C_Fisher + Disequilibrium, 
    data = analysis_data,
    maxit = 500,  # Maximum iterations
    trace = FALSE  # Suppress iteration output
  )
  
  cat("\n✅ Model fitted successfully!\n")
  
  # --------------------------------------------------------------------------
  # MODEL SUMMARY AND STATISTICS
  # --------------------------------------------------------------------------
  
  # Capture model summary
  sink(file.path(n_dir, paste0("Model_Summary_n", n_val, ".txt")))
  cat("================================================================================\n")
  cat(sprintf("MULTINOMIAL LOGISTIC REGRESSION MODEL SUMMARY (n = %d)\n", n_val))
  cat("================================================================================\n\n")
  print(summary(multinom_model))
  cat("\n\n")
  sink()
  
  # Get AIC
  model_aic <- AIC(multinom_model)
  cat(sprintf("\nModel AIC: %.2f\n", model_aic))
  
  # --------------------------------------------------------------------------
  # MODEL PREDICTIONS AND ACCURACY
  # --------------------------------------------------------------------------
  
  # Make predictions on the same data (training accuracy)
  predictions <- predict(multinom_model, newdata = analysis_data, type = "class")
  
  # Create confusion matrix
  conf_matrix <- confusionMatrix(predictions, analysis_data$Model)
  
  # Extract overall accuracy
  overall_accuracy <- conf_matrix$overall['Accuracy']
  cat(sprintf("Overall Model Accuracy: %.4f (%.2f%%)\n", 
              overall_accuracy, overall_accuracy * 100))
  
  # Per-class accuracy (sensitivity)
  per_class_accuracy <- conf_matrix$byClass[, 'Sensitivity']
  
  cat("\nPer-Class Accuracy (Sensitivity):\n")
  print(round(per_class_accuracy, 4))
  
  # --------------------------------------------------------------------------
  # SAVE CONFUSION MATRIX
  # --------------------------------------------------------------------------
  
  # Save confusion matrix
  sink(file.path(n_dir, paste0("Confusion_Matrix_n", n_val, ".txt")))
  cat("================================================================================\n")
  cat(sprintf("CONFUSION MATRIX (n = %d)\n", n_val))
  cat("================================================================================\n\n")
  print(conf_matrix)
  sink()
  
  # Save confusion matrix table as CSV
  conf_table <- as.data.frame.matrix(table(Predicted = predictions, Actual = analysis_data$Model))
  write.csv(conf_table,
            file = file.path(n_dir, paste0("Confusion_Matrix_Table_n", n_val, ".csv")),
            row.names = TRUE)
  
  # --------------------------------------------------------------------------
  # SAVE MODEL COEFFICIENTS
  # --------------------------------------------------------------------------
  
  # Extract coefficients
  coefficients <- summary(multinom_model)$coefficients
  std_errors <- summary(multinom_model)$standard.errors
  
  # Calculate z-values and p-values
  z_values <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # Save coefficients
  write.csv(coefficients,
            file = file.path(n_dir, paste0("Model_Coefficients_n", n_val, ".csv")),
            row.names = TRUE)
  
  write.csv(std_errors,
            file = file.path(n_dir, paste0("Model_StdErrors_n", n_val, ".csv")),
            row.names = TRUE)
  
  write.csv(p_values,
            file = file.path(n_dir, paste0("Model_PValues_n", n_val, ".csv")),
            row.names = TRUE)
  
  # --------------------------------------------------------------------------
  # SAVE PERFORMANCE METRICS
  # --------------------------------------------------------------------------
  
  # Create performance summary
  performance_summary <- data.frame(
    Metric = c("Sample Size", "Number of Observations", "Number of Classes", 
               "AIC", "Overall Accuracy", "Kappa"),
    Value = c(n_val, 
              nrow(analysis_data),
              length(unique(analysis_data$Model)),
              round(model_aic, 2),
              round(overall_accuracy, 4),
              round(conf_matrix$overall['Kappa'], 4))
  )
  
  write.csv(performance_summary,
            file = file.path(n_dir, paste0("Performance_Summary_n", n_val, ".csv")),
            row.names = FALSE)
  
  cat("\n### PERFORMANCE SUMMARY ###\n")
  print(performance_summary)
  
  # --------------------------------------------------------------------------
  # SAVE PER-CLASS METRICS
  # --------------------------------------------------------------------------
  
  # Extract per-class metrics
  per_class_metrics <- data.frame(
    Model = rownames(conf_matrix$byClass),
    Sensitivity = conf_matrix$byClass[, 'Sensitivity'],
    Specificity = conf_matrix$byClass[, 'Specificity'],
    Pos_Pred_Value = conf_matrix$byClass[, 'Pos Pred Value'],
    Neg_Pred_Value = conf_matrix$byClass[, 'Neg Pred Value'],
    Precision = conf_matrix$byClass[, 'Precision'],
    Recall = conf_matrix$byClass[, 'Recall'],
    F1 = conf_matrix$byClass[, 'F1'],
    Balanced_Accuracy = conf_matrix$byClass[, 'Balanced Accuracy']
  )
  
  write.csv(per_class_metrics,
            file = file.path(n_dir, paste0("Per_Class_Metrics_n", n_val, ".csv")),
            row.names = FALSE)
  
  cat("\n✅ All results saved for n =", n_val, "\n")
  
  # --------------------------------------------------------------------------
  # SAVE R MODEL OBJECT
  # --------------------------------------------------------------------------
  
  saveRDS(multinom_model,
          file = file.path(n_dir, paste0("Multinom_Model_n", n_val, ".rds")))
  
  cat(sprintf("\n✅ Model object saved as RDS file\n"))
  
}

#-----------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
# This script use time series clustering to analyze the similarity between the models 
# based on their entropy and complexity
#-------------------------------------------------------------------------------------
library(cluster)
library(factoextra)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(gridExtra)
library(reshape2)
library(RColorBrewer)
library(viridis)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension

# Paths
hc_data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Hybrid ARMA+LogisticMap_data/D4_Data/Combined_All_Models_HC_Results_D4.xlsx"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Hybrid_Analysis/Multinomial_Logistic_Regression"

# Create output directory for clustering
clustering_dir <- file.path(output_dir, "Time_Series_Clustering")
dir.create(clustering_dir, recursive = TRUE, showWarnings = FALSE)

# Sample sizes to process
sample_sizes <- c(1000, 5000)

# Maximum number of clusters to test
max_k <- 10

# ==============================================================================
# DEFINE MODEL NAMES
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data...\n")
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Check for Disequilibrium column
if ("Disequilibrium" %in% names(hc_data_all)) {
  disequil_col <- "Disequilibrium"
} else if ("Disequlibrium" %in% names(hc_data_all)) {
  disequil_col <- "Disequlibrium"
} else {
  stop("Disequilibrium column not found in data!")
}

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon,
         H_Renyi, C_Renyi,
         H_Tsallis, C_Tsallis,
         H_Fisher, C_Fisher,
         all_of(disequil_col))

# Rename for consistency
names(hc_data)[names(hc_data) == disequil_col] <- "Disequilibrium"

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

cat("\nTotal data points loaded:", nrow(hc_data), "\n")

# ==============================================================================
# DEFINE CLUSTERING FEATURES
# ==============================================================================
clustering_features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
                         "C_Shannon", "Disequilibrium", "C_Renyi", 
                         "C_Tsallis", "C_Fisher")

cat("\nFeatures for clustering:\n")
print(clustering_features)

# ==============================================================================
# CLUSCO ALGORITHM IMPLEMENTATION
# ==============================================================================

# Function to compute cluster center
compute_center <- function(data) {
  colMeans(data, na.rm = TRUE)
}

# Function to compute distance from point to center
compute_distance <- function(point, center) {
  sqrt(sum((point - center)^2))
}

# Function to compute cluster function f_l (based on silhouette-inspired metric)
compute_cluster_function <- function(data, centers) {
  n <- nrow(data)
  k <- nrow(centers)
  
  # Assign each point to nearest center
  distances <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    distances[, i] <- apply(data, 1, function(x) compute_distance(x, centers[i, ]))
  }
  
  assignments <- apply(distances, 1, which.min)
  
  # Compute average within-cluster distance (compactness)
  within_cluster_dist <- 0
  for (j in 1:k) {
    cluster_points <- data[assignments == j, , drop = FALSE]
    if (nrow(cluster_points) > 0) {
      within_cluster_dist <- within_cluster_dist + 
        mean(apply(cluster_points, 1, function(x) compute_distance(x, centers[j, ])))
    }
  }
  
  # Return negative (we want to minimize this)
  return(-within_cluster_dist / k)
}

# Discrete gradient method for center optimization
optimize_centers <- function(data, initial_centers, max_iter = 100, tol = 1e-6) {
  centers <- initial_centers
  k <- nrow(centers)
  
  for (iter in 1:max_iter) {
    old_centers <- centers
    
    # Assign points to nearest center
    n <- nrow(data)
    distances <- matrix(0, nrow = n, ncol = k)
    for (i in 1:k) {
      distances[, i] <- apply(data, 1, function(x) compute_distance(x, centers[i, ]))
    }
    assignments <- apply(distances, 1, which.min)
    
    # Update centers as mean of assigned points
    for (j in 1:k) {
      cluster_points <- data[assignments == j, , drop = FALSE]
      if (nrow(cluster_points) > 0) {
        centers[j, ] <- colMeans(cluster_points)
      }
    }
    
    # Check convergence
    if (max(abs(centers - old_centers)) < tol) {
      break
    }
  }
  
  return(centers)
}

# CLUSCO Algorithm Implementation
clusco_algorithm <- function(data, max_k) {
  n <- nrow(data)
  p <- ncol(data)
  
  results <- list()
  
  # Step 1: Initialize with k=1 (single cluster)
  cat("\n  Computing initial center (k=1)...\n")
  x1 <- compute_center(data)
  centers <- matrix(x1, nrow = 1)
  f_values <- c(compute_cluster_function(data, centers))
  
  results[[1]] <- list(
    k = 1,
    centers = centers,
    f_value = f_values[1]
  )
  
  # Step 2-3: Iteratively add clusters
  for (l in 2:max_k) {
    cat(sprintf("  Computing solution for k=%d...\n", l))
    
    # Find point that maximizes the improvement
    current_centers <- results[[l-1]]$centers
    f_prev <- results[[l-1]]$f_value
    
    best_point <- NULL
    best_improvement <- -Inf
    
    # Sample points to test (for efficiency, test every point)
    for (i in 1:min(n, 500)) {  # Limit to 500 random points for large datasets
      idx <- sample(1:n, 1)
      test_point <- as.numeric(data[idx, ])
      
      # Create temporary centers with this new point
      temp_centers <- rbind(current_centers, test_point)
      
      # Compute improvement
      f_new <- compute_cluster_function(data, temp_centers)
      improvement <- f_prev - f_new
      
      if (improvement > best_improvement) {
        best_improvement <- improvement
        best_point <- test_point
      }
    }
    
    # Apply discrete gradient method
    initial_centers <- rbind(current_centers, best_point)
    optimized_centers <- optimize_centers(data, initial_centers)
    
    f_l <- compute_cluster_function(data, optimized_centers)
    
    results[[l]] <- list(
      k = l,
      centers = optimized_centers,
      f_value = f_l
    )
  }
  
  return(results)
}

# ==============================================================================
# ANALYSIS FOR EACH SAMPLE SIZE
# ==============================================================================

for (n_val in sample_sizes) {
  
  cat(sprintf("\n\n================================================================================\n"))
  cat(sprintf("TIME SERIES CLUSTERING FOR SAMPLE SIZE n = %d\n", n_val))
  cat(sprintf("================================================================================\n\n"))
  
  # Create subdirectory for this sample size
  n_dir <- file.path(clustering_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # PREPARE DATA
  # --------------------------------------------------------------------------
  
  clustering_data <- hc_data %>%
    filter(n == n_val) %>%
    select(Model, all_of(clustering_features)) %>%
    na.omit()
  
  cat(sprintf("Data points: %d\n", nrow(clustering_data)))
  cat(sprintf("Models: %d\n", length(unique(clustering_data$Model))))
  
  # Extract feature matrix (without Model column)
  feature_matrix <- as.matrix(clustering_data[, clustering_features])
  
  # Standardize features (important for distance-based clustering)
  feature_matrix_scaled <- scale(feature_matrix)
  
  # --------------------------------------------------------------------------
  # APPLY CLUSCO ALGORITHM
  # --------------------------------------------------------------------------
  
  cat("\n### APPLYING CLUSCO ALGORITHM ###\n")
  
  clusco_results <- clusco_algorithm(feature_matrix_scaled, max_k)
  
  cat("\n✅ CLUSCO algorithm completed!\n")
  
  # --------------------------------------------------------------------------
  # COMPUTE SILHOUETTE COEFFICIENTS FOR EACH k
  # --------------------------------------------------------------------------
  
  cat("\n### COMPUTING SILHOUETTE COEFFICIENTS ###\n")
  
  silhouette_results <- data.frame(
    k = integer(),
    avg_silhouette = numeric(),
    f_value = numeric()
  )
  
  cluster_assignments_all <- list()
  silhouette_data_all <- list()
  
  for (l in 2:max_k) {  # Silhouette needs at least 2 clusters
    centers <- clusco_results[[l]]$centers
    
    # Assign points to nearest center
    distances <- matrix(0, nrow = nrow(feature_matrix_scaled), ncol = l)
    for (i in 1:l) {
      distances[, i] <- apply(feature_matrix_scaled, 1, 
                              function(x) compute_distance(x, centers[i, ]))
    }
    assignments <- apply(distances, 1, which.min)
    
    # Compute silhouette
    if (length(unique(assignments)) > 1) {
      sil <- silhouette(assignments, dist(feature_matrix_scaled))
      avg_sil <- mean(sil[, 3])
      
      silhouette_results <- rbind(silhouette_results, data.frame(
        k = l,
        avg_silhouette = avg_sil,
        f_value = clusco_results[[l]]$f_value
      ))
      
      cluster_assignments_all[[as.character(l)]] <- assignments
      silhouette_data_all[[as.character(l)]] <- sil
      
      cat(sprintf("  k=%d: Avg Silhouette = %.4f, f_value = %.4f\n", 
                  l, avg_sil, clusco_results[[l]]$f_value))
    }
  }
  
  # Save silhouette results
  write.csv(silhouette_results,
            file = file.path(n_dir, paste0("Silhouette_Results_n", n_val, ".csv")),
            row.names = FALSE)
  
  # --------------------------------------------------------------------------
  # FIND OPTIMAL NUMBER OF CLUSTERS
  # --------------------------------------------------------------------------
  
  optimal_k <- silhouette_results$k[which.max(silhouette_results$avg_silhouette)]
  cat(sprintf("\n✅ Optimal number of clusters: k = %d (max avg silhouette = %.4f)\n", 
              optimal_k, max(silhouette_results$avg_silhouette)))
  
  # --------------------------------------------------------------------------
  # FIGURE 1: SILHOUETTE COEFFICIENT VS NUMBER OF CLUSTERS
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 1: Silhouette Coefficient vs k ###\n")
  
  fig1 <- ggplot(silhouette_results, aes(x = k, y = avg_silhouette)) +
    geom_line(color = "#2E86AB", size = 1.2, linetype = "solid") +
    geom_point(color = "#2E86AB", size = 4, shape = 19) +
    geom_point(data = silhouette_results[silhouette_results$k == optimal_k, ],
               aes(x = k, y = avg_silhouette), 
               color = "#A23B72", size = 5, shape = 18) +
    geom_vline(xintercept = optimal_k, linetype = "dashed", 
               color = "#A23B72", size = 0.8, alpha = 0.7) +
    annotate("text", x = optimal_k + 0.3, 
             y = max(silhouette_results$avg_silhouette),
             label = paste0("Optimal k = ", optimal_k), 
             color = "#A23B72", hjust = 0, size = 5, fontface = "bold") +
    scale_x_continuous(breaks = 2:max_k) +
    labs(title = paste0("Average Silhouette Coefficient (n = ", n_val, ")"),
         x = "Number of Clusters (k)",
         y = "Average Silhouette Width") +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1)
    )
  
  print(fig1)
  
  ggsave(file.path(n_dir, paste0("Figure_1_Silhouette_vs_k_n", n_val, ".pdf")),
         plot = fig1, width = 10, height = 6)
  
  cat("✅ Figure 1 saved\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 2: CLUSTER FUNCTION VALUE VS k
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 2: Cluster Function Value vs k ###\n")
  
  f_values_df <- data.frame(
    k = 1:max_k,
    f_value = sapply(clusco_results, function(x) x$f_value)
  )
  
  fig2 <- ggplot(f_values_df, aes(x = k, y = f_value)) +
    geom_line(color = "#06A77D", size = 1.2) +
    geom_point(color = "#06A77D", size = 4, shape = 17) +
    geom_point(data = f_values_df[f_values_df$k == optimal_k, ],
               aes(x = k, y = f_value),
               color = "#D62246", size = 5, shape = 18) +
    geom_vline(xintercept = optimal_k, linetype = "dashed", 
               color = "#D62246", size = 0.8, alpha = 0.7) +
    scale_x_continuous(breaks = 1:max_k) +
    labs(title = paste0("CLUSCO Objective Function (n = ", n_val, ")"),
         x = "Number of Clusters (k)",
         y = "Cluster Function Value f(k)") +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1)
    )
  
  print(fig2)
  
  ggsave(file.path(n_dir, paste0("Figure_2_Cluster_Function_vs_k_n", n_val, ".pdf")),
         plot = fig2, width = 10, height = 6)
  
  cat("✅ Figure 2 saved\n")
  
  # --------------------------------------------------------------------------
  # ANALYZE OPTIMAL CLUSTERING
  # --------------------------------------------------------------------------
  
  cat(sprintf("\n### ANALYZING OPTIMAL CLUSTERING (k = %d) ###\n", optimal_k))
  
  optimal_assignments <- cluster_assignments_all[[as.character(optimal_k)]]
  
  # Add cluster assignments to data
  clustering_data$Cluster <- optimal_assignments
  
  # Cluster size summary
  cluster_sizes <- table(optimal_assignments)
  cat("\nCluster sizes:\n")
  print(cluster_sizes)
  
  # Model distribution across clusters
  model_cluster_table <- table(clustering_data$Model, clustering_data$Cluster)
  cat("\nModel distribution across clusters:\n")
  print(model_cluster_table)
  
  # Save cluster assignments
  cluster_assignment_df <- data.frame(
    Model = clustering_data$Model,
    Cluster = optimal_assignments
  )
  
  write.csv(cluster_assignment_df,
            file = file.path(n_dir, paste0("Cluster_Assignments_k", optimal_k, "_n", n_val, ".csv")),
            row.names = FALSE)
  
  # Save model-cluster cross-tabulation
  write.csv(model_cluster_table,
            file = file.path(n_dir, paste0("Model_Cluster_CrossTab_k", optimal_k, "_n", n_val, ".csv")),
            row.names = TRUE)
  
  # --------------------------------------------------------------------------
  # FIGURE 3: SILHOUETTE PLOT FOR OPTIMAL k
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 3: Silhouette Plot for Optimal k ###\n")
  
  sil_optimal <- silhouette_data_all[[as.character(optimal_k)]]
  
  pdf(file.path(n_dir, paste0("Figure_3_Silhouette_Plot_k", optimal_k, "_n", n_val, ".pdf")),
      width = 10, height = 7)
  
  plot(sil_optimal, 
       main = paste0("Silhouette Plot (k = ", optimal_k, ", n = ", n_val, ")"),
       col = brewer.pal(min(optimal_k, 9), "Set1"),
       border = NA,
       cex.names = 0.8)
  
  dev.off()
  
  cat("✅ Figure 3 saved (view in PDF viewer)\n")
  
  # --------------------------------------------------------------------------
  # CLUSTER CENTERS ANALYSIS
  # --------------------------------------------------------------------------
  
  optimal_centers <- clusco_results[[optimal_k]]$centers
  rownames(optimal_centers) <- paste0("Cluster_", 1:optimal_k)
  colnames(optimal_centers) <- clustering_features
  
  # Inverse transform to original scale
  optimal_centers_original <- optimal_centers
  for (i in 1:ncol(optimal_centers)) {
    optimal_centers_original[, i] <- optimal_centers[, i] * 
      attr(feature_matrix_scaled, "scaled:scale")[i] + 
      attr(feature_matrix_scaled, "scaled:center")[i]
  }
  
  # Save cluster centers
  write.csv(optimal_centers_original,
            file = file.path(n_dir, paste0("Cluster_Centers_k", optimal_k, "_n", n_val, ".csv")),
            row.names = TRUE)
  
  cat("\nCluster centers (original scale):\n")
  print(round(optimal_centers_original, 4))
  
  # --------------------------------------------------------------------------
  # FIGURE 4: PCA BIPLOT WITH CLUSTERS
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 4: PCA Biplot with Clusters ###\n")
  
  # Perform PCA for visualization
  pca_result <- prcomp(feature_matrix_scaled, center = FALSE, scale. = FALSE)
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  pca_data$Cluster <- factor(optimal_assignments)
  pca_data$Model <- clustering_data$Model
  
  # Project cluster centers onto PC space
  centers_pca <- as.data.frame(predict(pca_result, optimal_centers)[, 1:2])
  centers_pca$Cluster <- factor(1:optimal_k)
  
  # Variance explained
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # Get loadings for biplot arrows
  loadings <- pca_result$rotation[, 1:2]
  loadings_df <- as.data.frame(loadings * 4)  # Scale for visibility
  loadings_df$feature <- rownames(loadings)
  
  # Plot
  fig4 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_point(data = centers_pca, aes(x = PC1, y = PC2, color = Cluster),
               size = 10, shape = 8, stroke = 2.5, show.legend = FALSE) +
    # Add arrows for loadings
    geom_segment(data = loadings_df,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 color = "gray30", size = 0.7, alpha = 0.5,
                 inherit.aes = FALSE) +
    geom_text(data = loadings_df,
              aes(x = PC1 * 1.1, y = PC2 * 1.1, label = feature),
              color = "gray20", size = 3, fontface = "italic",
              inherit.aes = FALSE) +
    scale_color_brewer(palette = "Set1") +
    labs(title = paste0("PCA Biplot (k = ", optimal_k, ", n = ", n_val, ")"),
         subtitle = "Cluster centers marked with asterisks (*)",
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.title = element_text(face = "bold", size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  print(fig4)
  
  ggsave(file.path(n_dir, paste0("Figure_4_PCA_Biplot_k", optimal_k, "_n", n_val, ".pdf")),
         plot = fig4, width = 12, height = 8)
  
  cat("✅ Figure 4 saved\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 5: CLUSTER PROFILES HEATMAP
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 5: Cluster Profiles Heatmap ###\n")
  
  # Compute mean feature values per cluster
  cluster_profiles <- clustering_data %>%
    group_by(Cluster) %>%
    summarise(across(all_of(clustering_features), mean, .names = "{.col}")) %>%
    as.data.frame()
  
  write.csv(cluster_profiles,
            file = file.path(n_dir, paste0("Cluster_Profiles_k", optimal_k, "_n", n_val, ".csv")),
            row.names = FALSE)
  
  # Prepare for heatmap
  cluster_profiles_mat <- as.matrix(cluster_profiles[, -1])
  rownames(cluster_profiles_mat) <- paste0("Cluster ", cluster_profiles$Cluster)
  
  # Standardize for heatmap
  cluster_profiles_scaled <- scale(t(cluster_profiles_mat))
  
  pdf(file.path(n_dir, paste0("Figure_5_Cluster_Profiles_Heatmap_k", optimal_k, "_n", n_val, ".pdf")),
      width = 10, height = 8)
  
  par(mar = c(5, 10, 4, 2))
  heatmap(cluster_profiles_scaled,
          main = paste0("Cluster Profiles (k = ", optimal_k, ", n = ", n_val, ")"),
          xlab = "Clusters",
          ylab = "",
          scale = "none",
          col = colorRampPalette(c("#313695", "#4575B4", "#ABD9E9", 
                                   "#FFFFBF", "#FDAE61", "#F46D43", "#A50026"))(100),
          margins = c(6, 12),
          cexRow = 1.0,
          cexCol = 1.2)
  
  dev.off()
  
  cat("✅ Figure 5 saved (view in PDF viewer)\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 6: CLUSTER SIZES BAR PLOT
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 6: Cluster Sizes ###\n")
  
  cluster_size_df <- data.frame(
    Cluster = factor(names(cluster_sizes)),
    Size = as.numeric(cluster_sizes)
  )
  
  fig6 <- ggplot(cluster_size_df, aes(x = Cluster, y = Size, fill = Cluster)) +
    geom_bar(stat = "identity", color = "black", size = 0.5) +
    geom_text(aes(label = Size), vjust = -0.5, size = 5, fontface = "bold") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = paste0("Cluster Sizes (k = ", optimal_k, ", n = ", n_val, ")"),
         x = "Cluster",
         y = "Number of Time Series") +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(fig6)
  
  ggsave(file.path(n_dir, paste0("Figure_6_Cluster_Sizes_k", optimal_k, "_n", n_val, ".pdf")),
         plot = fig6, width = 10, height = 6)
  
  cat("✅ Figure 6 saved\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 7: MODEL DISTRIBUTION ACROSS CLUSTERS
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 7: Model Distribution Across Clusters ###\n")
  
  # Prepare data for stacked bar chart
  model_cluster_df <- as.data.frame(model_cluster_table)
  names(model_cluster_df) <- c("Model", "Cluster", "Count")
  model_cluster_df$Cluster <- factor(model_cluster_df$Cluster)
  
  
    fig7 <- ggplot(model_cluster_df, aes(x = Cluster, y = Count, fill = Model)) +
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.3) +
    scale_fill_manual(values = c(
      "ARMA(2,2)" = "#E41A1C", "ARMA(1,1)" = "blue", "AR(2)" = "#4DAF4A",
      "AR(1)" = "#984EA3", "MA(2)" = "#FF7F00", "MA(1)" = "#FFFF33",
      "Logistic" = "#000000",
      "Hybrid_ARMA(2,2)" = "#F781BF", "Hybrid_ARMA(1,1)" = "#A65628",
      "Hybrid_AR(2)" = "#999999", "Hybrid_AR(1)" = "#66C2A5",
      "Hybrid_MA(2)" = "#FC8D62", "Hybrid_MA(1)" = "#8DA0CB"
    )) +
    labs(title = paste0("Model Distribution Across Clusters (k = ", optimal_k, ", n = ", n_val, ")"),
         x = "Cluster",
         y = "Count",
         fill = "Model") +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major.x = element_blank()
    )
  
  print(fig7)
  
  ggsave(file.path(n_dir, paste0("Figure_7_Model_Distribution_k", optimal_k, "_n", n_val, ".pdf")),
         plot = fig7, width = 12, height = 7)
  
  cat("✅ Figure 7 saved\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 8: WITHIN-CLUSTER SUM OF SQUARES (ELBOW METHOD)
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 8: Within-Cluster Sum of Squares ###\n")
  
  # Calculate WSS for each k
  wss_values <- numeric(max_k)
  
  for (k in 1:max_k) {
    centers <- clusco_results[[k]]$centers
    
    # Assign points
    if (k == 1) {
      distances <- apply(feature_matrix_scaled, 1, 
                         function(x) compute_distance(x, centers))
      wss_values[k] <- sum(distances^2)
    } else {
      distances <- matrix(0, nrow = nrow(feature_matrix_scaled), ncol = k)
      for (i in 1:k) {
        distances[, i] <- apply(feature_matrix_scaled, 1, 
                                function(x) compute_distance(x, centers[i, ]))
      }
      assignments <- apply(distances, 1, which.min)
      
      wss <- 0
      for (j in 1:k) {
        cluster_points <- feature_matrix_scaled[assignments == j, , drop = FALSE]
        if (nrow(cluster_points) > 0) {
          wss <- wss + sum(apply(cluster_points, 1, 
                                 function(x) compute_distance(x, centers[j, ])^2))
        }
      }
      wss_values[k] <- wss
    }
  }
  
  wss_df <- data.frame(k = 1:max_k, WSS = wss_values)
  
  fig8 <- ggplot(wss_df, aes(x = k, y = WSS)) +
    geom_line(color = "#E07A5F", size = 1.2) +
    geom_point(color = "#E07A5F", size = 4, shape = 19) +
    geom_point(data = wss_df[wss_df$k == optimal_k, ],
               aes(x = k, y = WSS),
               color = "#3D405B", size = 5, shape = 18) +
    geom_vline(xintercept = optimal_k, linetype = "dashed", 
               color = "#3D405B", size = 0.8, alpha = 0.7) +
    scale_x_continuous(breaks = 1:max_k) +
    labs(title = paste0("Within-Cluster Sum of Squares (n = ", n_val, ")"),
         x = "Number of Clusters (k)",
         y = "Within-Cluster SS") +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1)
    )
  
  print(fig8)
  
  ggsave(file.path(n_dir, paste0("Figure_8_WSS_Elbow_n", n_val, ".pdf")),
         plot = fig8, width = 10, height = 6)
  
  cat("✅ Figure 8 saved\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 9: PAIRWISE FEATURE SCATTER PLOTS (SELECTED FEATURES)
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 9: Pairwise Feature Scatter Plots ###\n")
  
  # Select key features for visualization
  key_features <- c("H_Shannon", "C_Shannon", "H_Fisher", "C_Fisher")
  
  plot_data_pairs <- clustering_data %>%
    select(all_of(key_features), Cluster, Model)
  
  plot_data_pairs$Cluster <- factor(plot_data_pairs$Cluster)
  
  # Create pairs plot
  pdf(file.path(n_dir, paste0("Figure_9_Pairwise_Scatter_k", optimal_k, "_n", n_val, ".pdf")),
      width = 12, height = 12)
  
  pairs(plot_data_pairs[, key_features],
        col = brewer.pal(min(optimal_k, 9), "Set1")[plot_data_pairs$Cluster],
        pch = 19,
        cex = 1.2,
        main = paste0("Pairwise Feature Scatter Plots (k = ", optimal_k, ", n = ", n_val, ")"),
        upper.panel = NULL)
  
  dev.off()
  
  cat("✅ Figure 9 saved (view in PDF viewer)\n")
  
  # --------------------------------------------------------------------------
  # FIGURE 10: COMPARISON OF CLUSTERING QUALITY METRICS
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING FIGURE 10: Clustering Quality Metrics Comparison ###\n")
  
  # Normalize metrics to [0, 1] for comparison
  metrics_df <- data.frame(
    k = silhouette_results$k,
    Silhouette = (silhouette_results$avg_silhouette - min(silhouette_results$avg_silhouette)) /
      (max(silhouette_results$avg_silhouette) - min(silhouette_results$avg_silhouette)),
    WSS_normalized = 1 - (wss_values[2:max_k] - min(wss_values[2:max_k])) /
      (max(wss_values[2:max_k]) - min(wss_values[2:max_k])),
    F_value_normalized = (silhouette_results$f_value - min(silhouette_results$f_value)) /
      (max(silhouette_results$f_value) - min(silhouette_results$f_value))
  )
  
  # Reshape for plotting
  metrics_long <- metrics_df %>%
    pivot_longer(cols = c(Silhouette, WSS_normalized, F_value_normalized),
                 names_to = "Metric", values_to = "Value")
  
  metrics_long$Metric <- factor(metrics_long$Metric,
                                levels = c("Silhouette", "WSS_normalized", "F_value_normalized"),
                                labels = c("Avg Silhouette", "1 - Normalized WSS", "Normalized f(k)"))
  
  fig10 <- ggplot(metrics_long, aes(x = k, y = Value, color = Metric, group = Metric)) +
    geom_line(size = 1.2) +
    geom_point(size = 3.5, shape = 19) +
    geom_vline(xintercept = optimal_k, linetype = "dashed", 
               color = "black", size = 0.8, alpha = 0.5) +
    scale_color_manual(values = c("#2E86AB", "#A23B72", "#F18F01")) +
    scale_x_continuous(breaks = 2:max_k) +
    labs(title = paste0("Clustering Quality Metrics Comparison (n = ", n_val, ")"),
         x = "Number of Clusters (k)",
         y = "Normalized Metric Value",
         color = "Metric") +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank()
    )
  
  print(fig10)
  
  ggsave(file.path(n_dir, paste0("Figure_10_Quality_Metrics_Comparison_n", n_val, ".pdf")),
         plot = fig10, width = 12, height = 7)
  
  cat("✅ Figure 10 saved\n")
  
  # --------------------------------------------------------------------------
  # SAVE ALL CLUSCO RESULTS
  # --------------------------------------------------------------------------
  
  # Save all cluster solutions
  all_solutions <- data.frame(
    k = 1:max_k,
    f_value = sapply(clusco_results, function(x) x$f_value),
    avg_silhouette = c(NA, silhouette_results$avg_silhouette),
    WSS = wss_values
  )
  
  write.csv(all_solutions,
            file = file.path(n_dir, paste0("CLUSCO_All_Solutions_n", n_val, ".csv")),
            row.names = FALSE)
  
  # Save summary report
  sink(file.path(n_dir, paste0("Clustering_Summary_n", n_val, ".txt")))
  cat("================================================================================\n")
  cat(sprintf("TIME SERIES CLUSTERING SUMMARY (n = %d)\n", n_val))
  cat("================================================================================\n\n")
  cat("Algorithm: CLUSCO (Finding compact and well-separated clusters)\n")
  cat(sprintf("Maximum k tested: %d\n", max_k))
  cat(sprintf("Optimal number of clusters: %d\n", optimal_k))
  cat(sprintf("Average silhouette coefficient at optimal k: %.4f\n", 
              max(silhouette_results$avg_silhouette)))
  cat(sprintf("Within-cluster SS at optimal k: %.4f\n", wss_values[optimal_k]))
  cat("\n--- Cluster Sizes ---\n")
  print(cluster_sizes)
  cat("\n--- Model Distribution Across Clusters ---\n")
  print(model_cluster_table)
  cat("\n--- Cluster Centers (Original Scale) ---\n")
  print(round(optimal_centers_original, 4))
  sink()
  
  cat(sprintf("\n✅ All clustering results saved for n = %d\n", n_val))
}

#----------------------------------------------------------------
#--------------------------------------------------------------------
# This csript generate iteration wise cluster centers and cluster function 
# values for each k.
# The results are saved in the "clusco_results" list, which contains the
# final cluster centers and function values for each k.
#--------------------------------------------------------------------
library(cluster)
library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
library(RColorBrewer)
library(scales)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension

# Paths
hc_data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Hybrid ARMA+LogisticMap_data/D4_Data/Combined_All_Models_HC_Results_D4.xlsx"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Hybrid_Analysis/Multinomial_Logistic_Regression"

# Create output directory for detailed visualizations
detailed_viz_dir <- file.path(output_dir, "Time_Series_Clustering", "Detailed_Algorithm_Visualization")
dir.create(detailed_viz_dir, recursive = TRUE, showWarnings = FALSE)

# Sample sizes to process
sample_sizes <- c(1000, 5000)

# Maximum number of clusters to test
max_k <- 10

# ==============================================================================
# DEFINE MODEL NAMES
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data...\n")
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Check for Disequilibrium column
if ("Disequilibrium" %in% names(hc_data_all)) {
  disequil_col <- "Disequilibrium"
} else if ("Disequlibrium" %in% names(hc_data_all)) {
  disequil_col <- "Disequlibrium"
} else {
  stop("Disequilibrium column not found in data!")
}

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon,
         H_Renyi, C_Renyi,
         H_Tsallis, C_Tsallis,
         H_Fisher, C_Fisher,
         all_of(disequil_col))

# Rename for consistency
names(hc_data)[names(hc_data) == disequil_col] <- "Disequilibrium"
hc_data$Model <- factor(hc_data$Model, levels = model_names)

# ==============================================================================
# DEFINE CLUSTERING FEATURES
# ==============================================================================
clustering_features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
                         "C_Shannon", "Disequilibrium", "C_Renyi", 
                         "C_Tsallis", "C_Fisher")

# ==============================================================================
# HELPER FUNCTIONS FOR CLUSCO ALGORITHM
# ==============================================================================

compute_center <- function(data) {
  colMeans(data, na.rm = TRUE)
}

compute_distance <- function(point, center) {
  sqrt(sum((point - center)^2))
}

compute_cluster_function <- function(data, centers) {
  n <- nrow(data)
  k <- nrow(centers)
  
  distances <- matrix(0, nrow = n, ncol = k)
  for (i in 1:k) {
    distances[, i] <- apply(data, 1, function(x) compute_distance(x, centers[i, ]))
  }
  
  assignments <- apply(distances, 1, which.min)
  
  within_cluster_dist <- 0
  for (j in 1:k) {
    cluster_points <- data[assignments == j, , drop = FALSE]
    if (nrow(cluster_points) > 0) {
      within_cluster_dist <- within_cluster_dist + 
        mean(apply(cluster_points, 1, function(x) compute_distance(x, centers[j, ])))
    }
  }
  
  return(-within_cluster_dist / k)
}

optimize_centers_with_history <- function(data, initial_centers, max_iter = 100, tol = 1e-6) {
  centers <- initial_centers
  k <- nrow(centers)
  history <- list()
  
  for (iter in 1:max_iter) {
    old_centers <- centers
    
    # Assign points to nearest center
    n <- nrow(data)
    distances <- matrix(0, nrow = n, ncol = k)
    for (i in 1:k) {
      distances[, i] <- apply(data, 1, function(x) compute_distance(x, centers[i, ]))
    }
    assignments <- apply(distances, 1, which.min)
    
    # Store history
    history[[iter]] <- list(
      iteration = iter,
      centers = centers,
      assignments = assignments
    )
    
    # Update centers
    for (j in 1:k) {
      cluster_points <- data[assignments == j, , drop = FALSE]
      if (nrow(cluster_points) > 0) {
        centers[j, ] <- colMeans(cluster_points)
      }
    }
    
    # Check convergence
    if (max(abs(centers - old_centers)) < tol) {
      break
    }
  }
  
  return(list(centers = centers, history = history))
}

# Enhanced CLUSCO with iteration tracking
clusco_algorithm_detailed <- function(data, max_k) {
  n <- nrow(data)
  p <- ncol(data)
  
  results <- list()
  
  # Step 1: k=1
  cat("\n  Computing initial center (k=1)...\n")
  x1 <- compute_center(data)
  centers <- matrix(x1, nrow = 1)
  f_values <- c(compute_cluster_function(data, centers))
  
  results[[1]] <- list(
    k = 1,
    centers = centers,
    f_value = f_values[1],
    iteration_history = NULL
  )
  
  # Step 2-3: Iteratively add clusters
  for (l in 2:max_k) {
    cat(sprintf("  Computing solution for k=%d...\n", l))
    
    current_centers <- results[[l-1]]$centers
    f_prev <- results[[l-1]]$f_value
    
    best_point <- NULL
    best_improvement <- -Inf
    
    for (i in 1:min(n, 500)) {
      idx <- sample(1:n, 1)
      test_point <- as.numeric(data[idx, ])
      
      temp_centers <- rbind(current_centers, test_point)
      f_new <- compute_cluster_function(data, temp_centers)
      improvement <- f_prev - f_new
      
      if (improvement > best_improvement) {
        best_improvement <- improvement
        best_point <- test_point
      }
    }
    
    # Apply discrete gradient method with history tracking
    initial_centers <- rbind(current_centers, best_point)
    optimization_result <- optimize_centers_with_history(data, initial_centers)
    
    f_l <- compute_cluster_function(data, optimization_result$centers)
    
    results[[l]] <- list(
      k = l,
      centers = optimization_result$centers,
      f_value = f_l,
      iteration_history = optimization_result$history,
      initial_centers = initial_centers
    )
  }
  
  return(results)
}

# ==============================================================================
# ANALYSIS FOR EACH SAMPLE SIZE
# ==============================================================================

for (n_val in sample_sizes) {
  
  cat(sprintf("\n\n================================================================================\n"))
  cat(sprintf("DETAILED VISUALIZATION FOR SAMPLE SIZE n = %d\n", n_val))
  cat(sprintf("================================================================================\n\n"))
  
  # Create subdirectory for this sample size
  n_dir <- file.path(detailed_viz_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # PREPARE DATA
  # --------------------------------------------------------------------------
  
  clustering_data <- hc_data %>%
    filter(n == n_val) %>%
    select(Model, all_of(clustering_features)) %>%
    na.omit()
  
  cat(sprintf("Data points: %d\n", nrow(clustering_data)))
  
  # Extract feature matrix (WITHOUT Model column - only features)
  feature_matrix <- as.matrix(clustering_data[, clustering_features])
  feature_matrix_scaled <- scale(feature_matrix)
  
  # PCA for 2D visualization
  pca_result <- prcomp(feature_matrix_scaled, center = FALSE, scale. = FALSE)
  pca_data <- as.data.frame(pca_result$x[, 1:2])
  pca_data$Model <- clustering_data$Model
  var_explained <- summary(pca_result)$importance[2, 1:2] * 100
  
  # --------------------------------------------------------------------------
  # ADDITIONAL FIGURE 1: ALL DATA POINTS (GEOMETRIC, NO LABELS)
  # --------------------------------------------------------------------------
  
  cat("\n### GENERATING ADDITIONAL FIGURE 1: All Data Points (Geometric) ###\n")
  
  fig_add1 <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 2.5, alpha = 0.6, color = "#2C3E50", shape = 19) +
    labs(title = paste0("All Time Series Data Points (n = ", n_val, ")"),
         #subtitle = "PCA Projection - No Model Labels",
         x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      #plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 1)
    )
  
  print(fig_add1)
  
  ggsave(file.path(n_dir, paste0("Additional_Fig1_All_Points_Geometric_n", n_val, ".pdf")),
         plot = fig_add1, width = 10, height = 8)
  
  cat("✅ Additional Figure 1 saved\n")
  
  # --------------------------------------------------------------------------
  # RUN CLUSCO ALGORITHM WITH DETAILED TRACKING
  # --------------------------------------------------------------------------
  
  cat("\n### RUNNING CLUSCO ALGORITHM WITH ITERATION TRACKING ###\n")
  
  clusco_results_detailed <- clusco_algorithm_detailed(feature_matrix_scaled, max_k)
  
  # --------------------------------------------------------------------------
  # SELECT THREE COMPOSITIONS (k values) FOR DETAILED VISUALIZATION
  # We'll use k=3, k=5, k=7 as three compositions
  # --------------------------------------------------------------------------
  
  compositions <- c(3, 5, 7)
  comp_names <- c("First", "Second", "Third")
  
  for (comp_idx in 1:length(compositions)) {
    
    k_comp <- compositions[comp_idx]
    comp_name <- comp_names[comp_idx]
    
    if (k_comp > max_k) next
    
    cat(sprintf("\n### PROCESSING %s COMPOSITION (k = %d) ###\n", toupper(comp_name), k_comp))
    
    comp_result <- clusco_results_detailed[[k_comp]]
    
    # ------------------------------------------------------------------------
    # FIGURE: UNCONSTRAINED vs CONSTRAINED CLUSTERING
    # ------------------------------------------------------------------------
    
    cat(sprintf("  Generating unconstrained vs constrained plot for k=%d...\n", k_comp))
    
    # Unconstrained: initial centers before optimization
    initial_centers <- comp_result$initial_centers
    colnames(initial_centers) <- clustering_features  # Add column names
    initial_pca <- as.data.frame(predict(pca_result, initial_centers)[, 1:2])
    initial_pca$Type <- "Initial (Unconstrained)"
    
    # Assign points to initial centers
    distances_init <- matrix(0, nrow = nrow(feature_matrix_scaled), ncol = k_comp)
    for (i in 1:k_comp) {
      distances_init[, i] <- apply(feature_matrix_scaled, 1, 
                                   function(x) compute_distance(x, initial_centers[i, ]))
    }
    assignments_init <- apply(distances_init, 1, which.min)
    
    # Constrained: optimized centers
    final_centers <- comp_result$centers
    colnames(final_centers) <- clustering_features  # Add column names
    final_pca <- as.data.frame(predict(pca_result, final_centers)[, 1:2])
    final_pca$Type <- "Optimized (Constrained)"
    
    # Assign points to final centers
    distances_final <- matrix(0, nrow = nrow(feature_matrix_scaled), ncol = k_comp)
    for (i in 1:k_comp) {
      distances_final[, i] <- apply(feature_matrix_scaled, 1, 
                                    function(x) compute_distance(x, final_centers[i, ]))
    }
    assignments_final <- apply(distances_final, 1, which.min)
    
    # Plot unconstrained
    pca_data_unc <- pca_data
    pca_data_unc$Cluster <- factor(assignments_init)
    
    fig_unc <- ggplot(pca_data_unc, aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(size = 2, alpha = 0.5) +
      geom_point(data = initial_pca, aes(x = PC1, y = PC2), 
                 color = "black", size = 2, shape = 4, stroke = 2,
                 inherit.aes = FALSE) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "Unconstrained (Initial Centers)",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_bw(base_size = 12, base_family = "serif") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
    
    # Plot constrained
    pca_data_con <- pca_data
    pca_data_con$Cluster <- factor(assignments_final)
    
    fig_con <- ggplot(pca_data_con, aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(size = 2, alpha = 0.5) +
      geom_point(data = final_pca, aes(x = PC1, y = PC2), 
                 color = "black", size = 2, shape = 8, stroke = 2,
                 inherit.aes = FALSE) +
      scale_color_brewer(palette = "Set1") +
      labs(title = "Constrained (Optimized Centers)",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_bw(base_size = 12, base_family = "serif") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
    
    # Combine both plots
    fig_unc_con <- grid.arrange(fig_unc, fig_con, ncol = 2,
                                top = paste0(comp_name, " Composition (k = ", k_comp, 
                                             ", n = ", n_val, 
                                             "): Unconstrained vs Constrained"))
    
    print(fig_unc_con)
    
    ggsave(file.path(n_dir, paste0("Additional_Fig_", comp_name, "_Comp_Unconstrained_Constrained_k", k_comp, "_n", n_val, ".pdf")),
           plot = fig_unc_con, width = 16, height = 7)
    
    cat(sprintf("  ✅ Unconstrained vs Constrained plot saved for k=%d\n", k_comp))
    
    # ------------------------------------------------------------------------
    # FIGURE: COMPOSITION PLOT (FINAL STATE)
    # ------------------------------------------------------------------------
    
    cat(sprintf("  Generating %s composition plot...\n", comp_name))
    
    fig_comp <- ggplot(pca_data_con, aes(x = PC1, y = PC2, color = Cluster, shape = Cluster)) +
      geom_point(size = 3, alpha = 0.6) +
      geom_point(data = final_pca, aes(x = PC1, y = PC2), 
                 color = "black", size = 2, shape = 8, stroke = 2.5,
                 inherit.aes = FALSE) +
      scale_color_brewer(palette = "Set1") +
      scale_shape_manual(values = 15:(15 + k_comp - 1)) +
      labs(title = paste0(comp_name, " Composition (k = ", k_comp, ", n = ", n_val, ")"),
           #subtitle = "Final clustering result after optimization",
           x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
           y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
      theme_bw(base_size = 13, base_family = "serif") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 11),
        axis.title = element_text(face = "bold", size = 14),
        legend.position = "right"
      )
    
    print(fig_comp)
    
    ggsave(file.path(n_dir, paste0("Additional_Fig_", comp_name, "_Composition_k", k_comp, "_n", n_val, ".pdf")),
           plot = fig_comp, width = 12, height = 8)
    
    cat(sprintf("  ✅ %s composition plot saved\n", comp_name))
    
    # ------------------------------------------------------------------------
    # FIGURE: ITERATIONS 2, 3, 4 OF THE ALGORITHM
    # ------------------------------------------------------------------------
    
    cat(sprintf("  Generating iteration plots (2, 3, 4) for k=%d...\n", k_comp))
    
    iteration_history <- comp_result$iteration_history
    
    if (!is.null(iteration_history) && length(iteration_history) >= 4) {
      
      # Select iterations 2, 3, 4
      iterations_to_plot <- c(2, 3, 4)
      plot_list <- list()
      
      for (iter_idx in iterations_to_plot) {
        
        if (iter_idx <= length(iteration_history)) {
          
          iter_data <- iteration_history[[iter_idx]]
          iter_centers <- iter_data$centers
          colnames(iter_centers) <- clustering_features  # Add column names
          iter_assignments <- iter_data$assignments
          
          # Project centers to PCA space
          iter_centers_pca <- as.data.frame(predict(pca_result, iter_centers)[, 1:2])
          
          # Create plot data
          pca_data_iter <- pca_data
          pca_data_iter$Cluster <- factor(iter_assignments)
          
          # Create plot
          p_iter <- ggplot(pca_data_iter, aes(x = PC1, y = PC2, color = Cluster)) +
            geom_point(size = 2, alpha = 0.5) +
            geom_point(data = iter_centers_pca, aes(x = PC1, y = PC2), 
                       color = "black", size = 2, shape = 8, stroke = 2,
                       inherit.aes = FALSE) +
            scale_color_brewer(palette = "Set1") +
            labs(title = paste0("Iteration ", iter_idx),
                 x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
                 y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
            theme_bw(base_size = 11, base_family = "serif") +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
              legend.position = "none",
              axis.title = element_text(size = 11),
              axis.text = element_text(size = 9)
            )
          
          plot_list[[length(plot_list) + 1]] <- p_iter
        }
      }
      
      # Combine iteration plots
      if (length(plot_list) == 3) {
        fig_iterations <- grid.arrange(
          plot_list[[1]], plot_list[[2]], plot_list[[3]], 
          ncol = 3,
          top = paste0(comp_name, " Composition (k = ", k_comp, ", n = ", n_val, 
                       "): Algorithm Iterations 2, 3, 4")
        )
        
        print(fig_iterations)
        
        ggsave(file.path(n_dir, paste0("Additional_Fig_", comp_name, "_Comp_Iterations_234_k", k_comp, "_n", n_val, ".pdf")),
               plot = fig_iterations, width = 18, height = 6)
        
        cat(sprintf("  ✅ Iteration plots (2, 3, 4) saved for k=%d\n", k_comp))
      } else {
        cat(sprintf("  ⚠️ Not enough iterations recorded for k=%d\n", k_comp))
      }
      
    } else {
      cat(sprintf("  ⚠️ Iteration history not available for k=%d\n", k_comp))
    }
    
    cat(sprintf("✅ All visualizations completed for %s composition (k=%d)\n\n", comp_name, k_comp))
  }
  
  
  
  for (comp_idx in 1:length(compositions)) {
    k_comp <- compositions[comp_idx]
    comp_name <- comp_names[comp_idx]
      }
  
  }

