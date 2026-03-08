library(StatOrdPattHxC)
library(dplyr)
library(writexl)
library(here)

# ==============================================================================
# PARAMETERS
# ==============================================================================
set.seed(1234567890, kind = "Mersenne-Twister")
BETA <- 1.5               # Parameter for Rényi / Tsallis
#N    <- c(1000, 5000)     # Sample sizes
N    <- c(50000)           # Sample sizes (for testing)
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
#base_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Hybrid model_data/D3_Data"
base_path <- here("Data", "Hybrid model_data", paste0("D", D, "_Data"))

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

