# ======================================================
# 📊 Shannon Entropy–Complexity Simulation for Selected Models
# ======================================================

library(StatOrdPattHxC)
library(dplyr)
library(openxlsx)

set.seed(1234567890, kind = "Mersenne-Twister")

# --- Parameters ---
D <- 3
N_vals <- c(5000, 10000)
R <- 100

# --- Model Definitions ---
arma11_models <- list(ARMA11_M1 = list(ar = c(0.8), ma = c(0.8)))
arma22_models    <- list(ARMA22_M1    = list(ar = c(0.1, 0.8), ma = c(0.1, 0.8)))
ma1_models    <- list(MA1_M2   = list(ma = c(0.1)))

# --- Function to simulate and compute Shannon metrics ---
simulate_model <- function(ar = NULL, ma = NULL, model_name, N_vals, R, D) {
  hc_results <- data.frame()
  ts_store <- list()
  
  for (n in N_vals) {
    for (r in 1:R) {
      ts_data <- arima.sim(model = list(ar = ar, ma = ma), n = n)
      
      # Compute Shannon metrics
      ProbTS <- OPprob(ts_data, emb = D)
      Hs <- HShannon(ProbTS)
      C_Shannon <- StatComplexity(ProbTS)
      Var_H <- sigma2q(ts_data, emb = D, ent = "S")
      
      # Store entropy-complexity results
      hc_results <- rbind(
        hc_results,
        data.frame(
          H_Shannon = Hs,
          C_Shannon = C_Shannon,
          Var_Shannon = Var_H,
          n = n,
          Model = model_name,
          stringsAsFactors = FALSE
        )
      )
      
      # Store time series
      ts_store[[paste0(model_name, "_N", n, "_R", r)]] <- data.frame(
        Index = seq_along(ts_data),
        Value = ts_data,
        Model = model_name,
        N = n,
        Rep = r
      )
    }
  }
  
  return(list(hc = hc_results, ts = ts_store))
}

# --- Run simulations for all three models ---
results <- list()
results$ARMA11_M1 <- simulate_model(ar = c(0.8), ma = c(0.8), model_name = "ARMA11_M1", N_vals, R, D)
results$ARMA22_M1    <- simulate_model(ar = c(0.1, 0.8), ma = c(0.1, 0.8), model_name = "ARMA22_M1", N_vals, R, D)
results$MA1_M2    <- simulate_model(ar = NULL, ma = c(0.1), model_name = "MA1_M2", N_vals, R, D)

# --- Output directory ---
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/New ARMA Time series results/New_Case 4"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Save HC (Entropy–Complexity) results ---
hc_path <- file.path(output_dir, "HC_Shannon_Results_New_Case 4.xlsx")
wb_hc <- createWorkbook()
for (model_name in names(results)) {
  addWorksheet(wb_hc, model_name)
  writeData(wb_hc, model_name, results[[model_name]]$hc)
}
saveWorkbook(wb_hc, hc_path, overwrite = TRUE)
cat("✅ HC results saved to:", hc_path, "\n")

# --- Save Time Series data ---
ts_path <- file.path(output_dir, "Time_Series_Data_New_Case 4.xlsx")
wb_ts <- createWorkbook()
for (model_name in names(results)) {
  addWorksheet(wb_ts, model_name)
  ts_combined <- bind_rows(results[[model_name]]$ts)
  writeData(wb_ts, model_name, ts_combined)
}
saveWorkbook(wb_ts, ts_path, overwrite = TRUE)
cat("✅ Time series data saved to:", ts_path, "\n")

# --- Summary message ---
cat("\n🎉 Simulation complete!")
cat("\n🧩 Models simulated: ARMA11_M1, ARMA22_M1, MA1_M2")
cat("\n📏 Lengths: N = 5000, 10000 | Replications = 100")
cat("\n💾 Results stored in:", output_dir, "\n")

# ======================================================
# 📊 Shannon Entropy–Complexity Simulation for Selected Models
# ======================================================

library(StatOrdPattHxC)
library(dplyr)
library(openxlsx)

set.seed(1234567890, kind = "Mersenne-Twister")

# --- Parameters ---
D <- 3
N_vals <- c(5000, 10000)
R <- 100

# --- Model Definitions ---
arma11_models <- list(ARMA11_M3 = list(ar = c(-0.8), ma = c(-0.8)))
arma22_models    <- list(ARMA22_M4    = list(ar = c(-0.8, -0.1), ma = c(-0.8, -0.1)))
ma1_models    <- list(MA1_M4   = list(ma = c(-0.1)))

# --- Function to simulate and compute Shannon metrics ---
simulate_model <- function(ar = NULL, ma = NULL, model_name, N_vals, R, D) {
  hc_results <- data.frame()
  ts_store <- list()
  
  for (n in N_vals) {
    for (r in 1:R) {
      ts_data <- arima.sim(model = list(ar = ar, ma = ma), n = n)
      
      # Compute Shannon metrics
      ProbTS <- OPprob(ts_data, emb = D)
      Hs <- HShannon(ProbTS)
      C_Shannon <- StatComplexity(ProbTS)
      Var_H <- sigma2q(ts_data, emb = D, ent = "S")
      
      # Store entropy-complexity results
      hc_results <- rbind(
        hc_results,
        data.frame(
          H_Shannon = Hs,
          C_Shannon = C_Shannon,
          Var_Shannon = Var_H,
          n = n,
          Model = model_name,
          stringsAsFactors = FALSE
        )
      )
      
      # Store time series
      ts_store[[paste0(model_name, "_N", n, "_R", r)]] <- data.frame(
        Index = seq_along(ts_data),
        Value = ts_data,
        Model = model_name,
        N = n,
        Rep = r
      )
    }
  }
  
  return(list(hc = hc_results, ts = ts_store))
}

# --- Run simulations for all three models ---
results <- list()
results$ARMA11_M3 <- simulate_model(ar = c(-0.8), ma = c(-0.8), model_name = "ARMA11_M3", N_vals, R, D)
results$ARMA22_M4    <- simulate_model(ar = c(-0.8, -0.1), ma = c(-0.8, -0.1), model_name = "ARMA22_M4", N_vals, R, D)
results$MA1_M4    <- simulate_model(ar = NULL, ma = c(-0.1), model_name = "MA1_M4", N_vals, R, D)

# --- Output directory ---
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/New ARMA Time series results/New_Case 5"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Save HC (Entropy–Complexity) results ---
hc_path <- file.path(output_dir, "HC_Shannon_Results_New_Case 5.xlsx")
wb_hc <- createWorkbook()
for (model_name in names(results)) {
  addWorksheet(wb_hc, model_name)
  writeData(wb_hc, model_name, results[[model_name]]$hc)
}
saveWorkbook(wb_hc, hc_path, overwrite = TRUE)
cat("✅ HC results saved to:", hc_path, "\n")

# --- Save Time Series data ---
ts_path <- file.path(output_dir, "Time_Series_Data_New_Case 5.xlsx")
wb_ts <- createWorkbook()
for (model_name in names(results)) {
  addWorksheet(wb_ts, model_name)
  ts_combined <- bind_rows(results[[model_name]]$ts)
  writeData(wb_ts, model_name, ts_combined)
}
saveWorkbook(wb_ts, ts_path, overwrite = TRUE)
cat("✅ Time series data saved to:", ts_path, "\n")

# --- Summary message ---
cat("\n🎉 Simulation complete!")
cat("\n🧩 Models simulated: ARMA11_M3, ARMA22_M4, MA1_M4")
cat("\n📏 Lengths: N = 5000, 10000 | Replications = 100")
cat("\n💾 Results stored in:", output_dir, "\n")
# ======================================================


# ======================================================
# 📊 Shannon Entropy–Complexity Simulation for Selected Models
# ======================================================

library(StatOrdPattHxC)
library(dplyr)
library(openxlsx)

set.seed(1234567890, kind = "Mersenne-Twister")

# --- Parameters ---
D <- 3
N_vals <- c(5000, 10000)
R <- 100

# --- Model Definitions ---
arma22_m2_models <- list(ARMA22_M2 = list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1)))
arma22_m3_models    <- list(ARMA22_M3  = list(ar = c(0.1, -0.8), ma = c(0.1, -0.8)))
ma2_models    <- list(MA2_M2   = list(ma = c(-0.8 ,0.1)))

# --- Function to simulate and compute Shannon metrics ---
simulate_model <- function(ar = NULL, ma = NULL, model_name, N_vals, R, D) {
  hc_results <- data.frame()
  ts_store <- list()
  
  for (n in N_vals) {
    for (r in 1:R) {
      ts_data <- arima.sim(model = list(ar = ar, ma = ma), n = n)
      
      # Compute Shannon metrics
      ProbTS <- OPprob(ts_data, emb = D)
      Hs <- HShannon(ProbTS)
      C_Shannon <- StatComplexity(ProbTS)
      Var_H <- sigma2q(ts_data, emb = D, ent = "S")
      
      # Store entropy-complexity results
      hc_results <- rbind(
        hc_results,
        data.frame(
          H_Shannon = Hs,
          C_Shannon = C_Shannon,
          Var_Shannon = Var_H,
          n = n,
          Model = model_name,
          stringsAsFactors = FALSE
        )
      )
      
      # Store time series
      ts_store[[paste0(model_name, "_N", n, "_R", r)]] <- data.frame(
        Index = seq_along(ts_data),
        Value = ts_data,
        Model = model_name,
        N = n,
        Rep = r
      )
    }
  }
  
  return(list(hc = hc_results, ts = ts_store))
}

# --- Run simulations for all three models ---
results <- list()
results$ARMA22_M2 <- simulate_model(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), model_name = "ARMA22_M2", N_vals, R, D)
results$ARMA22_M3    <- simulate_model(ar = c(0.1, -0.8), ma = c(0.1, -0.8), model_name = "ARMA22_M3", N_vals, R, D)
results$MA2_M2    <- simulate_model(ar = NULL, ma = c(-0.8, 0.1), model_name = "MA2_M2", N_vals, R, D)

# --- Output directory ---
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/New ARMA Time series results/New_Case 6"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Save HC (Entropy–Complexity) results ---
hc_path <- file.path(output_dir, "HC_Shannon_Results_New_Case 6.xlsx")
wb_hc <- createWorkbook()
for (model_name in names(results)) {
  addWorksheet(wb_hc, model_name)
  writeData(wb_hc, model_name, results[[model_name]]$hc)
}
saveWorkbook(wb_hc, hc_path, overwrite = TRUE)
cat("✅ HC results saved to:", hc_path, "\n")

# --- Save Time Series data ---
ts_path <- file.path(output_dir, "Time_Series_Data_New_Case 6.xlsx")
wb_ts <- createWorkbook()
for (model_name in names(results)) {
  addWorksheet(wb_ts, model_name)
  ts_combined <- bind_rows(results[[model_name]]$ts)
  writeData(wb_ts, model_name, ts_combined)
}
saveWorkbook(wb_ts, ts_path, overwrite = TRUE)
cat("✅ Time series data saved to:", ts_path, "\n")

# --- Summary message ---
cat("\n🎉 Simulation complete!")
cat("\n🧩 Models simulated: ARMA22_M2, ARMA22_M3, MA2_M2")
cat("\n📏 Lengths: N = 5000, 10000 | Replications = 100")
cat("\n💾 Results stored in:", output_dir, "\n")
# ======================================================
