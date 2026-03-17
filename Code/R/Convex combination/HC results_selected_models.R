#--------------------------------------------------------------------------
# ══════════════════════════════════════════════════════════════════════════════
#  Simulation: All Entropies, Complexities, Variances, Semi-lengths
#  Output:
#    HC_Results.xlsx          — 2 sheets (n1000, n5000)
#    timeseries/n1000/Rep_1.xlsx ... Rep_R.xlsx
#    timeseries/n5000/Rep_1.xlsx ... Rep_R.xlsx
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(StatOrdPattHxC)
library(writexl)
library(here)

# ── Global parameters ─────────────────────────────────────────────────────────
D    <- 4
R    <- 100
N    <- c(1000, 5000)
r    <- 3.8
f    <- 0.04
BETA <- 1.5          # parameter for Renyi and Tsallis entropies
set.seed(1234567890, kind = "Mersenne-Twister")

# ── Output folders ────────────────────────────────────────────────────────────
base_path <- here("Data", "Convex_combination", "D3_Data")

dir.create(file.path(base_path, "timeseries/n1000"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_path, "timeseries/n5000"), recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  Model Specifications
# ══════════════════════════════════════════════════════════════════════════════
arma22_model <- list(ar = c(-0.8,  0.1), ma = c(-0.8, 0.1))
ar2_model    <- list(ar = c(-0.8,  0.1), ma = NULL)
ma2_model    <- list(ar = NULL,          ma = c(-0.8, 0.1))

# ── Generator functions ───────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

arma22   <- function(n) as.numeric(normalize(arima.sim(model = arma22_model, n = n)))
ar2      <- function(n) as.numeric(normalize(arima.sim(model = ar2_model,    n = n)))
ma2      <- function(n) as.numeric(normalize(arima.sim(model = ma2_model,    n = n)))

logistic <- function(n, r) {
  x <- numeric(n); x[1] <- 0.5
  for (i in 2:n) x[i] <- r * x[i-1] * (1 - x[i-1])
  return(x)
}

sine <- function(n, f) as.numeric(sin(2 * pi * f * 1:n))

# ── Normalized Jensen-Shannon divergence ──────────────────────────────────────
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12) / (m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12) / (m + 1e-12))))
  return(js / log(2))   # normalized to [0, 1]
}

# ── Complexity functions ──────────────────────────────────────────────────────
GeneralizedComplexity <- function(prob, entropy_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * entropy_value)
}

FisherBasedComplexity <- function(prob, fisher_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * fisher_value)
}

# ══════════════════════════════════════════════════════════════════════════════
#  Core function: compute all measures for one series
# ══════════════════════════════════════════════════════════════════════════════
compute_measures <- function(ts_data, label, n_eff) {
  # n_eff = effective sample size = n - (D - 1)
  z_alpha <- qnorm(0.975)   # 95% CI
  prob    <- OPprob(ts_data, emb = D)
  
  # ── Entropies ──────────────────────────────────────────────────────────────
  Hs <- HShannon(prob)
  Hr <- HRenyi(prob,   beta = BETA)
  Ht <- HTsallis(prob, beta = BETA)
  Hf <- HFisher(prob)
  
  # ── Complexities ───────────────────────────────────────────────────────────
  Cs <- StatComplexity(prob)               # Shannon-based
  Cr <- GeneralizedComplexity(prob, Hr)    # Renyi-based
  Ct <- GeneralizedComplexity(prob, Ht)    # Tsallis-based
  Cf <- FisherBasedComplexity(prob, Hf)    # Fisher-based
  
  # ── Variances of entropies ─────────────────────────────────────────────────
  Var_Hs <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
  Var_Hr <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "R", beta = BETA))
  Var_Ht <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "T", beta = BETA))
  Var_Hf <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "F"))
  
  # ── Variance of Shannon complexity ────────────────────────────────────────
  Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(prob, n_eff))
  Var_CI <- suppressWarnings(varC(prob, n_eff))
  a_ratio <- ifelse(Var_HI > 0, Var_Hs / Var_HI, NA)
  Var_Cs  <- a_ratio * Var_CI
  
  # ── Semi-lengths (half-width of 95% CI) ───────────────────────────────────
  semi <- function(v) {
    ifelse(!is.finite(v) | v <= 0, NA, sqrt(v) / sqrt(n_eff) * z_alpha)
  }
  
  data.frame(
    Model = label,
    
    # Entropies
    H_Shannon = Hs,
    H_Renyi   = Hr,
    H_Tsallis = Ht,
    H_Fisher  = Hf,
    
    # Complexities
    C_Shannon = Cs,
    C_Renyi   = Cr,
    C_Tsallis = Ct,
    C_Fisher  = Cf,
    
    # Variances of entropies
    Var_H_Shannon = Var_Hs,
    Var_H_Renyi   = Var_Hr,
    Var_H_Tsallis = Var_Ht,
    Var_H_Fisher  = Var_Hf,
    
    # Variance of Shannon complexity
    Var_C_Shannon = Var_Cs,
    
    # Semi-lengths of entropies
    Semi_H_Shannon = semi(Var_Hs),
    Semi_H_Renyi   = semi(Var_Hr),
    Semi_H_Tsallis = semi(Var_Ht),
    Semi_H_Fisher  = semi(Var_Hf),
    
    # Semi-length of Shannon complexity
    Semi_C_Shannon = semi(Var_Cs),
    
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  Main Simulation Loop
# ══════════════════════════════════════════════════════════════════════════════
hc_sheets <- list()

for (n_val in N) {
  
  message("\n▶ Simulating n = ", n_val)
  n_eff       <- n_val - (D - 1)
  hc_all_reps <- vector("list", R)
  
  # Deterministic signals — identical across all reps
  x_log <- logistic(n_val, r)
  x_sin <- sine(n_val, f)
  
  for (rep_i in seq_len(R)) {
    
    message("  Rep ", rep_i)
    
    # ── Stochastic pure signals ──────────────────────────────────────────────
    x_arma <- arma22(n_val)
    x_ar2  <- ar2(n_val)
    x_ma2  <- ma2(n_val)
    
    # ── Full series list ─────────────────────────────────────────────────────
    series_list <- list(
      
      # Pure models
      "ARMA(2,2)" = x_arma,
      "AR(2)"     = x_ar2,
      "MA(2)"     = x_ma2,
      "Logistic"  = x_log,
      "Sine"      = x_sin,
      
      # ── ARMA + Sine ───────────────────────────────────────────────────────
      "ARMA+Sine(w=0.1)" = 0.1*x_arma + 0.9*x_sin,
      "ARMA+Sine(w=0.2)" = 0.2*x_arma + 0.8*x_sin,
      "ARMA+Sine(w=0.3)" = 0.3*x_arma + 0.7*x_sin,
      
      # ── AR2 + Logistic ────────────────────────────────────────────────────
      "AR2+Logistic(w=0.1)" = 0.1*x_ar2 + 0.9*x_log,
      
      # ── AR2 + Sine ────────────────────────────────────────────────────────
      "AR2+Sine(w=0.8)" = 0.8*x_ar2 + 0.2*x_sin,
      
      # ── MA2 + Logistic ────────────────────────────────────────────────────
      "MA2+Logistic(w=0.2)" = 0.2*x_ma2 + 0.8*x_log,
      "MA2+Logistic(w=0.7)" = 0.7*x_ma2 + 0.3*x_log,
      
      # ── MA2 + Sine ────────────────────────────────────────────────────────
      "MA2+Sine(w=0.4)" = 0.4*x_ma2 + 0.6*x_sin,
      "MA2+Sine(w=0.6)" = 0.6*x_ma2 + 0.4*x_sin,
      "MA2+Sine(w=0.8)" = 0.8*x_ma2 + 0.2*x_sin
    )
    
    # ── Compute all measures for each series ─────────────────────────────────
    hc_rep <- bind_rows(mapply(
      function(s, nm) compute_measures(s, nm, n_eff),
      series_list, names(series_list),
      SIMPLIFY = FALSE
    ))
    hc_rep$Rep <- rep_i
    hc_rep$n   <- n_val
    hc_all_reps[[rep_i]] <- hc_rep
    
    # ── Save wide-format time series for this rep ─────────────────────────────
    ts_wide   <- as.data.frame(series_list)
    ts_wide$t <- seq_len(n_val)
    ts_wide   <- ts_wide[, c("t", setdiff(names(ts_wide), "t"))]
    
    ts_fname <- file.path(base_path, paste0("timeseries/n", n_val),
                          paste0("Rep_", rep_i, ".xlsx"))
    write_xlsx(ts_wide, path = ts_fname)
  }
  
  # ── Collect HC results for this sample size ───────────────────────────────
  hc_df <- bind_rows(hc_all_reps) %>%
    select(n, Rep, Model, everything())
  
  hc_sheets[[paste0("n", n_val)]] <- hc_df
  message("  ✓ Done n = ", n_val)
}

# ══════════════════════════════════════════════════════════════════════════════
#  Save HC Results
# ══════════════════════════════════════════════════════════════════════════════
# Define base path
base_path <- here("Data", "Convex_combination", "D3_Data")

# Create directory if it doesn't exist
if (!dir.exists(base_path)) {
  dir.create(base_path, recursive = TRUE)
  cat(sprintf("Created directory: %s\n", base_path))
}

summary_path <- file.path(base_path, "HC_Results.xlsx")
write_xlsx(hc_sheets, path = summary_path)

