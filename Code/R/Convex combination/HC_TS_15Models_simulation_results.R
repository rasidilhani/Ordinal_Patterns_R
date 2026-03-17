# ══════════════════════════════════════════════════════════════════════════════
#  Simulation: All Entropies, Complexities, Variances, Semi-lengths
#  D = 3  |  n = 1000, 5000  |  R = 50 replications
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(StatOrdPattHxC)
library(writexl)
library(here)

# ══════════════════════════════════════════════════════════════════════════════
#  GLOBAL PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
D    <- 4
R    <- 50
N    <- c(10000)
r    <- 3.8       # logistic map parameter
f    <- 0.04      # sine frequency
BETA <- 1.5         # Renyi / Tsallis parameter
set.seed(1234567890, kind = "Mersenne-Twister")

# ══════════════════════════════════════════════════════════════════════════════
#  OUTPUT PATHS
# ══════════════════════════════════════════════════════════════════════════════
base_path <- here("Data", "Convex_combination", "D4_Data")

dir.create(file.path(base_path, "timeseries", "n10000"), recursive = TRUE, showWarnings = FALSE)
#dir.create(file.path(base_path, "timeseries", "n5000"), recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL SPECIFICATIONS
# ══════════════════════════════════════════════════════════════════════════════
arma22_model <- list(ar = c(-0.8,  0.1), ma = c(-0.8, 0.1))
ar2_model    <- list(ar = c(-0.8,  0.1), ma = NULL)
ma2_model    <- list(ar = NULL,          ma = c(-0.8, 0.1))

# ══════════════════════════════════════════════════════════════════════════════
#  GENERATOR FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════
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

# ══════════════════════════════════════════════════════════════════════════════
#  COMPLEXITY FUNCTIONS (normalized Jensen-Shannon)
# ══════════════════════════════════════════════════════════════════════════════

# Normalized JS divergence — returns value in [0, 1]
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12) / (m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12) / (m + 1e-12))))
  return(js / log(2))
}

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
#  CORE MEASURE FUNCTION
# ══════════════════════════════════════════════════════════════════════════════
compute_measures <- function(ts_data, label, n_eff) {
  
  z_alpha <- qnorm(0.975)
  prob    <- OPprob(ts_data, emb = D)
  
  # Entropies
  Hs <- HShannon(prob)
  Hr <- HRenyi(prob,   beta = BETA)
  Ht <- HTsallis(prob, beta = BETA)
  Hf <- HFisher(prob)
  
  # Complexities
  Cs <- StatComplexity(prob)
  Cr <- GeneralizedComplexity(prob, Hr)
  Ct <- GeneralizedComplexity(prob, Ht)
  Cf <- FisherBasedComplexity(prob, Hf)
  
  # Variances of entropies
  Var_Hs <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
  Var_Hr <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "R", beta = BETA))
  Var_Ht <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "T", beta = BETA))
  Var_Hf <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "F"))
  
  # Variance of Shannon complexity
  Var_HI  <- suppressWarnings(asymptoticVarHShannonMultinomial(prob, n_eff))
  Var_CI  <- suppressWarnings(varC(prob, n_eff))
  a_ratio <- ifelse(Var_HI > 0, Var_Hs / Var_HI, NA)
  Var_Cs  <- a_ratio * Var_CI
  
  # Semi-lengths (half-width of 95% CI)
  semi <- function(v) ifelse(!is.finite(v) | v <= 0, NA,
                             sqrt(v) / sqrt(n_eff) * z_alpha)
  
  data.frame(
    Model          = label,
    H_Shannon      = Hs,
    H_Renyi        = Hr,
    H_Tsallis      = Ht,
    H_Fisher       = Hf,
    C_Shannon      = Cs,
    C_Renyi        = Cr,
    C_Tsallis      = Ct,
    C_Fisher       = Cf,
    Var_H_Shannon  = Var_Hs,
    Var_H_Renyi    = Var_Hr,
    Var_H_Tsallis  = Var_Ht,
    Var_H_Fisher   = Var_Hf,
    Var_C_Shannon  = Var_Cs,
    Semi_H_Shannon = semi(Var_Hs),
    Semi_H_Renyi   = semi(Var_Hr),
    Semi_H_Tsallis = semi(Var_Ht),
    Semi_H_Fisher  = semi(Var_Hf),
    Semi_C_Shannon = semi(Var_Cs),
    stringsAsFactors = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN SIMULATION LOOP
# ══════════════════════════════════════════════════════════════════════════════
hc_sheets <- list()   # stores HC results per sample size

for (n_val in N) {
  
  message("\n▶ n = ", n_val)
  n_eff       <- n_val - (D - 1)
  hc_all_reps <- vector("list", R)
  
  # Deterministic signals — same every rep
  x_log <- logistic(n_val, r)
  x_sin <- sine(n_val, f)
  
  for (rep_i in seq_len(R)) {
    
    message("  Rep ", rep_i)
    
    # Stochastic signals — new draw each rep
    x_arma <- arma22(n_val)
    x_ar2  <- ar2(n_val)
    x_ma2  <- ma2(n_val)
    
    # ── Series list ────────────────────────────────────────────────────────
    # To add a model:    uncomment the line
    # To remove a model: comment the line with #
    series_list <- list(
      
      # Pure models
      "ARMA(2,2)" = x_arma,
      "AR(2)"     = x_ar2,
      "MA(2)"     = x_ma2,
      "Logistic"  = x_log,
      "Sine"      = x_sin,
      
      # ARMA + Sine
      "ARMA+Sine(w=0.1)" = 0.1*x_arma + 0.9*x_sin,
      "ARMA+Sine(w=0.2)" = 0.2*x_arma + 0.8*x_sin,
      "ARMA+Sine(w=0.3)" = 0.3*x_arma + 0.7*x_sin,
      
      # AR2 + Logistic
      "AR2+Logistic(w=0.1)" = 0.1*x_ar2 + 0.9*x_log,
      
      # AR2 + Sine
      "AR2+Sine(w=0.8)" = 0.8*x_ar2 + 0.2*x_sin,
      
      # MA2 + Logistic
      "MA2+Logistic(w=0.2)" = 0.2*x_ma2 + 0.8*x_log,
      "MA2+Logistic(w=0.7)" = 0.7*x_ma2 + 0.3*x_log,
      
      # MA2 + Sine
      "MA2+Sine(w=0.4)" = 0.4*x_ma2 + 0.6*x_sin,
      "MA2+Sine(w=0.6)" = 0.6*x_ma2 + 0.4*x_sin,
      "MA2+Sine(w=0.8)" = 0.8*x_ma2 + 0.2*x_sin
      
    )
    
    # ── Compute all measures ───────────────────────────────────────────────
    hc_rep <- bind_rows(mapply(
      function(s, nm) compute_measures(s, nm, n_eff),
      series_list, names(series_list),
      SIMPLIFY = FALSE
    ))
    hc_rep$Rep <- rep_i
    hc_rep$n   <- n_val
    hc_all_reps[[rep_i]] <- hc_rep
    
    # ── Save time series (wide format, one file per rep) ───────────────────
    ts_wide     <- as.data.frame(series_list)
    ts_wide$t   <- seq_len(n_val)
    ts_wide     <- ts_wide[, c("t", setdiff(names(ts_wide), "t"))]
    
    ts_fname <- file.path(base_path, "timeseries",
                          paste0("n", n_val),
                          paste0("Rep_", rep_i, ".xlsx"))
    write_xlsx(ts_wide, path = ts_fname)
  }
  
  # ── Collect HC results for this sample size ──────────────────────────────
  hc_df <- bind_rows(hc_all_reps) %>%
    select(n, Rep, Model, everything())
  
  hc_sheets[[paste0("n", n_val)]] <- hc_df
  message("  ✓ Done n = ", n_val)
}

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE HC RESULTS
# ══════════════════════════════════════════════════════════════════════════════
hc_fname <- file.path(base_path, "HC_Results_D4_n10000.xlsx")
write_xlsx(hc_sheets, path = hc_fname)

