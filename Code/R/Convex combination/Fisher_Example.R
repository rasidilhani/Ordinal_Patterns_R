library(readxl)
library(StatOrdPattHxC)

Rep_1 <- read_excel("Data/Convex_combination/D4_Data/timeseries/n1000/Rep_1.xlsx")
#View(Rep_1)

df<- data.frame(Rep_1)

Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12)/(m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12)/(m + 1e-12))))
  return(js / log(2))
}

op_ARMA <- OPseq(df$ARMA.2.2., 4) 
prob_ARMA <- OPprob(df$ARMA.2.2., 4)
H_F <- HFisher(prob_ARMA)

Pe <- rep(1 / length(prob_ARMA), length(prob_ARMA))
JS <- Jensen_Shannon(prob_ARMA, Pe)
C_f <- JS * H_F

Var_F <- sigma2q(df$ARMA.2.2., 4, ent = "F")


getAnywhere(HFisher)


##################################################
# Ferri Fisher Formula

library(readxl)
library(StatOrdPattHxC)

# Load data
Rep_1 <- read_excel("Data/Convex_combination/D4_Data/timeseries/n1000/Rep_1.xlsx")
df <- data.frame(Rep_1)

# Jensen–Shannon divergence
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12)/(m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12)/(m + 1e-12))))
  return(js / log(2))
}

# ---------------------------------------------------------
# Fisher Information (Ferri et al. 2012 Physica A formula)
# ---------------------------------------------------------
Fisher_Ferri <- function(p) {
  eps <- 1e-12
  N <- length(p)
  total <- 0
  
  for (i in 1:(N-1)) {
    num <- (p[i+1] - p[i])^2
    den <- p[i+1] + p[i] + eps
    total <- total + num / den
  }
  
  return(0.5 * total)
}


# Ordinal patterns (canonical order)
op_ARMA <- OPseq(df$ARMA.2.2., 4)
prob_ARMA <- OPprob(df$ARMA.2.2., 4)

# Fisher using Ferri formula
F_Ferri <- Fisher_Ferri(prob_ARMA)

# Fisher using package (Hellinger form)
F_Hellinger <- HFisher(prob_ARMA)

# Uniform distribution
Pe <- rep(1 / length(prob_ARMA), length(prob_ARMA))

# JS divergence + complexity (using Hellinger Fisher and Ferrir Fisher)
JS <- Jensen_Shannon(prob_ARMA, Pe)
C_f_Hellinger <- JS * F_Hellinger
C_f_Ferri <- JS * F_Ferri


# Variance
Var_F <- sigma2q(df$ARMA.2.2., 4, ent = "F")

#========================================
# All model results for FIM saved 

library(readxl)
library(writexl)
library(StatOrdPattHxC)
library(here)

# ── Model column names (in data) → readable labels ────────────────────────────
model_labels <- c(
  "ARMA.2.2."           = "ARMA(2,2)",
  "AR.2."               = "AR(2)",
  "MA.2."               = "MA(2)",
  "Logistic"            = "Logistic",
  "Sine"                = "Sine",
  "ARMA.Sine.w.0.1."    = "ARMA+Sine(w=0.1)",
  "ARMA.Sine.w.0.2."    = "ARMA+Sine(w=0.2)",
  "ARMA.Sine.w.0.3."    = "ARMA+Sine(w=0.3)",
  "AR2.Logistic.w.0.1." = "AR2+Logistic(w=0.1)",
  "AR2.Sine.w.0.8."     = "AR2+Sine(w=0.8)",
  "MA2.Logistic.w.0.2." = "MA2+Logistic(w=0.2)",
  "MA2.Logistic.w.0.7." = "MA2+Logistic(w=0.7)",
  "MA2.Sine.w.0.4."     = "MA2+Sine(w=0.4)",
  "MA2.Sine.w.0.6."     = "MA2+Sine(w=0.6)",
  "MA2.Sine.w.0.8."     = "MA2+Sine(w=0.8)"
)

# ── Functions ──────────────────────────────────────────────────────────────────

Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12) / (m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12) / (m + 1e-12))))
  return(js / log(2))
}

Fisher_Ferri <- function(p) {
  total <- 0
  for (i in 1:(length(p) - 1)) {
    total <- total + (p[i+1] - p[i])^2 / (p[i+1] + p[i] + 1e-12)
  }
  return(0.5 * total)
}

# ── Main batch loop ────────────────────────────────────────────────────────────

base_dir <- here("Data", "Convex_combination", "D4_Data", "timeseries")
results  <- list()

for (n in c(1000, 10000)) {
  for (rep in 1:50) {
    
    df <- data.frame(read_excel(file.path(base_dir, paste0("n", n),
                                          paste0("Rep_", rep, ".xlsx"))))
    
    for (col in names(model_labels)) {
      prob <- OPprob(df[[col]], 4)
      Pe   <- rep(1 / length(prob), length(prob))   # uniform distribution
      
      F_Ferri       <- Fisher_Ferri(prob)
      F_Hellinger   <- HFisher(prob)
      JS            <- Jensen_Shannon(prob, Pe)
      
      results[[length(results) + 1]] <- data.frame(
        n             = n,
        Rep           = rep,
        Model         = model_labels[col],
        F_Ferri       = F_Ferri,
        F_Hellinger   = F_Hellinger,
        JS            = JS,
        C_f_Ferri     = JS * F_Ferri,
        C_f_Hellinger = JS * F_Hellinger
      )
    }
  }
  cat("Done: n =", n, "\n")
}

# ── Save results ───────────────────────────────────────────────────────────────

results_df <- do.call(rbind, results)
write_xlsx(results_df, here("Data", "Convex_combination", "D4_Data",
                            "Fisher_Complexity_Results.xlsx"))
cat("Results saved.\n")

