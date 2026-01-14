# --- Packages ---
library(StatOrdPattHxC)   
library(writexl)  
library(readxl)
library(ggplot2)          

# --- Settings ---
set.seed(1234567890, kind = "Mersenne-Twister")
D  <- 3                   # embedding dimension
n  <- 5000                 # length of time series
R  <- 10                 # number of replications

# ARMA(2,2) parameters
#ar_coefs <- c(-0.8, 0.1, 0.3, -0.2)
#ma_coefs <- c(-0.8, 0.1, -0.4, 0.2)
#ar_coefs <- c(-0.8897, -0.4858)  # AR coefficients
#ma_coefs <- c(-0.2279, 0.2488)   # MA coefficients
ar_coefs <- c(-0.8, 0.15)
ma_coefs <- c(-0.8, 0.15)

# --- Storage objects ---
ts_list      <- vector("list", R)      # store time series
results_list <- vector("list", R)      # store H, C, variances, etc.

# --- Main simulation loop ---
for (r in 1:R) {
  ## 1. Simulate ARMA(2,2) time series
  arma_ts <- arima.sim(
    n = n,
    #model = list(ar = ar_coefs, ma = ma_coefs), sd = sqrt(0.1796)
    model = list(ar = ar_coefs, ma = ma_coefs)
  )
  
  ts_list[[r]] <- data.frame(
    Replication = r,
    Time       = 1:n,
    X          = as.numeric(arma_ts)
  )
  
  ## 2. Ordinal pattern probability
  Prob <- OPprob(arma_ts, emb = D)
  
  ## 3. Shannon entropy and statistical complexity
  H <- HShannon(Prob)
  C <- StatComplexity(Prob)
  
  ## 4. Asymptotic variances and related quantities
  var_HD <- sigma2q(arma_ts, emb = D)               # variance of H (D-dependent) [web:19]
  var_HI <- asymptoticVarHShannonMultinomial(Prob, n - 2)  # variance of H (I-multinomial) [web:19]
  Var_C  <- varC(Prob, n - 2)                       # variance of complexity [web:22]
  
  a      <- var_HD / var_HI
  Var_CD <- a * Var_C
  
  ## 5. Semi-lengths for confidence intervals
  z_alpha <- qnorm(1 - 0.05 / 2)
  SemiLength_H <- sqrt(var_HD) * z_alpha / sqrt(n - 2)
  SemiLength_C <- sqrt(Var_CD) * z_alpha / sqrt(n - 2)
  
  ## 6. Store scalar results
  results_list[[r]] <- data.frame(
    Replication   = r,
    H             = H,
    C             = C,
    var_HD        = var_HD,
    var_HI        = var_HI,
    Var_C         = Var_C,
    a             = a,
    Var_CD        = Var_CD,
    SemiLength_H  = SemiLength_H,
    SemiLength_C  = SemiLength_C
    
  )
}

# --- Bind results ---
ts_all      <- do.call(rbind, ts_list)       
results_all <- do.call(rbind, results_list)  

# --- Save to Excel (in working directory) ---
#write_xlsx(ts_all,      path = "ARMA22_time_series_R30.xlsx")   
write_xlsx(results_all, path = "ARMA22_HC_results_R30.xlsx")    


ARMA22_results <- read_excel("ARMA22_HC_results_R100.xlsx")
#View(ARMA22_results)
ggplot(ARMA22_results, aes(x=H, y=C)) + 
  geom_point(alpha = 0.6, color = "blue", size = 2) +
                        theme_minimal() + 
  ggtitle("H-C Plane for ARMA(2,2) Time Series (R=30)") + xlab("Shannon Entropy (H)") + ylab("Statistical Complexity (C)")

collinnearity <- cor(ARMA22_results$H, ARMA22_results$C)
print(paste("Correlation between H and C:", collinnearity))
SD_H <- sd(ARMA22_results$H)
print(paste("Standard Deviation of H:", SD_H))
SD_C <- sd(ARMA22_results$C)
print(paste("Standard Deviation of C:", SD_C))

##-----------------------------------------------------------------------------------------------


# Include the tau parameter and modify the code accordingly
library(StatOrdPattHxC)
library(dplyr)
library(writexl)
library(ggplot2)

# --- Parameters ---
set.seed(1234567890, kind = "Mersenne-Twister")
D <- 5                     # Embedding dimension
tau <- 2                   # Delay parameter
N <- 100                   # Sample size
R <- 100                   # Number of replications

# --- Single Model Definition ---
model <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), type = "ARMA22_M2")

# --- Function to simulate data and compute entropy/complexity ---
generate_model_data <- function(model, n, D, tau, R) {
  results <- list()
  ts_store <- list()
  
  for (r in 1:R) {
    ts_data <- arima.sim(model = model, n = n)
    ProbTS <- OPprob(ts_data, emb = D)
    
    # Entropy and complexity
    Hs <- HShannon(ProbTS)
    Cs <- StatComplexity(ProbTS)
    
    # Variances (keep negative)
    Var_HD <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(ProbTS, n - (D - 1) * tau - 1))
    a_ratio <- Var_HD / Var_HI
    Var_CI <- suppressWarnings(varC(ProbTS, n - (D - 1) * tau - 1))
    Var_CD <- a_ratio * Var_CI
    
    # Semi-lengths: NA if variance <= 0
    SemiLengthH <- ifelse(!is.finite(Var_HD) | Var_HD <= 0, NA,
                          sqrt(Var_HD) / sqrt(n - (D - 1) * tau - 1) * qnorm(1 - 0.05 / 2))
    SemiLengthC <- ifelse(!is.finite(Var_CD) | Var_CD <= 0, NA,
                          sqrt(Var_CD) / sqrt(n - (D - 1) * tau - 1) * qnorm(1 - 0.05 / 2))
    
    results[[r]] <- data.frame(
      Model = model$type,
      n = n,
      D = D,
      tau = tau,
      Rep = r,
      H_Shannon = Hs,
      C_Shannon = Cs,
      Var_H = Var_HD,
      Var_C = Var_CD,
      SemiLength_H = SemiLengthH,
      SemiLength_C = SemiLengthC
    )
    
    ts_store[[r]] <- data.frame(
      Model = model$type,
      n = n,
      D = D,
      tau = tau,
      Rep = r,
      Value = as.numeric(ts_data)
    )
  }
  
  list(
    summary = bind_rows(results),
    timeseries = bind_rows(ts_store)
  )
}

# --- Run simulation ---
message("Simulating ARMA22_M2 model with n = ", N, ", D = ", D, ", tau = ", tau)
res <- generate_model_data(model, n = N, D = D, tau = tau, R = R)

# --- Output paths ---
summary_path <- "Entropy_Complexity_ARMA22_M2.xlsx"
ts_path <- "TimeSeries_Data_ARMA22_M2.xlsx"

# --- Save Excel files ---
write_xlsx(list(ARMA22_M2 = res$summary), path = summary_path)
write_xlsx(list(ARMA22_M2 = res$timeseries), path = ts_path)

hc_plot <- ggplot(res$summary, aes(x = H_Shannon, y = C_Shannon)) +
  geom_point(alpha = 0.6, color = "blue", size = 2) +
  labs(title = "Entropy-Complexity Plane (ARMA22_M2)",
       subtitle = paste0("n = ", N, ", D = ", D, ", tau = ", tau, ", R = ", R),
       x = "H",
       y = "C") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(hc_plot)


cor_HC <- cor(res$summary$H_Shannon, res$summary$C_Shannon)
cat("Correlation between H and C:", round(cor_HC, 4), "\n")

# 2. Correlation test with p-value
cor_test <- cor.test(res$summary$H_Shannon, res$summary$C_Shannon)

# --- End of file ---
#-----------------------------------------------------------------------------------
# I use the Op sequences function as I want to change the lag parameter (tau)

library(StatOrdPattHxC)
library(dplyr)
library(writexl)
library(ggplot2)

# --- Parameters ---
set.seed(1234567890, kind = "Mersenne-Twister")
D <- 6                     # Embedding dimension
tau <- 3                   # Delay parameter
N <- 1000                   # Sample size
R <- 50                   # Number of replications

# --- Single Model Definition ---
model <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), type = "ARMA22_M2")

# --- Function to simulate data and compute entropy/complexity ---
generate_model_data <- function(model, n, D, tau, R) {
  results <- list()
  ts_store <- list()
  
  for (r in 1:R) {
    ts_data <- arima.sim(model = model, n = n)
    
    # Get ordinal pattern sequences
    OP_sequences <- OPseq(ts_data, emb = D, lag = tau)
    
    # Calculate probabilities from sequences
    # Count frequency of each ordinal pattern
    pattern_counts <- table(OP_sequences)
    total_patterns <- length(OP_sequences)
    ProbTS <- as.numeric(pattern_counts / total_patterns)
    
    # Entropy and complexity
    Hs <- HShannon(ProbTS)
    Cs <- StatComplexity(ProbTS)
    
    # Variances (keep negative)
    Var_HD <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(ProbTS, n - (D - 1) * tau))
    a_ratio <- Var_HD / Var_HI
    Var_CI <- suppressWarnings(varC(ProbTS, n - (D - 1) * tau))
    Var_CD <- a_ratio * Var_CI
    
    # Semi-lengths: NA if variance <= 0
    SemiLengthH <- ifelse(!is.finite(Var_HD) | Var_HD <= 0, NA,
                          sqrt(Var_HD) / sqrt(n - (D - 1) * tau) * qnorm(1 - 0.05 / 2))
    SemiLengthC <- ifelse(!is.finite(Var_CD) | Var_CD <= 0, NA,
                          sqrt(Var_CD) / sqrt(n - (D - 1) * tau) * qnorm(1 - 0.05 / 2))
    
    results[[r]] <- data.frame(
      Model = model$type,
      n = n,
      D = D,
      tau = tau,
      Rep = r,
      H_Shannon = Hs,
      C_Shannon = Cs,
      Var_H = Var_HD,
      Var_C = Var_CD,
      SemiLength_H = SemiLengthH,
      SemiLength_C = SemiLengthC
    )
    
    ts_store[[r]] <- data.frame(
      Model = model$type,
      n = n,
      D = D,
      tau = tau,
      Rep = r,
      Value = as.numeric(ts_data)
    )
  }
  
  list(
    summary = bind_rows(results),
    timeseries = bind_rows(ts_store)
  )
}

# --- Run simulation ---
message("Simulating ARMA22_M2 model with n = ", N, ", D = ", D, ", tau = ", tau)
res <- generate_model_data(model, n = N, D = D, tau = tau, R = R)

# --- Output paths ---
summary_path <- "Entropy_Complexity_ARMA22_M2.xlsx"
ts_path <- "TimeSeries_Data_ARMA22_M2.xlsx"

# --- Save Excel files ---
write_xlsx(list(ARMA22_M2 = res$summary), path = summary_path)
write_xlsx(list(ARMA22_M2 = res$timeseries), path = ts_path)

# --- HC Point Plot ---
library(ggplot2)

# Create HC plot
hc_plot <- ggplot(res$summary, aes(x = H_Shannon, y = C_Shannon)) +
  geom_point(alpha = 0.6, color = "blue", size = 2) +
  #geom_errorbar(aes(ymin = C_Shannon - SemiLength_C, 
  #                  ymax = C_Shannon + SemiLength_C), 
  #              alpha = 0.3, width = 0.01) +
 # geom_errorbarh(aes(xmin = H_Shannon - SemiLength_H, 
 #                    xmax = H_Shannon + SemiLength_H), 
 #                alpha = 0.3, height = 0.01) +
  labs(title = "Entropy-Complexity Plane (ARMA22_M2)",
       subtitle = paste0("n = ", N, ", D = ", D, ", tau = ", tau, ", R = ", R),
       x = "Shannon Entropy (H)",
       y = "Statistical Complexity (C)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(hc_plot)
# Save plot
#ggsave("HC_Plot_ARMA22_M2.png", plot = hc_plot, width = 8, height = 6, dpi = 300)

# --- Collinearity Check ---
cat("\n=== Collinearity Analysis ===\n")

# 1. Correlation between H and C
cor_HC <- cor(res$summary$H_Shannon, res$summary$C_Shannon)
cat("Correlation between H and C:", round(cor_HC, 4), "\n")

# 2. Correlation test with p-value
cor_test <- cor.test(res$summary$H_Shannon, res$summary$C_Shannon)
cat("Correlation test p-value:", format.pval(cor_test$p.value), "\n")
cat("95% Confidence interval:", round(cor_test$conf.int, 4), "\n\n")

# 3. Collinearity interpretation
if (abs(cor_HC) > 0.9) {
  cat("⚠️  Very high collinearity detected (|r| > 0.9)\n")
} else if (abs(cor_HC) > 0.7) {
  cat("⚠️  High collinearity detected (|r| > 0.7)\n")
} else if (abs(cor_HC) > 0.5) {
  cat("ℹ️  Moderate correlation detected (|r| > 0.5)\n")
} else {
  cat("✅ Low collinearity (|r| < 0.5)\n")
}

##-----------------------------------------------------------------------------
#ARMA-GARCH model------------------------------------------------------------
library(StatOrdPattHxC)
library(dplyr)
library(writexl)
library(ggplot2)
library(rugarch)   # <-- ARMA–GARCH simulation [web:41]

# --- Parameters ---
set.seed(1234567890, kind = "Mersenne-Twister")
D   <- 5                     # Embedding dimension
tau <- 1                     # Delay parameter
N   <- 100                  # Sample size
R   <- 50                    # Number of replications

# --- ARMA–GARCH specification (ARMA(2,2) + sGARCH(1,1)) ---
spec_arma22_garch <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder    = c(2, 2),
    include.mean = TRUE
  ),
  distribution.model = "norm"
)

# (Optional) fix parameters for reproducible dynamics
setfixed(spec_arma22_garch) <- list(
  mu     = 0.0,
  ar1    = -0.8,
  ar2    =  0.15,
  ma1    = -0.8,
  ma2    =  0.15,
  omega  =  0.1,
  alpha1 =  0.15,
  beta1  =  0.8
)

model <- list(
  spec = spec_arma22_garch,
  type = "ARMA22_M2_GARCH"
)

# --- Function to simulate data and compute entropy/complexity ---
generate_model_data <- function(model, n, D, tau, R) {
  results  <- list()
  ts_store <- list()
  
  for (r in 1:R) {
    # --- ARMA–GARCH simulation ---
    sim <- ugarchpath(
      model$spec,
      n.sim   = n,
      n.start = 1000,   # burn-in to reach stationarity [web:41][web:61]
      m.sim   = 1
    )
    
    # Use conditional mean as ARMA–GARCH series (returns)
    ts_data <- as.numeric(fitted(sim)[, 1])
    
    # --- Ordinal patterns ---
    OP_sequences <- OPseq(ts_data, emb = D, lag = tau)
    
    # Probabilities from sequences
    pattern_counts <- table(OP_sequences)
    total_patterns <- length(OP_sequences)
    ProbTS <- as.numeric(pattern_counts / total_patterns)
    
    # Entropy and complexity
    Hs <- HShannon(ProbTS)
    Cs <- StatComplexity(ProbTS)
    
    # Variances
    Var_HD <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HI <- suppressWarnings(
      asymptoticVarHShannonMultinomial(ProbTS, n - (D - 1) * tau)
    )
    a_ratio <- Var_HD / Var_HI
    Var_CI  <- suppressWarnings(varC(ProbTS, n - (D - 1) * tau))
    Var_CD  <- a_ratio * Var_CI
    
    # Semi-lengths
    SemiLengthH <- ifelse(
      !is.finite(Var_HD) | Var_HD <= 0, NA,
      sqrt(Var_HD) / sqrt(n - (D - 1) * tau) * qnorm(1 - 0.05 / 2)
    )
    SemiLengthC <- ifelse(
      !is.finite(Var_CD) | Var_CD <= 0, NA,
      sqrt(Var_CD) / sqrt(n - (D - 1) * tau) * qnorm(1 - 0.05 / 2)
    )
    
    results[[r]] <- data.frame(
      Model         = model$type,
      n             = n,
      D             = D,
      tau           = tau,
      Rep           = r,
      H_Shannon     = Hs,
      C_Shannon     = Cs,
      Var_H         = Var_HD,
      Var_C         = Var_CD,
      SemiLength_H  = SemiLengthH,
      SemiLength_C  = SemiLengthC
    )
    
    ts_store[[r]] <- data.frame(
      Model = model$type,
      n     = n,
      D     = D,
      tau   = tau,
      Rep   = r,
      Value = ts_data
    )
  }
  
  list(
    summary    = bind_rows(results),
    timeseries = bind_rows(ts_store)
  )
}

# --- Run simulation ---
message("Simulating ", model$type, " with n = ", N, ", D = ", D, ", tau = ", tau)
res <- generate_model_data(model, n = N, D = D, tau = tau, R = R)

# --- Output paths ---
summary_path <- "Entropy_Complexity_ARMA22_M2_GARCH.xlsx"
ts_path      <- "TimeSeries_Data_ARMA22_M2_GARCH.xlsx"

# --- Save Excel files ---
write_xlsx(list(ARMA22_M2_GARCH = res$summary),    path = summary_path)
write_xlsx(list(ARMA22_M2_GARCH = res$timeseries), path = ts_path)

# --- HC Point Plot ---
hc_plot <- ggplot(res$summary, aes(x = H_Shannon, y = C_Shannon)) +
  geom_point(alpha = 0.6, color = "blue", size = 2) +
  labs(
    title    = "Entropy-Complexity Plane (ARMA22_M2 + GARCH)",
    subtitle = paste0("n = ", N, ", D = ", D, ", tau = ", tau, ", R = ", R),
    x        = "Shannon Entropy (H)",
    y        = "Statistical Complexity (C)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )
print(hc_plot)

# --- Collinearity Check ---
cat("\n=== Collinearity Analysis ===\n")
cor_HC   <- cor(res$summary$H_Shannon, res$summary$C_Shannon)
cat("Correlation between H and C:", round(cor_HC, 4), "\n")
cor_test <- cor.test(res$summary$H_Shannon, res$summary$C_Shannon)
cat("Correlation test p-value:", format.pval(cor_test$p.value), "\n")

