library(StatOrdPattHxC)
library(writexl)
library(readxl)
library(ggplot2)

# Logistic map function (chaotic time series)
logistic_map <- function(r, x0, n) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) {
    x[i] <- r * x[i-1] * (1 - x[i-1])
  }
  return(x)
}

set.seed(1234567890, kind = "Mersenne-Twister")
D  <- 5          # embedding dimension
tau <- 1         # delay parameter
n  <- 2000     # length of time series
R  <- 50         # number of replications

# Logistic map parameters (chaotic regime)
r_param <- 3.8   # chaos parameter (3.9 = fully chaotic)
x0      <- 0.1   # initial condition

# Storage
ts_list      <- vector("list", R)
results_list <- vector("list", R)

for (r in 1:R) {
  x0_r <- 0.1 + 1e-6 * r  # tiny different x0 for each rep
  ts_data <- logistic_map(r = r_param, x0 = x0_r, n = n)
  
  ts_list[[r]] <- data.frame(
    Replication = r,
    Time        = 1:n,
    X           = ts_data
  )
  
  ## 2A. Direct ordinal probability
  Prob_OPprob <- OPprob(ts_data, emb = D)
  
  ## 2B. OP sequences -> empirical probabilities
  OP_sequences <- OPseq(ts_data, emb = D, lag = tau)
  pattern_counts <- table(OP_sequences)
  total_patterns <- length(OP_sequences)
  Prob_OPseq <- as.numeric(pattern_counts / total_patterns)
  
  ## Use OPprob (change to Prob_OPseq if you prefer)
  Prob <- Prob_OPprob
  
  ## 3. Shannon entropy and statistical complexity
  H <- HShannon(Prob)
  C <- StatComplexity(Prob)
  
  ## 4. Asymptotic variances and related quantities
  var_HD <- sigma2q(ts_data, emb = D, ent = "S")
  var_HI <- asymptoticVarHShannonMultinomial(Prob, n - 2)
  Var_C  <- varC(Prob, n - 2)
  
  a      <- var_HD / var_HI
  Var_CD <- a * Var_C
  
  ## 5. Semi-lengths for 95% CIs
  z_alpha      <- qnorm(1 - 0.05 / 2)
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

# Bind results
ts_all      <- do.call(rbind, ts_list)
results_all <- do.call(rbind, results_list)

# Save to Excel
write_xlsx(ts_all,      path = "LogisticMap_ts_R50.xlsx")
write_xlsx(results_all, path = "LogisticMap_HC_results_R50.xlsx")

# Reload and plot
Logistic_results <- read_excel("LogisticMap_HC_results_R50.xlsx")

ggplot(Logistic_results, aes(x = H, y = C)) + 
  geom_point(alpha = 0.6, color = "blue", size = 2) + 
  theme_minimal() + 
  ggtitle("H-C Plane for Logistic Map (r = 3.9, R = 50)") +
  xlab("Shannon Entropy (H)") +
  ylab("Statistical Complexity (C)")

collinnearity <- cor(Logistic_results$H, Logistic_results$C)
print(paste("Correlation between H and C:", round(collinnearity, 4)))

