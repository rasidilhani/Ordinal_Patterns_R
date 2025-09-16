# Load required libraries
library(StatOrdPattHxC)
library(dplyr)
library(ggplot2)

# 1. Create grids of AR and MA values
ar1_vals <- c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8)
ma1_vals <- c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8)
param_grid <- expand.grid(ar = ar1_vals, ma = ma1_vals, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
param_grid$label <- with(param_grid, paste0("AR=", ar, ", MA=", ma))

# 2. Generate and save standard normal sample data
set.seed(1234567890, kind="Mersenne-Twister")
n <- 5000
sample_data <- rnorm(n)
save(sample_data, file = "sample_data_Renyi.RData")  

# 3. Simulation and feature calculation for ARMA combinations
D <- 5  
results <- data.frame()  

for (i in seq_len(nrow(param_grid))) {
  ar <- param_grid$ar[i]
  ma <- param_grid$ma[i]
  label <- param_grid$label[i]
  
  ts_data <- arima.sim(model = list(ar = ar, ma = ma), n = n)
  ProbARMA <- OPprob(ts_data, emb = D)
  REntropy <- HRenyi(ProbARMA, beta = 1.5)
  Complexity <- StatComplexity(ProbARMA)
  variancePD <- sigma2q(ts_data, emb = D)
  VarianceM <- asymptoticVarHShannonMultinomial(ProbARMA, n - D + 1)
  #VarianceM <- max(VarianceM, 1e-8)
  varCM <- varC(ProbARMA, n - D + 1)
  a <- variancePD / VarianceM
  VariancePDC <- a * varCM
  semiLengthPD <- sqrt(variancePD) / sqrt(n - D + 1) * qnorm(0.975)
  semiLengthPDC <- sqrt(VariancePDC) / sqrt(n - D + 1)
  
  results <- rbind(
    results,
    data.frame(
      ParameterSet = label,
      REntropy,
      Complexity,
      semiLengthPD,
      semiLengthPDC
    )
  )
}

save(results, file = "arma_results_Renyi.RData")  # Save results for flexible plotting

# 4. Plot results (can reload the results for future plots without re-running simulation)
load("arma_results_Renyi.RData") # Use if you've already run the simulation and saved results

#data("LinfLsup")  

ggplot(results, aes(x = REntropy, y = Complexity, color = ParameterSet)) +
  geom_point(size = 2) +
  #geom_errorbarh(aes(xmin = REntropy - semiLengthPD, xmax = REntropy + semiLengthPD), height = 0.02) +
  #geom_errorbar(aes(ymin = Complexity - semiLengthPDC, ymax = Complexity + semiLengthPDC), width = 0.02) +
  #geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
           # aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE) +
 # geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
  #          aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0.7, 1.0), ylim = c(0, 0.4)) +
  labs(
    title = paste("Renyi Entropy and Complexity for ARMA Time Series (D =", D, ", n =", n, ")"),
    x = expression(italic(H)),
    y = expression(italic(C)),
    color = "Parameter Set"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None")
