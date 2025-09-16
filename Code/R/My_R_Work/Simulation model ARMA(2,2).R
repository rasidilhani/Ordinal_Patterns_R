# Load required libraries
library(StatOrdPattHxC)
library(dplyr)
library(ggplot2)
library(writexl)

# 1. Create grids of AR and MA values
ar_vals1 <- c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8)
ar_vals2 <- c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8)
ma_vals1 <- c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8)
ma_vals2 <- c(-0.8, -0.5, -0.3, 0, 0.3, 0.5, 0.8)
#param_grid <- expand.grid(ar1 = ar_vals1, ma1 = ma_vals1, ma2 = ma_vals2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#param_grid$label <- with(param_grid, paste0("AR1=", ar1, "MA1=", ma1, ", MA2=", ma2))
param_grid <- expand.grid(ar1 = ar_vals1, ar2 = ar_vals2, ma1 = ma_vals1, ma2 = ma_vals2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
param_grid$label <- with(param_grid, paste0("AR1=", ar1, "AR2=", ar2, ", MA1=", ma1, ", MA2=", ma2))


# 2. Generate and save standard normal sample data
set.seed(1234567890, kind="Mersenne-Twister")
n <- 5000
sample_data <- rnorm(n)
save(sample_data, file = "sample_data.RData")  # Use load("sample_data.RData") later to reload

# 3. Simulation and feature calculation for ARMA combinations
D <- 5  
results <- data.frame()  # Store all results

for (i in seq_len(nrow(param_grid))) {
  ar1 <- param_grid$ar1[i]
  ar2 <- param_grid$ar2[i]
  ma1 <- param_grid$ma1[i]
  ma2 <- param_grid$ma2[i]
  label <- param_grid$label[i]
  
  if(ar1 + ar2 < 1 || ar2 - ar1 < 1) {
  print(label)
  ts_data <- arima.sim(model = list(ar1 = ar1, ar2 = ar2, ma1 = ma1, ma2 = ma2), n = n)
  ProbARMA <- OPprob(ts_data, emb = D)
  Entropy <- HShannon(ProbARMA)
  Complexity <- StatComplexity(ProbARMA)
  variancePD <- sigma2q(ts_data, emb = D)
  VarianceM <- asymptoticVarHShannonMultinomial(ProbARMA, n - D + 1)
  VarianceM <- max(VarianceM, 1e-8)
  varCM <- varC(ProbARMA, n - D + 1)
  a <- variancePD / VarianceM
  VariancePDC <- a * varCM
  semiLengthPD <- sqrt(variancePD) / sqrt(n - D + 1) * qnorm(0.975)
  semiLengthPDC <- sqrt(VariancePDC) / sqrt(n - D + 1)
  
  results <- rbind(
    results,
    data.frame(
      ParameterSet = label,
      Entropy,
      Complexity,
      semiLengthPD,
      semiLengthPDC
    )
  )
    }
  
}


save(results, file = "arma_grid_with_condition_results.RData")  # Save results for flexible plotting
write_xlsx(results, path = "arma(12)_grid_results.xlsx")

# 4. Plot results 
load("arma_grid_with_condition_results.RData") # Use if you've already run the simulation and saved results

data("LinfLsup")  

ggplot(results, aes(x = Entropy, y = Complexity, color = ParameterSet)) +
  geom_point(size = 2) +
  #geom_errorbarh(aes(xmin = Entropy - semiLengthPD, xmax = Entropy + semiLengthPD), height = 0.02) +
  #geom_errorbar(aes(ymin = Complexity - semiLengthPDC, ymax = Complexity + semiLengthPDC), width = 0.02) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0.7, 1.0), ylim = c(0, 0.4)) +
  labs(
    title = paste("Entropy and Complexity for ARMA(2,2) Time Series"),
    x = expression(italic(H)),
    y = expression(italic(C)),
    color = "Parameter Set"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None")
