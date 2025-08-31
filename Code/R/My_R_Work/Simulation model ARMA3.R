library(StatOrdPattHxC)
library(dplyr)
library(tidyverse)
library(ggplot2)

# 1. SETUP: Define the AR and MA parameter options I want to test
ar_list <- list(
  numeric(0),          # AR(0)
  0.6,                 # AR(1)
  c(0.8, -0.4),        # AR(2) with positive and negative
  c(0.5, 0.3)          # AR(2) with two positives
)
ma_list <- list(
  numeric(0),          # MA(0)
  0.3,                 # MA(1)
  c(0.7, 0.2),         # MA(2)
  0.4                  # MA(1)
)

# 2. GENERATE ALL PARAMETER COMBINATIONS
param_combinations <- list()
labels <- c()
for (ar in ar_list) {
  for (ma in ma_list) {
    # Build label for legend and tracking
    ar_label <- if(length(ar) == 0) "0" else paste(ar, collapse = ",")
    ma_label <- if(length(ma) == 0) "0" else paste(ma, collapse = ",")
    label <- paste0("AR=", ar_label, ", MA=", ma_label)
    param_combinations[[label]] <- list(ar = ar, ma = ma)
    labels <- c(labels, label)
  }
}
labels <- unique(labels)  # Ensure all labels are unique

# 3. MAIN PARAMETERS FOR THE SIMULATION
D <- 5      # Embedding dimension for ordinal patterns
n <- 1000   # Length of time series
set.seed(1234567890)  # For reproducibility

# 4. SIMULATION AND FEATURE CALCULATION
results <- data.frame()
for (label in labels) {
  params <- param_combinations[[label]]
  ts_data <- arima.sim(model = list(ar = params$ar, ma = params$ma), n = n)
  ProbARMA <- OPprob(ts_data, emb = D)
  Entropy <- HShannon(ProbARMA)
  Complexity <- StatComplexity(ProbARMA)
  variancePD <- sigma2q(ts_data, emb = D)
  VarianceM <- asymptoticVarHShannonMultinomial(ProbARMA, n - D + 1)
  VarianceM <- max(VarianceM, 1e-8)  # Avoid division by zero
  varCM <- varC(ProbARMA, n - D + 1)
  a <- variancePD / VarianceM
  VariancePDC <- a * varCM
  semiLengthPD <- sqrt(variancePD) / sqrt(n - D + 1) * qnorm(0.975)
  semiLengthPDC <- sqrt(VariancePDC) / sqrt(n - D + 1)
  results <- rbind(
    results,
    data.frame(
      ParameterSet = factor(label, levels = labels),
      Entropy,
      Complexity,
      semiLengthPD,
      semiLengthPDC
    )
  )
}

# 5. PLOTTING
data("LinfLsup") # Theoretical boundary reference
ggplot(results, aes(x = Entropy, y = Complexity, color = ParameterSet)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = Entropy - semiLengthPD, xmax = Entropy + semiLengthPD), height = 0.02) +
  geom_errorbar(aes(ymin = Complexity - semiLengthPDC, ymax = Complexity + semiLengthPDC), width = 0.02) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0.25, 1.0), ylim = c(0, 0.6)) +
  labs(
    title = paste("Entropy and Complexity for ARMA Time Series (D =", D, ", n =", n, ")"),
    x = expression(italic(H)),
    y = expression(italic(C)),
    color = "Parameter Set"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")
