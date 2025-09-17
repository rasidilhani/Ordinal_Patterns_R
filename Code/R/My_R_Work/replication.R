library(StatOrdPattHxC)
library(dplyr)
library(ggplot2)

D <- 3
n <- 1000

ts_data <- arima.sim(model = list(ar1 = 0.8, ar2 = 0.5, ma1 = 0.8, ma2 = 0.5), n = n)
ProbARMA <- OPprob(ts_data, emb = D)
Entropy <- HShannon(ProbARMA)
Complexity <- StatComplexity(ProbARMA)


#####################################################################################################
library(StatOrdPattHxC)
D <- 3
n <- 100

# Function for a single run
arma_entropy_complexity <- function() {
  ts_data <- arima.sim(model = list(ar1 = 0.8, ar2 = 0.5, ma1 = 0.8, ma2 = 0.5), n = n)
  ProbARMA <- OPprob(ts_data, emb = D)
  H <- HShannon(ProbARMA)
  C <- StatComplexity(ProbARMA)
  c(Entropy = H, Complexity = C)
}

set.seed(123)  # for reproducibility!
results_matrix <- replicate(300, arma_entropy_complexity())

# Convert to a data frame for easy use (transpose if needed)
results_df <- data.frame(t(results_matrix))
head(results_df)

