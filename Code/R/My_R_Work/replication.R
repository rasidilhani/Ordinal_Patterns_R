library(StatOrdPattHxC)
library(dplyr)
library(ggplot2)
D <- 3
n <- 500

ts_data <- arima.sim(model = list(ar1 = 0.8, ar2 = 0.5, ma1 = 0.8, ma2 = 0.5), n = n)
ProbARMA <- OPprob(ts_data, emb = D)
Entropy <- HShannon(ProbARMA)
Complexity <- StatComplexity(ProbARMA)
