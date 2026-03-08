library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

# ==============================================================================
# PARAMETERS
# ==============================================================================
D <- 4  # Embedding dimension
n <- 1000
alpha_1 <- 0.1  
alpha_2 <- 0.5
alpha_3 <- 0.9
burnin <- 500 
ARMA22 = list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), name = "ARMA(2,2)")

# ==============================================================================
generate_logistic_map <- function(n, r_log =3.8) {
  x <- numeric(n + burnin)
  x[1] <- 0.5  # Initial condition
  
  for (i in 2:(n + burnin)) {
    x[i] <- r_log * x[i-1] * (1 - x[i-1])
  }
  
  # Remove burnin period
  x <- x[(burnin + 1):(n + burnin)]
  return(x)
}

# ==============================================================================
# PURE ARMA MODEL
# ==============================================================================
generate_arma <- function(n, ar, ma) {
  arma_model <- list(ar = ar, ma = ma)
  arma_ts <- arima.sim(model = arma_model, n = n)
  return(as.numeric(arma_ts))
}

Hybrid_model_1 <- alpha_1 * generate_logistic_map(n) + (1 - alpha_1) * generate_arma(n, ar = ARMA22$ar, ma = ARMA22$ma)
Hybrid_model_2 <- alpha_2 * generate_logistic_map(n) + (1 - alpha_2) * generate_arma(n, ar = ARMA22$ar, ma = ARMA22$ma)
Hybrid_model_3 <- alpha_3 * generate_logistic_map(n) + (1 - alpha_3) * generate_arma(n, ar = ARMA22$ar, ma = ARMA22$ma)
plot(Hybrid_model_1, type = "l", main = "Hybrid Model Time Series_1", xlab = "Time", ylab = "Value")
plot(Hybrid_model_2, type = "l", main = "Hybrid Model Time Series_2", xlab = "Time", ylab = "Value")
plot(Hybrid_model_3, type = "l", main = "Hybrid Model Time Series_3", xlab = "Time", ylab = "Value")

ProbTS_1  <- OPprob(Hybrid_model_1, emb = D)
ProbTS_2  <- OPprob(Hybrid_model_2, emb = D)
ProbTS_3  <- OPprob(Hybrid_model_3, emb = D)

Hs_1 <- HShannon(ProbTS_1)
Hs_2 <- HShannon(ProbTS_2)
Hs_3 <- HShannon(ProbTS_3)



C_Shannon_1 <- StatComplexity(ProbTS_1)
C_Shannon_2 <- StatComplexity(ProbTS_2)
C_Shannon_3 <- StatComplexity(ProbTS_3)



data("LinfLsup")

# Full boundary curve
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# Separate upper and lower boundaries
boundary_lower <- LinfLsup_full %>% filter(Side == "Lower")
boundary_upper <- LinfLsup_full %>% filter(Side == "Upper")


boundary <- ggplot() + geom_line(data = boundary_lower, aes(x = H, y = C), 
            color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
  geom_line(data = boundary_upper, aes(x = H, y = C), 
            color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8)

a <- boundary +
  geom_point(data = data.frame(Hs_1, C_Shannon_1),
             aes(x = Hs_1, y = C_Shannon_1),
             color = "blue", size = 3, alpha = 0.7) +
  geom_point(data = data.frame(Hs_2, C_Shannon_2),
             aes(x = Hs_2, y = C_Shannon_2),
             color = "red", size = 3, alpha = 0.7) +
  geom_point(data = data.frame(Hs_3, C_Shannon_3),
             aes(x = Hs_3, y = C_Shannon_3),
             color = "green", size = 3, alpha = 0.7) +
  labs(title = "Complexity-Entropy Plane for Hybrid Model",
       x = "H",
       y = "C") +
  theme_minimal()

print(a)

#---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(StatOrdPattHxC)

D <- 4
n <- 1000
alpha <- 0.5
r <- 3.8

data("LinfLsup")

set.seed(1234567890, kind = "Mersenne-Twister")

arma_list <- list(ar = c(-0.8,0.1), ma = c(-0.8,0.1))

arma <- function(n){
  ts_data <- arima.sim(model = arma_list, n) 
  return(as.numeric(ts_data))
}

logistic_model <- function(r,n){
  x <- numeric(n)
  x[1] <- 0.5
  for (i in 2:n) {
   x[i] <- r * x[i-1]*(1-x[i-1]) 
  }
  return(x)
}

hybrid_model <- alpha * logistic_model(r,n) + (1-alpha) * arma(n)

plot(hybrid_model, type = "l")

op <- OPprob(hybrid_model, D)
Entropy <- HShannon(op)
Complexity <- StatComplexity(op)

varD <- sigma2q(hybrid_model, D)
varM <- asymptoticVarHShannonMultinomial(op, n-D)

a <- varD/varM

var_C <- varC(op, n-D)

Var_C_D <- a * var_C

semi_length_H <- sqrt(varD)/sqrt(n-D)*qnorm(1-0.05)/2
semi_length_C <- sqrt(Var_C_D)/sqrt(n-D)*qnorm(1-0.05)/2

upper_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Upper")
lower_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Lower")

p <- ggplot()+
  geom_line(data = upper_bound, aes(x= H, y= C), color="black")+
  geom_line(data = lower_bound, aes(x=H, y=C), color = "black")+
  geom_point(data = data.frame(Entropy, Complexity), aes(x= Entropy, y= Complexity), color = "blue")+
  geom_errorbarh(aes(xmin = Entropy - semi_length_H, xmax = Entropy + semi_length_H, y= Complexity))+
  geom_errorbar(aes(x= Entropy, ymin = Complexity - semi_length_C, ymax = Complexity+semi_length_C))+
  theme_bw(base_size = 11, base_family = "serif")+
  labs(x= expression(italic(H)),
       y= expression(italic(C)),
       title = "HC plane")
print(p)







