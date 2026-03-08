library(dplyr)
library(ggplot2)
library(StatOrdPattHxC)

D <- 4
n <- 1000
alpha_values <- c(0.1,0.5,0.8)
r1 <- 3.8
r2 <- 3.6
f <- 0.04

arma_list <- list(ar = c(-0.8,0.1), ma = c(-0.8,0.1))

arma_model <- function(n){
  ts_data <- arima.sim(model = arma_list, n)
  return(as.numeric(ts_data))
}


logistic_model <- function(r,n){
  x <- numeric(n)
  x[1] <- 0.5
  for (i in 2:n){
    x[i] <- r*x[i-1]*(1-x[i-1])
  }
  return(x)
}

sine_wave <- function(n,f){
  t <- 1:n
  x <- sin(2*pi*f*t)
  return(as.numeric(x))
}

results_1 <- data.frame()
for (i in seq_along(alpha_values)){
  alpha <- alpha_values[i]
  hybrid_arma <- alpha*logistic_model(r1,n) + (1-alpha)*arma_model(n)
  
  op <- OPprob(hybrid_arma, D)
  Entropy <- HShannon(op)
  Complexity <- StatComplexity(op)
  
  H_variance_D <- sigma2q(hybrid_arma, D)
  H_variance_M <- asymptoticVarHShannonMultinomial(op, n-D+1)
  
  C_variance_M <- varC(op, n-D+1)
  ratio <- H_variance_D/H_variance_M
  
  C_variance_D <- ratio*C_variance_M
  
  semi_H <- sqrt(H_variance_D/(n-D+1))*qnorm(1-0.05)/2
  semi_C <- sqrt(C_variance_D/(n-D+1))*qnorm(1-0.05)/2
  
  results_1 <-rbind(results_1, data.frame(alpha,
                                          Entropy,
                                          Complexity,
                                          semi_H,
                                          semi_C))
}


results_2 <- data.frame()

for (i in seq_along(alpha_values)){
  alpha <- alpha_values[i]
  hybrid_sine <- alpha*logistic_model(r2,n) + (1-alpha)*sine_wave(n,f)
  op <- OPprob(hybrid_sine, D)
  Entropy <- HShannon(op)
  Complexity <- StatComplexity(op)
  
  H_variance_D <- sigma2q(hybrid_sine, D)
  H_variance_M <- asymptoticVarHShannonMultinomial(op, n-D+1)
  
  C_variance_M <- varC(op, n-D+1)
  ratio <- H_variance_D/H_variance_M
  
  C_variance_D <- ratio*C_variance_M
  
  semi_H <- sqrt(H_variance_D/(n-D+1))*qnorm(1-0.05)/2
  semi_C <- sqrt(C_variance_D/(n-D+1))*qnorm(1-0.05)/2
  
  results_2 <-rbind(results_2, data.frame(alpha,
                                          Entropy,
                                          Complexity,
                                          semi_H,
                                          semi_C))
}

data("LinfLsup")

lower_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Lower")
upper_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Upper")

results_1$model <- "Arma_Hybrid"
results_2$model <- "Sine_Hybrid"

p <- ggplot() +
  geom_line(data = lower_bound, aes(x = H, y = C), color = "grey", linetype = "dashed")+
  geom_line(data = upper_bound, aes(x = H, y = C), color = "grey", linetype = "dashed")+
  geom_point(data = results_1, aes(x = Entropy, y = Complexity, color = factor(alpha), shape = model))+
  geom_point(data = results_2, aes(x = Entropy, y = Complexity, color = factor(alpha), shape = model))+
  geom_errorbarh(data = results_1, aes(xmin = Entropy-semi_H, xmax = Entropy+semi_H, y = Complexity))+
  geom_errorbar(data = results_1, aes(x = Entropy, ymin = Complexity-semi_C, ymax = Complexity+semi_C))+
  geom_errorbarh(data = results_2, aes(xmin = Entropy-semi_H, xmax = Entropy+semi_H, y = Complexity))+
  geom_errorbar(data = results_2, aes(x = Entropy, ymin = Complexity-semi_C, ymax = Complexity+semi_C))+
  theme_bw(base_size = 11, base_family = "serif")+
  labs(x = expression(italic(H)),
       y = expression(italic(C)),
       title = "HC Plane")

print(p)
