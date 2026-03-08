library(tidyverse)
library(StatOrdPattHxC)

D <- 4
n <- 1000
r <- 3.8
f <- 0.04

set.seed(1234567890, kind = "Mersenne-Twister")
arma_list <- list(ar = c(-0.8,0.1), ma = c(-0.8,0.1))

arma <- function(n){
  ts_data <- arima.sim(model = arma_list, n)
  ts_norm <- (ts_data - min(ts_data))/(max(ts_data)-min(ts_data))
  return(as.numeric(ts_norm))
}

x <- arma(1000)

plot(x, type = "l")
 
op_arma <- OPprob(x, 4)
H_arma <- HShannon(op_arma)
C_arma <- StatComplexity(op_arma)

logistic <- function(r, n){
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) {
    a[i] <- r*a[i-1]*(1-a[i-1])
  }
  return(a)
}

y <- logistic(3.8,1000)

plot(y, type = "l")

op_logistic <- OPprob(y, 4)
H_logistic <- HShannon(op_logistic)
C_logistic <- StatComplexity(op_logistic)

sine <- function(n, f){
  t <- 1:n
  b <- sin(2*pi*f*t)
  return(as.numeric(b))
}
 z <- sine(1000, 0.04)

plot(z, type = "l")

op_sine <- OPprob(z, 4)
H_sine <- HShannon(op_sine)
C_sine <- StatComplexity(op_sine)

data("LinfLsup")

lower_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Lower")
upper_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Upper")

ArmaLogistic <- 0.1*x + 0.9*y
op_ArmaLogistic <- OPprob(ArmaLogistic, 4)
H_ArmaLogistic <- HShannon(op_ArmaLogistic)
C_ArmaLogistic <- StatComplexity(op_ArmaLogistic)

ArmaSine <- 0.1*x + 0.9* z
op_ArmaSine <- OPprob(ArmaSine, 4)
H_ArmaSine <- HShannon(op_ArmaSine)
C_ArmaSine <- StatComplexity(op_ArmaSine)

LogisticSine <- 0.1*y + 0.9*z
op_LogisticSine <- OPprob(LogisticSine, 4)
H_LogisticSine <- HShannon(op_LogisticSine)
C_LogisticSine <- StatComplexity(op_LogisticSine)

points_df <- data.frame(
  H = c(H_arma, H_logistic, H_sine, H_ArmaLogistic, H_ArmaSine, H_LogisticSine),
  C = c(C_arma, C_logistic, C_sine, C_ArmaLogistic, C_ArmaSine, C_LogisticSine),
  Model = c("ARMA(2,2)", "Logistic", "Sine Wave", 
            "ARMA+Logistic", "ARMA+Sine", "Logistic+Sine")
)

model_colors <- c(
  "ARMA(2,2)"     = "red",
  "Logistic"      = "blue",
  "Sine Wave"     = "green",
  "ARMA+Logistic" = "purple",
  "ARMA+Sine"     = "black",
  "Logistic+Sine" = "orange")

model_shapes <- c(
  "ARMA(2,2)"     = 16,
  "Logistic"      = 17,
  "Sine Wave"     = 15,
  "ARMA+Logistic" = 18,
  "ARMA+Sine"     = 8,
  "Logistic+Sine" = 10)

p <- ggplot()+
  geom_line(data = lower_bound, aes(x= H, y= C), color = "grey", linetype = "dashed")+
  geom_line(data = upper_bound, aes(x= H, y= C), color = "grey", linetype = "dashed")+
  geom_point(data = points_df, aes(x = H, y = C, color = Model, shape = Model), size = 2) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  theme_bw(base_size = 11, base_family = "serif")+
  labs(x = expression(italic(H)),
       y = expression(italic(C)),
       title = "HC plane",
       color = "Model",
       shape = "Model")

#p <- p + geom_point(aes(x = H_arma, y = C_arma), color = "red")+
#  geom_point(aes(x = H_logistic, y = C_logistic), color = "blue")+
#  geom_point(aes(x = H_sine, y = C_sine), color = "green")+
#  geom_point(aes(x = H_ArmaLogistic, y = C_ArmaLogistic), color = "purple")+
#  geom_point(aes(x = H_ArmaSine, y = C_ArmaSine), color = "black")+
#  geom_point(aes(x = H_LogisticSine, y = C_LogisticSine), color = "orange")

p

cor(points_df$H, points_df$C, method = "spearman")

#----------------------------------------------------------------
# End of the code
#------------------------------------------------------------

#--------------------------------------------------------------
# In this code, I improved the same model structure with 50 replications
#--------------------------------------------------------------
library(tidyverse)
library(StatOrdPattHxC)

D <- 4
n <- 1000
r <- 3.8
f <- 0.04
R <- 50

arma_list <- list(ar = c(-0.8,0.1), ma = c(-0.8,0.1))

arma <- function(n){
  ts_data <- arima.sim(model = arma_list, n)
  ts_norm <- (ts_data - min(ts_data))/(max(ts_data)-min(ts_data))
  return(as.numeric(ts_norm))
}

logistic <- function(r, n){
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) {
    a[i] <- r*a[i-1]*(1-a[i-1])
  }
  return(a)
}

sine <- function(n, f){
  t <- 1:n
  b <- sin(2*pi*f*t)
  return(as.numeric(b))
}

results_df <- data.frame()

for (i in 1:R) {
  x <- arma(n)
  y <- logistic(r, n)
  z <- sine(n, f)
  
  ArmaLogistic <- 0.1*x + 0.9*y
  ArmaSine <- 0.1*x + 0.9* z
  LogisticSine <- 0.1*y + 0.9*z
  
  H_arma   <- HShannon(OPprob(x, D))
  C_arma   <- StatComplexity(OPprob(x, D))
  
  H_logistic <- HShannon(OPprob(y, D))
  C_logistic <- StatComplexity(OPprob(y, D))
  
  H_sine <- HShannon(OPprob(z, D))
  C_sine <- StatComplexity(OPprob(z, D))
  
  H_ArmaLogistic <- HShannon(OPprob(ArmaLogistic, D))
  C_ArmaLogistic <- StatComplexity(OPprob(ArmaLogistic, D))
  
  H_ArmaSine <- HShannon(OPprob(ArmaSine, D))
  C_ArmaSine <- StatComplexity(OPprob(ArmaSine, D))
  
  H_LogisticSine <- HShannon(OPprob(LogisticSine, D))
  C_LogisticSine <- StatComplexity(OPprob(LogisticSine, D))
  
  points_df <- data.frame(
    H = c(H_arma, H_logistic, H_sine, H_ArmaLogistic, H_ArmaSine, H_LogisticSine),
    
    C = c(C_arma, C_logistic, C_sine, C_ArmaLogistic, C_ArmaSine, C_LogisticSine),
    
    Model = c("ARMA(2,2)", "Logistic", "Sine Wave", 
              "ARMA+Logistic", "ARMA+Sine", "Logistic+Sine")
  )
  
  results_df <- rbind(results_df, points_df)
  
}

model_colors <- c(
  "ARMA(2,2)"     = "red",
  "Logistic"      = "blue",
  "Sine Wave"     = "green",
  "ARMA+Logistic" = "purple",
  "ARMA+Sine"     = "black",
  "Logistic+Sine" = "orange")

model_shapes <- c(
  "ARMA(2,2)"     = 16,
  "Logistic"      = 17,
  "Sine Wave"     = 15,
  "ARMA+Logistic" = 18,
  "ARMA+Sine"     = 8,
  "Logistic+Sine" = 10)

data("LinfLsup")

lower_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Lower")
upper_bound <- LinfLsup%>% filter(Dimension == as.character(D), Side == "Upper")


p <- ggplot()+
  geom_line(data = lower_bound, aes(x= H, y= C), color = "grey", linetype = "dashed")+
  geom_line(data = upper_bound, aes(x= H, y= C), color = "grey", linetype = "dashed")+
  geom_point(data = results_df, aes(x = H, y = C, color = Model, shape = Model), size = 2, alpha = 0.6) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  theme_bw(base_size = 11, base_family = "serif")+
  labs(x = expression(italic(H)),
       y = expression(italic(C)),
       title = paste0("HC Plane (50 replications), D = ", D),
       color = "Model",
       shape = "Model")

print(p)

#-----------------------------------------------------
# End of the code
#----------------------------------------------------

#----------------------------------------------------------
# In this code, I have added AR(2) and MA(2) models and combined them with the previous models. 
# The purpose of this code is to include these additional cases in order to extend my study design.
#----------------------------------------------------------------
library(tidyverse)
library(StatOrdPattHxC)
library(here)

D <- 4
n <- 1000
r <- 3.8
f <- 0.04

arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_list  <- list(ar = NULL,          ma = c(-0.8, 0.1))

arma <- function(n) {
  ts_data <- arima.sim(model = arma_list, n)
  ts_norm <- (ts_data - min(ts_data)) / (max(ts_data) - min(ts_data))
  return(as.numeric(ts_norm))
}

ar2 <- function(n) {
  ts_data <- arima.sim(model = ar2_list, n)
  ts_norm <- (ts_data - min(ts_data)) / (max(ts_data) - min(ts_data))
  return(as.numeric(ts_norm))
}

ma2 <- function(n) {
  ts_data <- arima.sim(model = ma2_list, n)
  ts_norm <- (ts_data - min(ts_data)) / (max(ts_data) - min(ts_data))
  return(as.numeric(ts_norm))
}

logistic <- function(r, n) {
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) {
    a[i] <- r * a[i-1] * (1 - a[i-1])
  }
  return(a)
}

sine <- function(n, f) {
  t <- 1:n
  b <- sin(2 * pi * f * t)
  return(as.numeric(b))
}

x    <- arma(n)
x_ar <- ar2(n)
x_ma <- ma2(n)
y    <- logistic(r, n)
z    <- sine(n, f)

ArmaLogistic <- 0.1 * x    + 0.9 * y
ArmaSine     <- 0.1 * x    + 0.9 * z
AR2Logistic  <- 0.1 * x_ar + 0.9 * y
AR2Sine      <- 0.1 * x_ar + 0.9 * z
MA2Logistic  <- 0.1 * x_ma + 0.9 * y
MA2Sine      <- 0.1 * x_ma + 0.9 * z
LogisticSine <- 0.1 * y    + 0.9 * z

H_arma         <- HShannon(OPprob(x, D))
C_arma         <- StatComplexity(OPprob(x, D))

H_ar2          <- HShannon(OPprob(x_ar, D))
C_ar2          <- StatComplexity(OPprob(x_ar, D))

H_ma2          <- HShannon(OPprob(x_ma, D))
C_ma2          <- StatComplexity(OPprob(x_ma, D))

H_logistic     <- HShannon(OPprob(y, D))
C_logistic     <- StatComplexity(OPprob(y, D))

H_sine         <- HShannon(OPprob(z, D))
C_sine         <- StatComplexity(OPprob(z, D))

H_ArmaLogistic <- HShannon(OPprob(ArmaLogistic, D))
C_ArmaLogistic <- StatComplexity(OPprob(ArmaLogistic, D))

H_ArmaSine     <- HShannon(OPprob(ArmaSine, D))
C_ArmaSine     <- StatComplexity(OPprob(ArmaSine, D))

H_AR2Logistic  <- HShannon(OPprob(AR2Logistic, D))
C_AR2Logistic  <- StatComplexity(OPprob(AR2Logistic, D))

H_AR2Sine      <- HShannon(OPprob(AR2Sine, D))
C_AR2Sine      <- StatComplexity(OPprob(AR2Sine, D))

H_MA2Logistic  <- HShannon(OPprob(MA2Logistic, D))
C_MA2Logistic  <- StatComplexity(OPprob(MA2Logistic, D))

H_MA2Sine      <- HShannon(OPprob(MA2Sine, D))
C_MA2Sine      <- StatComplexity(OPprob(MA2Sine, D))

H_LogisticSine <- HShannon(OPprob(LogisticSine, D))
C_LogisticSine <- StatComplexity(OPprob(LogisticSine, D))

results_df <- data.frame(
  H = c(H_arma, H_ar2, H_ma2, H_logistic, H_sine,
        H_ArmaLogistic, H_ArmaSine,
        H_AR2Logistic,  H_AR2Sine,
        H_MA2Logistic,  H_MA2Sine,
        H_LogisticSine),
  C = c(C_arma, C_ar2, C_ma2, C_logistic, C_sine,
        C_ArmaLogistic, C_ArmaSine,
        C_AR2Logistic,  C_AR2Sine,
        C_MA2Logistic,  C_MA2Sine,
        C_LogisticSine),
  Model = c("ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
            "ARMA+Logistic", "ARMA+Sine",
            "AR2+Logistic",  "AR2+Sine",
            "MA2+Logistic",  "MA2+Sine",
            "Logistic+Sine")
)

#write_xlsx(results_df, "HC_results.xlsx")

model_colors <- c(
  "ARMA(2,2)"     = "red",
  "AR(2)"         = "darkred",
  "MA(2)"         = "tomato",
  "Logistic"      = "blue",
  "Sine"          = "green",
  "ARMA+Logistic" = "purple",
  "ARMA+Sine"     = "black",
  "AR2+Logistic"  = "darkorange",
  "AR2+Sine"      = "brown",
  "MA2+Logistic"  = "deeppink",
  "MA2+Sine"      = "cyan4",
  "Logistic+Sine" = "orange"
)

model_shapes <- c(
  "ARMA(2,2)"     = 16,
  "AR(2)"         = 17,
  "MA(2)"         = 15,
  "Logistic"      = 18,
  "Sine"          = 8,
  "ARMA+Logistic" = 10,
  "ARMA+Sine"     = 11,
  "AR2+Logistic"  = 12,
  "AR2+Sine"      = 13,
  "MA2+Logistic"  = 14,
  "MA2+Sine"      = 3,
  "Logistic+Sine" = 4
)

data("LinfLsup")
lower_bound <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Lower")
upper_bound <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Upper")

p <- ggplot() +
  geom_line(data = lower_bound, aes(x = H, y = C), color = "grey", linetype = "dashed") +
  geom_line(data = upper_bound, aes(x = H, y = C), color = "grey", linetype = "dashed") +
  geom_point(data = results_df, aes(x = H, y = C, color = Model, shape = Model), size = 2, alpha = 0.8) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  theme_bw(base_size = 11, base_family = "serif") +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane, D = ", D),
    color = "Model",
    shape = "Model"
  )
print(p)

# End of the code --------------------------------------


#------------------------------------------------------------
# With the 50 replications of the above model.
#-----------------------------------------------------------------------

library(tidyverse)
library(StatOrdPattHxC)
library(writexl)

D <- 4
n <- 1000
r <- 3.8
f <- 0.04
R <- 50

arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_list  <- list(ar = NULL,          ma = c(-0.8, 0.1))

arma <- function(n) {
  ts_data <- arima.sim(model = arma_list, n)
  ts_norm <- (ts_data - min(ts_data)) / (max(ts_data) - min(ts_data))
  return(as.numeric(ts_norm))
}

ar2 <- function(n) {
  ts_data <- arima.sim(model = ar2_list, n)
  ts_norm <- (ts_data - min(ts_data)) / (max(ts_data) - min(ts_data))
  return(as.numeric(ts_norm))
}

ma2 <- function(n) {
  ts_data <- arima.sim(model = ma2_list, n)
  ts_norm <- (ts_data - min(ts_data)) / (max(ts_data) - min(ts_data))
  return(as.numeric(ts_norm))
}

logistic <- function(r, n) {
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) {
    a[i] <- r * a[i-1] * (1 - a[i-1])
  }
  return(a)
}

sine <- function(n, f) {
  t <- 1:n
  b <- sin(2 * pi * f * t)
  return(as.numeric(b))
}

results_df <- data.frame()

for (i in 1:R) {
  
  x    <- arma(n)
  x_ar <- ar2(n)
  x_ma <- ma2(n)
  y    <- logistic(r, n)
  z    <- sine(n, f)
  
  ArmaLogistic   <- 0.1 * x    + 0.9 * y
  ArmaSine       <- 0.1 * x    + 0.9 * z
  AR2Logistic    <- 0.1 * x_ar + 0.9 * y
  AR2Sine        <- 0.1 * x_ar + 0.9 * z
  MA2Logistic    <- 0.1 * x_ma + 0.9 * y
  MA2Sine        <- 0.1 * x_ma + 0.9 * z
  LogisticSine   <- 0.1 * y    + 0.9 * z
  
  H_arma         <- HShannon(OPprob(x, D))
  C_arma         <- StatComplexity(OPprob(x, D))
  
  H_ar2          <- HShannon(OPprob(x_ar, D))
  C_ar2          <- StatComplexity(OPprob(x_ar, D))
  
  H_ma2          <- HShannon(OPprob(x_ma, D))
  C_ma2          <- StatComplexity(OPprob(x_ma, D))
  
  H_logistic     <- HShannon(OPprob(y, D))
  C_logistic     <- StatComplexity(OPprob(y, D))
  
  H_sine         <- HShannon(OPprob(z, D))
  C_sine         <- StatComplexity(OPprob(z, D))
  
  H_ArmaLogistic <- HShannon(OPprob(ArmaLogistic, D))
  C_ArmaLogistic <- StatComplexity(OPprob(ArmaLogistic, D))
  
  H_ArmaSine     <- HShannon(OPprob(ArmaSine, D))
  C_ArmaSine     <- StatComplexity(OPprob(ArmaSine, D))
  
  H_AR2Logistic  <- HShannon(OPprob(AR2Logistic, D))
  C_AR2Logistic  <- StatComplexity(OPprob(AR2Logistic, D))
  
  H_AR2Sine      <- HShannon(OPprob(AR2Sine, D))
  C_AR2Sine      <- StatComplexity(OPprob(AR2Sine, D))
  
  H_MA2Logistic  <- HShannon(OPprob(MA2Logistic, D))
  C_MA2Logistic  <- StatComplexity(OPprob(MA2Logistic, D))
  
  H_MA2Sine      <- HShannon(OPprob(MA2Sine, D))
  C_MA2Sine      <- StatComplexity(OPprob(MA2Sine, D))
  
  H_LogisticSine <- HShannon(OPprob(LogisticSine, D))
  C_LogisticSine <- StatComplexity(OPprob(LogisticSine, D))
  
  points_df <- data.frame(
    H = c(H_arma, H_ar2, H_ma2, H_logistic, H_sine,
          H_ArmaLogistic, H_ArmaSine,
          H_AR2Logistic,  H_AR2Sine,
          H_MA2Logistic,  H_MA2Sine,
          H_LogisticSine),
    C = c(C_arma, C_ar2, C_ma2, C_logistic, C_sine,
          C_ArmaLogistic, C_ArmaSine,
          C_AR2Logistic,  C_AR2Sine,
          C_MA2Logistic,  C_MA2Sine,
          C_LogisticSine),
    Model = c("ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
              "ARMA+Logistic", "ARMA+Sine",
              "AR2+Logistic",  "AR2+Sine",
              "MA2+Logistic",  "MA2+Sine",
              "Logistic+Sine"),
    rep = i
  )
  
  results_df <- rbind(results_df, points_df)
}

#write_xlsx(results_df, "HC_results.xlsx")

model_colors <- c(
  "ARMA(2,2)"     = "red",
  "AR(2)"         = "darkred",
  "MA(2)"         = "tomato",
  "Logistic"      = "blue",
  "Sine"          = "green",
  "ARMA+Logistic" = "purple",
  "ARMA+Sine"     = "black",
  "AR2+Logistic"  = "darkorange",
  "AR2+Sine"      = "brown",
  "MA2+Logistic"  = "deeppink",
  "MA2+Sine"      = "cyan4",
  "Logistic+Sine" = "orange"
)

model_shapes <- c(
  "ARMA(2,2)"     = 16,
  "AR(2)"         = 17,
  "MA(2)"         = 15,
  "Logistic"      = 18,
  "Sine"          = 8,
  "ARMA+Logistic" = 10,
  "ARMA+Sine"     = 11,
  "AR2+Logistic"  = 12,
  "AR2+Sine"      = 13,
  "MA2+Logistic"  = 14,
  "MA2+Sine"      = 3,
  "Logistic+Sine" = 4
)

data("LinfLsup")
lower_bound <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Lower")
upper_bound <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Upper")

p <- ggplot() +
  geom_line(data = lower_bound, aes(x = H, y = C), color = "grey", linetype = "dashed") +
  geom_line(data = upper_bound, aes(x = H, y = C), color = "grey", linetype = "dashed") +
  geom_point(data = results_df, aes(x = H, y = C, color = Model, shape = Model), size = 2, alpha = 0.6) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  theme_bw(base_size = 11, base_family = "serif") +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane (", R, " replications), D = ", D),
    color = "Model",
    shape = "Model"
  )
print(p)

#cor_matrix <- cor(points_df[, c("H", "C")])
#print(cor_matrix)