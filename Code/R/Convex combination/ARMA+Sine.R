library(tidyverse)
library(StatOrdPattHxC)

D <- 4
n <- 1000
r <- 3.8
f <- 0.04

set.seed(1234567890, kind = "Mersenne-Twister")
arma_list <- list(ar = c(-0.8,0.1), ma = c(-0.8,0.1))

data("LinfLsup")
lower_bound <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Lower")
upper_bound <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Upper")


arma <- function(n){
  ts_data <- arima.sim(model = arma_list, n)
  ts_norm <- (ts_data - min(ts_data))/(max(ts_data)-min(ts_data))
  return(as.numeric(ts_norm))
}

x <- arma(1000)

Entropy_x <- HShannon(OPprob(x,4))
Complexity_x <- StatComplexity(OPprob(x, 4))

sine <- function(n, f){
  t <- 1:n
  b <- sin(2*pi*f*t)
  return(as.numeric(b))
}
z <- sine(1000, 0.04)

Entropy_z <- HShannon(OPprob(z,4))
Complexity_z <- StatComplexity(OPprob(z, 4))

L1 <- 0.1*x + 0.9*z 

Entropy_L1 <- HShannon(OPprob(L1,4))
Complexity_L1 <- StatComplexity(OPprob(L1, 4))

L2 <- 0.2*x + 0.8*z 

Entropy_L2 <- HShannon(OPprob(L2,4))
Complexity_L2 <- StatComplexity(OPprob(L2, 4))

L3 <- 0.3*x + 0.7*z 

Entropy_L3 <- HShannon(OPprob(L3,4))
Complexity_L3 <- StatComplexity(OPprob(L3, 4))

L4 <- 0.4*x + 0.6*z 

Entropy_L4 <- HShannon(OPprob(L4,4))
Complexity_L4 <- StatComplexity(OPprob(L4, 4))

L5 <- 0.5*x + 0.5*z 

Entropy_L5 <- HShannon(OPprob(L5,4))
Complexity_L5 <- StatComplexity(OPprob(L5, 4))

L6 <- 0.6*x + 0.4*z 

Entropy_L6 <- HShannon(OPprob(L6,4))
Complexity_L6 <- StatComplexity(OPprob(L6, 4))

L7 <- 0.7*x + 0.3*z 

Entropy_L7 <- HShannon(OPprob(L7,4))
Complexity_L7 <- StatComplexity(OPprob(L7, 4))

L8 <- 0.8*x + 0.2*z 

Entropy_L8 <- HShannon(OPprob(L8,4))
Complexity_L8 <- StatComplexity(OPprob(L8, 4))

L9 <- 0.9*x + 0.1*z 

Entropy_L9 <- HShannon(OPprob(L9,4))
Complexity_L9 <- StatComplexity(OPprob(L9, 4))


p <- ggplot()+
  geom_line(data = lower_bound, aes(x= H, y= C), color = "grey", linetype = "dashed")+
  geom_line(data = upper_bound, aes(x= H, y= C), color = "grey", linetype = "dashed")+
  geom_point(aes(x = Entropy_x, y = Complexity_x)) +
  geom_point(aes(x = Entropy_z, y = Complexity_z)) +
  geom_point(aes(x = Entropy_L1, y = Complexity_L1), color = "red") +
  geom_point(aes(x = Entropy_L2, y = Complexity_L2), color = "blue") +
  geom_point(aes(x = Entropy_L3, y = Complexity_L3), color = "purple") +
  geom_point(aes(x = Entropy_L4, y = Complexity_L4), color = "green") +
  geom_point(aes(x = Entropy_L5, y = Complexity_L5), color = "orange") +
  geom_point(aes(x = Entropy_L6, y = Complexity_L6), color = "yellow") +
  geom_point(aes(x = Entropy_L7, y = Complexity_L7), color = "brown") +
  geom_point(aes(x = Entropy_L8, y = Complexity_L8), color = "pink") +
  geom_point(aes(x = Entropy_L9, y = Complexity_L9), color = "grey") +
  
  theme_bw(base_size = 11, base_family = "serif")+
  labs(x = expression(italic(H)),
       y = expression(italic(C)),
       title = "HC plane"
       )
p


#-------------------------------------------------------------

library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
f         <- 0.04
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
arma <- function(n) {
  ts_data <- arima.sim(model = arma_list, n)
  as.numeric((ts_data - min(ts_data)) / (max(ts_data) - min(ts_data)))
}

sine <- function(n, f) {
  as.numeric(sin(2 * pi * f * 1:n))
}

get_HC <- function(series) {
  prob <- OPprob(series, D)
  tibble(H = HShannon(prob), C = StatComplexity(prob))
}

# ── Signals ───────────────────────────────────────────────────────────────────
x <- arma(n)
z <- sine(n, f)

# ── HC values ─────────────────────────────────────────────────────────────────
weights <- seq(0.1, 0.9, by = 0.1)

HC_all <- bind_rows(
  mutate(get_HC(x), Model = "ARMA(2,2)"),
  mutate(get_HC(z), Model = "Sine"),
  map(weights, function(w) {
    mix <- w * x + (1 - w) * z
    hc  <- get_HC(mix)
    mutate(hc, Model = paste0("Mix ", w))
  })
)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- LinfLsup |> filter(Dimension == as.character(D))

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "dashed") +
  geom_point(data = HC_all, aes(x = H, y = C, color = Model), size = 2) +
  scale_color_manual(values = c(
    "ARMA(2,2)"    = "red",
    "Sine"         = "blue",
    setNames(viridis::viridis(9), paste0("Mix ", weights))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane for ARMA(2,2), Sine and Mixed model  (D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

# End of the code
#------------------------------------------------------------------------------


#__________________________________________________________________________________
# In this code I replicated the 11 models to see how it appear in the HC plane
#__________________________________________________________________________________
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
f         <- 0.04
reps      <- 50
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
arma <- function(n) {
  ts_data <- arima.sim(model = arma_list, n)
  as.numeric((ts_data - min(ts_data)) / (max(ts_data) - min(ts_data)))
}

sine <- function(n, f) {
  as.numeric(sin(2 * pi * f * 1:n))
}

get_HC <- function(series) {
  prob <- OPprob(series, D)
  tibble(H = HShannon(prob), C = StatComplexity(prob))
}

# ── Signals ───────────────────────────────────────────────────────────────────
z       <- sine(n, f)
weights <- seq(0.1, 0.9, by = 0.1)

# ── 50 Replications ───────────────────────────────────────────────────────────
HC_all <- map(1:reps, function(rep) {
  x      <- arma(n)
  HC_x   <- mutate(get_HC(x), Model = "ARMA(2,2)")
  HC_z   <- mutate(get_HC(z), Model = "Sine")
  
  HC_mix <- map(weights, function(w) {
    mix   <- w * x + (1 - w) * z
    hc    <- get_HC(mix)
    mutate(hc, Model = paste0("Mix ", w))
  })
  
  hc_rep <- bind_rows(HC_x, HC_z, HC_mix)
  mutate(hc_rep, Rep = rep)
})

HC_all <- bind_rows(HC_all)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "dashed") +
  geom_point(data = HC_all, aes(x = H, y = C, color = Model),
             size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c(
    "ARMA(2,2)" = "red",
    "Sine"      = "blue",
    setNames(viridis::viridis(9), paste0("Mix ", weights))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane — ", reps, " replications (D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

