library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
r        <- 3.8
ma2_list <- list(ar = NULL, ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ma2       <- function(n) as.numeric(normalize(arima.sim(model = ma2_list, n)))

logistic <- function(r, n) {
  a    <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  return(a)
}

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x_ma    <- ma2(n)
y       <- logistic(r, n)
weights <- seq(0.1, 0.9, by = 0.1)

# ── HC values ─────────────────────────────────────────────────────────────────
pure <- bind_rows(get_HC(x_ma, "MA(2)"), get_HC(y, "Logistic"))

mix  <- bind_rows(lapply(weights, function(w) {
  get_HC(w * x_ma + (1 - w) * y, paste0("Mix (w=", w, ")"))
}))

results_df <- rbind(pure, mix)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "dashed") +
  geom_point(data = results_df, aes(x = H, y = C, color = Model),
             size = 2) +
  scale_color_manual(values = c(
    "MA(2)"    = "tomato",
    "Logistic" = "blue",
    setNames(viridis::viridis(9), paste0("Mix (w=", weights, ")"))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane — MA(2) + Logistic (D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

# End of the code 
#----------------------------------------------------------------------------

# MA2 + Logistic with 50 replications

library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
r        <- 3.8
R        <- 50
ma2_list <- list(ar = NULL, ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ma2       <- function(n) as.numeric(normalize(arima.sim(model = ma2_list, n)))

logistic <- function(r, n) {
  a    <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  return(a)
}

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Fixed signal ──────────────────────────────────────────────────────────────
y       <- logistic(r, n)
weights <- seq(0.1, 0.9, by = 0.1)

# ── 50 Replications ───────────────────────────────────────────────────────────
results_df <- bind_rows(lapply(1:R, function(i) {
  x_ma <- ma2(n)
  pure <- bind_rows(get_HC(x_ma, "MA(2)"), get_HC(y, "Logistic"))
  mix  <- bind_rows(lapply(weights, function(w) {
    get_HC(w * x_ma + (1 - w) * y, paste0("Mix (w=", w, ")"))
  }))
  mutate(rbind(pure, mix), Rep = i)
}))

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "dashed") +
  geom_point(data = results_df, aes(x = H, y = C, color = Model),
             size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c(
    "MA(2)"    = "tomato",
    "Logistic" = "blue",
    setNames(viridis::viridis(9), paste0("Mix (w=", weights, ")"))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane — MA(2) + Logistic (", R, " reps, D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))