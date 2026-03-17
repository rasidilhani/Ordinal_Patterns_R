library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
f        <- 0.04
ar2_list <- list(ar = c(-0.8, 0.1), ma = NULL)

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ar2       <- function(n) as.numeric(normalize(arima.sim(model = ar2_list, n)))
sine      <- function(n, f) as.numeric(sin(2 * pi * f * 1:n))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x_ar    <- ar2(n)
z       <- sine(n, f)
weights <- seq(0.1, 0.9, by = 0.1)

# ── HC values ─────────────────────────────────────────────────────────────────
pure <- bind_rows(get_HC(x_ar, "AR(2)"), get_HC(z, "Sine"))

mix  <- bind_rows(lapply(weights, function(w) {
  get_HC(w * x_ar + (1 - w) * z, paste0("AR2+Sine (w=", w, ")"))
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
    "AR(2)" = "darkred",
    "Sine"  = "green",
    setNames(viridis::viridis(9), paste0("AR2+Sine (w=", weights, ")"))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane — AR(2) + Sine (D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

# End of the code
#------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# AR2 + Sine with 50 replications

library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
f        <- 0.04
R        <- 50
ar2_list <- list(ar = c(-0.8, 0.1), ma = NULL)

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ar2       <- function(n) as.numeric(normalize(arima.sim(model = ar2_list, n)))
sine      <- function(n, f) as.numeric(sin(2 * pi * f * 1:n))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Fixed signal ──────────────────────────────────────────────────────────────
z       <- sine(n, f)
weights <- seq(0.1, 0.9, by = 0.1)

# ── 50 Replications ───────────────────────────────────────────────────────────
results_df <- bind_rows(lapply(1:R, function(i) {
  x_ar <- ar2(n)
  pure <- bind_rows(get_HC(x_ar, "AR(2)"), get_HC(z, "Sine"))
  mix  <- bind_rows(lapply(weights, function(w) {
    get_HC(w * x_ar + (1 - w) * z, paste0("AR2+Sine (w=", w, ")"))
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
    "AR(2)" = "darkred",
    "Sine"  = "green",
    setNames(viridis::viridis(9), paste0("AR2+Sine (w=", weights, ")"))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane — AR(2) + Sine (", R, " reps, D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))