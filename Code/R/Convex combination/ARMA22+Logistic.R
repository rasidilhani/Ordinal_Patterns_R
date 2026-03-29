library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

arma <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))

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
x       <- arma(n)
y       <- logistic(r, n)
weights <- seq(0.1, 0.9, by = 0.1)

# ── HC values ─────────────────────────────────────────────────────────────────
pure <- bind_rows(
  get_HC(x, "ARMA(2,2)"),
  get_HC(y, "Logistic")
)

mix <- bind_rows(lapply(weights, function(w) {
  get_HC(w * x + (1 - w) * y, paste0("ARMA+Logistic(w=", w, ")"))
}))

results_df <- rbind(pure, mix)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "ARMA(2,2)" = expression(italic(ARMA)(2*","*2)),
  "Logistic"  = "Logistic",
  setNames(
    lapply(weights, function(w)
      bquote(italic(ARMA) + Logistic ~ (italic(w) == .(w)))
    ),
    paste0("ARMA+Logistic(w=", weights, ")")
  )
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "solid") +
  geom_point(data = results_df, aes(x = H, y = C, color = Model),
             size = 2) +
  scale_color_manual(
    values = c(
      "ARMA(2,2)" = "red",
      "Logistic"  = "blue",
      setNames(viridis::viridis(9), paste0("ARMA+Logistic(w=", weights, ")"))
    ),
    labels = legend_labels
  ) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = bquote(HC - Plane ~ italic(ARMA)(2*","*2) + Logistic ~ (italic(D) == .(D))),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))




# ------------------End of the code----------------------
#------------------------------------------------------------------



# ARMA22+Sine with 50 replications

library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
R         <- 50
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
arma      <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))

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
  x    <- arma(n)
  pure <- bind_rows(get_HC(x, "ARMA(2,2)"), get_HC(y, "Logistic"))
  mix  <- bind_rows(lapply(weights, function(w) {
    get_HC(w * x + (1 - w) * y, paste0("ARMA+Logistic(w=", w, ")"))
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
    "ARMA(2,2)" = "red",
    "Logistic"  = "blue",
    setNames(viridis::viridis(9), paste0("ARMA+Logistic(w=", weights, ")"))
  )) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane — ARMA(2,2) + Logistic (", R, " reps, D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

# End of the code--------------------------------------------
#================================================================================



# This code is make based on the comments from supervisor to make simple and understandable replication
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D   <- 4
n   <- 1000
r   <- 3.8
R   <- 50
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890)

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
arma      <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))

logistic <- function(r, n) {
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  a
}

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Only selected weights ─────────────────────────────────────────────────────
weights <- c(0.1, 0.3, 0.6)

# ── Fixed logistic series ─────────────────────────────────────────────────────
y <- logistic(r, n)

# ── 50 replications ───────────────────────────────────────────────────────────
results_df <- bind_rows(lapply(1:R, function(i) {
  x <- arma(n)
  
  pure <- bind_rows(
    get_HC(x, "ARMA(2,2)"),
    get_HC(y, "Logistic")
  )
  
  mix <- bind_rows(lapply(weights, function(w) {
    get_HC(
      w * x + (1 - w) * y,
      paste0("ARMA+Logistic(w=", w, ")")
    )
  }))
  
  mutate(bind_rows(pure, mix), Rep = i)
}))

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Legend labels ─────────────────────────────────────────────────────────────  # <-- CHANGED BLOCK
legend_labels <- c(
  "ARMA(2,2)" = expression(italic(ARMA)(2*","*2)),
  "Logistic"  = "Logistic",
  setNames(
    lapply(weights, function(w)
      bquote(italic(ARMA) + Logistic ~ (italic(w) == .(w)))
    ),
    paste0("ARMA+Logistic(w=", weights, ")")
  )
)
# ── Colors ────────────────────────────────────────────────────────────────────
legend_cols <- c(
  "ARMA(2,2)" = "red",
  "Logistic"  = "blue",
  setNames(viridis::viridis(length(weights)),
           paste0("ARMA+Logistic(w=", weights, ")"))
)
# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "solid") +
  geom_point(data = results_df, aes(x = H, y = C, color = Model),
             size = 1.5, alpha = 0.5) +
  scale_color_manual(values = legend_cols, labels = legend_labels) +  # <-- CHANGED: labels = legend_labels
  labs(
    x = expression(italic(H)),
    y = expression(italic(C)),
    title = bquote(HC ~ Plane ~ italic(ARMA)(2*","*2) + Logistic ~ (.(R) ~ reps*"," ~ italic(D) == .(D))),  # <-- CHANGED
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))
