library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
f         <- 0.04
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_list  <- list(ar = NULL,          ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

arma <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))
ar2  <- function(n) as.numeric(normalize(arima.sim(model = ar2_list,  n)))
ma2  <- function(n) as.numeric(normalize(arima.sim(model = ma2_list,  n)))

logistic <- function(r, n) {
  a    <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  return(a)
}

sine <- function(n, f) as.numeric(sin(2 * pi * f * 1:n))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x    <- arma(n)
x_ar <- ar2(n)
x_ma <- ma2(n)
y    <- logistic(r, n)
z    <- sine(n, f)

weights <- seq(0.1, 0.9, by = 0.1)

# ── Pairs to mix ──────────────────────────────────────────────────────────────
mix_pairs <- list(
  list(s1 = x,    s2 = y, name = "ARMA+Logistic"),
  list(s1 = x,    s2 = z, name = "ARMA+Sine"),
  list(s1 = x_ar, s2 = y, name = "AR2+Logistic"),
  list(s1 = x_ar, s2 = z, name = "AR2+Sine"),
  list(s1 = x_ma, s2 = y, name = "MA2+Logistic"),
  list(s1 = x_ma, s2 = z, name = "MA2+Sine")
  #list(s1 = y,    s2 = z, name = "Logistic+Sine")
)

# ── HC values: 5 pure + 6 pairs x 9 weights = 5 + 54 = 59 points ─────────────
pure_df <- bind_rows(
  get_HC(x,    "ARMA(2,2)"),
  get_HC(x_ar, "AR(2)"),
  get_HC(x_ma, "MA(2)"),
  get_HC(y,    "Logistic"),
  get_HC(z,    "Sine")
)

mix_df <- bind_rows(lapply(mix_pairs, function(pair) {
  bind_rows(lapply(weights, function(w) {
    mix <- w * pair$s1 + (1 - w) * pair$s2
    get_HC(mix, paste0(pair$name, " (w=", w, ")"))
  }))
}))

results_df <- rbind(pure_df, mix_df)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Colors ────────────────────────────────────────────────────────────────────
pure_colors <- c(
  "ARMA(2,2)" = "red",
  "AR(2)"     = "darkred",
  "MA(2)"     = "tomato",
  "Logistic"  = "blue",
  "Sine"      = "green"
)

pair_colors <- c(
  "ARMA+Logistic" = "purple",
  "ARMA+Sine"     = "black",
  "AR2+Logistic"  = "darkorange",
  "AR2+Sine"      = "brown",
  "MA2+Logistic"  = "deeppink",
  "MA2+Sine"      = "cyan4"
  #"Logistic+Sine" = "orange"
)

# assign each mixture point the color of its pair
results_df <- mutate(results_df,
                     PairGroup = case_when(
                       Model %in% names(pure_colors) ~ Model,
                       TRUE ~ sub(" \\(w=.*\\)", "", Model)   # strip weight label
                     )
)

all_colors <- c(pure_colors, pair_colors)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "dashed") +
  geom_point(data = results_df, aes(x = H, y = C, color = PairGroup),
             size = 2, alpha = 0.7) +
  scale_color_manual(values = all_colors) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane (D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

#End of the code

#_________________________________________________
# This code has 50 replications for the above code.

library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
f         <- 0.04
R         <- 50
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_list  <- list(ar = NULL,          ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

arma <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))
ar2  <- function(n) as.numeric(normalize(arima.sim(model = ar2_list,  n)))
ma2  <- function(n) as.numeric(normalize(arima.sim(model = ma2_list,  n)))

logistic <- function(r, n) {
  a    <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  return(a)
}

sine <- function(n, f) as.numeric(sin(2 * pi * f * 1:n))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Fixed signals (deterministic — no need to replicate) ──────────────────────
y <- logistic(r, n)
z <- sine(n, f)

weights <- seq(0.1, 0.9, by = 0.1)

# ── 50 Replications ───────────────────────────────────────────────────────────
results_df <- bind_rows(lapply(1:R, function(i) {
  
  x    <- arma(n)
  x_ar <- ar2(n)
  x_ma <- ma2(n)
  
  mix_pairs <- list(
    list(s1 = x,    s2 = y, name = "ARMA+Logistic"),
    list(s1 = x,    s2 = z, name = "ARMA+Sine"),
    list(s1 = x_ar, s2 = y, name = "AR2+Logistic"),
    list(s1 = x_ar, s2 = z, name = "AR2+Sine"),
    list(s1 = x_ma, s2 = y, name = "MA2+Logistic"),
    list(s1 = x_ma, s2 = z, name = "MA2+Sine"),
    list(s1 = y,    s2 = z, name = "Logistic+Sine")
  )
  
  pure_df <- bind_rows(
    get_HC(x,    "ARMA(2,2)"),
    get_HC(x_ar, "AR(2)"),
    get_HC(x_ma, "MA(2)"),
    get_HC(y,    "Logistic"),
    get_HC(z,    "Sine")
  )
  
  mix_df <- bind_rows(lapply(mix_pairs, function(pair) {
    bind_rows(lapply(weights, function(w) {
      mix <- w * pair$s1 + (1 - w) * pair$s2
      get_HC(mix, paste0(pair$name, " (w=", w, ")"))
    }))
  }))
  
  rep_df <- rbind(pure_df, mix_df)
  mutate(rep_df, Rep = i)
}))

# strip weight label for grouping color
results_df <- mutate(results_df,
                     PairGroup = sub(" \\(w=.*\\)", "", Model)
)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Colors ────────────────────────────────────────────────────────────────────
all_colors <- c(
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

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "dashed") +
  geom_point(data = results_df, aes(x = H, y = C, color = PairGroup),
             size = 1.5, alpha = 0.5) +
  scale_color_manual(values = all_colors) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane (", R, " replications, D = ", D, ")"),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))