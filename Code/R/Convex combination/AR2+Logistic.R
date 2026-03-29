library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
r        <- 3.8
ar2_list <- list(ar = c(-0.8, 0.1), ma = NULL)

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ar2       <- function(n) as.numeric(normalize(arima.sim(model = ar2_list, n)))

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
x_ar    <- ar2(n)
y       <- logistic(r, n)
weights <- seq(0.1, 0.9, by = 0.1)

# ── HC values ─────────────────────────────────────────────────────────────────
pure <- bind_rows(get_HC(x_ar, "AR(2)"), get_HC(y, "Logistic"))

mix  <- bind_rows(lapply(weights, function(w) {
  #get_HC(w * x_ar + (1 - w) * y, paste0("AR2+Logistic (w=", w, ")"))
  get_HC(w * x_ar + (1 - w) * y,
         paste0("AR2+Logistic(w=", w, ")"))
}))

results_df <- rbind(pure, mix)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "AR(2)" = expression(italic(AR)(2)),
  "Logistic" = expression(Logistic),
  setNames(
    lapply(weights, function(w)
      bquote(italic(AR)(2) + Logistic ~ "(" * italic(w) == .(w) * ")")
    ),
    paste0("AR2+Logistic(w=", weights, ")")
  )
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(
    data = bounds,
    aes(x = H, y = C, group = Side),
    color = "grey60",
    linetype = "solid"
  ) +
  geom_point(
    data = results_df,
    aes(x = H, y = C, color = Model),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "AR(2)"    = "darkred",
      "Logistic" = "blue",
      setNames(viridis::viridis(length(weights)),
               paste0("AR2+Logistic(w=", weights, ")"))
    ),
    breaks = names(legend_labels),
    labels = legend_labels
  ) +
  labs(
    x = expression(italic(H)),
    y = expression(italic(C)),
    title = bquote(
      "HC Plane — " * italic(AR)(2) *
        " + Logistic (" * italic(D) == .(D) * ")"
    ),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))

#End of the code 
#------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# AR2 +Logistic with 50 replications

library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
r        <- 3.8
R        <- 50
ar2_list <- list(ar = c(-0.8, 0.1), ma = NULL)

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ar2       <- function(n) as.numeric(normalize(arima.sim(model = ar2_list, n)))

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

# ── Fixed signal ──────────────────────────────────────────────────────────────
y       <- logistic(r, n)
weights <- seq(0.1, 0.9, by = 0.1)

# ── 50 Replications ───────────────────────────────────────────────────────────
results_df <- bind_rows(lapply(1:R, function(i) {
  x_ar <- ar2(n)
  
  pure <- bind_rows(
    get_HC(x_ar, "AR(2)"),
    get_HC(y, "Logistic")
  )
  
  mix <- bind_rows(lapply(weights, function(w) {
    # ✅ CHANGED HERE: consistent model name
    get_HC(
      w * x_ar + (1 - w) * y,
      paste0("AR2+Logistic(w=", w, ")")
    )
  }))
  
  mutate(bind_rows(pure, mix), Rep = i)
}))

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

legend_labels <- c(
  "AR(2)" = expression(italic(AR)(2)),
  "Logistic" = expression(Logistic),
  setNames(
    lapply(weights, function(w)
      bquote(italic(AR)(2) + Logistic ~ "(" * italic(w) == .(w) * ")")
    ),
    paste0("AR2+Logistic(w=", weights, ")")
  )
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(
    data = bounds,
    aes(x = H, y = C, group = Side),
    color = "grey60",
    linetype = "solid"
  ) +
  geom_point(
    data = results_df,
    aes(x = H, y = C, color = Model),
    size = 1.5,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "AR(2)"    = "darkred",
      "Logistic" = "blue",
      setNames(
        viridis::viridis(length(weights)),
        paste0("AR2+Logistic(w=", weights, ")")
      )
    ),
    breaks = names(legend_labels),
    labels = legend_labels
  ) +
  labs(
    x = expression(italic(H)),
    y = expression(italic(C)),
    title = bquote(
      "HC Plane — " * italic(AR)(2) *
        " + Logistic (" * .(R) * " reps, " * italic(D) == .(D) * ")"
    ),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5))
