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
  get_HC(w * x_ma + (1 - w) * y, paste0("MA2+Logistic(w=", w, ")"))
}))

results_df <- rbind(pure, mix)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "MA(2)" = "MA(2)",
  "Logistic" = "Logistic",
  setNames(
    lapply(weights, function(w)
      bquote(MA(2) + Logistic ~ "(" * italic(w) == .(w) * ")")
    ),
    paste0("MA2+Logistic(w=", weights, ")")
  )
)

# ── Output path ───────────────────────────────────────────────────────────────
output_dir  <- file.path("Results", "Convex_combination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "MA2+Logistic.pdf")

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(
    data = bounds,
    aes(x = H, y = C, group = Side),
    color = "grey60", linetype = "solid"
  ) +
  geom_point(
    data = results_df,
    aes(x = H, y = C, color = Model),
    size = 2
  ) +
  scale_color_manual(
    values = c(
      "MA(2)" = "tomato",
      "Logistic" = "blue",
      setNames(
        viridis::viridis(length(weights)),
        paste0("MA2+Logistic(w=", weights, ")")
      )
    ),
    breaks = names(legend_labels),
    labels = legend_labels
  ) +
  labs(
    x = expression(italic(H)),
    y = expression(italic(C)),
   # title = bquote(
  #     italic(H) %*% italic(C) ~ "Plane," ~ MA(2) + Logistic ~
  #    (italic(D) == .(D))
  #    ),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave(output_file, width = 8, height = 5, dpi = 300)

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
    get_HC(w * x_ma + (1 - w) * y, paste0("MA2+Logistic(w=", w, ")"))
  }))
  mutate(rbind(pure, mix), Rep = i)
}))

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "MA(2)" = "MA(2)",
  "Logistic" = "Logistic",
  setNames(
    lapply(weights, function(w)
      bquote(MA(2) + Logistic ~ "(" * italic(w) == .(w) * ")")
    ),
    paste0("MA2+Logistic(w=", weights, ")")
  )
)

# ── Output path ───────────────────────────────────────────────────────────────
output_dir  <- file.path("Results", "Convex_combination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "MA2+Logistic_50Rep.pdf")

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot() +
  geom_line(
    data = bounds,
    aes(x = H, y = C, group = Side),
    color = "grey60", linetype = "solid"
  ) +
  geom_point(
    data = results_df,
    aes(x = H, y = C, color = Model),
    size = 1.5, alpha = 0.5
  ) +
  scale_color_manual(
    values = c(
      "MA(2)" = "tomato",
      "Logistic" = "blue",
      setNames(
        viridis::viridis(length(weights)),
        paste0("MA2+Logistic(w=", weights, ")")
      )
    ),
    breaks = names(legend_labels),
    labels = legend_labels
  ) +
  labs(
    x = expression(italic(H)),
    y = expression(italic(C)),
    # title = bquote(
    #    italic(H) %*% italic(C) ~ "Plane," ~ MA(2) + Logistic ~
    #      (.(R) ~ "reps," ~ italic(D) == .(D))
    #  ),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
ggsave(output_file, width = 8, height = 5, dpi = 300)