# This code generate to show the couple of examples stochastic components riding sine function
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
f         <- 0.04
ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

ar2  <- function(n) as.numeric(normalize(arima.sim(model = ar2_list,  n)))

sine <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

# ── Weights ──────────────────────────────────────────────
weights <- seq(0.1, 0.9, by = 0.1)

# ── Signals from your base code ───────────────────────────────────────────────
x_ar <- ar2(n)
z    <- sine(n, f)
t    <- 1:n


x_ar <- normalize(x_ar)
z    <- normalize(z)

# ── AR(2) riding the sine  ─────────────────────────────────────
example1_df <- bind_rows(lapply(weights, function(w) {
  tibble(
    t = t,
    value = w * x_ar + (1 - w) * z,
    w = w
  )
}))

w_labels <- setNames(
  lapply(weights, function(w)
    bquote(italic(w) == .(w))
  ),
  as.character(weights)
)

# ── Plot ──────────────────────────────────────────────────────────────────────
      ggplot(example1_df, aes(t, value, color = factor(w))) +
        geom_line(alpha = 0.7) +
        theme_bw(base_size = 11, base_family = "serif") +
        scale_color_discrete(
          labels = w_labels
        ) +
        labs(
          title = bquote("AR(2) stochastic component riding the sine"),
          x = "time",
          y = "",
          color = expression(italic(w))
        )
      
      
#--------------------------------------------------------------------------

# This code generate to show the MA2 stochastic components riding sine function
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
f         <- 0.04
ma2_list  <- list(ar = NULL,          ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

ma2  <- function(n) as.numeric(normalize(arima.sim(model = ma2_list,  n)))

sine <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

# ── Weights ──────────────────────────────────────────────
weights <- seq(0.1, 0.9, by = 0.1)

# ── Signals from base code ───────────────────────────────────────────────
x_ma <- ma2(n)
z    <- sine(n, f)
t    <- 1:n


x_ma <- normalize(x_ma)
z    <- normalize(z)

# ── ARMA riding the sine  ─────────────────────────────────────
example2_df <- bind_rows(lapply(weights, function(w) {
  tibble(
    t = t,
    value = w * x_ma + (1 - w) * z,
    w = w
  )
}))

w_labels <- setNames(
  lapply(weights, function(w)
    bquote(italic(w) == .(w))
  ),
  as.character(weights)
)

# ── Plot ──────────────────────────────────────────────────────────────────────
      ggplot(example1_df, aes(t, value, color = factor(w))) +
        geom_line(alpha = 0.7) +
        theme_bw(base_size = 11, base_family = "serif") +
        scale_color_discrete(
          labels = w_labels
        ) +
        labs(
          title = bquote("MA(2) stochastic component riding the sine"),
          x = "time",
          y = "",
          color = expression(italic(w))
        )

#==========================================================================================
# This code generate to show the ARMA components riding sine function
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
r         <- 3.8
f         <- 0.04
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

arma <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))

logistic <- function(r, n) {
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  a
}

sine <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

# ── Weights ──────────────────────────────────────────────
weights <- seq(0.1, 0.9, by = 0.1)

# ── Signals from your base code ───────────────────────────────────────────────
x_arma <- arma(n)
z    <- sine(n, f)
t    <- 1:n


x_arma <- normalize(x_arma)
z    <- normalize(z)

# ── AR(2) riding the sine  ─────────────────────────────────────
example1_df <- bind_rows(lapply(weights, function(w) {
  tibble(
    t = t,
    value = w * x_arma + (1 - w) * z,
    w = w
  )
}))

w_labels <- setNames(
  lapply(weights, function(w)
    bquote(italic(w) == .(w))
  ),
  as.character(weights)
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(example1_df, aes(t, value, color = factor(w))) +
  geom_line(alpha = 0.7) +
  theme_bw(base_size = 11, base_family = "serif") +
  scale_color_discrete(
  labels = w_labels
  ) +
  labs(
  title = bquote("ARMA(2,2) stochastic component riding the sine"),
  x = "time",
  y = "",
  color = expression(italic(w))
)

#=======================================================================
# This code shows Stochastic function riding to the Logistic

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
  a <- numeric(n)
  a[1] <- 0.5
  for (i in 2:n) a[i] <- r * a[i-1] * (1 - a[i-1])
  a
}

sine <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

# ── Weights ──────────────────────────────────────────────
weights <- seq(0.1, 0.9, by = 0.1)

# ── Signals from your base code ───────────────────────────────────────────────
x_arma <- arma(n)
z    <- sine(n, f)
t    <- 1:n


x_arma <- normalize(x_arma)
z    <- normalize(z)

# ── AR(2) riding the sine  ─────────────────────────────────────
example1_df <- bind_rows(lapply(weights, function(w) {
  tibble(
    t = t,
    value = w * x_arma + (1 - w) * z,
    w = w
  )
}))

w_labels <- setNames(
  lapply(weights, function(w)
    bquote(italic(w) == .(w))
  ),
  as.character(weights)
)

# ── Plot ──────────────────────────────────────────────────────────────────────
ggplot(example1_df, aes(t, value, color = factor(w))) +
  geom_line(alpha = 0.7) +
  theme_bw(base_size = 11, base_family = "serif") +
  scale_color_discrete(
    labels = w_labels
  ) +
  labs(
    title = bquote("ARMA(2,2) stochastic component riding the sine"),
    x = "time",
    y = "",
    color = expression(italic(w))
  )

#End of the code
#======================================================================
library(tidyverse)

# ── Parameters ───────────────────────────────────────────────
n <- 1000
f <- 0.04
weights <- seq(0.1, 0.9, by = 0.1)

ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_list  <- list(ar = NULL, ma = c(-0.8, 0.1))
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890)

# ── Functions ───────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

ar2  <- function(n) normalize(arima.sim(model = ar2_list,  n))
ma2  <- function(n) normalize(arima.sim(model = ma2_list,  n))
arma <- function(n) normalize(arima.sim(model = arma_list, n))

sine <- function(n, f) as.numeric(normalize(sin(2 * pi * f * (1:n))))

t <- 1:n
z <- sine(n, f)

# ── Helper to build riding-sine data ───────────────────────
build_df <- function(x, label) {
  bind_rows(lapply(weights, function(w) {
    tibble(
      t = t,
      value = w * x + (1 - w) * z,
      w = factor(w),
      ModelType = label
    )
  }))
}

df_all <- bind_rows(
  build_df(ar2(n),  "AR(2) riding sine"),
  build_df(ma2(n),  "MA(2) riding sine"),
  build_df(arma(n), "ARMA(2,2) riding sine")
)

# ── Math-mode legend labels for w ───────────────────────────
w_labels <- setNames(
  lapply(weights, function(w) bquote(italic(w) == .(w))),
  as.character(weights)
)

# ── Combined plot ──────────────────────────────────────────
ggplot(df_all, aes(t, value, color = w)) +
  geom_line(alpha = 0.6) +
  facet_wrap(~ ModelType, ncol = 1) +
  scale_color_discrete(labels = w_labels) +
  theme_bw(base_size = 11, base_family = "serif") +
  labs(
    title = "Stochastic components riding a sinusoidal signal",
    x = "time",
    y = "Amplitude",
    color = expression(italic(w))
  )


#=====================================================
# In this code, we calculate entropy and complexity based on the previous work
#====================================================================
#================================
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ───────────────────────────────────────────────
D <- 4
n <- 1000
f <- 0.04
weights <- seq(0, 1, by = 0.1)   # includes w = 0 (sine) and w = 1 (pure stochastic)

set.seed(1234567890)

# ── Model definitions ────────────────────────────────────────
ar2_list  <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_list  <- list(ar = NULL, ma = c(-0.8, 0.1))
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

# ── Utility functions ────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

ar2   <- function(n) normalize(arima.sim(model = ar2_list,  n))
ma2   <- function(n) normalize(arima.sim(model = ma2_list,  n))
arma <- function(n) normalize(arima.sim(model = arma_list, n))

sine <- function(n, f) normalize(sin(2 * pi * f * (1:n)))

get_HC <- function(series, family, w) {
  p <- OPprob(series, D)
  tibble(
    H = HShannon(p),
    C = StatComplexity(p),
    Family = family,
    w = w
  )
}

# ── Base sine ────────────────────────────────────────────────
z <- sine(n, f)

# ── Build HC points for one family ───────────────────────────
build_family <- function(x, family_name) {
  bind_rows(lapply(weights, function(w) {
    series <- if (w == 0) z else if (w == 1) x else w * x + (1 - w) * z
    get_HC(series, family_name, w)
  }))
}

# ── Compute HC data for all families ─────────────────────────
HC_df <- bind_rows(
  build_family(ar2(n),  "AR(2)"),
  build_family(ma2(n),  "MA(2)"),
  build_family(arma(n), "ARMA(2,2)")
)

# ── HC plane bounds ──────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Combined HC‑plane plot ───────────────────────────────────
ggplot() +
  geom_line(
    data = bounds,
    aes(H, C, group = Side),
    color = "grey60",
    linetype = "solid"
  ) +
  geom_path(
    data = HC_df %>% filter(w > 0 & w < 1),
    aes(H, C, group = Family),
    color = "black"
  ) +
  geom_point(
    data = HC_df,
    aes(H, C, color = w, shape = Family),
    size = 3
  ) +
  facet_wrap(~ Family, ncol = 3) +
  scale_color_viridis_c(
    name = expression(italic(w))   
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  labs(
    title = bquote(
      "Entropy–complexity plane: stochastic processes riding a sine (" *
        italic(D) == .(D) * ")"
    ),
    x = expression(italic(H)),
    y = expression(italic(C)),
    shape = "Model"
  )

ggplot() +
  geom_line(
    data = bounds,
    aes(H, C, group = Side),
    color = "grey60",
    linetype = "solid"
  ) +
  geom_path(
    data = HC_df %>% filter(w > 0 & w < 1),
    aes(H, C, group = Family),
    color = "black"
  ) +
  geom_point(
    data = HC_df,
    aes(H, C, color = w, shape = Family),
    size = 3
  ) +
  facet_wrap(~ Family, ncol = 3) +
  scale_color_viridis_c(
    name = expression(italic(w))   
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  labs(
    title = bquote(
      "Entropy–complexity plane: stochastic processes riding a sine (" *
        italic(D) == .(D) * ")"
    ),
    x = expression(italic(H)),
    y = expression(italic(C)),
    shape = "Model"
  )
