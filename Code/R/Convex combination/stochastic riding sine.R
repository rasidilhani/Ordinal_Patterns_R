library(tidyverse)
library(StatOrdPattHxC)
library(gridExtra)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
f         <- 0.04
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
arma      <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))
sine      <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x <- arma(n)
z <- sine(n, f)

# w = 0 is pure Sine, w = 1 is pure ARMA
weights     <- seq(0.1, 0.9, by = 0.1)
all_weights <- c(0, weights, 1)

# ── Time series data ──────────────────────────────────────────────────────────
ts_df <- bind_rows(
  tibble(t = 1:n, value = z,    Model = "Sine (w=0)"),
  tibble(t = 1:n, value = x,    Model = "ARMA(2,2) (w=1)"),
  bind_rows(lapply(weights, function(w) {
    tibble(t = 1:n, value = w * x + (1 - w) * z,
           Model = paste0("w=", w))
  }))
)

# ── HC values ─────────────────────────────────────────────────────────────────
HC_df <- bind_rows(
  get_HC(z, "Sine (w=0)"),
  get_HC(x, "ARMA(2,2) (w=1)"),
  bind_rows(lapply(weights, function(w) {
    get_HC(w * x + (1 - w) * z, paste0("w=", w))
  }))
)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Colors ────────────────────────────────────────────────────────────────────
mix_colors  <- setNames(viridis::viridis(length(weights)), paste0("w=", weights))
model_colors <- c(
  "Sine (w=0)"     = "blue",
  "ARMA(2,2) (w=1)" = "red",
  mix_colors
)

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "Sine (w=0)"      = "Sine (w=0)",
  "ARMA(2,2) (w=1)" = "ARMA(2,2) (w=1)",
  setNames(
    lapply(weights, function(w) bquote(italic(w) == .(w))),
    paste0("w=", weights)
  )
)

# ── Output path ───────────────────────────────────────────────────────────────
output_dir  <- file.path("Results", "Convex_combination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "HC_plane_ARMA_Sine_surrounding.pdf")

# ── HC plane (center plot, no legend) ────────────────────────────────────────
p_hc <- ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "solid") +
  geom_point(data = HC_df, aes(x = H, y = C, color = Model), size = 2) +
  scale_color_manual(values = model_colors, labels = legend_labels) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = bquote(italic(H) %*% italic(C) ~ "Plane," ~
                     ARMA(2*","*2) + Sine ~ (italic(D) == .(D))),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p_hc_no_legend <- p_hc + theme(legend.position = "none")

# ── Time series thumbnail function ───────────────────────────────────────────
make_ts_plot <- function(model_label) {
  df     <- filter(ts_df, Model == model_label)
  n_show <- min(500, nrow(df))
  col    <- model_colors[model_label]
  
  ggplot(df[1:n_show, ], aes(x = t, y = value)) +
    geom_line(color = col, linewidth = 0.5) +
    labs(title = model_label, x = NULL, y = NULL) +
    theme_bw(base_family = "serif", base_size = 8) +
    theme(
      plot.title       = element_text(size = 7, face = "bold", hjust = 0.5,
                                      color = col),
      axis.text        = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
      panel.border     = element_rect(color = col, linewidth = 1.2),
      plot.margin      = unit(c(3, 3, 3, 3), "pt")
    )
}

# ── Build all thumbnail plots ─────────────────────────────────────────────────
model_names <- c(
  "Sine (w=0)",
  paste0("w=", weights),
  "ARMA(2,2) (w=1)"
)

ts_plots <- setNames(
  lapply(model_names, make_ts_plot),
  model_names
)

# ── Layout: thumbnails surrounding the HC plane ───────────────────────────────
# Top row:    w=0, w=0.1, w=0.2, w=0.3
# Left col:   w=0.4, w=0.5
# Right col:  w=0.6, w=0.7
# Bottom row: w=0.8, w=0.9, w=1
# Center:     HC plane

all_plots        <- ts_plots
all_plots[["HC"]] <- p_hc_no_legend

top_row    <- c("Sine (w=0)", "w=0.1", "w=0.2", "w=0.3")
left_col   <- c("w=0.4", "w=0.5")
right_col  <- c("w=0.6", "w=0.7")
bottom_row <- c("w=0.8", "w=0.9", "ARMA(2,2) (w=1)")

layout_matrix <- rbind(
  c( NA,  1,  2,  3, 4),
  c( 5, 12, 12, 12,  10),
  c( 6, 12, 12, 12,  11),
  c(NA, 7, 8, 9, NA)
)

plot_order <- c(top_row, left_col, right_col, bottom_row, "HC")
grob_list  <- lapply(plot_order, function(nm) all_plots[[nm]])

combined_plot <- grid.arrange(
  grobs         = grob_list,
  layout_matrix = layout_matrix,
  widths        = c(1.2, 1.2, 1.2, 1.2, 1.2),
  heights       = c(1.2, 1.2, 1.2, 1.2)
)

ggsave(output_file, plot = combined_plot,
       width = 7, height = 6, dpi = 300, device = cairo_pdf)

# -----------End of code------------------------------------------------
#-----------------------------------------------------------------------


#-------------------------------------------------
# AR2+Sine riding
#------------------------------------

library(tidyverse)
library(StatOrdPattHxC)
library(gridExtra)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
f         <- 0.04
ar2_list <- list(ar = c(-0.8, 0.1), ma = NULL)

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ar2      <- function(n) as.numeric(normalize(arima.sim(model = ar2_list, n)))
sine      <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x <- ar2(n)
z <- sine(n, f)

# w = 0 is pure Sine, w = 1 is pure ARMA
weights     <- seq(0.1, 0.9, by = 0.1)
all_weights <- c(0, weights, 1)

# ── Time series data ──────────────────────────────────────────────────────────
ts_df <- bind_rows(
  tibble(t = 1:n, value = z,    Model = "Sine (w=0)"),
  tibble(t = 1:n, value = x,    Model = "AR(2) (w=1)"),
  bind_rows(lapply(weights, function(w) {
    tibble(t = 1:n, value = w * x + (1 - w) * z,
           Model = paste0("w=", w))
  }))
)

# ── HC values ─────────────────────────────────────────────────────────────────
HC_df <- bind_rows(
  get_HC(z, "Sine (w=0)"),
  get_HC(x, "AR(2) (w=1)"),
  bind_rows(lapply(weights, function(w) {
    get_HC(w * x + (1 - w) * z, paste0("w=", w))
  }))
)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Colors ────────────────────────────────────────────────────────────────────
mix_colors  <- setNames(viridis::viridis(length(weights)), paste0("w=", weights))
model_colors <- c(
  "Sine (w=0)"     = "blue",
  "AR(2) (w=1)" = "darkred",
  mix_colors
)

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "Sine (w=0)"      = "Sine (w=0)",
  "AR(2) (w=1)" = "AR(2) (w=1)",
  setNames(
    lapply(weights, function(w) bquote(italic(w) == .(w))),
    paste0("w=", weights)
  )
)

# ── Output path ───────────────────────────────────────────────────────────────
output_dir  <- file.path("Results", "Convex_combination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "HC_plane_AR2_Sine_surrounding.pdf")

# ── HC plane (center plot, no legend) ────────────────────────────────────────
p_hc <- ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "solid") +
  geom_point(data = HC_df, aes(x = H, y = C, color = Model), size = 2) +
  scale_color_manual(values = model_colors, labels = legend_labels) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = bquote(italic(H) %*% italic(C) ~ "Plane," ~
                     AR(2) + Sine ~ (italic(D) == .(D))),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p_hc_no_legend <- p_hc + theme(legend.position = "none")

# ── Time series thumbnail function ───────────────────────────────────────────
make_ts_plot <- function(model_label) {
  df     <- filter(ts_df, Model == model_label)
  n_show <- min(500, nrow(df))
  col    <- model_colors[model_label]
  
  ggplot(df[1:n_show, ], aes(x = t, y = value)) +
    geom_line(color = col, linewidth = 0.5) +
    labs(title = model_label, x = NULL, y = NULL) +
    theme_bw(base_family = "serif", base_size = 8) +
    theme(
      plot.title       = element_text(size = 7, face = "bold", hjust = 0.5,
                                      color = col),
      axis.text        = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
      panel.border     = element_rect(color = col, linewidth = 1.2),
      plot.margin      = unit(c(3, 3, 3, 3), "pt")
    )
}

# ── Build all thumbnail plots ─────────────────────────────────────────────────
model_names <- c(
  "Sine (w=0)",
  paste0("w=", weights),
  "AR(2) (w=1)"
)

ts_plots <- setNames(
  lapply(model_names, make_ts_plot),
  model_names
)

# ── Layout: thumbnails surrounding the HC plane ───────────────────────────────
# Top row:    w=0, w=0.1, w=0.2, w=0.3
# Left col:   w=0.4, w=0.5
# Right col:  w=0.6, w=0.7
# Bottom row: w=0.8, w=0.9, w=1
# Center:     HC plane

all_plots        <- ts_plots
all_plots[["HC"]] <- p_hc_no_legend

top_row    <- c("Sine (w=0)", "w=0.1", "w=0.2", "w=0.3")
left_col   <- c("w=0.4", "w=0.5")
right_col  <- c("w=0.6", "w=0.7")
bottom_row <- c("w=0.8", "w=0.9", "AR(2) (w=1)")

layout_matrix <- rbind(
  c( NA,  1,  2,  3, 4),
  c( 5, 12, 12, 12,  11),
  c( 6, 12, 12, 12,  10),
  c(NA, 7, 8, 9, NA)
)

plot_order <- c(top_row, left_col, right_col, bottom_row, "HC")
grob_list  <- lapply(plot_order, function(nm) all_plots[[nm]])

combined_plot <- grid.arrange(
  grobs         = grob_list,
  layout_matrix = layout_matrix,
  widths        = c(1.2, 1.2, 1.2, 1.2, 1.2),
  heights       = c(1.2, 1.2, 1.2, 1.2)
)

ggsave(output_file, plot = combined_plot,
       width = 7, height = 6, dpi = 300, device = "pdf")
# -----------End of code------------------------------------------------
#-----------------------------------------------------------------------


#-------------------------------------------------
# MA2+Sine riding
#------------------------------------

library(tidyverse)
library(StatOrdPattHxC)
library(gridExtra)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
f         <- 0.04
ma2_list <- list(ar = NULL, ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
ma2      <- function(n) as.numeric(normalize(arima.sim(model = ma2_list, n)))
sine      <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x <- ma2(n)
z <- sine(n, f)

# w = 0 is pure Sine, w = 1 is pure MA
weights     <- seq(0.1, 0.9, by = 0.1)
all_weights <- c(0, weights, 1)

# ── Time series data ──────────────────────────────────────────────────────────
ts_df <- bind_rows(
  tibble(t = 1:n, value = z,    Model = "Sine (w=0)"),
  tibble(t = 1:n, value = x,    Model = "MA(2) (w=1)"),
  bind_rows(lapply(weights, function(w) {
    tibble(t = 1:n, value = w * x + (1 - w) * z,
           Model = paste0("w=", w))
  }))
)

# ── HC values ─────────────────────────────────────────────────────────────────
HC_df <- bind_rows(
  get_HC(z, "Sine (w=0)"),
  get_HC(x, "MA(2) (w=1)"),
  bind_rows(lapply(weights, function(w) {
    get_HC(w * x + (1 - w) * z, paste0("w=", w))
  }))
)

# ── Bounds ────────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Colors ────────────────────────────────────────────────────────────────────
mix_colors  <- setNames(viridis::viridis(length(weights)), paste0("w=", weights))
model_colors <- c(
  "Sine (w=0)"     = "blue",
  "MA(2) (w=1)" = "tomato",
  mix_colors
)

# ── Legend labels ─────────────────────────────────────────────────────────────
legend_labels <- c(
  "Sine (w=0)"      = "Sine (w=0)",
  "MA(2) (w=1)" = "MA(2) (w=1)",
  setNames(
    lapply(weights, function(w) bquote(italic(w) == .(w))),
    paste0("w=", weights)
  )
)

# ── Output path ───────────────────────────────────────────────────────────────
output_dir  <- file.path("Results", "Convex_combination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "HC_plane_MA2_Sine_surrounding.pdf")

# ── HC plane (center plot, no legend) ────────────────────────────────────────
p_hc <- ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey60", linetype = "solid") +
  geom_point(data = HC_df, aes(x = H, y = C, color = Model), size = 2) +
  scale_color_manual(values = model_colors, labels = legend_labels) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = bquote(italic(H) %*% italic(C) ~ "Plane," ~
                     MA(2) + Sine ~ (italic(D) == .(D))),
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

p_hc_no_legend <- p_hc + theme(legend.position = "none")

# ── Time series thumbnail function ───────────────────────────────────────────
make_ts_plot <- function(model_label) {
  df     <- filter(ts_df, Model == model_label)
  n_show <- min(500, nrow(df))
  col    <- model_colors[model_label]
  
  ggplot(df[1:n_show, ], aes(x = t, y = value)) +
    geom_line(color = col, linewidth = 0.5) +
    labs(title = model_label, x = NULL, y = NULL) +
    theme_bw(base_family = "serif", base_size = 8) +
    theme(
      plot.title       = element_text(size = 7, face = "bold", hjust = 0.5,
                                      color = col),
      axis.text        = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "gray90"),
      panel.border     = element_rect(color = col, linewidth = 1.2),
      plot.margin      = unit(c(3, 3, 3, 3), "pt")
    )
}

# ── Build all thumbnail plots ─────────────────────────────────────────────────
model_names <- c(
  "Sine (w=0)",
  paste0("w=", weights),
  "MA(2) (w=1)"
)

ts_plots <- setNames(
  lapply(model_names, make_ts_plot),
  model_names
)

# ── Layout: thumbnails surrounding the HC plane ───────────────────────────────
# Top row:    w=0, w=0.1, w=0.2, w=0.3
# Left col:   w=0.4, w=0.5
# Right col:  w=0.6, w=0.7
# Bottom row: w=0.8, w=0.9, w=1
# Center:     HC plane

all_plots        <- ts_plots
all_plots[["HC"]] <- p_hc_no_legend

top_row    <- c("Sine (w=0)", "w=0.1", "w=0.2", "w=0.3")
left_col   <- c("w=0.4", "w=0.5")
right_col  <- c("w=0.6", "w=0.7")
bottom_row <- c("w=0.8", "w=0.9", "MA(2) (w=1)")

layout_matrix <- rbind(
  c( NA,  1,  2,  3, 4),
  c( 5, 12, 12, 12,  11),
  c( 6, 12, 12, 12,  10),
  c(NA, 7, 8, 9, NA)
)

plot_order <- c(top_row, left_col, right_col, bottom_row, "HC")
grob_list  <- lapply(plot_order, function(nm) all_plots[[nm]])

combined_plot <- grid.arrange(
  grobs         = grob_list,
  layout_matrix = layout_matrix,
  widths        = c(1.2, 1.2, 1.2, 1.2, 1.2),
  heights       = c(1.2, 1.2, 1.2, 1.2)
)

ggsave(output_file, plot = combined_plot,
       width = 7, height = 6, dpi = 300, device = "pdf")


# End of the code
#-----------------------------------------------------------------
#---------------------------------------------------------

library(tidyverse)
library(gridExtra)

# ── Parameters ────────────────────────────────────────────────────────────────
D         <- 4
n         <- 1000
f         <- 0.04
arma_list <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_list  <- list(ar = c(-0.8, 0.1))
ma2_list  <- list(ma = c(-0.8, 0.1))

set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))
arma      <- function(n) as.numeric(normalize(arima.sim(model = arma_list, n)))
ar2       <- function(n) as.numeric(normalize(arima.sim(model = ar2_list,  n)))
ma2       <- function(n) as.numeric(normalize(arima.sim(model = ma2_list,  n)))
sine      <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

# ── Signals ───────────────────────────────────────────────────────────────────
x_arma <- arma(n)
x_ar2  <- ar2(n)
x_ma2  <- ma2(n)
z      <- sine(n, f)

weights     <- seq(0.1, 0.9, by = 0.1)
all_weights <- c(0, weights, 1)

# ── Colors ────────────────────────────────────────────────────────────────────
mix_colors <- setNames(viridis::viridis(length(weights)), paste0("w=", weights))

colors_arma <- c("Sine" = "blue", mix_colors, "ARMA(2,2)" = "red")
colors_ar2  <- c("Sine" = "blue", mix_colors, "AR(2)"     = "darkred")
colors_ma2  <- c("Sine" = "blue", mix_colors, "MA(2)"     = "tomato")

# ── Time series data builder ──────────────────────────────────────────────────
make_ts_df <- function(x_stoch, stoch_label) {
  bind_rows(
    tibble(t = 1:n, value = z,       Model = "Sine"),
    bind_rows(lapply(weights, function(w) {
      tibble(t     = 1:n,
             value = w * x_stoch + (1 - w) * z,
             Model = paste0("w=", w))
    })),
    tibble(t = 1:n, value = x_stoch, Model = stoch_label)
  )
}

ts_arma <- make_ts_df(x_arma, "ARMA(2,2)")
ts_ar2  <- make_ts_df(x_ar2,  "AR(2)")
ts_ma2  <- make_ts_df(x_ma2,  "MA(2)")

# ── Combined TS plot function ─────────────────────────────────────────────────
make_combined_ts_plot <- function(ts_df, model_colors, stoch_label, title_expr) {
  n_show  <- min(500, n)
  plot_df <- filter(ts_df, t <= n_show)
  
  # Fix factor order: Sine at bottom, stoch at top of legend
  w_levels <- c("Sine", paste0("w=", weights), stoch_label)
  plot_df  <- mutate(plot_df, Model = factor(Model, levels = w_levels))
  
  ggplot(plot_df, aes(x = t, y = value, color = Model)) +
    geom_line(linewidth = 0.4, alpha = 0.8) +
    scale_color_manual(
      values = model_colors,
      labels = c(
        "Sine" = "Sine",
        setNames(
          lapply(weights, function(w) bquote(italic(w) == .(w))),
          paste0("w=", weights)
        ),
        setNames(list(stoch_label), stoch_label)
      )
    ) +
    labs(
      title = title_expr,
      x     = expression(italic(t)),
      y     = "Value",
      color = expression(italic(w) == "weight")
    ) +
    theme_bw(base_size = 10, base_family = "serif") +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 11),
      legend.position = "right"
    )
}

# ── Build three plots ─────────────────────────────────────────────────────────
p_arma <- make_combined_ts_plot(
  ts_arma, colors_arma, "ARMA(2,2)",
  bquote(ARMA(2*","*2) + Sine ~ "convex combination")
)

p_ar2 <- make_combined_ts_plot(
  ts_ar2, colors_ar2, "AR(2)",
  bquote(AR(2) + Sine ~ "convex combination")
)

p_ma2 <- make_combined_ts_plot(
  ts_ma2, colors_ma2, "MA(2)",
  bquote(MA(2) + Sine ~ "convex combination")
)

# ── Output path ───────────────────────────────────────────────────────────────
output_dir <- file.path("Results", "Convex_combination")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Save individually ─────────────────────────────────────────────────────────
ggsave(file.path(output_dir, "TS_ARMA_Sine.pdf"),
       plot = p_arma, width = 10, height = 4, dpi = 300)

ggsave(file.path(output_dir, "TS_AR2_Sine.pdf"),
       plot = p_ar2, width = 10, height = 4, dpi = 300)

ggsave(file.path(output_dir, "TS_MA2_Sine.pdf"),
       plot = p_ma2, width = 10, height = 4, dpi = 300)

# ── Save as stacked combined figure ──────────────────────────────────────────
ts_combined <- grid.arrange(p_arma, p_ar2, p_ma2, ncol = 1)
ggsave(file.path(output_dir, "TS_all_models_combined.pdf"),
       plot = ts_combined, width = 10, height = 11, dpi = 300)