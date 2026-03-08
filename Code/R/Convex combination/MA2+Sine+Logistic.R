# ══════════════════════════════════════════════════════════════════════════════
#  HC Plane — MA(2), Sine, Logistic, MA(2)+Logistic, MA(2)+Sine
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(StatOrdPattHxC)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
f        <- 0.04          # Sine frequency
r        <- 3.8          # Logistic map parameter (chaos)
ma2_list <- list(ar = NULL, ma = c(-0.8, 0.1))
set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ──────────────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

# MA(2) signal
ma2 <- function(n) as.numeric(normalize(arima.sim(model = ma2_list, n = n)))

# Sine signal (already in [-1,1]; normalize to [0,1] for consistency)
sine <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))

# Logistic map signal
logistic <- function(n, r, x0 = 0.5) {
  x <- numeric(n)
  x[1] <- x0
  for (i in 2:n) x[i] <- r * x[i - 1] * (1 - x[i - 1])
  as.numeric(x)          # already in [0,1]
}

# Compute H and C for a labelled series
get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}

# ── Generate pure signals ─────────────────────────────────────────────────────
x_ma  <- ma2(n)
x_sin <- sine(n, f)
x_log <- logistic(n, r)

# ── Mixing weights ────────────────────────────────────────────────────────────
weights <- seq(0.1, 0.9, by = 0.1)   # w applied to MA(2); (1-w) to other

# ── HC for pure models ────────────────────────────────────────────────────────
pure <- bind_rows(
  get_HC(x_ma,  "MA(2)"),
  get_HC(x_sin, "Sine"),
  get_HC(x_log, "Logistic")
)

# ── HC for MA(2)+Sine mixtures ────────────────────────────────────────────────
mix_sine <- bind_rows(lapply(weights, function(w) {
  get_HC(w * x_ma + (1 - w) * x_sin,
         paste0("MA2+Sine (w=", w, ")"))
}))

# ── HC for MA(2)+Logistic mixtures ────────────────────────────────────────────
mix_log <- bind_rows(lapply(weights, function(w) {
  get_HC(w * x_ma + (1 - w) * x_log,
         paste0("MA2+Log (w=", w, ")"))
}))

# ── Combine all results ───────────────────────────────────────────────────────
results_df <- bind_rows(pure, mix_sine, mix_log)

# ── HC bounds ─────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Color palette ─────────────────────────────────────────────────────────────
#   Pure models    : bold named colours
#   MA2+Sine mix   : viridis gradient (blues/greens)
#   MA2+Logistic   : magma gradient  (oranges/reds)
sine_colors <- setNames(viridis::viridis(9),          paste0("MA2+Sine (w=", weights, ")"))
log_colors  <- setNames(viridis::magma(9, begin = 0.2, end = 0.85),
                        paste0("MA2+Log (w=",  weights, ")"))

all_colors <- c(
  "MA(2)"      = "tomato",
  "Sine"       = "forestgreen",
  "Logistic"   = "dodgerblue",
  sine_colors,
  log_colors
)

# ── Group/shape mapping for legend clarity ─────────────────────────────────────
results_df <- results_df %>%
  mutate(
    Group = case_when(
      Model == "MA(2)"                      ~ "Pure",
      Model == "Sine"                       ~ "Pure",
      Model == "Logistic"                   ~ "Pure",
      str_starts(Model, "MA2\\+Sine")       ~ "MA2+Sine",
      str_starts(Model, "MA2\\+Log")        ~ "MA2+Log"
    ),
    Shape = case_when(
      Group == "Pure"     ~ 18,   # diamond
      Group == "MA2+Sine" ~ 16,   # circle
      Group == "MA2+Log"  ~ 17    # triangle
    )
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  # HC boundary curves
  geom_line(data = bounds,
            aes(x = H, y = C, group = Side),
            color = "grey55", linetype = "dashed", linewidth = 0.5) +
  # mixture paths (lines connecting mix points for each pair)
  geom_path(data = filter(results_df, Group == "MA2+Sine") %>%
              arrange(as.numeric(str_extract(Model, "0\\.\\d+"))),
            aes(x = H, y = C), color = "steelblue", linewidth = 0.4,
            linetype = "dotted") +
  geom_path(data = filter(results_df, Group == "MA2+Log") %>%
              arrange(as.numeric(str_extract(Model, "0\\.\\d+"))),
            aes(x = H, y = C), color = "darkorange", linewidth = 0.4,
            linetype = "dotted") +
  # all points
  geom_point(data = results_df,
             aes(x = H, y = C, color = Model, shape = Group),
             size = 2) +
  # pure-model labels
  ggrepel::geom_label_repel(
    data = filter(results_df, Group == "Pure"),
    aes(x = H, y = C, label = Model, color = Model),
    size = 2, family = "serif", show.legend = FALSE,
    box.padding = 0.4, label.padding = 0.2
  ) +
  scale_color_manual(values = all_colors) +
  scale_shape_manual(values = c("Pure" = 18, "MA2+Sine" = 16, "MA2+Log" = 17)) +
  guides(
    color = guide_legend(nrow =  3, override.aes = list(size = 3)),
   shape = guide_legend(title = "Series type", override.aes = list(size = 2))
  ) +
  labs(
    x     = expression(italic(H)),
    y     = expression(italic(C)),
    title = paste0("HC Plane D = ", D),
    subtitle = "MA(2),  Sine,  Logistic,  MA(2)+Sine, MA(2)+Logistic",
    color = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position  = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold")
  )

print(p)
# ggsave("HC_plane_all_models.pdf", p, width = 9, height = 6)

# End of the code

#-----------------------------------------------------------------------
# ══════════════════════════════════════════════════════════════════════════════
#  HC Plane + Ordinal Patterns + Time Series Thumbnails
#  Models: MA(2), Sine, Logistic, MA(2)+Sine, MA(2)+Logistic
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(StatOrdPattHxC)
library(ggpubr)      # for ggarrange
library(ggrepel)
library(patchwork)   # for combining plots

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
f        <- 0.04
r        <- 3.8
ma2_list <- list(ar = NULL, ma = c(-0.8, 0.1))
set.seed(1234567890, kind = "Mersenne-Twister")

# ── Functions ─────────────────────────────────────────────────────────────────
normalize  <- function(x) (x - min(x)) / (max(x) - min(x))
ma2        <- function(n) as.numeric(normalize(arima.sim(model = ma2_list, n = n)))
sine       <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))
logistic   <- function(n, r, x0 = 0.5) {
  x <- numeric(n); x[1] <- x0
  for (i in 2:n) x[i] <- r * x[i-1] * (1 - x[i-1])
  as.numeric(x)
}
get_HC     <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label)
}
get_probs  <- function(series, label) {
  prob <- as.numeric(OPprob(series, D))
  data.frame(Pattern = seq_along(prob), Prob = prob, Model = label)
}

# ── Signals ───────────────────────────────────────────────────────────────────
x_ma  <- ma2(n)
x_sin <- sine(n, f)
x_log <- logistic(n, r)
weights <- seq(0.1, 0.9, by = 0.1)

# ── Build named list of ALL 21 series ─────────────────────────────────────────
series_list <- c(
  list(
    "MA(2)"     = x_ma,
    "Sine"      = x_sin,
    "Logistic"  = x_log
  ),
  setNames(
    lapply(weights, function(w) w * x_ma + (1 - w) * x_sin),
    paste0("MA2+Sine\nw=", weights)
  ),
  setNames(
    lapply(weights, function(w) w * x_ma + (1 - w) * x_log),
    paste0("MA2+Log\nw=", weights)
  )
)

# ── HC results ────────────────────────────────────────────────────────────────
results_df <- bind_rows(mapply(get_HC, series_list, names(series_list),
                               SIMPLIFY = FALSE))

# Add group column for coloring
results_df <- results_df %>%
  mutate(Group = case_when(
    Model %in% c("MA(2)", "Sine", "Logistic") ~ "Pure",
    str_starts(Model, "MA2\\+Sine")            ~ "MA2+Sine",
    str_starts(Model, "MA2\\+Log")             ~ "MA2+Log"
  ))

# ── Ordinal pattern probabilities ─────────────────────────────────────────────
# Only pure + representative mixes for clarity; or all 5 model types
op_models <- list(
  "MA(2)"    = x_ma,
  "Sine"     = x_sin,
  "Logistic" = x_log,
  "MA2+Sine\n(w=0.5)"  = 0.5 * x_ma + 0.5 * x_sin,
  "MA2+Log\n(w=0.5)"   = 0.5 * x_ma + 0.5 * x_log
)
op_df <- bind_rows(mapply(get_probs, op_models, names(op_models),
                          SIMPLIFY = FALSE)) %>%
  mutate(Model = factor(Model, levels = names(op_models)),
         Pattern = factor(Pattern))

# ── Color palettes ────────────────────────────────────────────────────────────
sine_colors <- setNames(viridis::viridis(9),
                        paste0("MA2+Sine\nw=", weights))
log_colors  <- setNames(viridis::magma(9, begin = 0.2, end = 0.85),
                        paste0("MA2+Log\nw=", weights))
all_colors  <- c("MA(2)" = "tomato", "Sine" = "forestgreen",
                 "Logistic" = "dodgerblue", sine_colors, log_colors)

pure_5_colors <- c(
  "MA(2)"             = "tomato",
  "Sine"              = "forestgreen",
  "Logistic"          = "dodgerblue",
  "MA2+Sine\n(w=0.5)" = "#35B779",
  "MA2+Log\n(w=0.5)"  = "#E07B39"
)

# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 1 — HC Plane with time series thumbnails
# ══════════════════════════════════════════════════════════════════════════════
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

hc_main <- ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey55", linetype = "dashed", linewidth = 0.5) +
  geom_path(data = filter(results_df, Group == "MA2+Sine") %>%
              arrange(as.numeric(str_extract(Model, "0\\.\\d+"))),
            aes(x = H, y = C), color = "steelblue",
            linewidth = 0.4, linetype = "dotted") +
  geom_path(data = filter(results_df, Group == "MA2+Log") %>%
              arrange(as.numeric(str_extract(Model, "0\\.\\d+"))),
            aes(x = H, y = C), color = "darkorange",
            linewidth = 0.4, linetype = "dotted") +
  geom_point(data = results_df,
             aes(x = H, y = C, color = Model), size = 2.5) +
  geom_label_repel(
    data  = filter(results_df, Group == "Pure"),
    aes(x = H, y = C, label = Model, color = Model),
    size  = 3, family = "serif", show.legend = FALSE,
    box.padding = 0.4, label.padding = 0.2
  ) +
  scale_color_manual(values = all_colors) +
  guides(color = guide_legend(nrow = 4, override.aes = list(size = 3))) +
  labs(
    x        = expression(italic(H)),
    y        = expression(italic(C)),
    title    = paste0("HC Plane  —  D = ", D),
    subtitle = "MA(2),  Sine,  Logistic,  MA(2)+Sine,  MA(2)+Logistic",
    color    = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title     = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle  = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.text    = element_text(size = 7),
    legend.title   = element_text(size = 8, face = "bold")
  )

# ── Thumbnail helper ───────────────────────────────────────────────────────────
make_thumb <- function(series, label, col) {
  df <- data.frame(t = 1:200, y = series[1:200])
  ggplot(df, aes(t, y)) +
    geom_line(color = col, linewidth = 0.3) +
    labs(title = label) +
    theme_void(base_size = 6, base_family = "serif") +
    theme(
      plot.title      = element_text(size = 5, hjust = 0.5,
                                     margin = margin(b = 1)),
      plot.background = element_rect(fill = "white", color = "grey80",
                                     linewidth = 0.3),
      plot.margin     = unit(c(1,1,1,1), "pt")
    )
}

# Assign a color to every series
thumb_colors <- all_colors
thumb_colors["Logistic"] <- "dodgerblue"

# Build all 21 thumbnails
thumb_plots <- mapply(
  function(s, nm) make_thumb(s, nm, thumb_colors[nm]),
  series_list, names(series_list), SIMPLIFY = FALSE
)

# ── Arrange: HC centre + thumbnails around edges ───────────────────────────────
# Layout: 3 pure on top, 9 MA2+Sine on left, 9 MA2+Log on right, HC in centre
# We use patchwork with a manual layout

# Top row: 3 pure model thumbnails + spacers
pure_names    <- c("MA(2)", "Sine", "Logistic")
sine_names    <- paste0("MA2+Sine\nw=", weights)
log_names     <- paste0("MA2+Log\nw=", weights)

# Build composite using patchwork inset_element on a base HC plot
# Strategy: wrap all thumbs into a grid around the HC plot using cowplot

library(cowplot)

# Place HC plot as base canvas
base <- ggdraw() +
  draw_plot(hc_main, x = 0.18, y = 0.08, width = 0.64, height = 0.82)

# ── Thumbnail grid positions ───────────────────────────────────────────────────
# Top (pure):     evenly spaced across top
# Left  (Sine mix w=0.1..0.9): stacked vertically on left
# Right (Log mix w=0.1..0.9):  stacked vertically on right

tw <- 0.13; th <- 0.08   # thumb width / height

# Top: 3 pure models
top_x  <- c(0.25, 0.44, 0.63)
top_y  <- 0.92
for (i in seq_along(pure_names)) {
  base <- base + draw_plot(thumb_plots[[pure_names[i]]],
                           x = top_x[i], y = top_y, width = tw, height = th)
}

# Left: MA2+Sine w=0.1 (top) → w=0.9 (bottom)
left_x <- 0.01
left_y <- seq(0.83, 0.08, length.out = 9)
for (i in 1:9) {
  base <- base + draw_plot(thumb_plots[[sine_names[i]]],
                           x = left_x, y = left_y[i], width = tw, height = th)
}

# Right: MA2+Log w=0.1 (top) → w=0.9 (bottom)
right_x <- 0.86
right_y <- seq(0.83, 0.08, length.out = 9)
for (i in 1:9) {
  base <- base + draw_plot(thumb_plots[[log_names[i]]],
                           x = right_x, y = right_y[i], width = tw, height = th)
}

# ── Save Plot 1 ───────────────────────────────────────────────────────────────
ggsave("HC_plane_with_thumbnails.pdf", base,
       width = 14, height = 10, device = cairo_pdf)
ggsave("HC_plane_with_thumbnails.png", base,
       width = 14, height = 10, dpi = 200)
message("✓ Saved HC_plane_with_thumbnails.pdf/.png")

# ══════════════════════════════════════════════════════════════════════════════
#  PLOT 2 — Ordinal Patterns Bar Chart (all 5 models, 24 patterns)
# ══════════════════════════════════════════════════════════════════════════════
op_plot <- ggplot(op_df, aes(x = Pattern, y = Prob, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.75) +
  geom_hline(yintercept = 1/factorial(D), linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = 24.5, y = 1/factorial(D) + 0.002,
           label = "Uniform", size = 2.8, hjust = 1,
           color = "grey40", family = "serif") +
  scale_fill_manual(values = pure_5_colors) +
  scale_x_discrete(breaks = as.character(seq(1, 24, by = 2))) +
  labs(
    x     = "Ordinal Pattern Index",
    y     = "Probability",
    title = paste0("Ordinal Pattern Probabilities  —  D = ", D),
    subtitle = "Dashed line = uniform distribution (1/4! = 1/24)",
    fill  = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "lines"),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    axis.text.x  = element_text(size = 8)
  )

ggsave("ordinal_patterns_barchart.pdf", op_plot,
       width = 10, height = 5, device = cairo_pdf)
ggsave("ordinal_patterns_barchart.png", op_plot,
       width = 10, height = 5, dpi = 200)
message("✓ Saved ordinal_patterns_barchart.pdf/.png")