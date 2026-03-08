# ══════════════════════════════════════════════════════════════════════════════
#  HC Plane + Time Series Thumbnails
#  37 Models: Pure (ARMA22, AR2, MA2, Logistic, Sine) + Mixtures
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(StatOrdPattHxC)
library(cowplot)
library(ggrepel)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
f        <- 0.04
r        <- 3.8
arma22_model <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_model    <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_model    <- list(ar = NULL,         ma = c(-0.8, 0.1))
set.seed(1234567890, kind = "Mersenne-Twister")

# ── Generator functions ───────────────────────────────────────────────────────
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

gen_arma22  <- function(n) as.numeric(normalize(arima.sim(model = arma22_model, n = n)))
gen_ar2     <- function(n) as.numeric(normalize(arima.sim(model = ar2_model,    n = n)))
gen_ma2     <- function(n) as.numeric(normalize(arima.sim(model = ma2_model,    n = n)))
gen_sine    <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))
gen_logistic <- function(n, r, x0 = 0.5) {
  x <- numeric(n); x[1] <- x0
  for (i in 2:n) x[i] <- r * x[i-1] * (1 - x[i-1])
  as.numeric(x)
}

# ── Pure signals ──────────────────────────────────────────────────────────────
x_arma <- gen_arma22(n)
x_ar2  <- gen_ar2(n)
x_ma2  <- gen_ma2(n)
x_sin  <- gen_sine(n, f)
x_log  <- gen_logistic(n, r)

# ── Build full series list (37 models) ────────────────────────────────────────
series_list <- list(
  # Pure
  "ARMA(2,2)"            = x_arma,
  "AR(2)"                = x_ar2,
  "MA(2)"                = x_ma2,
  "Logistic"             = x_log,
  "Sine"                 = x_sin,
  # ARMA + Logistic
  "ARMA+Logistic(w=0.1)" = 0.1*x_arma + 0.9*x_log,
  "ARMA+Logistic(w=0.2)" = 0.2*x_arma + 0.8*x_log,
  "ARMA+Logistic(w=0.4)" = 0.4*x_arma + 0.6*x_log,
  "ARMA+Logistic(w=0.7)" = 0.7*x_arma + 0.3*x_log,
  "ARMA+Logistic(w=0.9)" = 0.9*x_arma + 0.1*x_log,
  # ARMA + Sine
  "ARMA+Sine(w=0.1)"     = 0.1*x_arma + 0.9*x_sin,
  "ARMA+Sine(w=0.2)"     = 0.2*x_arma + 0.8*x_sin,
  "ARMA+Sine(w=0.3)"     = 0.3*x_arma + 0.7*x_sin,
  "ARMA+Sine(w=0.4)"     = 0.4*x_arma + 0.6*x_sin,
  "ARMA+Sine(w=0.7)"     = 0.7*x_arma + 0.3*x_sin,
  "ARMA+Sine(w=0.9)"     = 0.9*x_arma + 0.1*x_sin,
  # AR2 + Logistic
  "AR2+Logistic(w=0.2)"  = 0.2*x_ar2  + 0.8*x_log,
  "AR2+Logistic(w=0.3)"  = 0.3*x_ar2  + 0.7*x_log,
  "AR2+Logistic(w=0.4)"  = 0.4*x_ar2  + 0.6*x_log,
  "AR2+Logistic(w=0.7)"  = 0.7*x_ar2  + 0.3*x_log,
  "AR2+Logistic(w=0.8)"  = 0.8*x_ar2  + 0.2*x_log,
  "AR2+Logistic(w=0.9)"  = 0.9*x_ar2  + 0.1*x_log,
  # AR2 + Sine
  "AR2+Sine(w=0.2)"      = 0.2*x_ar2  + 0.8*x_sin,
  "AR2+Sine(w=0.3)"      = 0.3*x_ar2  + 0.7*x_sin,
  "AR2+Sine(w=0.4)"      = 0.4*x_ar2  + 0.6*x_sin,
  "AR2+Sine(w=0.7)"      = 0.7*x_ar2  + 0.3*x_sin,
  # MA2 + Logistic
  "MA2+Logistic(w=0.1)"  = 0.1*x_ma2  + 0.9*x_log,
  "MA2+Logistic(w=0.2)"  = 0.2*x_ma2  + 0.8*x_log,
  "MA2+Logistic(w=0.3)"  = 0.3*x_ma2  + 0.7*x_log,
  "MA2+Logistic(w=0.4)"  = 0.4*x_ma2  + 0.6*x_log,
  "MA2+Logistic(w=0.5)"  = 0.5*x_ma2  + 0.5*x_log,
  "MA2+Logistic(w=0.7)"  = 0.7*x_ma2  + 0.3*x_log,
  # MA2 + Sine
  "MA2+Sine(w=0.1)"      = 0.1*x_ma2  + 0.9*x_sin,
  "MA2+Sine(w=0.2)"      = 0.2*x_ma2  + 0.8*x_sin,
  "MA2+Sine(w=0.3)"      = 0.3*x_ma2  + 0.7*x_sin,
  "MA2+Sine(w=0.4)"      = 0.4*x_ma2  + 0.6*x_sin,
  "MA2+Sine(w=0.5)"      = 0.5*x_ma2  + 0.5*x_sin
)

# ── HC computation ────────────────────────────────────────────────────────────
get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob), Model = label,
             stringsAsFactors = FALSE)
}

results_df <- bind_rows(mapply(get_HC, series_list, names(series_list),
                               SIMPLIFY = FALSE))

# ── Group classification ──────────────────────────────────────────────────────
results_df <- results_df %>%
  mutate(Group = case_when(
    Model %in% c("ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine") ~ "Pure",
    str_detect(Model, "ARMA\\+Logistic") ~ "ARMA+Logistic",
    str_detect(Model, "ARMA\\+Sine")     ~ "ARMA+Sine",
    str_detect(Model, "AR2\\+Logistic")  ~ "AR2+Logistic",
    str_detect(Model, "AR2\\+Sine")      ~ "AR2+Sine",
    str_detect(Model, "MA2\\+Logistic")  ~ "MA2+Logistic",
    str_detect(Model, "MA2\\+Sine")      ~ "MA2+Sine"
  ))

# ── Color palette ─────────────────────────────────────────────────────────────
group_colors <- c(
  "Pure"           = NA,   # handled individually below
  "ARMA+Logistic"  = "#E07B39",
  "ARMA+Sine"      = "#3B9AB2",
  "AR2+Logistic"   = "#D4A017",
  "AR2+Sine"       = "#7EC8A4",
  "MA2+Logistic"   = "#C0392B",
  "MA2+Sine"       = "#8E44AD"
)

pure_colors <- c(
  "ARMA(2,2)" = "black",
  "AR(2)"     = "tomato",
  "MA(2)"     = "navy",
  "Logistic"  = "dodgerblue",
  "Sine"      = "forestgreen"
)

# Assign per-model color based on group (gradient within group by w)
make_group_gradient <- function(models, base_color, n) {
  cols <- colorRampPalette(c(lighten(base_color, 0.5), base_color))(n)
  setNames(cols, models)
}

lighten <- function(col, factor = 0.5) {
  rgb_vals <- col2rgb(col) / 255
  rgb((rgb_vals + factor) / (1 + factor),
      (rgb_vals[2,] + factor) / (1 + factor),
      (rgb_vals[3,] + factor) / (1 + factor))
}

# Build full color vector
model_colors <- c(
  pure_colors,
  setNames(
    colorRampPalette(c("#F5C48A","#E07B39"))(5),
    c("ARMA+Logistic(w=0.1)","ARMA+Logistic(w=0.2)","ARMA+Logistic(w=0.4)",
      "ARMA+Logistic(w=0.7)","ARMA+Logistic(w=0.9)")
  ),
  setNames(
    colorRampPalette(c("#A8D8EA","#3B9AB2"))(6),
    c("ARMA+Sine(w=0.1)","ARMA+Sine(w=0.2)","ARMA+Sine(w=0.3)",
      "ARMA+Sine(w=0.4)","ARMA+Sine(w=0.7)","ARMA+Sine(w=0.9)")
  ),
  setNames(
    colorRampPalette(c("#F7E08A","#D4A017"))(6),
    c("AR2+Logistic(w=0.2)","AR2+Logistic(w=0.3)","AR2+Logistic(w=0.4)",
      "AR2+Logistic(w=0.7)","AR2+Logistic(w=0.8)","AR2+Logistic(w=0.9)")
  ),
  setNames(
    colorRampPalette(c("#C8EDD8","#7EC8A4"))(4),
    c("AR2+Sine(w=0.2)","AR2+Sine(w=0.3)","AR2+Sine(w=0.4)","AR2+Sine(w=0.7)")
  ),
  setNames(
    colorRampPalette(c("#F5AAAA","#C0392B"))(6),
    c("MA2+Logistic(w=0.1)","MA2+Logistic(w=0.2)","MA2+Logistic(w=0.3)",
      "MA2+Logistic(w=0.4)","MA2+Logistic(w=0.5)","MA2+Logistic(w=0.7)")
  ),
  setNames(
    colorRampPalette(c("#D9B8E8","#8E44AD"))(5),
    c("MA2+Sine(w=0.1)","MA2+Sine(w=0.2)","MA2+Sine(w=0.3)",
      "MA2+Sine(w=0.4)","MA2+Sine(w=0.5)")
  )
)

# ── HC Bounds ─────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Path data per mixture group ───────────────────────────────────────────────
path_df <- results_df %>%
  filter(Group != "Pure") %>%
  mutate(w = as.numeric(str_extract(Model, "0\\.\\d+"))) %>%
  arrange(Group, w)

# ── HC Main Plot ──────────────────────────────────────────────────────────────
hc_main <- ggplot() +
  geom_line(data = bounds, aes(x = H, y = C, group = Side),
            color = "grey55", linetype = "dashed", linewidth = 0.5) +
  # Mixture paths
  geom_path(data = path_df,
            aes(x = H, y = C, group = Group, color = Group),
            linewidth = 0.4, linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = c(
    "ARMA+Logistic" = "#E07B39", "ARMA+Sine"    = "#3B9AB2",
    "AR2+Logistic"  = "#D4A017", "AR2+Sine"     = "#7EC8A4",
    "MA2+Logistic"  = "#C0392B", "MA2+Sine"     = "#8E44AD"
  ), guide = "none") +
  # All points
  ggnewscale::new_scale_color() +
  geom_point(data = results_df,
             aes(x = H, y = C, color = Model,
                 size  = ifelse(Group == "Pure", 4, 3.5),
                 shape = ifelse(Group == "Pure", 18, 16)),
             show.legend = TRUE) +
  scale_color_manual(values = model_colors) +
  scale_size_identity() +
  scale_shape_identity() +
  # Pure model labels
  geom_label_repel(
    data = filter(results_df, Group == "Pure"),
    aes(x = H, y = C, label = Model, color = Model),
    size = 3, family = "serif", show.legend = FALSE,
    box.padding = 0.5, label.padding = 0.2, fontface = "bold"
  ) +
  guides(color = guide_legend(nrow = 4, override.aes = list(size = 3))) +
  labs(
    x        = expression(italic(H)),
    y        = expression(italic(C)),
    title    = paste0("HC Plane  —  D = ", D),
    subtitle = "Convex Combination model",
    color    = "Model"
  ) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "lines"),
    legend.text     = element_text(size = 6),
    legend.title    = element_text(size = 7, face = "bold")
  )

# ══════════════════════════════════════════════════════════════════════════════
#  TIME SERIES THUMBNAILS around HC plane
# ══════════════════════════════════════════════════════════════════════════════

make_thumb <- function(nm) {
  s   <- series_list[[nm]]
  col <- model_colors[nm]
  df  <- data.frame(t = 1:200, y = s[1:200])
  # short display label
  lbl <- nm %>%
    str_replace("ARMA\\(2,2\\)", "ARMA22") %>%
    str_replace("\\(w=", "\nw=") %>%
    str_replace("\\)", "")
  ggplot(df, aes(t, y)) +
    geom_line(color = col, linewidth = 0.35) +
    labs(title = lbl) +
    theme_void(base_size = 5, base_family = "serif") +
    theme(
      plot.title      = element_text(size = 4.5, hjust = 0.5,
                                     margin = margin(b = 1),
                                     color = col, face = "bold"),
      plot.background = element_rect(fill = "white", color = col,
                                     linewidth = 0.4),
      plot.margin     = unit(c(2,2,2,2), "pt")
    )
}

# ── Layout strategy ───────────────────────────────────────────────────────────
# Top row    (8):  5 pure + first 3 of ARMA+Logistic
# Bottom row (8):  last 8 mixtures (MA2+Sine w=0.1..0.5 + AR2+Sine w=0.7 + 2 more)
# Left col   (11): ARMA+Logistic(w=0.4..0.9)[2] + ARMA+Sine(6) + AR2+Logistic(5 top)
# Right col  (11): AR2+Logistic(w=0.8,0.9) + AR2+Sine(3) + MA2+Logistic(6)
# ─────────────────────────────────────────────────────────────────────────────

all_names <- names(series_list)

pure_names     <- c("ARMA(2,2)","AR(2)","MA(2)","Logistic","Sine")
arma_log_names <- c("ARMA+Logistic(w=0.1)","ARMA+Logistic(w=0.2)",
                    "ARMA+Logistic(w=0.4)","ARMA+Logistic(w=0.7)",
                    "ARMA+Logistic(w=0.9)")
arma_sin_names <- c("ARMA+Sine(w=0.1)","ARMA+Sine(w=0.2)","ARMA+Sine(w=0.3)",
                    "ARMA+Sine(w=0.4)","ARMA+Sine(w=0.7)","ARMA+Sine(w=0.9)")
ar2_log_names  <- c("AR2+Logistic(w=0.2)","AR2+Logistic(w=0.3)",
                    "AR2+Logistic(w=0.4)","AR2+Logistic(w=0.7)",
                    "AR2+Logistic(w=0.8)","AR2+Logistic(w=0.9)")
ar2_sin_names  <- c("AR2+Sine(w=0.2)","AR2+Sine(w=0.3)",
                    "AR2+Sine(w=0.4)","AR2+Sine(w=0.7)")
ma2_log_names  <- c("MA2+Logistic(w=0.1)","MA2+Logistic(w=0.2)",
                    "MA2+Logistic(w=0.3)","MA2+Logistic(w=0.4)",
                    "MA2+Logistic(w=0.5)","MA2+Logistic(w=0.7)")
ma2_sin_names  <- c("MA2+Sine(w=0.1)","MA2+Sine(w=0.2)","MA2+Sine(w=0.3)",
                    "MA2+Sine(w=0.4)","MA2+Sine(w=0.5)")

# ── Row / column assignments ──────────────────────────────────────────────────
# Top row (8):    5 pure + first 3 ARMA+Logistic
top_names    <- c(pure_names, arma_log_names[1:3])

# Bottom row (8): last 8 of all remaining after top + left + right
# Left (11):  ARMA+Logistic[4:5] + ARMA+Sine[1:6] + AR2+Logistic[1:3]
left_names   <- c(arma_log_names[4:5], arma_sin_names, ar2_log_names[1:3])  # 2+6+3 = 11

# Right (11): AR2+Logistic[4:6] + AR2+Sine[1:4] + MA2+Logistic[1:4]
right_names  <- c(ar2_log_names[4:6], ar2_sin_names, ma2_log_names[1:4])    # 3+4+4 = 11

# Bottom (8): MA2+Logistic[5:6] + MA2+Sine[1:5] + 1 spare = MA2+Log(w=0.5,0.7) + all MA2+Sine
bottom_names <- c(ma2_log_names[5:6], ma2_sin_names)                        # 2+5 = 7
# pad to 8 with AR2+Sine(w=0.7) moved here if needed — keep as 7, evenly spaced
bottom_names <- c(ar2_sin_names[4], ma2_log_names[5:6], ma2_sin_names)      # 1+2+5 = 8

# Pre-build all thumbnails
thumbs <- setNames(lapply(all_names, make_thumb), all_names)

# ── Canvas geometry ───────────────────────────────────────────────────────────
# Total canvas: width=24, height=16 (inches) — set in ggsave
# HC plot occupies centre region in normalised [0,1] coords:
#   x: [0.14, 0.86]   y: [0.13, 0.88]
# Thumb size:  tw=0.092 (width), th=0.068 (height) in normalised units

tw <- 0.092; th <- 0.068
hc_x0 <- 0.14; hc_x1 <- 0.86
hc_y0 <- 0.13; hc_y1 <- 0.88

base <- ggdraw() +
  draw_plot(hc_main,
            x = hc_x0, y = hc_y0,
            width  = hc_x1 - hc_x0,
            height = hc_y1 - hc_y0)

# ── Top row: 8 thumbnails evenly across top ───────────────────────────────────
top_xs <- seq(hc_x0, hc_x1 - tw, length.out = 8)
top_y  <- hc_y1 + 0.004
for (i in seq_along(top_names)) {
  base <- base + draw_plot(thumbs[[top_names[i]]],
                           x = top_xs[i], y = top_y,
                           width = tw, height = th)
}

# ── Bottom row: 8 thumbnails evenly across bottom ────────────────────────────
bot_xs <- seq(hc_x0, hc_x1 - tw, length.out = 8)
bot_y  <- hc_y0 - th - 0.004
for (i in seq_along(bottom_names)) {
  base <- base + draw_plot(thumbs[[bottom_names[i]]],
                           x = bot_xs[i], y = bot_y,
                           width = tw, height = th)
}

# ── Left col: 11 thumbnails stacked top→bottom ───────────────────────────────
left_x  <- hc_x0 - tw - 0.004
left_ys <- seq(hc_y1 - th, hc_y0, length.out = 11)
for (i in seq_along(left_names)) {
  base <- base + draw_plot(thumbs[[left_names[i]]],
                           x = left_x, y = left_ys[i],
                           width = tw, height = th)
}

# ── Right col: 11 thumbnails stacked top→bottom ──────────────────────────────
right_x  <- hc_x1 + 0.004
right_ys <- seq(hc_y1 - th, hc_y0, length.out = 11)
for (i in seq_along(right_names)) {
  base <- base + draw_plot(thumbs[[right_names[i]]],
                           x = right_x, y = right_ys[i],
                           width = tw, height = th)
}

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave("HC_37models_thumbnails.pdf", base,
       width = 24, height = 16, device = cairo_pdf)
ggsave("HC_37models_thumbnails.png", base,
       width = 24, height = 16, dpi = 180)
message("✓ Saved HC_37models_thumbnails.pdf/.png")

# End of the code
#=============================================================================================
#=============================================================================================
# ══════════════════════════════════════════════════════════════════════════════
#  HC Plane — Scatter Plot Only
#  37 Models: ARMA(2,2), AR(2), MA(2), Logistic, Sine + Mixtures
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(StatOrdPattHxC)
library(ggrepel)

# ── Parameters ────────────────────────────────────────────────────────────────
D        <- 4
n        <- 1000
f        <- 0.04
r        <- 3.8
arma22_model <- list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1))
ar2_model    <- list(ar = c(-0.8, 0.1), ma = NULL)
ma2_model    <- list(ar = NULL,         ma = c(-0.8, 0.1))
set.seed(1234567890, kind = "Mersenne-Twister")

# ── Generator functions ───────────────────────────────────────────────────────
normalize    <- function(x) (x - min(x)) / (max(x) - min(x))
gen_arma22   <- function(n) as.numeric(normalize(arima.sim(model = arma22_model, n = n)))
gen_ar2      <- function(n) as.numeric(normalize(arima.sim(model = ar2_model,    n = n)))
gen_ma2      <- function(n) as.numeric(normalize(arima.sim(model = ma2_model,    n = n)))
gen_sine     <- function(n, f) as.numeric(normalize(sin(2 * pi * f * 1:n)))
gen_logistic <- function(n, r, x0 = 0.5) {
  x <- numeric(n); x[1] <- x0
  for (i in 2:n) x[i] <- r * x[i-1] * (1 - x[i-1])
  as.numeric(x)
}

# ── Pure signals ──────────────────────────────────────────────────────────────
x_arma <- gen_arma22(n)
x_ar2  <- gen_ar2(n)
x_ma2  <- gen_ma2(n)
x_sin  <- gen_sine(n, f)
x_log  <- gen_logistic(n, r)

# ── Full series list (37 models) ──────────────────────────────────────────────
series_list <- list(
  # Pure
  "ARMA(2,2)"            = x_arma,
  "AR(2)"                = x_ar2,
  "MA(2)"                = x_ma2,
  "Logistic"             = x_log,
  "Sine"                 = x_sin,
  # ARMA + Logistic
  "ARMA+Logistic(w=0.1)" = 0.1*x_arma + 0.9*x_log,
  "ARMA+Logistic(w=0.2)" = 0.2*x_arma + 0.8*x_log,
  "ARMA+Logistic(w=0.4)" = 0.4*x_arma + 0.6*x_log,
  "ARMA+Logistic(w=0.7)" = 0.7*x_arma + 0.3*x_log,
  "ARMA+Logistic(w=0.9)" = 0.9*x_arma + 0.1*x_log,
  # ARMA + Sine
  "ARMA+Sine(w=0.1)"     = 0.1*x_arma + 0.9*x_sin,
  "ARMA+Sine(w=0.2)"     = 0.2*x_arma + 0.8*x_sin,
  "ARMA+Sine(w=0.3)"     = 0.3*x_arma + 0.7*x_sin,
  "ARMA+Sine(w=0.4)"     = 0.4*x_arma + 0.6*x_sin,
  "ARMA+Sine(w=0.7)"     = 0.7*x_arma + 0.3*x_sin,
  "ARMA+Sine(w=0.9)"     = 0.9*x_arma + 0.1*x_sin,
  # AR2 + Logistic
  "AR2+Logistic(w=0.2)"  = 0.2*x_ar2  + 0.8*x_log,
  "AR2+Logistic(w=0.3)"  = 0.3*x_ar2  + 0.7*x_log,
  "AR2+Logistic(w=0.4)"  = 0.4*x_ar2  + 0.6*x_log,
  "AR2+Logistic(w=0.7)"  = 0.7*x_ar2  + 0.3*x_log,
  "AR2+Logistic(w=0.8)"  = 0.8*x_ar2  + 0.2*x_log,
  "AR2+Logistic(w=0.9)"  = 0.9*x_ar2  + 0.1*x_log,
  # AR2 + Sine
  "AR2+Sine(w=0.2)"      = 0.2*x_ar2  + 0.8*x_sin,
  "AR2+Sine(w=0.3)"      = 0.3*x_ar2  + 0.7*x_sin,
  "AR2+Sine(w=0.4)"      = 0.4*x_ar2  + 0.6*x_sin,
  "AR2+Sine(w=0.7)"      = 0.7*x_ar2  + 0.3*x_sin,
  # MA2 + Logistic
  "MA2+Logistic(w=0.1)"  = 0.1*x_ma2  + 0.9*x_log,
  "MA2+Logistic(w=0.2)"  = 0.2*x_ma2  + 0.8*x_log,
  "MA2+Logistic(w=0.3)"  = 0.3*x_ma2  + 0.7*x_log,
  "MA2+Logistic(w=0.4)"  = 0.4*x_ma2  + 0.6*x_log,
  "MA2+Logistic(w=0.5)"  = 0.5*x_ma2  + 0.5*x_log,
  "MA2+Logistic(w=0.7)"  = 0.7*x_ma2  + 0.3*x_log,
  # MA2 + Sine
  "MA2+Sine(w=0.1)"      = 0.1*x_ma2  + 0.9*x_sin,
  "MA2+Sine(w=0.2)"      = 0.2*x_ma2  + 0.8*x_sin,
  "MA2+Sine(w=0.3)"      = 0.3*x_ma2  + 0.7*x_sin,
  "MA2+Sine(w=0.4)"      = 0.4*x_ma2  + 0.6*x_sin,
  "MA2+Sine(w=0.5)"      = 0.5*x_ma2  + 0.5*x_sin
)

# ── HC computation ────────────────────────────────────────────────────────────
get_HC <- function(series, label) {
  prob <- OPprob(series, D)
  data.frame(H = HShannon(prob), C = StatComplexity(prob),
             Model = label, stringsAsFactors = FALSE)
}

results_df <- bind_rows(mapply(get_HC, series_list, names(series_list),
                               SIMPLIFY = FALSE))

# ── Group classification ──────────────────────────────────────────────────────
results_df <- results_df %>%
  mutate(Group = case_when(
    Model %in% c("ARMA(2,2)","AR(2)","MA(2)","Logistic","Sine") ~ "Pure",
    str_detect(Model, "ARMA\\+Logistic") ~ "ARMA+Logistic",
    str_detect(Model, "ARMA\\+Sine")     ~ "ARMA+Sine",
    str_detect(Model, "AR2\\+Logistic")  ~ "AR2+Logistic",
    str_detect(Model, "AR2\\+Sine")      ~ "AR2+Sine",
    str_detect(Model, "MA2\\+Logistic")  ~ "MA2+Logistic",
    str_detect(Model, "MA2\\+Sine")      ~ "MA2+Sine"
  ))

# ── Color palette ─────────────────────────────────────────────────────────────
model_colors <- c(
  # Pure
  "ARMA(2,2)" = "black",
  "AR(2)"     = "tomato",
  "MA(2)"     = "navy",
  "Logistic"  = "dodgerblue",
  "Sine"      = "forestgreen",
  # ARMA+Logistic: orange gradient
  setNames(colorRampPalette(c("#F5C48A","#E07B39"))(5),
           c("ARMA+Logistic(w=0.1)","ARMA+Logistic(w=0.2)",
             "ARMA+Logistic(w=0.4)","ARMA+Logistic(w=0.7)",
             "ARMA+Logistic(w=0.9)")),
  # ARMA+Sine: blue gradient
  setNames(colorRampPalette(c("#A8D8EA","#3B9AB2"))(6),
           c("ARMA+Sine(w=0.1)","ARMA+Sine(w=0.2)","ARMA+Sine(w=0.3)",
             "ARMA+Sine(w=0.4)","ARMA+Sine(w=0.7)","ARMA+Sine(w=0.9)")),
  # AR2+Logistic: gold gradient
  setNames(colorRampPalette(c("#F7E08A","#D4A017"))(6),
           c("AR2+Logistic(w=0.2)","AR2+Logistic(w=0.3)","AR2+Logistic(w=0.4)",
             "AR2+Logistic(w=0.7)","AR2+Logistic(w=0.8)","AR2+Logistic(w=0.9)")),
  # AR2+Sine: teal gradient
  setNames(colorRampPalette(c("#C8EDD8","#7EC8A4"))(4),
           c("AR2+Sine(w=0.2)","AR2+Sine(w=0.3)",
             "AR2+Sine(w=0.4)","AR2+Sine(w=0.7)")),
  # MA2+Logistic: red gradient
  setNames(colorRampPalette(c("#F5AAAA","#C0392B"))(6),
           c("MA2+Logistic(w=0.1)","MA2+Logistic(w=0.2)","MA2+Logistic(w=0.3)",
             "MA2+Logistic(w=0.4)","MA2+Logistic(w=0.5)","MA2+Logistic(w=0.7)")),
  # MA2+Sine: purple gradient
  setNames(colorRampPalette(c("#D9B8E8","#8E44AD"))(5),
           c("MA2+Sine(w=0.1)","MA2+Sine(w=0.2)","MA2+Sine(w=0.3)",
             "MA2+Sine(w=0.4)","MA2+Sine(w=0.5)"))
)

# ── HC Bounds ─────────────────────────────────────────────────────────────────
data("LinfLsup")
bounds <- filter(LinfLsup, Dimension == as.character(D))

# ── Path data (dotted lines connecting mixture points per group) ───────────────
path_df <- results_df %>%
  filter(Group != "Pure") %>%
  mutate(w = as.numeric(str_extract(Model, "0\\.\\d+"))) %>%
  arrange(Group, w)

group_path_colors <- c(
  "ARMA+Logistic" = "#E07B39",
  "ARMA+Sine"     = "#3B9AB2",
  "AR2+Logistic"  = "#D4A017",
  "AR2+Sine"      = "#7EC8A4",
  "MA2+Logistic"  = "#C0392B",
  "MA2+Sine"      = "#8E44AD"
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot() +
  # HC boundary
  geom_line(data = bounds,
            aes(x = H, y = C, group = Side),
            color = "grey55", linetype = "dashed", linewidth = 0.5) +
  # Mixture path lines
  geom_path(data = path_df,
            aes(x = H, y = C, group = Group, color = Group),
            linewidth = 0.45, linetype = "dotted") +
  scale_color_manual(name = "Mixture path", values = group_path_colors,
                     guide = guide_legend(order = 2,
                                          override.aes = list(linewidth = 1))) +
  # All 37 points
  ggnewscale::new_scale_color() +
  geom_point(data = filter(results_df, Group != "Pure"),
             aes(x = H, y = C, color = Model),
             shape = 16, size = 2.2) +
  geom_point(data = filter(results_df, Group == "Pure"),
             aes(x = H, y = C, color = Model),
             shape = 18, size = 4) +
  scale_color_manual(name = "Model", values = model_colors,
                     guide = guide_legend(order = 1, nrow = 5,
                                          override.aes = list(size = 3,
                                                              shape = 16))) +
  # Pure model labels
  geom_label_repel(
    data      = filter(results_df, Group == "Pure"),
    aes(x = H, y = C, label = Model),
    color     = "grey20",
    size      = 3,
    family    = "serif",
    fontface  = "bold",
    box.padding   = 0.5,
    label.padding = 0.2,
    show.legend   = FALSE
  ) +
  labs(
    x        = expression(italic(H)),
    y        = expression(italic(C)),
    title    = paste0("HC Plane  —  D = ", D),
    subtitle = "Convex combination model"
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 10),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.key.size  = unit(0.5, "lines"),
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave("HC_scatter_only.pdf", p, width = 10, height = 8, device = cairo_pdf)
ggsave("HC_scatter_only.png", p, width = 10, height = 8, dpi = 200)
message("✓ Saved HC_scatter_only.pdf/.png")