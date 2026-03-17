# ══════════════════════════════════════════════════════════════════════════════
#  HC Plane with Time Series Thumbnails — Shannon Only
#  Reads Rep_1.xlsx (wide format) and HC_Results.xlsx (sheet n10000)
#  Uses first replication only
# ══════════════════════════════════════════════════════════════════════════════
library(StatOrdPattHxC)
library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
library(here)
library(ggrepel)

# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
D      <- 4
n_val  <- 10000    # sample size
rep_id <- 1       # first replication only

# ══════════════════════════════════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════════════════════════════════
ts_dir       <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     "timeseries", paste0("n", n_val))
hc_data_path <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))
output_dir   <- here("Results", "Convex_combination", paste0("n", n_val))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL DEFINITIONS
# ══════════════════════════════════════════════════════════════════════════════
model_names <- c(
  "ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
  "ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)",
  "AR2+Logistic(w=0.1)", "AR2+Sine(w=0.8)",
  "MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)",
  "MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"
)

pure_models <- c("ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine")

# Excel renames special characters to dots
col_to_model <- c(
  "ARMA.2.2."           = "ARMA(2,2)",
  "AR.2."               = "AR(2)",
  "MA.2."               = "MA(2)",
  "Logistic"            = "Logistic",
  "Sine"                = "Sine",
  "ARMA.Sine.w.0.1."    = "ARMA+Sine(w=0.1)",
  "ARMA.Sine.w.0.2."    = "ARMA+Sine(w=0.2)",
  "ARMA.Sine.w.0.3."    = "ARMA+Sine(w=0.3)",
  "AR2.Logistic.w.0.1." = "AR2+Logistic(w=0.1)",
  "AR2.Sine.w.0.8."     = "AR2+Sine(w=0.8)",
  "MA2.Logistic.w.0.2." = "MA2+Logistic(w=0.2)",
  "MA2.Logistic.w.0.7." = "MA2+Logistic(w=0.7)",
  "MA2.Sine.w.0.4."     = "MA2+Sine(w=0.4)",
  "MA2.Sine.w.0.6."     = "MA2+Sine(w=0.6)",
  "MA2.Sine.w.0.8."     = "MA2+Sine(w=0.8)"
)

# ══════════════════════════════════════════════════════════════════════════════
#  COLORS AND SHAPES
# ══════════════════════════════════════════════════════════════════════════════
model_colors <- c(
  "ARMA(2,2)" = "black",
  "AR(2)"     = "tomato",
  "MA(2)"     = "navy",
  "Logistic"  = "dodgerblue",
  "Sine"      = "forestgreen",
  setNames(colorRampPalette(c("#A8D8EA", "#3B9AB2"))(3),
           c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)")),
  "AR2+Logistic(w=0.1)" = "#D4A017",
  "AR2+Sine(w=0.8)"     = "#7EC8A4",
  setNames(colorRampPalette(c("#F5AAAA", "#C0392B"))(2),
           c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  setNames(colorRampPalette(c("#D9B8E8", "#8E44AD"))(3),
           c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"))
)

model_shapes <- c(
  "ARMA(2,2)" = 16, "AR(2)" = 17, "MA(2)" = 15,
  "Logistic"  = 18, "Sine"  = 8,
  setNames(rep(17, 3), c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)",
                         "ARMA+Sine(w=0.3)")),
  "AR2+Logistic(w=0.1)" = 15,
  "AR2+Sine(w=0.8)"     = 18,
  setNames(rep(1, 2), c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  setNames(rep(2, 3), c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)",
                        "MA2+Sine(w=0.8)"))
)

# ══════════════════════════════════════════════════════════════════════════════
#  READ TIME SERIES — Rep_1.xlsx
# ══════════════════════════════════════════════════════════════════════════════
ts_file <- file.path(ts_dir, paste0("Rep_", rep_id, ".xlsx"))
if (!file.exists(ts_file)) stop(sprintf("File not found: %s", ts_file))

ts_wide <- read_xlsx(ts_file)
cat(sprintf("✅ Read time series: %s  (%d rows)\n", basename(ts_file), nrow(ts_wide)))

# ══════════════════════════════════════════════════════════════════════════════
#  READ HC RESULTS — sheet n10000, Rep 1, Shannon only
# ══════════════════════════════════════════════════════════════════════════════
hc_raw <- read_xlsx(hc_data_path, sheet = paste0("n", n_val))

hc_data <- hc_raw %>%
  filter(Rep == rep_id) %>%
  select(Model, H_Shannon, C_Shannon, Semi_H_Shannon, Semi_C_Shannon) %>%
  rename(H_val  = H_Shannon,
         C_val  = C_Shannon,
         Semi_H = Semi_H_Shannon,
         Semi_C = Semi_C_Shannon) %>%
  filter(!is.na(H_val), !is.na(C_val),
         is.finite(H_val), is.finite(C_val))

hc_data$Model <- factor(hc_data$Model, levels = model_names)
hc_data <- hc_data %>% arrange(Model)

cat(sprintf("✅ HC data loaded: %d models for Rep %d\n", nrow(hc_data), rep_id))
print(hc_data)

# ══════════════════════════════════════════════════════════════════════════════
#  HC BOUNDS (Shannon only)
# ══════════════════════════════════════════════════════════════════════════════
data("LinfLsup")
boundary_lower <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Lower")
boundary_upper <- LinfLsup %>% filter(Dimension == as.character(D), Side == "Upper")

# ══════════════════════════════════════════════════════════════════════════════
#  HC PLANE PLOT (Shannon, with CI, no legend — for combining)
# ══════════════════════════════════════════════════════════════════════════════
hc_data_ci <- hc_data %>%
  mutate(
    has_H_error = !is.na(Semi_H) & Semi_H > 0 & is.finite(Semi_H),
    has_C_error = !is.na(Semi_C) & Semi_C > 0 & is.finite(Semi_C)
  )

p_hc <- ggplot() +
  geom_line(data = boundary_lower, aes(x = H, y = C),
            color = "gray50", linetype = "dashed", linewidth = 1.0, alpha = 0.8) +
  geom_line(data = boundary_upper, aes(x = H, y = C),
            color = "gray50", linetype = "dashed", linewidth = 1.0, alpha = 0.8) +
  geom_errorbarh(data = hc_data_ci %>% filter(has_H_error),
                 aes(y = C_val, xmin = H_val - Semi_H, xmax = H_val + Semi_H,
                     color = Model),
                 height = 0.005, linewidth = 0.8, alpha = 0.8,
                 show.legend = FALSE) +
  geom_errorbar(data = hc_data_ci %>% filter(has_C_error),
                aes(x = H_val, ymin = C_val - Semi_C, ymax = C_val + Semi_C,
                    color = Model),
                width = 0.01, linewidth = 0.8, alpha = 0.8,
                show.legend = FALSE) +
  geom_point(data = hc_data,
             aes(x = H_val, y = C_val, color = Model, shape = Model),
             size = 3, stroke = 1.2) +
  geom_label_repel(data = hc_data %>% filter(Model %in% pure_models),
                   aes(x = H_val, y = C_val, label = Model, color = Model),
                   size = 2.5, family = "serif", fontface = "bold",
                   box.padding = 0.4, label.padding = 0.15,
                   show.legend = FALSE) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  scale_x_continuous(limits = c(0, 1),   breaks = seq(0, 1,   0.25)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title    = paste0("Shannon HC Plane  (D = ", D, ",  n = ", n_val,
                      ",  Rep = ", rep_id, ")"),
    x        = expression(italic(H)[Shannon]),
    y        = expression(italic(C)[Shannon]),
    color    = "Model",
    shape    = "Model"
  ) +
  guides(color = guide_legend(nrow = 3, override.aes = list(size = 3)),
         shape = guide_legend(nrow = 3)) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text        = element_text(size = 10),
    legend.position  = "bottom",
    legend.text      = element_text(size = 7),
    legend.title     = element_text(size = 8, face = "bold"),
    legend.key.size  = unit(0.5, "cm"),
    panel.grid.minor = element_blank()
  )

# Save standalone HC plot
ggsave(file.path(output_dir,
                 paste0("HC_Shannon_D", D, "_n", n_val, "_Rep", rep_id, ".pdf")),
       p_hc, width = 10, height = 8, device = cairo_pdf)
cat(sprintf("✅ Saved: HC_Shannon_D%d_n%d_Rep%d.pdf\n", D, n_val, rep_id))
print(p_hc)

# Version without legend (for combined layout)
p_hc_no_legend <- p_hc + theme(legend.position = "none")

# ══════════════════════════════════════════════════════════════════════════════
#  TIME SERIES THUMBNAIL HELPER
# ══════════════════════════════════════════════════════════════════════════════
make_ts_plot <- function(col_name, model_label) {
  
  if (!col_name %in% names(ts_wide)) {
    cat(sprintf("⚠️  Column not found: %s\n", col_name))
    return(ggplot() + theme_void() +
             labs(title = paste(model_label, "- NOT FOUND")))
  }
  
  ts_vals <- as.numeric(ts_wide[[col_name]])
  ts_vals <- ts_vals[!is.na(ts_vals)]
  n_show  <- min(500, length(ts_vals))
  col     <- model_colors[model_label]
  
  df <- data.frame(t = seq_len(n_show), y = ts_vals[seq_len(n_show)])
  
  ggplot(df, aes(x = t, y = y)) +
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
      plot.margin      = unit(c(3, 3, 3, 3), "pt"),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white")
    )
}

# Build all 15 thumbnail plots
ts_plots <- mapply(make_ts_plot,
                   names(col_to_model),
                   col_to_model,
                   SIMPLIFY = FALSE)
names(ts_plots) <- col_to_model   # key by model label

cat(sprintf("✅ Created %d time series thumbnails\n", length(ts_plots)))

# ══════════════════════════════════════════════════════════════════════════════
#  COMBINED LAYOUT: HC PLANE CENTRE, TIME SERIES AROUND EDGES
#  Layout (5 cols x 6 rows):
#    Row 1: [NA] [TS1] [TS2] [TS3] [TS4]       ← top: 4 models
#    Row 2: [TS5] [HC ] [HC ] [HC ] [TS9]       ─┐
#    Row 3: [TS6] [HC ] [HC ] [HC ] [TS10]       │ sides: 4 left, 4 right
#    Row 4: [TS7] [HC ] [HC ] [HC ] [TS11]       │
#    Row 5: [TS8] [HC ] [HC ] [HC ] [TS12]      ─┘
#    Row 6: [TS13][TS14][TS15][ NA][ NA]        ← bottom: 3 models
# ══════════════════════════════════════════════════════════════════════════════

# Assign positions (15 thumbnails + 1 HC = 16 slots)
# Use model_names order: 1-15 for thumbnails, 16 = HC plane
top_models    <- model_names[1:4]      # row 1: models 1-4
left_models   <- model_names[5:8]      # col 1: models 5-8
right_models  <- model_names[9:12]     # col 5: models 9-12
bottom_models <- model_names[13:15]    # row 6: models 13-15

# Build indexed grob list: slots 1-15 = thumbnails, slot 16 = HC
all_plots <- vector("list", 16)
for (i in seq_along(model_names)) {
  nm <- model_names[i]
  all_plots[[i]] <- if (nm %in% names(ts_plots)) ts_plots[[nm]] else
    ggplot() + theme_void() + labs(title = paste(nm, "- missing"))
}
all_plots[[16]] <- p_hc_no_legend

layout_matrix <- rbind(
  c(NA,  4,  9, 11,  1),
  c( 5, 16, 16, 16,  2),
  c( 6, 16, 16, 16, 10),
  c( 7, 16, 16, 16, 12),
  c( 8, 16, 16, 16, 3),
  c(NA, 13, 14, 15, NA)
)

combined_plot <- gridExtra::grid.arrange(
  grobs       = all_plots,
  layout_matrix = layout_matrix,
  widths      = c(1.2, 1.2, 1.2, 1.2, 1.2),
  heights     = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2)
)

ggsave(
  filename = file.path(output_dir,
                       paste0("HC_TimeSeries_D", D, "_n", n_val,
                              "_Rep", rep_id, ".pdf")),
  plot   = combined_plot,
  width  = 20,
  height = 20,
  bg     = "white",
  device = cairo_pdf
)
