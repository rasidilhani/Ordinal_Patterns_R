# ══════════════════════════════════════════════════════════════════════════════
#  Ordinal Pattern Distributions — 15 Models
#  Reads: timeseries/n10000/Rep_1.xlsx (wide format, one col per model)
#  Saves: faceted plot + combined overlay plot
# ══════════════════════════════════════════════════════════════════════════════
library(StatOrdPattHxC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(here)

# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
D      <- 4       # embedding dimension
n_val  <- 10000    
rep_id <- 1       

ts_dir     <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                   "timeseries", paste0("n", n_val))
output_dir <- here("Results", "Convex_combination", paste0("n", n_val))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL NAMES (as they appear in your model definitions)
# ══════════════════════════════════════════════════════════════════════════════
model_names <- c(
  "ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
  "ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)",
  "AR2+Logistic(w=0.1)", "AR2+Sine(w=0.8)",
  "MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)",
  "MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"
)

# ── Map Excel column names → model names ──────────────────────────────────────
col_to_model <- c(
  "ARMA.2.2."          = "ARMA(2,2)",
  "AR.2."              = "AR(2)",
  "MA.2."              = "MA(2)",
  "Logistic"           = "Logistic",
  "Sine"               = "Sine",
  "ARMA.Sine.w.0.1."   = "ARMA+Sine(w=0.1)",
  "ARMA.Sine.w.0.2."   = "ARMA+Sine(w=0.2)",
  "ARMA.Sine.w.0.3."   = "ARMA+Sine(w=0.3)",
  "AR2.Logistic.w.0.1."= "AR2+Logistic(w=0.1)",
  "AR2.Sine.w.0.8."    = "AR2+Sine(w=0.8)",
  "MA2.Logistic.w.0.2."= "MA2+Logistic(w=0.2)",
  "MA2.Logistic.w.0.7."= "MA2+Logistic(w=0.7)",
  "MA2.Sine.w.0.4."    = "MA2+Sine(w=0.4)",
  "MA2.Sine.w.0.6."    = "MA2+Sine(w=0.6)",
  "MA2.Sine.w.0.8."    = "MA2+Sine(w=0.8)"
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
#  READ TIME SERIES FROM EXCEL
# ══════════════════════════════════════════════════════════════════════════════
fname <- file.path(ts_dir, paste0("Rep_", rep_id, ".xlsx"))

if (!file.exists(fname)) {
  stop(sprintf("File not found: %s", fname))
}

ts_wide <- read_xlsx(fname)

# ══════════════════════════════════════════════════════════════════════════════
#  COMPUTE ORDINAL PATTERN PROBABILITIES
# ══════════════════════════════════════════════════════════════════════════════
all_op_probs <- list()

for (col_name in names(col_to_model)) {
  
  model_label <- col_to_model[[col_name]]
  
  # Check column exists
  if (!col_name %in% names(ts_wide)) {
    next
  }
  
  ts_values <- as.numeric(ts_wide[[col_name]])
  ts_values <- ts_values[!is.na(ts_values)]
  
  # Compute ordinal pattern probabilities
  op_prob <- as.numeric(OPprob(ts_values, emb = D))
  
  all_op_probs[[model_label]] <- data.frame(
    Model       = model_label,
    Pattern     = seq_along(op_prob),
    Probability = op_prob
  )
 
}

# Combine
combined_op <- bind_rows(all_op_probs)
combined_op$Model <- factor(combined_op$Model, levels = model_names)

n_patterns <- factorial(D)


# ══════════════════════════════════════════════════════════════════════════════
#  PLOT — Faceted (one panel per model)
# ══════════════════════════════════════════════════════════════════════════════
p_faceted <- ggplot(combined_op,
                    aes(x = Pattern, y = Probability, fill = Model)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ Model, ncol = 3, scales = "free_y") +
  geom_hline(yintercept = 1 / n_patterns,
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(1, n_patterns, by = 1)) +
  scale_fill_manual(values = model_colors) +
  labs(
    title    = paste0("Ordinal Pattern Distributions  (D = ", D,
                      ",  n = ", n_val, ",  Rep = ", rep_id, ")"),
    #subtitle = "Dashed line = uniform distribution",
    x        = "Ordinal Pattern",
    y        = "Probability"
  ) +
  theme_bw(base_family = "serif", base_size = 12) +
  theme(
    plot.title        = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "grey40"),
    strip.text        = element_text(size = 9,  face = "bold"),
    strip.background  = element_rect(fill = "white"),
    axis.text.x       = element_text(size = 8),
    axis.text.y       = element_text(size = 8),
    axis.title        = element_text(size = 11, face = "bold"),
    legend.position   = "none",
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(output_dir,
                 paste0("OP_Faceted_D", D, "_n", n_val, "_Rep", rep_id, ".pdf")),
       p_faceted, width = 14, height = 12, device = cairo_pdf)
print(p_faceted)

cat("\n✅ All ordinal pattern plots complete!\n")