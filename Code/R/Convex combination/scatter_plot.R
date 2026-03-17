# ══════════════════════════════════════════════════════════════════════════════
#  HC Plane Scatter Plots — Shannon, Tsallis, Renyi, Fisher
#  Reads HC_Results.xlsx (sheets: n1000, n510000)
# ══════════════════════════════════════════════════════════════════════════════
library(tidyverse)
library(readxl)
library(ggrepel)
library(here)
library(StatOrdPattHxC)

# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS AND PATHS
# ══════════════════════════════════════════════════════════════════════════════
D            <- 4
sample_sizes <- c(1000,10000)

hc_data_path <- here("Data", "Convex_combination", "D4_Data", "HC_Results_D4.xlsx")
output_dir   <- here("Results", "Convex_combination")

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL NAMES, COLORS, SHAPES
# ══════════════════════════════════════════════════════════════════════════════

model_names <- c(
  "ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
  "ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)",
  "AR2+Logistic(w=0.1)", "AR2+Sine(w=0.8)",
  "MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)",
  "MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"
)

model_colors <- c(
  # Pure models
  "ARMA(2,2)" = "black",
  "AR(2)"     = "tomato",
  "MA(2)"     = "navy",
  "Logistic"  = "dodgerblue",
  "Sine"      = "forestgreen",
  
  # ARMA+Sine — 3 models (blue gradient)
  setNames(colorRampPalette(c("#A8D8EA", "#3B9AB2"))(3),
           c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)")),
  
  # AR2+Logistic — 1 model (gold)
  "AR2+Logistic(w=0.1)" = "#D4A017",
  
  # AR2+Sine — 1 model (teal)
  "AR2+Sine(w=0.8)" = "#7EC8A4",
  
  # MA2+Logistic — 2 models (red gradient)
  setNames(colorRampPalette(c("#F5AAAA", "#C0392B"))(2),
           c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  
  # MA2+Sine — 3 models (purple gradient)
  setNames(colorRampPalette(c("#D9B8E8", "#8E44AD"))(3),
           c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"))
)

model_shapes <- c(
  # Pure models
  "ARMA(2,2)" = 16,   # filled circle
  "AR(2)"     = 17,   # filled triangle
  "MA(2)"     = 15,   # filled square
  "Logistic"  = 18,   # filled diamond
  "Sine"      = 8,    # asterisk
  
  # ARMA+Sine — 3 models
  setNames(rep(17, 3),
           c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)")),
  
  # AR2+Logistic — 1 model
  "AR2+Logistic(w=0.1)" = 15,
  
  # AR2+Sine — 1 model
  "AR2+Sine(w=0.8)" = 18,
  
  # MA2+Logistic — 2 models
  setNames(rep(1, 2),
           c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  
  # MA2+Sine — 3 models
  setNames(rep(2, 3),
           c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"))
)

# ══════════════════════════════════════════════════════════════════════════════
#  LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
hc_data_n1000 <- read_xlsx(hc_data_path, sheet = "n1000")
hc_data_n10000 <- read_xlsx(hc_data_path, sheet = "n10000")

hc_data_all <- bind_rows(hc_data_n1000, hc_data_n10000)

# Set model factor order
hc_data_all$Model <- factor(hc_data_all$Model, levels = model_names)

# ══════════════════════════════════════════════════════════════════════════════
#  HC BOUNDS (Shannon only)
# ══════════════════════════════════════════════════════════════════════════════
data("LinfLsup")
bounds <- LinfLsup %>% filter(Dimension == as.character(D))
boundary_lower <- bounds %>% filter(Side == "Lower")
boundary_upper <- bounds %>% filter(Side == "Upper")

# ══════════════════════════════════════════════════════════════════════════════
#  ENTROPY CONFIGURATIONS
# ══════════════════════════════════════════════════════════════════════════════
entropy_configs <- list(
  
  Shannon = list(
    H_col      = "H_Shannon",
    C_col      = "C_Shannon",
    # Semi_H_col = "Semi_H_Shannon",   
    # Semi_C_col = "Semi_C_Shannon",   
    add_bounds = TRUE,
    title      = "Shannon Entropy-Complexity Plane",
    xlab       = expression(italic(H)[Shannon]),
    ylab       = expression(italic(C)[Shannon])
  ),
  
  Renyi = list(
    H_col      = "H_Renyi",
    C_col      = "C_Renyi",
    # Semi_H_col = "Semi_H_Renyi",    
    # Semi_C_col = NA,
    add_bounds = FALSE,
    title      = "Rényi Entropy-Complexity Plane",
    xlab       = expression(italic(H)[Rényi]),
    ylab       = expression(italic(C)[Rényi])
  ),
  
  Tsallis = list(
    H_col      = "H_Tsallis",
    C_col      = "C_Tsallis",
    # Semi_H_col = "Semi_H_Tsallis",   
    # Semi_C_col = NA,
    add_bounds = FALSE,
    title      = "Tsallis Entropy-Complexity Plane",
    xlab       = expression(italic(H)[Tsallis]),
    ylab       = expression(italic(C)[Tsallis])
  ),
  
  Fisher = list(
    H_col      = "H_Fisher",
    C_col      = "C_Fisher",
    # Semi_H_col = "Semi_H_Fisher",    
    # Semi_C_col = NA,
    add_bounds = FALSE,
    title      = "Fisher Information-Complexity Plane",
    xlab       = expression(italic(H)[Fisher]),
    ylab       = expression(italic(C)[Fisher])
  )
)

# ══════════════════════════════════════════════════════════════════════════════
#  PLOT LOOP
# ══════════════════════════════════════════════════════════════════════════════
for (n_val in sample_sizes) {
  
  # Output subfolder per sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  hc_data_n <- hc_data_all %>% filter(n == n_val)
  
  for (ent_name in names(entropy_configs)) {
    
    cfg   <- entropy_configs[[ent_name]]
    H_col <- cfg$H_col
    C_col <- cfg$C_col
    
    # Check columns exist
    if (!all(c(H_col, C_col) %in% names(hc_data_n))) {
      cat(sprintf("⚠️  Skipping %s — columns not found\n", ent_name))
      next
    }
    
    # Prepare plot data
    plot_data <- hc_data_n %>%
      select(Model, Rep, all_of(c(H_col, C_col))) %>%
      rename(H_val = !!H_col, C_val = !!C_col) %>%
      filter(!is.na(H_val), !is.na(C_val),
             is.finite(H_val), is.finite(C_val))
    
    if (nrow(plot_data) == 0) {
      cat(sprintf("⚠️  Skipping %s — no valid data\n", ent_name))
      next
    }
    
    # ── Semi-length columns (only used when uncommented in entropy_configs) ──
    # Safely read Semi_H_col and Semi_C_col — default to NULL if not defined
    Semi_H_col <- if (!is.null(cfg$Semi_H_col)) cfg$Semi_H_col else NULL
    Semi_C_col <- if (!is.null(cfg$Semi_C_col)) cfg$Semi_C_col else NULL
    
    # Add semi-length columns to plot_data if available
    plot_data$Semi_H <- if (!is.null(Semi_H_col) && Semi_H_col %in% names(hc_data_n))
      hc_data_n[[Semi_H_col]][match(paste(plot_data$Model, plot_data$Rep),
                                    paste(hc_data_n$Model, hc_data_n$Rep))]
    else NA
    
    plot_data$Semi_C <- if (!is.null(Semi_C_col) && Semi_C_col %in% names(hc_data_n))
      hc_data_n[[Semi_C_col]][match(paste(plot_data$Model, plot_data$Rep),
                                    paste(hc_data_n$Model, hc_data_n$Rep))]
    else NA
    
    # Error bar data from Rep 1 only
    plot_data_ci <- plot_data %>%
      filter(Rep == 1) %>%
      mutate(
        has_H_error = !is.na(Semi_H) & Semi_H > 0 & is.finite(Semi_H),
        has_C_error = !is.na(Semi_C) & Semi_C > 0 & is.finite(Semi_C)
      )
    
    # Pure model names for labels
    pure_models <- c("ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine")
    
    # ── Build plot ────────────────────────────────────────────────────────────
    p <- ggplot()
    
    # Bounds (Shannon only)
    if (cfg$add_bounds) {
      p <- p +
        geom_line(data = boundary_lower, aes(x = H, y = C),
                  color = "gray50", linetype = "dashed",
                  linewidth = 1.0, alpha = 0.8) +
        geom_line(data = boundary_upper, aes(x = H, y = C),
                  color = "gray50", linetype = "dashed",
                  linewidth = 1.0, alpha = 0.8)
    }
    
    # Horizontal error bars — only drawn when Semi_H_col is uncommented
    if (any(plot_data_ci$has_H_error, na.rm = TRUE)) {
      p <- p +
        geom_errorbarh(
          data = plot_data_ci %>% filter(has_H_error),
          aes(y = C_val, xmin = H_val - Semi_H, xmax = H_val + Semi_H,
              color = Model),
          height = 0.002, linewidth = 0.6, alpha = 0.5,
          show.legend = FALSE
        )
    }
    
    # Vertical error bars — only drawn when Semi_C_col is uncommented
    if (any(plot_data_ci$has_C_error, na.rm = TRUE)) {
      p <- p +
        geom_errorbar(
          data = plot_data_ci %>% filter(has_C_error),
          aes(x = H_val, ymin = C_val - Semi_C, ymax = C_val + Semi_C,
              color = Model),
          width = 0.01, linewidth = 0.6, alpha = 0.5,
          show.legend = FALSE
        )
    }
    
    # All points + pure model labels
    p <- p +
      geom_point(
        data = plot_data,
        aes(x = H_val, y = C_val, color = Model, shape = Model),
        size = 3, stroke = 1.2, alpha = 0.7
      ) +
      geom_label_repel(
        data          = plot_data %>% filter(Rep == 1, Model %in% pure_models),
        aes(x = H_val, y = C_val, label = Model, color = Model),
        size          = 3,
        family        = "serif",
        fontface      = "bold",
        box.padding   = 0.4,
        label.padding = 0.2,
        show.legend   = FALSE
      ) +
      scale_color_manual(
        values = model_colors,
        name   = "Model",
        guide  = guide_legend(override.aes = list(size = 4, alpha = 1,
                                                  stroke = 1.5))
      ) +
      scale_shape_manual(values = model_shapes, guide = "none") +
      labs(
        title    = cfg$title,
        subtitle = paste0("D = ", D, ",  n = ", n_val,
                          "  (", max(plot_data$Rep, na.rm = TRUE),
                          " replications)"),
        x = cfg$xlab,
        y = cfg$ylab
      ) +
      theme_bw(base_family = "serif", base_size = 13) +
      theme(
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = 11, hjust = 0.5),
        axis.title       = element_text(size = 14, face = "bold"),
        axis.text        = element_text(size = 11),
        legend.position  = "right",
        legend.title     = element_text(size = 11, face = "bold"),
        legend.text      = element_text(size = 9),
        legend.key.size  = unit(0.6, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.background = element_rect(fill = "white"),
        plot.background  = element_rect(fill = "white")
      )
    
    # Save
    fname <- file.path(n_dir,
                       paste0("HC_", ent_name, "_D", D, "_n", n_val, ".pdf"))
    ggsave(fname, p, width = 12, height = 8, bg = "white")
    print(p)
  }
}

cat("\n✅ All plots complete!\n")
