# =====================================================
# 📊 Combined HC Results: CI Scatter, Heatmap & Time Series
# =====================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(here)
library(viridis)
library(StatOrdPattHxC)

# ---- 1. Define Paths ----
data_path <- here("GitHub", "Ordinal_Patterns_R", "Data",
                  "New ARMA Time series results", "All_HC_Central_SemiLength.xlsx")

plot_dir <- here("GitHub", "Ordinal_Patterns_R", "Data",
                 "New ARMA Time series results", "All_HC_Combined_Plots")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---- 2. Load Data ----
df_central <- read_excel(data_path, sheet = "Central_HC_SemiLength")
df_ts <- read_excel(data_path, sheet = "time series_n5000")

# ---- 3. Load LinfLsup (for boundary reference) ----
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# ---- 4. Prepare Central Data ----
df_central <- df_central %>%
  rename(
    H_data = H_star,
    C_data = C_star,
    SemiLength = SemiLength
  ) %>%
  mutate(
    n = case_when(
      grepl("5000", Entropy_Type) ~ 5000,
      grepl("10000", Entropy_Type) ~ 10000,
      TRUE ~ NA_real_
    ),
    n = as.factor(n),
    Model = as.factor(Model)
  )

# ---- 5. Define Consistent Color and Shape Palette ----
model_colors <- c(
  "ARMA11_M1" = "#1b9e77",
  "AR2_M1"    = "#d95f02",
  "MA1_M2"    = "#7570b3",
  "ARMA11_M3" = "#e7298a",
  "AR2_M4"    = "#66a61e",
  "MA1_M4"    = "#e6ab02",
  "ARMA22_M2" = "#a6761d",
  "AR2_M3"    = "#666666",
  "MA2_M3"    = "#1f78b4",
  "ARMA22_M1" = "#b15928",
  "ARMA22_M4" = "#6a3d9a",
  "MA2_M2"    = "#33a02c",
  "ARMA22_M3" = "#fb9a99"
)

model_shapes <- c(
  "ARMA11_M1" = 21, "AR2_M1" = 22, "MA1_M2" = 23,
  "ARMA11_M3" = 24, "AR2_M4" = 25, "MA1_M4" = 8,
  "ARMA22_M2" = 15, "AR2_M3" = 16, "MA2_M3" = 17,
  "ARMA22_M1" = 18, "ARMA22_M4" = 19, "MA2_M2" = 4, "ARMA22_M3" = 3
)

# ---- 6. Generate Plots by Sample Size ----
for (sample_size in levels(df_central$n)) {
  df_n <- df_central %>% filter(n == sample_size)
  models_to_plot <- intersect(unique(df_n$Model), names(model_colors))
  
  # Adjusted LinfLsup boundaries near data range
  h_range <- range(df_n$H_data, na.rm = TRUE)
  c_range <- range(df_n$C_data, na.rm = TRUE)
  LinfLsup_focus <- LinfLsup_subset %>%
    filter(H >= min(h_range) - 0.05, H <= max(h_range) + 0.05,
           C >= min(c_range) - 0.05, C <= max(c_range) + 0.05)
  
  # ---------------- Scatter Plot with CI ----------------
  p_scatter <- ggplot() +
    # Add trimmed LinfLsup boundaries (no legend)
    geom_line(data = subset(LinfLsup_focus, Side == "Lower"),
              aes(x = H, y = C),
              linetype = "dashed", color = "gray40", linewidth = 0.6) +
    geom_line(data = subset(LinfLsup_focus, Side == "Upper"),
              aes(x = H, y = C),
              linetype = "dashed", color = "gray40", linewidth = 0.6) +
    
    # Confidence intervals
    geom_errorbarh(data = df_n,
                   aes(y = C_data, xmin = H_data - SemiLength, xmax = H_data + SemiLength),
                   height = 0.002, linewidth = 0.7, color = "black") +
    
    # Central points
    geom_point(data = df_n,
               aes(x = H_data, y = C_data, fill = Model, shape = Model),
               size = 4, color = "black") +
    
    geom_text(data = df_n, aes(x = H_data, y = C_data, label = Model),
              vjust = -1, fontface = "bold", size = 3, show.legend = FALSE) +
    
    scale_fill_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    labs(
      title = paste("Entropy–Complexity Central Points with CI (n =", sample_size, ")"),
      x = expression(italic(H)), y = expression(italic(C))
    ) +
    theme_minimal(base_size = 13, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(file.path(plot_dir, paste0("Combined_HC_CI_Scatter_n", sample_size, ".pdf")),
         p_scatter, width = 8, height = 6)
  print(p_scatter)
  
  # ---------------- Heatmap ----------------
  p_heat <- ggplot(df_n, aes(x = H_data, y = C_data)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", bins = 60, alpha = 0.6) +
    scale_fill_viridis_c(option = "plasma") +
    geom_point(aes(color = Model, shape = Model), size = 3) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    labs(
      title = paste("Heatmap of Entropy–Complexity Plane (n =", sample_size, ")"),
      x = expression(italic(H)), y = expression(italic(C))
    ) +
    theme_minimal(base_size = 13, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(file.path(plot_dir, paste0("Combined_HC_Heatmap_n", sample_size, ".pdf")),
         p_heat, width = 8, height = 6)
  print(p_heat)
  
  # ---------------- Combine with Time Series (only n = 5000) ----------------
  if (sample_size == 5000) {
    ts_plots <- lapply(models_to_plot, function(m) {
      df_m <- df_ts %>% filter(Model == m)
      ggplot(df_m, aes(x = Index, y = Value)) +
        geom_line(color = model_colors[m], size = 0.7) +
        labs(title = paste("Time Series -", m), x = "Time", y = NULL) +
        theme_minimal(base_size = 12, base_family = "serif") +
        theme(axis.text = element_blank(), axis.ticks = element_blank())
    })
    
    # Layout (3 top, 3 bottom, middle scatter)
    layout_matrix <- rbind(
      c(1,  2,  3,  4),
      c(5, 10, 10,  6),
      c(7, 10, 10,  8),
      c(9,  NA, NA, NA)
    )
    
    plots_all <- c(ts_plots, list(p_scatter))
    combined_plot <- grid.arrange(grobs = plots_all, layout_matrix = layout_matrix,
                                  top = textGrob(paste("Combined Entropy–Complexity and Time Series (n =", sample_size, ")"),
                                                 gp = gpar(fontsize = 15, fontface = "bold")))
    
    ggsave(file.path(plot_dir, paste0("Combined_HC_TimeSeries_n", sample_size, ".pdf")),
           combined_plot, width = 14, height = 10)
  }
}

cat("\n✅ All combined plots saved in:\n", plot_dir, "\n")



