
###############################################################################################################
############################################################################################################

# ===============================
# 📊 HC Central Points + CI + Time Series + Heatmaps (Case General)
# ===============================

library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(here)
library(viridis)
library(StatOrdPattHxC)

# ---- 1. Define file paths ----
case_name <- "New_Case6"

data_path <- here("GitHub", "Ordinal_Patterns_R", "Data",
                  "New ARMA Time series results", case_name,
                  paste0("HC_Central_SemiLength_", case_name, ".xlsx"))

plot_dir <- here("GitHub", "Ordinal_Patterns_R", "Data",
                 "New ARMA Time series results", paste0(case_name, " Plots"))

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ---- 2. Read Excel Sheets ----
df_central <- read_excel(data_path, sheet = 2)
df_ts <- read_excel(data_path, sheet = 3) %>%
  rename(Index = 1, Value = 2, Model = 3)  # standardize columns

# ---- 3. LinfLsup boundaries ----
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# ---- 4. Prepare central data ----
df_central <- df_central %>%
  rename(H_data = H_star, C_data = C_star, SemiLength = SemiLength) %>%
  mutate(
    n = case_when(
      grepl("5000", Entropy_Type) ~ 5000,
      grepl("10000", Entropy_Type) ~ 10000,
      TRUE ~ NA_real_
    ),
    n = factor(n),
    Model = factor(Model)
  )

# ---- 5. Define consistent color & shape palette for all models ----
model_list <- c("ARMA11_M1", "AR2_M1", "MA1_M2",
                "ARMA11_M3", "AR2_M4", "MA1_M4",
                "ARMA22_M2", "AR2_M3", "MA2_M3",
                "ARMA22_M1", "ARMA22_M4", "MA2_M2", "ARMA22_M3")

# consistent color palette
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

# filled shapes for each model
model_shapes <- c(
  "ARMA11_M1" = 21, "AR2_M1" = 22, "MA1_M2" = 23,
  "ARMA11_M3" = 24, "AR2_M4" = 25, "MA1_M4" = 8,
  "ARMA22_M2" = 15, "AR2_M3" = 16, "MA2_M3" = 17,
  "ARMA22_M1" = 18, "ARMA22_M4" = 19, "MA2_M2" = 4, "ARMA22_M3" = 3
)

# ---- 6. Generate plots per sample size ----
for (sample_size in levels(df_central$n)) {
  
  df_n <- df_central %>% filter(n == sample_size)
  models_to_plot <- intersect(unique(df_n$Model), names(model_colors))
  
  # Focus LinfLsup region near the data
  h_range <- range(df_n$H_data, na.rm = TRUE)
  c_range <- range(df_n$C_data, na.rm = TRUE)
  LinfLsup_focus <- LinfLsup_subset %>%
    filter(H >= min(h_range) - 0.05, H <= max(h_range) + 0.05,
           C >= min(c_range) - 0.05, C <= max(c_range) + 0.05)
  
  # ===================================================
  # 📍 PLOT 1 — Scatter with Confidence Intervals (H*C*)
  # ===================================================
  p_scatter <- ggplot() +
    # Focused boundary lines (not in legend)
    geom_line(
      data = subset(LinfLsup_focus, Side == "Lower"),
      aes(x = H, y = C),
      linetype = "dashed", color = "gray40", linewidth = 0.6
    ) +
    geom_line(
      data = subset(LinfLsup_focus, Side == "Upper"),
      aes(x = H, y = C),
      linetype = "dashed", color = "gray40", linewidth = 0.6
    ) +
    
    # Confidence intervals (horizontal)
    geom_errorbarh(
      data = df_n,
      aes(y = C_data,
          xmin = H_data - SemiLength,
          xmax = H_data + SemiLength,
          color = Model),
      height = 0.002, linewidth = 0.7
    ) +
    
    # Central points
    geom_point(
      data = df_n,
      aes(x = H_data, y = C_data, fill = Model, shape = Model),
      size = 4, color = "black"
    ) +
    
    geom_text(
      data = df_n,
      aes(x = H_data, y = C_data, label = Model),
      vjust = -1, fontface = "bold", size = 3
    ) +
    
    scale_fill_manual(values = model_colors) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    labs(
      title = paste0("Entropy–Complexity Plane with Confidence Intervals (n = ", sample_size, ")"),
      x = expression(italic(H)), y = expression(italic(C))
    ) +
    theme_minimal(base_size = 13, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(file.path(plot_dir, paste0("HC_Scatter_with_CI_", case_name, "_n", sample_size, ".pdf")),
         p_scatter, width = 8, height = 6)
  print(p_scatter)
  
  # ===================================================
  # 🔥 PLOT 2 — Heatmap of Central Points
  # ===================================================
  p_heat <- ggplot() +
    stat_density_2d(
      data = df_n,
      aes(x = H_data, y = C_data, fill = after_stat(level)),
      geom = "polygon", bins = 60, alpha = 0.6
    ) +
    geom_point(
      data = df_n,
      aes(x = H_data, y = C_data, shape = Model, color = Model),
      size = 3
    ) +
    scale_fill_viridis_c(option = "plasma") +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    labs(
      title = paste("Heatmap of Central Points (n =", sample_size, ")"),
      x = expression(italic(H)), y = expression(italic(C))
    ) +
    theme_minimal(base_size = 13, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(file.path(plot_dir, paste0("HC_Heatmap_", case_name, "_n", sample_size, ".pdf")),
         p_heat, width = 8, height = 6)
  print(p_heat)
  
  
  # ===================================================
  # 🕒 PLOT 3 — Time Series + HC Combined
  # ===================================================
  ts_plots <- lapply(models_to_plot, function(m) {
    df_m <- df_ts %>% filter(Model == m)
    ggplot(df_m, aes(x = Index, y = Value)) +
      geom_line(color = model_colors[m], size = 0.7) +
      labs(title = paste("Time Series -", m), x = "Time", y = NULL) +
      theme_minimal(base_size = 12, base_family = "serif") +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  })
  
  layout_matrix <- rbind(
    c(1, 2, 3),
    c(4, 7, 5),
    c(6, 7, NA)
  )
  
  plots_all <- c(ts_plots, list(p_scatter))
  
  combined_plot <- grid.arrange(
    grobs = plots_all,
    layout_matrix = layout_matrix,
    top = textGrob(
      paste("Entropy–Complexity and Time Series (n =", sample_size, ")"),
      gp = gpar(fontsize = 15, fontface = "bold")
    )
  )
  
  ggsave(file.path(plot_dir, paste0("HC_Central_with_TimeSeries_", case_name, "_n", sample_size, ".pdf")),
         combined_plot, width = 14, height = 10)
  
  message("✅ Saved all plots for ", case_name, " (n = ", sample_size, ")")
}

cat("\n🎯 All case plots saved in:", plot_dir, "\n")
###############################################################################################################
