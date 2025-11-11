##########################################################################
# ===========================================================
# Scatter plots for H_Shannon vs C_Shannon (HxC plane)
# From HC_seperate_results.xlsx — single label per model
# ===========================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(StatOrdPattHxC)
library(tidyr)

# -----------------------------------------------------------
# Paths and setup
# -----------------------------------------------------------
base_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results"
excel_file <- file.path(base_dir, "HC_seperate_results.xlsx")
sheets <- c("case1", "case2", "case3")

# -----------------------------------------------------------
# Load LinfLsup boundaries
# -----------------------------------------------------------
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# -----------------------------------------------------------
# Define colors & shapes
# -----------------------------------------------------------
model_colors <- c(
  "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
  "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
  "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
  "ARMA22_M3"="#fb9a99",
  "AR1_M1"="#8dd3c7","AR1_M2"="#ffffb3","ARMA11_M2"="#bebada","MA1_M1"="#fb8072",
  "MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69","AMA11_M4"="#fccde5",
  "MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
)

model_shapes <- c(
  "ARMA11_M1"=21,"AR2_M1"=22,"MA1_M2"=23,"ARMA11_M3"=24,"AR2_M4"=25,"MA1_M4"=8,
  "ARMA22_M2"=15,"AR2_M3"=16,"MA2_M3"=17,"ARMA22_M1"=18,"ARMA22_M4"=19,"MA2_M2"=4,
  "ARMA22_M3"=3,
  "AR1_M1"=1,"AR1_M2"=2,"ARMA11_M2"=5,"MA1_M1"=6,"MA2_M1"=7,"AR1_M3"=9,
  "AR1_M4"=10,"AMA11_M4"=11,"MA1_M3"=12,"MA2_M4"=13,"AR2_M2"=14
)

# -----------------------------------------------------------
# Function to plot each case
# -----------------------------------------------------------
plot_case <- function(sheet_name) {
  message("Processing ", sheet_name, "...")
  
  # Load Excel sheet
  df <- read_excel(excel_file, sheet = sheet_name) %>%
    mutate(
      H_Shannon = as.numeric(H_Shannon),
      C_Shannon = as.numeric(C_Shannon)
    ) %>%
    drop_na(H_Shannon, C_Shannon)
  
  # Get one label per model
  df_label <- df %>%
    group_by(Model) %>%
    summarise(
      H_Shannon = mean(H_Shannon, na.rm = TRUE),
      C_Shannon = mean(C_Shannon, na.rm = TRUE)
    )
  
  # Output folder
  case_folder <- paste0("Case ", gsub("case", "", sheet_name))
  plot_dir <- file.path(base_dir, case_folder, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # Linf focus
  Linf_focus <- LinfLsup_subset %>%
    filter(H >= min(df$H_Shannon) - 0.05,
           H <= max(df$H_Shannon) + 0.05,
           C >= min(df$C_Shannon) - 0.05,
           C <= max(df$C_Shannon) + 0.05)
  
  # ------------------------------------------------------
  # 📊 Scatter Plot (HxC Shannon Plane)
  # ------------------------------------------------------
  p <- ggplot() +
    geom_line(data = subset(Linf_focus, Side == "Lower"),
              aes(x = H, y = C),
              color = "gray40", linetype = "dashed") +
    geom_line(data = subset(Linf_focus, Side == "Upper"),
              aes(x = H, y = C),
              color = "gray40", linetype = "dashed") +
    geom_point(data = df,
               aes(x = H_Shannon, y = C_Shannon,
                   color = Model, shape = Model),
               size = 3) +
    # Add only one label per model (mean location)
    geom_text(data = df_label,
              aes(x = H_Shannon, y = C_Shannon, label = Model),
              vjust = -1, size = 3, fontface = "bold") +
    scale_color_manual(values = model_colors, na.value = "gray60") +
    scale_shape_manual(values = model_shapes, na.value = 16) +
    labs(
      title = paste0("Entropy–Complexity Plane (", case_folder, " — Shannon)"),
      x = expression(italic(H)[Shannon]),
      y = expression(italic(C)[Shannon])
    ) +
    theme_minimal(base_size = 13, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  # Save PDF
  plot_path <- file.path(plot_dir, paste0("HC_Scatter_", case_folder, "_Shannon.pdf"))
  ggsave(plot_path, p, width = 8, height = 6)
  
  print(p)
  message("✅ Saved clean scatter plot for ", case_folder, " at:\n", plot_path)
}

# -----------------------------------------------------------
# Run all cases
# -----------------------------------------------------------
for (sheet in sheets) {
  plot_case(sheet)
}

message("🎯 All scatter plots generated successfully!")

