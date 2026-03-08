#----------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# This script generates time series plots versus HC plane. It reads the time series data for each model, 
# read the corresponding HC values, and creates scatter plots of entropy vs. complexity, colored by 
# model type. 

library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(gridExtra)
library(grid)
library(gtable)
library(here)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension

# Paths
ts_dir <- here("Data", "Hybrid model_data", "D4_Data", "All_Models_TimeSeries_D4")
hc_data_path <- here("Data", "Hybrid model_data", "D4_Data", "Combined_All_Models_HC_Results_D4.xlsx")
output_dir <- here("Plots", "Hybrid Model plots", "D4 plots")

# ==============================================================================
# DEFINE MODEL ORDER, COLORS, AND SHAPES
# ==============================================================================

# Original names (for data matching)
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", 
  "Logistic_r3_6", "Sine_Wave", "Logistic_Sine_Combined"
)

# DISPLAY NAMES (for plots only - nice and readable)
model_display_names <- c(
  "ARMA(2,2)" = "ARMA(2,2)",
  "ARMA(1,1)" = "ARMA(1,1)",
  "AR(2)" = "AR(2)",
  "AR(1)" = "AR(1)",
  "MA(2)" = "MA(2)",
  "MA(1)" = "MA(1)",
  "Logistic" = "Logistic (r=3.8)",
  "Hybrid_ARMA(2,2)" = "Hybrid ARMA(2,2)",
  "Hybrid_ARMA(1,1)" = "Hybrid ARMA(1,1)",
  "Hybrid_AR(2)" = "Hybrid AR(2)",
  "Hybrid_AR(1)" = "Hybrid AR(1)",
  "Hybrid_MA(2)" = "Hybrid MA(2)",
  "Hybrid_MA(1)" = "Hybrid MA(1)",
  "Logistic_r3_6" = "Logistic (r=3.6)",
  "Sine_Wave" = "Sine Wave",
  "Logistic_Sine_Combined" = "Logistic + Sine"
)

model_colors <- c(
  "ARMA(2,2)" = "#D55E00",   # Vermillion
  "ARMA(1,1)" = "#0072B2",   # Blue
  "AR(2)" = "#009E73",       # Bluish green
  "AR(1)" = "#CC79A7",       # Reddish purple
  "MA(2)" = "#E69F00",       # Orange
  "MA(1)" = "#56B4E9",       # Sky blue
  "Logistic" = "#000000",    # Black
  
  "Hybrid_ARMA(2,2)" = "#8B4513", # Dark brown
  "Hybrid_ARMA(1,1)" = "#4B0082", # Indigo
  "Hybrid_AR(2)" = "#696969",     # Dark gray
  "Hybrid_AR(1)" = "#20B2AA",     # Light sea green
  "Hybrid_MA(2)" = "#B22222",     # Firebrick
  "Hybrid_MA(1)" = "#4169E1",     # Royal blue
  
  "Logistic_r3_6" = "#2F4F4F",    # Slate gray
  "Sine_Wave" = "#228B22",        # Forest green
  "Logistic_Sine_Combined" = "#8A2BE2" # Blue violet
)

model_shapes <- c(
  "ARMA(2,2)" = 16,
  "ARMA(1,1)" = 17,
  "AR(2)" = 15,
  "AR(1)" = 18,
  "MA(2)" = 3,
  "MA(1)" = 7,
  "Logistic" = 8,
  
  "Hybrid_ARMA(2,2)" = 1,
  "Hybrid_ARMA(1,1)" = 2,
  "Hybrid_AR(2)" = 0,
  "Hybrid_AR(1)" = 5,
  "Hybrid_MA(2)" = 6,
  "Hybrid_MA(1)" = 4,
  
  "Logistic_r3_6" = 9,
  "Sine_Wave" = 10,
  "Logistic_Sine_Combined" = 11
)

# ==============================================================================
# LOAD HC DATA
# ==============================================================================
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Filter for n=5000 and first repetition
hc_data <- hc_data_all %>%
  filter(n == 5000, rep == 1) %>%
  select(Model, H_Shannon, C_Shannon, Semi_HS, Semi_CS)

# Rename columns for convenience
names(hc_data) <- c("Model", "H_val", "C_val", "Semi_H", "Semi_C")

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

# Sort by model order
hc_data <- hc_data %>% arrange(Model)

print("HC Data loaded:")
print(hc_data)

# ==============================================================================
# LOAD LinfLsup BOUNDARIES
# ==============================================================================
data("LinfLsup")

# Filter by Dimension
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# Choose which boundary to use
boundary_data <- LinfLsup_full

# Separate upper and lower boundaries
boundary_lower <- boundary_data %>% filter(Side == "Lower")
boundary_upper <- boundary_data %>% filter(Side == "Upper")

# ==============================================================================
# CREATE HC PLANE PLOT WITH CONFIDENCE INTERVALS (WITH LEGEND)
# ==============================================================================

# Prepare data for error bars
hc_data_ci <- hc_data %>%
  mutate(
    has_H_error = !is.na(Semi_H) & Semi_H > 0 & is.finite(Semi_H),
    has_C_error = !is.na(Semi_C) & Semi_C > 0 & is.finite(Semi_C)
  )

# Create the plot with display names in legend
p_hc <- ggplot() +
  # Add LOWER boundary (dashed gray)
  geom_line(data = boundary_lower, aes(x = H, y = C), 
            color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
  
  # Add UPPER boundary (dashed gray)
  geom_line(data = boundary_upper, aes(x = H, y = C), 
            color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
  
  # Add horizontal error bars (H_Shannon confidence intervals)
  geom_errorbarh(data = hc_data_ci %>% filter(has_H_error),
                 aes(y = C_val, xmin = H_val - Semi_H, xmax = H_val + Semi_H, 
                     color = Model),
                 height = 0.005, size = 0.8, alpha = 0.8) +
  
  # Add vertical error bars (C_Shannon confidence intervals)
  geom_errorbar(data = hc_data_ci %>% filter(has_C_error),
                aes(x = H_val, ymin = C_val - Semi_C, ymax = C_val + Semi_C, 
                    color = Model),
                width = 0.01, size = 0.8, alpha = 0.8) +
  
  # Add points (colors only)
  geom_point(data = hc_data, aes(x = H_val, y = C_val, color = Model),
             size = 3, stroke = 1) + 
  
  # Styling
  theme_bw(base_family = "serif", base_size = 13) +
  labs(
    x = expression(italic(H)),
    y = expression(italic(C))
  ) +
  scale_color_manual(values = model_colors, labels = model_display_names) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = guide_legend(nrow = 3, override.aes = list(size = 3), byrow = TRUE)
  )

print(p_hc)

ggsave(
  filename = file.path(output_dir, paste0("HC_Plane_with_CI_D", D, ".pdf")),
  plot = p_hc,
  width = 10,
  height = 10,
  bg = "white"
)

# ==============================================================================
# CREATE HC PLANE PLOT WITHOUT LEGEND (FOR COMBINED PLOT)
# ==============================================================================

p_hc_no_legend <- p_hc + theme(legend.position = "none")

# ==============================================================================
# READ TIME SERIES DATA
# ==============================================================================
all_ts_data <- list()

for (model in model_names) {
  filename <- file.path(ts_dir, paste0(model, "_n5000_rep001.csv"))
  
  if (!file.exists(filename)) {
    cat(sprintf("Warning: File not found - %s\n", filename))
    next
  }
  
  ts_data <- read.csv(filename)
  all_ts_data[[model]] <- ts_data
  
  cat(sprintf("Loaded time series: %s\n", model))
}

# ==============================================================================
# CREATE INDIVIDUAL TIME SERIES PLOTS WITH DISPLAY NAMES
# ==============================================================================

create_small_ts_plot <- function(model_name, ts_data, color) {
  # Limit to first 500 points for visibility
  plot_data <- ts_data[1:min(500, nrow(ts_data)), ]
  
  # Get display name
  display_name <- model_display_names[model_name]
  
  p <- ggplot(plot_data, aes(x = time, y = value)) +
    geom_line(color = color, size = 0.6) +
    theme_bw(base_family = "serif", base_size = 9) +
    labs(title = display_name, x = NULL, y = NULL) +  # Use display name
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, 
                                color = color),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 7),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.2, color = "gray90"),
      plot.margin = unit(c(4, 4, 4, 4), "pt"),
      panel.border = element_rect(color = color, size = 1.5),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  
  return(p)
}

# Create all time series plots
ts_plots <- list()
for (model in model_names) {
  if (model %in% names(all_ts_data)) {
    ts_plots[[model]] <- create_small_ts_plot(
      model, 
      all_ts_data[[model]], 
      model_colors[model]
    )
  }
}

# ==============================================================================
# CREATE COMBINED LAYOUT: HC PLANE SURROUNDED BY TIME SERIES
# ==============================================================================

cat("\n=== CHECKING PLOTS ===\n")
cat("Number of time series plots created:", length(ts_plots), "\n")
cat("Models with time series:", names(ts_plots), "\n")

# Create the plot list - ensure we have exactly 16 time series plots
all_plots <- list()
for (i in 1:16) {
  model <- model_names[i]
  if (model %in% names(ts_plots)) {
    all_plots[[i]] <- ts_plots[[model]]
  } else {
    # Create empty plot if time series is missing
    all_plots[[i]] <- ggplot() + 
      theme_void() + 
      labs(title = paste(model_display_names[model], "- NOT FOUND"))
  }
}
# Add HC plane WITHOUT LEGEND as plot 17
all_plots[[17]] <- p_hc_no_legend  # USE THE VERSION WITHOUT LEGEND

cat("Final plot list length:", length(all_plots), "\n")

# Define layout matrix
layout_matrix <- rbind(
  c(NA,  1,  2,  3,  4),     # Row 1: Top 4 time series
  c(5,  17, 17, 17, 6),      # Row 2: TS left, HC center, TS right
  c(7,  17, 17, 17, 8),      # Row 3: TS left, HC center, TS right
  c(9,  17, 17, 17, 10),     # Row 4: TS left, HC center, TS right
  c(11, 17, 17, 17, 12),     # Row 5: TS left, HC center, TS right
  c(13, 14, 15, 16, NA)      # Row 6: Bottom 4 time series
)

cat("Layout matrix dimensions:", nrow(layout_matrix), "x", ncol(layout_matrix), "\n")
cat("Max plot index in layout:", max(layout_matrix, na.rm = TRUE), "\n")

# Create combined plot with MATCHING dimensions (6 HEIGHTS FOR 6 ROWS!)
combined_plot <- grid.arrange(
  grobs = all_plots,
  layout_matrix = layout_matrix,
  widths = c(1.2, 1.2, 1.2, 1.2, 1.2),      # 5 columns
  heights = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2) # 6 HEIGHTS FOR 6 ROWS!
)

ggsave(
  filename = file.path(output_dir, paste0("HC_Plane_with_TimeSeries_Surrounded_D", D, ".pdf")),
  plot = combined_plot,
  width = 20,
  height = 20,  # Increased height for 6 rows
  bg = "white"
)

cat("\n✅ Combined plot saved!\n")

# ==============================================================================
# CREATE SUMMARY INFORMATION WITH DISPLAY NAMES
# ==============================================================================

summary_info <- hc_data %>%
  mutate(
    Model_Display = model_display_names[as.character(Model)],  # Add display name
    Color = model_colors[as.character(Model)],
    Has_H_CI = !is.na(Semi_H) & Semi_H > 0,
    Has_C_CI = !is.na(Semi_C) & Semi_C > 0,
    H_CI_lower = ifelse(Has_H_CI, H_val - Semi_H, NA),
    H_CI_upper = ifelse(Has_H_CI, H_val + Semi_H, NA),
    C_CI_lower = ifelse(Has_C_CI, C_val - Semi_C, NA),
    C_CI_upper = ifelse(Has_C_CI, C_val + Semi_C, NA)
  ) %>%
  select(Model_Display, Color, H_val, C_val, Has_H_CI, H_CI_lower, H_CI_upper, 
         Has_C_CI, C_CI_lower, C_CI_upper) %>%
  rename(Model = Model_Display)  # Rename for cleaner output

write.csv(summary_info,
          file = file.path(output_dir, paste0("HC_Plane_Summary_with_CI_D", D, ".csv")),
          row.names = FALSE)

print(summary_info)

cat("\n✅ Analysis complete!\n")
