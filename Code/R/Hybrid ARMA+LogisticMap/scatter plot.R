#--------------------------------------------------------------------------------------------
#This script generates scatter plots of entropy vs. complexity for various ARMA models,
#using extended entropy measures (Shannon, Rényi, Tsallis, Fisher).
#======================================================================
#----------------------------------------------------------------------------------------------
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
hc_data_path <- here("Data", "Hybrid model_data", "D4_Data", "Combined_All_Models_HC_Results_D4.xlsx")
output_dir <- here("Plots", "Hybrid Model plots", "D4 plots")

# Sample sizes to process
sample_sizes <- c(1000, 5000)

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
# LOAD HC DATA (ALL REPETITIONS)
# ==============================================================================
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Check available sample sizes
available_n <- unique(hc_data_all$n)
cat("Available sample sizes in data:", available_n, "\n")

# Verify we have the expected sample sizes
if (!all(sample_sizes %in% available_n)) {
  cat("WARNING: Not all expected sample sizes found in data!\n")
  cat("Expected:", sample_sizes, "\n")
  cat("Found:", available_n, "\n")
  sample_sizes <- intersect(sample_sizes, available_n)
  cat("Will process:", sample_sizes, "\n")
}

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon, Semi_HS, Semi_CS,
         H_Renyi, C_Renyi, Semi_HR, 
         H_Tsallis, C_Tsallis, Semi_HT,
         H_Fisher, C_Fisher, Semi_HF)

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

cat("\nTotal data points loaded:", nrow(hc_data), "\n")
cat("Sample sizes found:", unique(hc_data$n), "\n")
cat("Repetitions per sample size:\n")
print(table(hc_data$n))

# ==============================================================================
# LOAD LinfLsup BOUNDARIES (for Shannon only)
# ==============================================================================
data("LinfLsup")

# Full boundary curve
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# Separate upper and lower boundaries
boundary_lower <- LinfLsup_full %>% filter(Side == "Lower")
boundary_upper <- LinfLsup_full %>% filter(Side == "Upper")

cat("\nBoundary curves prepared for Shannon entropy\n")
cat("Boundary points - Lower:", nrow(boundary_lower), "Upper:", nrow(boundary_upper), "\n")

# ==============================================================================
# ENTROPY CONFIGURATIONS
# ==============================================================================
entropy_configs <- list(
  Shannon = list(
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    Semi_H_col = "Semi_HS",
    Semi_C_col = "Semi_CS",
    add_bounds = TRUE,
    title = "Shannon Entropy-Complexity Plane",
    xlab = expression(italic(H)[Shannon]),
    ylab = expression(italic(C)[Shannon]),
    xlim = c(0, 1),
    ylim = c(0, 0.5)
  ),
  Renyi = list(
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    Semi_H_col = "Semi_HR",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Rényi Entropy-Complexity Plane",
    xlab = expression(italic(H)[Rényi]),
    ylab = expression(italic(C)[Rényi]),
    xlim = NULL,
    ylim = NULL
  ),
  Tsallis = list(
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    Semi_H_col = "Semi_HT",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Tsallis Entropy-Complexity Plane",
    xlab = expression(italic(H)[Tsallis]),
    ylab = expression(italic(C)[Tsallis]),
    xlim = NULL,
    ylim = NULL
  ),
  Fisher = list(
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    Semi_H_col = "Semi_HF",
    Semi_C_col = NA,
    add_bounds = FALSE,
    title = "Fisher Entropy-Complexity Plane",
    xlab = expression(italic(H)[Fisher]),
    ylab = expression(italic(C)[Fisher]),
    xlim = NULL,
    ylim = NULL
  )
)

# ==============================================================================
# CREATE HC PLANE PLOTS FOR EACH SAMPLE SIZE AND ENTROPY TYPE
# ==============================================================================

for (n_val in sample_sizes) {
  
  # Create subdirectory for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n========== Processing Sample Size n = %d ==========\n", n_val))
  
  # Filter data for this sample size
  hc_data_n <- hc_data %>%
    filter(n == n_val)
  
  cat(sprintf("Data points for n=%d: %d\n", n_val, nrow(hc_data_n)))
  cat(sprintf("Models: %d, Max repetitions: %d\n", 
              length(unique(hc_data_n$Model)),
              max(hc_data_n$rep, na.rm = TRUE)))
  
  for (ent_name in names(entropy_configs)) {
    
    cfg <- entropy_configs[[ent_name]]
    
    cat(sprintf("\n=== Creating %s HC plane plot for n=%d ===\n", ent_name, n_val))
    
    # Extract relevant columns
    H_col <- cfg$H_col
    C_col <- cfg$C_col
    Semi_H_col <- cfg$Semi_H_col
    Semi_C_col <- cfg$Semi_C_col
    
    # Check if columns exist
    if (!all(c(H_col, C_col) %in% names(hc_data_n))) {
      cat(sprintf("⚠️ Skipping %s - columns not found\n", ent_name))
      next
    }
    
    # Prepare data with renamed columns for easier plotting
    plot_data <- hc_data_n %>%
      select(Model, rep, all_of(c(H_col, C_col))) %>%
      rename(H_val = !!H_col, C_val = !!C_col)
    
    # Add semi-length columns if they exist
    if (!is.na(Semi_H_col) && Semi_H_col %in% names(hc_data_n)) {
      plot_data$Semi_H <- hc_data_n[[Semi_H_col]]
    } else {
      plot_data$Semi_H <- NA
    }
    
    if (!is.na(Semi_C_col) && Semi_C_col %in% names(hc_data_n)) {
      plot_data$Semi_C <- hc_data_n[[Semi_C_col]]
    } else {
      plot_data$Semi_C <- NA
    }
    
    # Remove rows with NA in H or C
    plot_data <- plot_data %>%
      filter(!is.na(H_val), !is.na(C_val), is.finite(H_val), is.finite(C_val))
    
    if (nrow(plot_data) == 0) {
      cat(sprintf("⚠️ Skipping %s - no valid data\n", ent_name))
      next
    }
    
    cat(sprintf("Valid data points: %d\n", nrow(plot_data)))
    
    # Prepare error bar data (only for first repetition to avoid clutter)
    plot_data_ci <- plot_data %>%
      filter(rep == 1) %>%
      mutate(
        has_H_error = !is.na(Semi_H) & Semi_H > 0 & is.finite(Semi_H),
        has_C_error = !is.na(Semi_C) & Semi_C > 0 & is.finite(Semi_C)
      )
    
    # Create the plot
    p <- ggplot()
    
    # Add boundaries (Shannon only)
    if (cfg$add_bounds) {
      p <- p +
        geom_line(data = boundary_lower, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
        geom_line(data = boundary_upper, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8)
    }
    
    # Add error bars (only for rep 1 to avoid visual clutter)
    # Horizontal error bars (Entropy confidence intervals)
    if (any(plot_data_ci$has_H_error, na.rm = TRUE)) {
      p <- p +
        geom_errorbarh(data = plot_data_ci %>% filter(has_H_error),
                       aes(y = C_val, xmin = H_val - Semi_H, xmax = H_val + Semi_H, 
                           color = Model),
                       height = 0.002, size = 0.6, alpha = 0.5, show.legend = FALSE)
      cat(sprintf("  Added H error bars for %d models\n", 
                  sum(plot_data_ci$has_H_error, na.rm = TRUE)))
    } else {
      cat("  No H error bars (Semi_H not available)\n")
    }
    
    # Vertical error bars (Complexity confidence intervals)
    if (any(plot_data_ci$has_C_error, na.rm = TRUE)) {
      p <- p +
        geom_errorbar(data = plot_data_ci %>% filter(has_C_error),
                      aes(x = H_val, ymin = C_val - Semi_C, ymax = C_val + Semi_C, 
                          color = Model),
                      width = 0.01, size = 0.6, alpha = 0.5, show.legend = FALSE)
      cat(sprintf("  Added C error bars for %d models\n", 
                  sum(plot_data_ci$has_C_error, na.rm = TRUE)))
    } else {
      cat("  No C error bars (Semi_C not available)\n")
    }
    
    # Add points (all repetitions)
    p <- p +
      geom_point(data = plot_data, 
                 aes(x = H_val, y = C_val, 
                     color = Model, 
                     shape = Model,
                     fill = Model),
                 size = 3, stroke = 1.2, alpha = 0.7) +
      
      # Styling
      theme_bw(base_family = "serif", base_size = 13) +
      labs(
        title = cfg$title,
        subtitle = paste0("D = ", D, ", n = ", n_val, 
                          " (", max(plot_data$rep, na.rm = TRUE), " repetitions)"),
        x = cfg$xlab,
        y = cfg$ylab
      ) +
      
      # Use display names in legend
      scale_color_manual(
        values = model_colors,
        labels = model_display_names,
        name = "Model",
        guide = guide_legend(
          override.aes = list(
            size = 4,
            alpha = 1,
            stroke = 1.5
          )
        )
      ) +
      scale_shape_manual(values = model_shapes, labels = model_display_names, guide = "none") +
      scale_fill_manual(values = model_colors, labels = model_display_names, guide = "none") +
      
      theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.7, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.3),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      )
    
    # Set axis limits if specified
    if (!is.null(cfg$xlim)) {
      p <- p + scale_x_continuous(limits = cfg$xlim, 
                                  breaks = seq(cfg$xlim[1], cfg$xlim[2], length.out = 5))
    }
    if (!is.null(cfg$ylim)) {
      p <- p + scale_y_continuous(limits = cfg$ylim, 
                                  breaks = seq(cfg$ylim[1], cfg$ylim[2], length.out = 5))
    }
    
    # Save plot
    filename_pdf <- file.path(n_dir, paste0("HC_Plane_", ent_name, "_D", D, "_n", n_val, "_with_CI.pdf"))
    
    ggsave(
      filename = filename_pdf,
      plot = p,
      width = 12,
      height = 8,
      bg = "white"
    )
    
    cat(sprintf("✅ %s HC plane plot saved for n=%d\n", ent_name, n_val))
    
    # Print the plot
    print(p)
  }
} 
#----------------------------------------------------------------------------------------------
# This script has no CI in all scatter plots. It generates scatter plots of entropy vs. complexity 
#for various ARMA models, using extended entropy measures (Shannon, Rényi, Tsallis, Fisher). 
#----------------------------------------------------------------------------------------------

library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readxl)
library(here)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension

# Paths
hc_data_path <- here("Data", "Hybrid model_data", "D4_Data", "Combined_All_Models_HC_Results_D4.xlsx")
output_dir <- here("Plots", "Hybrid Model plots", "D4 plots")

# Sample sizes to process
sample_sizes <- c(1000, 5000)

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
# LOAD HC DATA (ALL REPETITIONS)
# ==============================================================================
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Check available sample sizes
available_n <- unique(hc_data_all$n)
cat("Available sample sizes in data:", available_n, "\n")

# Verify we have the expected sample sizes
if (!all(sample_sizes %in% available_n)) {
  cat("WARNING: Not all expected sample sizes found in data!\n")
  cat("Expected:", sample_sizes, "\n")
  cat("Found:", available_n, "\n")
  sample_sizes <- intersect(sample_sizes, available_n)
  cat("Will process:", sample_sizes, "\n")
}

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon,
         H_Renyi, C_Renyi,
         H_Tsallis, C_Tsallis,
         H_Fisher, C_Fisher)

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)

cat("\nTotal data points loaded:", nrow(hc_data), "\n")
cat("Sample sizes found:", unique(hc_data$n), "\n")
cat("Repetitions per sample size:\n")
print(table(hc_data$n))

# ==============================================================================
# LOAD LinfLsup BOUNDARIES (for Shannon only)
# ==============================================================================
data("LinfLsup")

# Full boundary curve
LinfLsup_full <- LinfLsup %>%
  filter(Dimension == as.character(D))

# Separate upper and lower boundaries
boundary_lower <- LinfLsup_full %>% filter(Side == "Lower")
boundary_upper <- LinfLsup_full %>% filter(Side == "Upper")

cat("\nBoundary curves prepared for Shannon entropy\n")
cat("Boundary points - Lower:", nrow(boundary_lower), "Upper:", nrow(boundary_upper), "\n")

# ==============================================================================
# ENTROPY CONFIGURATIONS
# ==============================================================================
entropy_configs <- list(
  Shannon = list(
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    add_bounds = TRUE,
    title = "Shannon Entropy-Complexity Plane",
    xlab = expression(italic(H)[Shannon]),
    ylab = expression(italic(C)[Shannon]),
    xlim = c(0, 1),
    ylim = c(0, 0.5)
  ),
  Renyi = list(
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    add_bounds = FALSE,
    title = "Rényi Entropy-Complexity Plane",
    xlab = expression(italic(H)[Rényi]),
    ylab = expression(italic(C)[Rényi]),
    xlim = NULL,
    ylim = NULL
  ),
  Tsallis = list(
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    add_bounds = FALSE,
    title = "Tsallis Entropy-Complexity Plane",
    xlab = expression(italic(H)[Tsallis]),
    ylab = expression(italic(C)[Tsallis]),
    xlim = NULL,
    ylim = NULL
  ),
  Fisher = list(
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    add_bounds = FALSE,
    title = "Fisher Entropy-Complexity Plane",
    xlab = expression(italic(H)[Fisher]),
    ylab = expression(italic(C)[Fisher]),
    xlim = NULL,
    ylim = NULL
  )
)

# ==============================================================================
# CREATE HC PLANE PLOTS FOR EACH SAMPLE SIZE AND ENTROPY TYPE
# ==============================================================================

for (n_val in sample_sizes) {
  
  # Create subdirectory for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n========== Processing Sample Size n = %d ==========\n", n_val))
  
  # Filter data for this sample size
  hc_data_n <- hc_data %>%
    filter(n == n_val)
  
  cat(sprintf("Data points for n=%d: %d\n", n_val, nrow(hc_data_n)))
  cat(sprintf("Models: %d, Max repetitions: %d\n", 
              length(unique(hc_data_n$Model)),
              max(hc_data_n$rep, na.rm = TRUE)))
  
  for (ent_name in names(entropy_configs)) {
    
    cfg <- entropy_configs[[ent_name]]
    
    cat(sprintf("\n=== Creating %s HC plane plot for n=%d ===\n", ent_name, n_val))
    
    # Extract relevant columns
    H_col <- cfg$H_col
    C_col <- cfg$C_col
    
    # Check if columns exist
    if (!all(c(H_col, C_col) %in% names(hc_data_n))) {
      cat(sprintf("⚠️ Skipping %s - columns not found\n", ent_name))
      next
    }
    
    # Prepare data with renamed columns for easier plotting
    plot_data <- hc_data_n %>%
      select(Model, rep, all_of(c(H_col, C_col))) %>%
      rename(H_val = !!H_col, C_val = !!C_col)
    
    # Remove rows with NA in H or C
    plot_data <- plot_data %>%
      filter(!is.na(H_val), !is.na(C_val), is.finite(H_val), is.finite(C_val))
    
    if (nrow(plot_data) == 0) {
      cat(sprintf("⚠️ Skipping %s - no valid data\n", ent_name))
      next
    }
    
    cat(sprintf("Valid data points: %d\n", nrow(plot_data)))
    
    # Create the plot
    p <- ggplot()
    
    # Add boundaries (Shannon only)
    if (cfg$add_bounds) {
      p <- p +
        geom_line(data = boundary_lower, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8) +
        geom_line(data = boundary_upper, aes(x = H, y = C), 
                  color = "gray50", linetype = "dashed", size = 1.0, alpha = 0.8)
    }
    
    # Add points (all repetitions) - NO CONFIDENCE INTERVALS
    p <- p +
      geom_point(data = plot_data, 
                 aes(x = H_val, y = C_val, 
                     color = Model, 
                     shape = Model,
                     fill = Model),
                 size = 3, stroke = 1.2, alpha = 0.7) +
      
      # Styling
      theme_bw(base_family = "serif", base_size = 13) +
      labs(
        subtitle = paste0("D = ", D, ", n = ", n_val),
        x = cfg$xlab,
        y = cfg$ylab
      ) +
      
      # Use display names in legend
      scale_color_manual(
        values = model_colors,
        labels = model_display_names,
        name = "Model",
        guide = guide_legend(
          override.aes = list(
            size = 4,
            alpha = 1,
            stroke = 1.5
          )
        )
      ) +
      scale_shape_manual(values = model_shapes, labels = model_display_names, guide = "none") +
      scale_fill_manual(values = model_colors, labels = model_display_names, guide = "none") +
      
      theme(
        plot.subtitle = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.7, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90", size = 0.3),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white")
      )
    
    # Set axis limits if specified
    if (!is.null(cfg$xlim)) {
      p <- p + scale_x_continuous(limits = cfg$xlim, 
                                  breaks = seq(cfg$xlim[1], cfg$xlim[2], length.out = 5))
    }
    if (!is.null(cfg$ylim)) {
      p <- p + scale_y_continuous(limits = cfg$ylim, 
                                  breaks = seq(cfg$ylim[1], cfg$ylim[2], length.out = 5))
    }
    
    # Save plot
    filename_pdf <- file.path(n_dir, paste0("HC_Plane_", ent_name, "_D", D, "_n", n_val, "_no_CI.pdf"))
    
    ggsave(
      filename = filename_pdf,
      plot = p,
      width = 12,
      height = 8,
      bg = "white"
    )
    
    cat(sprintf("✅ %s HC plane plot saved for n=%d\n", ent_name, n_val))
    
    # Print the plot
    print(p)
  }
}