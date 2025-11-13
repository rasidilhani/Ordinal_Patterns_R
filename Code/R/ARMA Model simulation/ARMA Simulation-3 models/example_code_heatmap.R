# ARMA Time Series Heatmap Analysis
# Generate density heatmaps for H-C plane by case and sample size

library(readxl)
library(ggplot2)
library(StatOrdPattHxC)
library(viridis)
library(MASS)
library(dplyr)
library(gridExtra)

# Load boundary data
data("LinfLsup")

# Set data path
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results"

# Read all three case files
case1_file <- file.path(data_path, "case1.xlsx")
case2_file <- file.path(data_path, "case2.xlsx")
case3_file <- file.path(data_path, "case3.xlsx")

# Read data
cat("Reading data files...\n")
case1_data <- read_excel(case1_file)
case2_data <- read_excel(case2_file)
case3_data <- read_excel(case3_file)

# Function to create density heatmap
create_heatmap <- function(data, sample_size, case_name, 
                           xlim_range = NULL, ylim_range = NULL) {
  
  # Filter by sample size
  data_filtered <- data %>% filter(n == sample_size)
  
  cat(paste("Processing", case_name, "- n =", sample_size, 
            "- Observations:", nrow(data_filtered), "\n"))
  
  # Use H_Star and C_Star columns (or Entropy/Complexity if different names)
  if("H_Star" %in% names(data_filtered)) {
    h_col <- "H_Star"
    c_col <- "C_Star"
  } else if("Emblematic_H_Star" %in% names(data_filtered)) {
    h_col <- "Emblematic_H_Star"
    c_col <- "C_Star"
  } else {
    h_col <- "Entropy"
    c_col <- "Complexity"
  }
  
  # Calculate 2D kernel density estimation
  h_vals <- data_filtered[[h_col]]
  c_vals <- data_filtered[[c_col]]
  
  # Remove any NA values
  valid_idx <- !is.na(h_vals) & !is.na(c_vals)
  h_vals <- h_vals[valid_idx]
  c_vals <- c_vals[valid_idx]
  
  if(length(h_vals) < 3) {
    warning(paste("Not enough data points for", case_name, "n =", sample_size))
    return(NULL)
  }
  
  # Determine limits
  if(is.null(xlim_range)) {
    xlim_range <- range(h_vals, na.rm = TRUE)
    xlim_range <- c(xlim_range[1] - 0.02, xlim_range[2] + 0.02)
  }
  if(is.null(ylim_range)) {
    ylim_range <- range(c_vals, na.rm = TRUE)
    ylim_range <- c(max(0, ylim_range[1] - 0.01), ylim_range[2] + 0.01)
  }
  
  # Create density estimation
  tryCatch({
    kde <- kde2d(h_vals, c_vals, n = 100,
                 lims = c(xlim_range, ylim_range))
    
    # Convert to data frame for ggplot
    kde_df <- expand.grid(H = kde$x, C = kde$y)
    kde_df$density <- as.vector(kde$z)
    
    # Create plot
    p <- ggplot() +
      # Add density heatmap
      geom_tile(data = kde_df, aes(x = H, y = C, fill = density)) +
      scale_fill_viridis(option = "plasma", name = "Density") +
      # Add boundary curves
      geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
                aes(x = H, y = C), col = "gray30", size = 1, alpha = 0.7) +
      geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
                aes(x = H, y = C), col = "gray30", size = 1, alpha = 0.7) +
      # Add actual data points (optional, semi-transparent)
      geom_point(data = data_filtered, 
                 aes_string(x = h_col, y = c_col),
                 color = "white", alpha = 0.3, size = 0.5) +
      xlim(xlim_range) +
      ylim(ylim_range) +
      xlab(expression(italic(H))) +
      ylab(expression(italic(C))) +
      ggtitle(paste(case_name, "- n =", sample_size)) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        legend.position = "right",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )
    
    return(p)
    
  }, error = function(e) {
    warning(paste("Error creating density for", case_name, "n =", sample_size, ":", e$message))
    return(NULL)
  })
}

# Function to create contour plot (alternative visualization)
create_contour_plot <- function(data, sample_size, case_name,
                                xlim_range = NULL, ylim_range = NULL) {
  
  # Filter by sample size
  data_filtered <- data %>% filter(n == sample_size)
  
  # Detect column names
  if("H_Star" %in% names(data_filtered)) {
    h_col <- "H_Star"
    c_col <- "C_Star"
  } else if("Emblematic_H_Star" %in% names(data_filtered)) {
    h_col <- "Emblematic_H_Star"
    c_col <- "C_Star"
  } else {
    h_col <- "Entropy"
    c_col <- "Complexity"
  }
  
  h_vals <- data_filtered[[h_col]]
  c_vals <- data_filtered[[c_col]]
  
  # Remove NA values
  valid_idx <- !is.na(h_vals) & !is.na(c_vals)
  h_vals <- h_vals[valid_idx]
  c_vals <- c_vals[valid_idx]
  
  if(length(h_vals) < 3) return(NULL)
  
  # Determine limits
  if(is.null(xlim_range)) {
    xlim_range <- range(h_vals, na.rm = TRUE)
    xlim_range <- c(xlim_range[1] - 0.02, xlim_range[2] + 0.02)
  }
  if(is.null(ylim_range)) {
    ylim_range <- range(c_vals, na.rm = TRUE)
    ylim_range <- c(max(0, ylim_range[1] - 0.01), ylim_range[2] + 0.01)
  }
  
  # Create plot with 2D density contours
  p <- ggplot(data_filtered, aes_string(x = h_col, y = c_col)) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.7) +
    scale_fill_viridis(option = "plasma", name = "Density") +
    # Add boundary curves
    geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
              aes(x = H, y = C), col = "gray30", size = 1, alpha = 0.7, inherit.aes = FALSE) +
    geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
              aes(x = H, y = C), col = "gray30", size = 1, alpha = 0.7, inherit.aes = FALSE) +
    geom_point(alpha = 0.4, size = 0.8, color = "white") +
    xlim(xlim_range) +
    ylim(ylim_range) +
    xlab(expression(italic(H))) +
    ylab(expression(italic(C))) +
    ggtitle(paste(case_name, "- n =", sample_size)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  return(p)
}

# Generate heatmaps for all cases and sample sizes
cat("\n=== GENERATING HEATMAPS ===\n\n")

# Case 1 heatmaps
cat("Case 1:\n")
case1_n500 <- create_heatmap(case1_data, 500, "Case 1")
case1_n1000 <- create_heatmap(case1_data, 1000, "Case 1")

# Case 2 heatmaps
cat("\nCase 2:\n")
case2_n500 <- create_heatmap(case2_data, 500, "Case 2")
case2_n1000 <- create_heatmap(case2_data, 1000, "Case 2")

# Case 3 heatmaps
cat("\nCase 3:\n")
case3_n500 <- create_heatmap(case3_data, 500, "Case 3")
case3_n1000 <- create_heatmap(case3_data, 1000, "Case 3")

# Save individual heatmaps
cat("\n=== SAVING INDIVIDUAL HEATMAPS ===\n")

if(!is.null(case1_n500)) {
  ggsave("Case1_n500_heatmap.pdf", case1_n500, width = 8, height = 6, dpi = 300)
  ggsave("Case1_n500_heatmap.png", case1_n500, width = 8, height = 6, dpi = 300)
}
if(!is.null(case1_n1000)) {
  ggsave("Case1_n1000_heatmap.pdf", case1_n1000, width = 8, height = 6, dpi = 300)
  ggsave("Case1_n1000_heatmap.png", case1_n1000, width = 8, height = 6, dpi = 300)
}

if(!is.null(case2_n500)) {
  ggsave("Case2_n500_heatmap.pdf", case2_n500, width = 8, height = 6, dpi = 300)
  ggsave("Case2_n500_heatmap.png", case2_n500, width = 8, height = 6, dpi = 300)
}
if(!is.null(case2_n1000)) {
  ggsave("Case2_n1000_heatmap.pdf", case2_n1000, width = 8, height = 6, dpi = 300)
  ggsave("Case2_n1000_heatmap.png", case2_n1000, width = 8, height = 6, dpi = 300)
}

if(!is.null(case3_n500)) {
  ggsave("Case3_n500_heatmap.pdf", case3_n500, width = 8, height = 6, dpi = 300)
  ggsave("Case3_n500_heatmap.png", case3_n500, width = 8, height = 6, dpi = 300)
}
if(!is.null(case3_n1000)) {
  ggsave("Case3_n1000_heatmap.pdf", case3_n1000, width = 8, height = 6, dpi = 300)
  ggsave("Case3_n1000_heatmap.png", case3_n1000, width = 8, height = 6, dpi = 300)
}

# Create combined plots for each case
cat("\n=== CREATING COMBINED PLOTS ===\n")

# Combined plot for Case 1
if(!is.null(case1_n500) && !is.null(case1_n1000)) {
  case1_combined <- grid.arrange(case1_n500, case1_n1000, ncol = 2,
                                 top = "Case 1: Entropy-Complexity Density Maps")
  ggsave("Case1_combined_heatmaps.pdf", case1_combined, width = 14, height = 6, dpi = 300)
  ggsave("Case1_combined_heatmaps.png", case1_combined, width = 14, height = 6, dpi = 300)
}

# Combined plot for Case 2
if(!is.null(case2_n500) && !is.null(case2_n1000)) {
  case2_combined <- grid.arrange(case2_n500, case2_n1000, ncol = 2,
                                 top = "Case 2: Entropy-Complexity Density Maps")
  ggsave("Case2_combined_heatmaps.pdf", case2_combined, width = 14, height = 6, dpi = 300)
  ggsave("Case2_combined_heatmaps.png", case2_combined, width = 14, height = 6, dpi = 300)
}

# Combined plot for Case 3
if(!is.null(case3_n500) && !is.null(case3_n1000)) {
  case3_combined <- grid.arrange(case3_n500, case3_n1000, ncol = 2,
                                 top = "Case 3: Entropy-Complexity Density Maps")
  ggsave("Case3_combined_heatmaps.pdf", case3_combined, width = 14, height = 6, dpi = 300)
  ggsave("Case3_combined_heatmaps.png", case3_combined, width = 14, height = 6, dpi = 300)
}

# Create mega combined plot (all 6 heatmaps)
all_plots <- list(case1_n500, case1_n1000, case2_n500, case2_n1000, case3_n500, case3_n1000)
all_plots <- all_plots[!sapply(all_plots, is.null)]

if(length(all_plots) == 6) {
  mega_combined <- grid.arrange(grobs = all_plots, ncol = 2, nrow = 3,
                                top = "ARMA Time Series: Entropy-Complexity Density Maps")
  ggsave("All_Cases_Combined_heatmaps.pdf", mega_combined, width = 14, height = 18, dpi = 300)
  ggsave("All_Cases_Combined_heatmaps.png", mega_combined, width = 14, height = 18, dpi = 300)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Generated files:\n")
cat("- Individual heatmaps: Case[1-3]_n[500/1000]_heatmap.pdf/png\n")
cat("- Combined by case: Case[1-3]_combined_heatmaps.pdf/png\n")
cat("- All combined: All_Cases_Combined_heatmaps.pdf/png\n")