library(StatOrdPattHxC)
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

# ==============================================================================
# PARAMETERS
# ==============================================================================
D <- 4  # Embedding dimension

ts_dir <- here("Data", "Hybrid model_data", "D4_Data", "All_Models_TimeSeries_D4")
output_dir <- here("Plots", "Hybrid Model plots", "D4 plots")

# ==============================================================================
# DEFINE MODEL ORDER, COLORS, AND SHAPES (CONSISTENT FOR ALL FUTURE WORK)
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

# FIXED COLOR PALETTE (use these colors consistently in all future plots)
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

# FIXED SHAPE PALETTE (use these shapes consistently in all future plots)
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
# READ TIME SERIES AND CALCULATE ORDINAL PATTERN PROBABILITIES
# ==============================================================================
all_op_probs <- list()

for (model in model_names) {
  # Construct filename
  filename <- file.path(ts_dir, paste0(model, "_n5000_rep001.csv"))
  
  # Check if file exists
  if (!file.exists(filename)) {
    cat(sprintf("Warning: File not found - %s\n", filename))
    next
  }
  
  # Read time series
  ts_data <- read.csv(filename)
  ts_values <- ts_data$value
  
  # Calculate ordinal pattern probabilities
  op_prob <- OPprob(ts_values, emb = D)
  
  # Store results with ORIGINAL model name (for data consistency)
  all_op_probs[[model]] <- data.frame(
    Model = model,
    Probability = op_prob,
    Pattern = 1:length(op_prob)
  )
  
  cat(sprintf("Processed: %s (n = %d)\n", model, length(ts_values)))
}

# Combine all results
combined_op <- bind_rows(all_op_probs)

# Factor the Model variable to control plotting order (use original names)
combined_op$Model <- factor(combined_op$Model, levels = model_names)

# Add display names as a separate column for plotting
combined_op$Model_Display <- model_display_names[as.character(combined_op$Model)]
combined_op$Model_Display <- factor(combined_op$Model_Display, 
                                    levels = model_display_names[model_names])

# ==============================================================================
# CREATE FACETED ORDINAL PATTERNS PLOT
# ==============================================================================

# Calculate number of possible patterns
n_patterns <- factorial(D)

# Create the faceted plot (using display names)
p_faceted <- ggplot(combined_op, aes(x = Pattern, y = Probability, fill = Model)) +
  geom_bar(stat = "identity", position = "identity", width = 0.8) +
  facet_wrap(~ Model_Display, ncol = 3, scales = "free_y") +  # Use display names
  theme_bw(base_family = "serif", base_size = 13) +
  labs(
    title = paste0("Ordinal Pattern Distributions (D = ", D, ", n = 5000)"),
    x = "Ordinal Pattern",
    y = "Probability"
  ) +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, n_patterns, by = 1)) +
  scale_fill_manual(values = model_colors)

# Save the faceted plot
ggsave(
  filename = file.path(output_dir, paste0("Ordinal_Patterns_Faceted_D", D,".png")),
  plot = p_faceted,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, paste0("Ordinal_Patterns_Faceted_D", D,".pdf")),
  plot = p_faceted,
  width = 14,
  height = 10
)

cat(sprintf("\nFaceted plot saved to: %s\n", 
            file.path(output_dir, paste0("Ordinal_Patterns_Faceted_D", D,".png"))))

# Display the plot
print(p_faceted)

# ==============================================================================
# CREATE COMBINED PLOT (All models overlaid with lines and points)
# ==============================================================================

p_combined <- ggplot(combined_op, aes(x = Pattern, y = Probability, 
                                      color = Model, shape = Model, group = Model)) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw(base_family = "serif", base_size = 13) +
  labs(
    title = paste0("Ordinal Pattern Distributions - All Models (D = ", D, ", n = 5000)"),
    x = "Ordinal Pattern",
    y = "Probability",
    color = "Model",
    shape = "Model"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.key.size = unit(0.8, "cm"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(1, n_patterns, by = 1)) +
  scale_color_manual(values = model_colors, labels = model_display_names) +  # Add labels
  scale_shape_manual(values = model_shapes, labels = model_display_names) +  # Add labels
  guides(
    color = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1))
  )

ggsave(
  filename = file.path(output_dir, paste0("Ordinal_Patterns_Hybrid_Models_D",D,".pdf")),
  plot = p_combined,
  width = 14,
  height = 8
)

print(p_combined)

# ==============================================================================
# CREATE SUMMARY TABLE OF ORDINAL PATTERNS
# ==============================================================================
summary_table <- combined_op %>%
  group_by(Model, Model_Display) %>%  # Include both for clarity
  summarise(
    N_Nonzero_Patterns = sum(Probability > 0),
    Total_Patterns = n_patterns,
    Coverage_Percent = round(100 * sum(Probability > 0) / n_patterns, 2),
    Max_Probability = round(max(Probability), 4),
    Min_Nonzero_Prob = round(min(Probability[Probability > 0]), 6),
    Shannon_Entropy = round(-sum(ifelse(Probability > 0, 
                                        Probability * log(Probability), 0)), 4),
    .groups = 'drop'
  ) %>%
  select(Model_Display, everything(), -Model) %>%  # Put display name first, remove original
  rename(Model = Model_Display)  # Rename for cleaner output

print(summary_table)

# Save summary table
write.csv(summary_table, 
          file = file.path(output_dir, paste0("Ordinal_Patterns_Summary_D",D,".csv")),
          row.names = FALSE)

cat(sprintf("\nSummary table saved to: %s\n", 
            file.path(output_dir, paste0("Ordinal_Patterns_Summary_D",D,".csv"))))