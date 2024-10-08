# Load required libraries
library(readxl)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(statcomp)

# Function to compute ordinal patterns for a given embedding dimension
compute_ordinal_patterns <- function(series, emb_dim) {
  n <- length(series)
  patterns <- vector("list", n - emb_dim + 1)
  
  for (i in 1:(n - emb_dim + 1)) {
    subseries <- series[i:(i + emb_dim - 1)]
    patterns[[i]] <- order(subseries)
  }
  return(patterns)
}

# Function to find the frequency of each ordinal pattern
ordinal_pattern_frequencies <- function(series, emb_dim) {
  patterns <- compute_ordinal_patterns(series, emb_dim)
  patterns_str <- sapply(patterns, paste, collapse = "")
  freq_table <- table(patterns_str)
  
  # Convert to data frame for better readability
  freq_df <- as.data.frame(freq_table)
  colnames(freq_df) <- c("Pattern", "Frequency")
  return(freq_df)
}

# Function to calculate probabilities from frequencies
calculate_probabilities <- function(freq_df) {
  total_patterns <- sum(freq_df$Frequency)
  probabilities <- freq_df$Frequency / total_patterns
  freq_df$Probability <- probabilities
  return(freq_df)
}

# Function to calculate permutation entropy (Shannon entropy)
calculate_permutation_entropy <- function(probabilities) {
  entropy <- -sum(probabilities * log2(probabilities), na.rm = TRUE)
  return(entropy / log2(length(probabilities))) # Normalized
}

# Function to calculate Jensen-Shannon divergence
js_divergence <- function(P, Q) {
  M <- 0.5 * (P + Q)
  kl_divergence <- function(P, Q) {
    sum(P * log2(P / Q), na.rm = TRUE)
  }
  jsd <- 0.5 * kl_divergence(P, M) + 0.5 * kl_divergence(Q, M)
  return(jsd)
}

# Function to calculate Jensen-Shannon entropy
calculate_js_entropy <- function(probabilities) {
  uniform_prob <- rep(1 / length(probabilities), length(probabilities))
  jsd <- js_divergence(probabilities, uniform_prob)
  js_entropy <- sqrt(jsd)
  return(js_entropy)
}

# Function to plot time series
plot_time_series <- function(time_series, title) {
  data <- data.frame(Time = seq_along(time_series), Value = time_series)
  ggplot(data, aes(x = Time, y = Value)) +
    geom_line() +
    ggtitle(title) +
    theme_minimal()
}

# Function to plot frequencies of the ordinal patterns (Histogram)
plot_frequencies <- function(freq_df, title_txt) {
  ggplot(freq_df, aes(x = factor(Pattern, levels = unique(Pattern)), y = Probability)) +
    geom_bar(stat = "identity", fill = "gray30", color = "blue", width = 0.7) +
    xlab("Pattern") +
    ylab("Probability") +
    ggtitle(title_txt) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylim(0, 1)
}

# Read the Excel file
file_path <- "Normal Baseline Data.xlsx"

# Define sheet names
sheet_names <- c("Normal 0", "Normal 1", "Normal 2", "Normal 3")

# Initialize lists to store time series data and results
DE_Time <- list()
FE_Time <- list()
results <- list()

# Read each sheet and extract the time series data
for (sheet in sheet_names) {
  data <- read_excel(file_path, sheet = sheet)
  
  # Print column names for debugging
  print(paste("Sheet:", sheet, "Columns:", paste(colnames(data), collapse = ", ")))
  
  if ("DE_Time" %in% colnames(data) & "FE_Time" %in% colnames(data)) {
    DE_Time[[sheet]] <- data$DE_Time
    FE_Time[[sheet]] <- data$FE_Time
    
    # Plot the time series
    print(plot_time_series(data$DE_Time, paste(sheet, "DE_Time")))
    print(plot_time_series(data$FE_Time, paste(sheet, "FE_Time")))
    
    for (emb_dim in 3:6) {
      # Calculate ordinal patterns and their frequencies
      pattern_frequencies_DE <- ordinal_pattern_frequencies(data$DE_Time, emb_dim)
      pattern_frequencies_FE <- ordinal_pattern_frequencies(data$FE_Time, emb_dim)
      
      # Calculate probabilities
      probabilities_df_DE <- calculate_probabilities(pattern_frequencies_DE)
      probabilities_df_FE <- calculate_probabilities(pattern_frequencies_FE)
      probabilities_DE <- probabilities_df_DE$Probability
      probabilities_FE <- probabilities_df_FE$Probability
      
      # Calculate normalized Shannon entropy
      shannon_entropy_DE <- calculate_permutation_entropy(probabilities_DE)
      shannon_entropy_FE <- calculate_permutation_entropy(probabilities_FE)
      
      # Calculate normalized Jensen-Shannon entropy
      js_entropy_DE <- calculate_js_entropy(probabilities_DE)
      js_entropy_FE <- calculate_js_entropy(probabilities_FE)
      
      # Store the results
      results[[paste(sheet, emb_dim, "DE_Time")]] <- list(
        Sheet = sheet,
        Emb_Dim = emb_dim,
        ShannonEntropy = shannon_entropy_DE,
        JSEntropy = js_entropy_DE
      )
      results[[paste(sheet, emb_dim, "FE_Time")]] <- list(
        Sheet = sheet,
        Emb_Dim = emb_dim,
        ShannonEntropy = shannon_entropy_FE,
        JSEntropy = js_entropy_FE
      )
      
      # Plot the frequencies
      print(plot_frequencies(probabilities_df_DE, paste(sheet, "DE_Time", "Embedding Dimension =", emb_dim)))
      print(plot_frequencies(probabilities_df_FE, paste(sheet, "FE_Time", "Embedding Dimension =", emb_dim)))
    }
  } else {
    warning(paste("Sheet", sheet, "does not contain DE_Time or FE_Time columns."))
  }
}

# Prepare data for complexity plane
complexity_data <- do.call(rbind, lapply(names(results), function(name) {
  result <- results[[name]]
  data.frame(
    Sheet = result$Sheet,
    Emb_Dim = result$Emb_Dim,
    ShannonEntropy = result$ShannonEntropy,
    JSEntropy = result$JSEntropy
  )
}))

# Function to plot the complexity plane
plot_complexity_plane <- function(data) {
  ggplot(data, aes(x = ShannonEntropy, y = JSEntropy, color = as.factor(Emb_Dim))) +
    geom_point(size = 3) +
    theme_minimal() +
    xlab("Normalized Shannon Entropy") +
    ylab("Normalized Jensen-Shannon Entropy") +
    ggtitle("Complexity Plane H*C") +
    scale_color_discrete(name = "Embedding Dimension")
}

# Convert results to a data frame for plotting
complexity_data_df <- as.data.frame(complexity_data)

# Plot the complexity plane
print(plot_complexity_plane(complexity_data_df))

# Create a summary table for each time series
summary_table <- do.call(rbind, lapply(names(results), function(name) {
  result <- results[[name]]
  data.frame(
    Sheet = result$Sheet,
    Emb_Dim = result$Emb_Dim,
    ShannonEntropy = result$ShannonEntropy,
    JSEntropy = result$JSEntropy
  )
}))

# Print the summary table
print(summary_table)

# # Create a confusion matrix based on the patterns
# confusion_matrix <- function(data) {
#   table(data$Pattern)
# }
# 
# # Print confusion matrix for each sheet
# for (sheet in sheet_names) {
#   if (sheet %in% names(DE_Time) & sheet %in% names(FE_Time)) {
#     print(confusion_matrix(ordinal_pattern_frequencies(DE_Time[[sheet]], 3)))
#     print(confusion_matrix(ordinal_pattern_frequencies(FE_Time[[sheet]], 3)))
#   }
# }

