# Required libraries
library(readxl)
library(ggplot2)

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
  patterns_str <- sapply(patterns, paste, collapse = "-")
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
  return(entropy)
}

# Function to calculate Kullback-Leibler divergence
kl_divergence <- function(P, Q) {
  divergence <- sum(P * log2(P / Q), na.rm = TRUE)
  return(divergence)
}

# Function to calculate Jensen-Shannon divergence
js_divergence <- function(P, Q) {
  M <- 0.5 * (P + Q)
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

# Function to plot the frequencies of the ordinal patterns (Histogram)
plot_frequencies <- function(freq_df, title_txt) {
  ggplot(freq_df, aes(x = reorder(Pattern, -Probability), y = Probability)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    xlab("Ordinal Patterns") +
    ylab(expression(P(π))) +
    ggtitle(title_txt) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Function to obtain a series from modified series
get_modified_series <- function(series, method) {
  if (method == "Complete") {
    return(series)
  } else if (method == "Time Ordered") {
    return(series)
  } else if (method == "Random") {
    set.seed(123)
    return(sample(series))
  } else if (method == "Data-Driven") {
    return(series)
  }
}

# Example usage
CAC <- read_excel("Load/CAC.xlsx")
View(CAC)
series <- CAC[[2]]  # Assuming the second column contains the series
timestamp <- CAC[[1]]

# Define the embedding dimensions and titles
embedding_dimensions <- 3:7
title <- c("Complete", "Time Ordered", "Random", "Data-Driven")

# Initialize list to store complexity results
complexity_results <- list()

# Calculate and plot for each embedding dimension and each method
for (emb_dim in embedding_dimensions) {
  for (method in title) {
    modified_series <- get_modified_series(series, method)
    pattern_frequencies <- ordinal_pattern_frequencies(modified_series, emb_dim)
    probabilities_df <- calculate_probabilities(pattern_frequencies)
    probabilities <- probabilities_df$Probability
    shannon_entropy <- calculate_permutation_entropy(probabilities)
    js_entropy <- calculate_js_entropy(probabilities)
    
    complexity_results[[paste(emb_dim, method, sep = "_")]] <- data.frame(
      EmbeddingDimension = emb_dim,
      Method = method,
      PermutationEntropy = shannon_entropy,
      JensenShannonEntropy = js_entropy
    )
    
    print(paste("Embedding Dimension:", emb_dim, "Method:", method))
    print(paste("Permutation Entropy:", shannon_entropy))
    print(paste("Jensen-Shannon Entropy:", js_entropy))
    
    # Plot frequencies for a specific embedding dimension and method
    print(plot_frequencies(probabilities_df, paste("Embedding Dimension:", emb_dim, "Method:", method)))
  }
}

# Combine all results into a single data frame
all_complexity_results <- do.call(rbind, complexity_results)

# Function to plot Complexity-Entropy Causality Plane
plot_cecp <- function(cecp_results) {
  ggplot(cecp_results, aes(x = PermutationEntropy, y = JensenShannonEntropy, color = Method)) +
    geom_point(size = 3) +
    geom_text(aes(label = EmbeddingDimension), vjust = -1, hjust = 1) +
    xlab("Permutation Entropy") +
    ylab("Jensen-Shannon Entropy") +
    ggtitle("Complexity-Entropy Causality Plane") +
    theme_minimal()
}

# Plot Complexity-Entropy Causality Plane
plot_cecp(all_complexity_results)

