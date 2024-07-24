# Load necessary libraries
library(ggplot2)
library(readxl)

# Function to calculate permutation entropy
calculate_permutation_entropy <- function(probabilities) {
  entropy <- -sum(probabilities * log2(probabilities), na.rm = TRUE)
  return(entropy)
}

# Function to calculate KL divergence
kl_divergence <- function(P, Q) {
  divergence <- sum(P * log2(P / Q), na.rm = TRUE)
  return(divergence)
}

# Function to calculate JS divergence
js_divergence <- function(P, Q) {
  M <- 0.5 * (P + Q)
  jsd <- 0.5 * kl_divergence(P, M) + 0.5 * kl_divergence(Q, M)
  return(jsd)
}

# Function to calculate JS entropy
calculate_js_entropy <- function(probabilities) {
  uniform_prob <- rep(1 / length(probabilities), length(probabilities))
  jsd <- js_divergence(probabilities, uniform_prob)
  js_entropy <- sqrt(jsd)
  return(js_entropy)
}

# Function to calculate normalized entropy
calculate_normalized_entropy <- function(entropy, emb_dim) {
  max_entropy <- log2(factorial(emb_dim))
  normalized_entropy <- entropy / max_entropy
  return(normalized_entropy)
}

# Function to plot the frequencies of the ordinal patterns (Histogram)
plot_frequencies <- function(freq_df, title_txt) {
  ggplot(freq_df, aes(x = factor(Pattern, levels = unique(freq_df$Pattern)), y = Probability)) +
    geom_bar(stat = "identity", fill = "gray30", color = "black", width = 0.7) +
    xlab("Pattern") +
    ylab(expression(P(π))) +
    ggtitle(title_txt) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ylim(0, 0.3)
}

# Function to compute ordinal patterns from a given embedding dimension
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

# Load the data
CAC <- read_excel("Load/CAC.xlsx")
series <- CAC[[2]]   # the second column contains the series
timestamp <- CAC[[1]]

# Define the embedding dimensions
embedding_dimensions <- 3:7

# Initialize a data frame to store the complexity measures
complexity_results <- data.frame()

# Loop through the embedding dimensions
for (emb_dim in embedding_dimensions) {
  # Calculate ordinal pattern frequencies and probabilities
  pattern_frequencies <- ordinal_pattern_frequencies(series, emb_dim)
  probabilities_df <- calculate_probabilities(pattern_frequencies)
  probabilities <- probabilities_df$Probability
  
  # Calculate Shannon entropy, Jensen-Shannon entropy, and normalized entropy
  shannon_entropy <- calculate_permutation_entropy(probabilities)
  js_entropy <- calculate_js_entropy(probabilities)
  normalized_entropy <- calculate_normalized_entropy(shannon_entropy, emb_dim)
  
  # Print the results
  print(paste("Embedding Dimension:", emb_dim))
  print(paste("Permutation Entropy:", shannon_entropy))
  print(paste("Jensen Shannon Entropy:", js_entropy))
  print(paste("Normalized Entropy:", normalized_entropy))
  
  # Plot frequencies for the embedding dimension
  print(plot_frequencies(probabilities_df, paste("Embedding Dimension", emb_dim)))
  
  # Append the results to the complexity_results data frame
  complexity_results <- rbind(complexity_results, data.frame(
    Emb_Dim = emb_dim,
    ShannonEntropy = shannon_entropy,
    JSEntropy = js_entropy,
    NormalizedEntropy = normalized_entropy
  ))
}

# Plot the complexity plane
ggplot(complexity_results, aes(x = ShannonEntropy, y = JSEntropy, color = factor(Emb_Dim))) +
  geom_point(size = 3) +
  xlab("Shannon Entropy") +
  ylab("Jensen-Shannon Entropy") +
  ggtitle("Complexity Plane for Embedding Dimensions 3 to 7") +
  scale_color_discrete(name = "Embedding Dimension") +
  theme_minimal()

# Plot normalized entropy for each embedding dimension
ggplot(complexity_results, aes(x = Emb_Dim, y = NormalizedEntropy, color = factor(Emb_Dim))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  xlab("Embedding Dimension") +
  ylab("Normalized Entropy") +
  ggtitle("Normalized Entropy for Embedding Dimensions 3 to 7") +
  scale_color_discrete(name = "Embedding Dimension") +
  theme_minimal()
