library(ggplot2)
library(gridExtra)
library(grid)
library(StatOrdPattHxC)
library(ggthemes)
library(ggExtra)
library(writexl)
library(dplyr)

data("LinfLsup")
set.seed(123456789)
D <- 3
N <- c(500, 1000)
R <- 100

ar2_models <- list(
  AR2_M1 = list(ar = c(0.1, 0.8), type = "AR2_M1"),
  AR2_M2 = list(ar = c(-0.8, 0.1), type = "AR2_M2"),
  AR2_M3 = list(ar = c(0.1, -0.8), type = "AR2_M3"),
  AR2_M4 = list(ar = c(-0.8, -0.1), type = "AR2_M4")
)

generate_ar2_data <- function(ar_coef, model_name, n, r_iterations = R) {
  Output <- NULL
  for(r in 1:r_iterations){
    ts_data <- arima.sim(model = list(ar = ar_coef), n = n)
    ProbTS <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbTS)
    Complexity <- StatComplexity(ProbTS)
    Output <- rbind(Output, c(Entropy, Complexity, n, model_name))
  }
  Output <- data.frame(Output, stringsAsFactors = FALSE)
  names(Output) <- c("Entropy", "Complexity", "n", "Model")
  Output$Entropy <- as.numeric(Output$Entropy)
  Output$Complexity <- as.numeric(Output$Complexity)
  Output$n <- as.factor(Output$n)
  return(Output)
}

get_axis_range <- function(x, pad_lo = 0.03, pad_hi = 0.03, x_axis = FALSE) {
  rng <- range(x, na.rm = TRUE)
  padding_lo <- pad_lo * diff(rng)
  padding_hi <- pad_hi * diff(rng)
  if (padding_lo == 0) padding_lo <- pad_lo * abs(rng[1])
  if (padding_hi == 0) padding_hi <- pad_hi * abs(rng[2])
  if (x_axis) {
    return(c(rng[1] - padding_lo, 1))
  } else {
    return(c(0, rng[2] + padding_hi))
  }
}

create_ar2_plot <- function(output_data, title) {
  x_rng <- get_axis_range(output_data$Entropy, pad_lo = 0.02, pad_hi = 0.02, x_axis = TRUE)
  y_rng <- get_axis_range(output_data$Complexity, pad_lo = 0.05, pad_hi = 0.05, x_axis = FALSE)
  p <- ggplot() +
    geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
              aes(x = H, y = C), col = "lightgray", size = 1) +
    geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
              aes(x = H, y = C), col = "lightgray", size = 1) +
    geom_point(data = output_data,
               aes(x = Entropy, y = Complexity, col = n), alpha = 0.7, size = 1.5) +
    xlab(expression(italic(H))) +
    ylab(expression(italic(C))) +
    ggtitle(title) +
    theme_pander() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    guides(color = guide_legend(title = "Sample Size")) +
    xlim(x_rng) +
    ylim(y_rng)
  p_marginal <- ggMarginal(p, groupColour = TRUE, groupFill = FALSE)
  return(p_marginal)
}

all_ar2_data <- data.frame()
individual_plots <- list()
for(i in 1:length(ar2_models)) {
  model_name <- names(ar2_models)[i]
  ar_coef <- ar2_models[[i]]$ar
  model_data <- data.frame()
  for(n in N) {
    temp_data <- generate_ar2_data(ar_coef, model_name, n)
    model_data <- rbind(model_data, temp_data)
  }
  individual_plots[[i]] <- create_ar2_plot(
    model_data,
    paste(model_name, "(AR =", paste(ar_coef, collapse = ", "), ")")
  )
  all_ar2_data <- rbind(all_ar2_data, model_data)
}

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

get_plot_bounds <- function(df) {
  x_rng <- get_axis_range(df$Entropy, pad_lo = 0.02, pad_hi = 0.02, x_axis = TRUE)
  y_rng <- get_axis_range(df$Complexity, pad_lo = 0.05, pad_hi = 0.05, x_axis = FALSE)
  list(x = x_rng, y = y_rng)
}
dummy_bounds <- get_plot_bounds(all_ar2_data)

dummy_plot <- ggplot() +
  geom_point(data = all_ar2_data,
             aes(x = Entropy, y = Complexity, col = n), alpha = 0.7, size = 1.5) +
  theme_pander() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Sample Size")) +
  xlim(dummy_bounds$x) +
  ylim(dummy_bounds$y)
shared_legend <- get_only_legend(dummy_plot)

combined_plot <- arrangeGrob(
  grobs = individual_plots, ncol = 2, nrow = 2,
  top = "AR(2) Models - Ordinal Pattern Analysis"
)
grid.arrange(combined_plot, shared_legend, nrow = 2, heights = c(10, 1))

ggsave(
  "AR2_combined_analysis.pdf",
  arrangeGrob(combined_plot, shared_legend, nrow = 2, heights = c(10, 1)),
  width = 12, height = 10, dpi = 300
)

################################################################################################
library(ggplot2)
library(gridExtra)
library(grid)
library(StatOrdPattHxC)
library(ggthemes)
library(ggExtra)
library(writexl)
library(dplyr)

# Load necessary data and set parameters
data("LinfLsup")
set.seed(1234567890, kind="Mersenne-Twister")

D <- 3          # Embedding dimension
N <- 1000       # Sample size
R <- 100        # Number of repetitions

# Define AR(2) model configurations
ar2_models <- list(
  AR2_M1 = list(ar = c(0.1, 0.8), type = "AR2_M1"),
  AR2_M2 = list(ar = c(-0.8, 0.1), type = "AR2_M2"),
  AR2_M3 = list(ar = c(0.1, -0.8), type = "AR2_M3"),
  AR2_M4 = list(ar = c(-0.8, -0.1), type = "AR2_M4")
)

# Function to generate AR(2) data and compute entropy & complexity
generate_ar2_data <- function(ar_coef, model_name, n, r_iterations = R) {
  Output <- NULL
  for (r in 1:r_iterations) {
    ts_data <- arima.sim(model = list(ar = ar_coef), n = n)
    ProbTS <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbTS)
    Complexity <- StatComplexity(ProbTS)
    Output <- rbind(Output, c(Entropy, Complexity, n, model_name))
  }
  Output <- data.frame(Output, stringsAsFactors = FALSE)
  names(Output) <- c("Entropy", "Complexity", "n", "Model")
  Output$EntropyH <- as.numeric(Output$Entropy)
  Output$Complexity <- as.numeric(Output$Complexity)
  Output$n <- as.factor(Output$n)
  return(Output)
}

# Generate results for all AR(2) models (N = 1000 only)
all_ar2_data <- data.frame()

for (i in 1:length(ar2_models)) {
  model_name <- names(ar2_models)[i]
  ar_coef <- ar2_models[[i]]$ar
  temp_data <- generate_ar2_data(ar_coef, model_name, N)
  all_ar2_data <- rbind(all_ar2_data, temp_data)
}

# Save final combined results to CSV
write.csv(all_ar2_data, "AR2_Entropy_Complexity_N1000.csv", row.names = FALSE)

cat("✅ Entropy–Complexity calculations for N = 1000 completed and saved as 'AR2_Entropy_Complexity_N1000.csv'\n")

                                                                                                             
##########################################################################################################
# ----------------------------------------------------------
# Ordinal Pattern Entropy–Complexity Framework: AR(2)
# ----------------------------------------------------------
library(StatOrdPattHxC)
library(dplyr)

# Load data and parameters
data("LinfLsup")
set.seed(1234567890, kind = "Mersenne-Twister")

D <- 3          # Embedding dimension (fixed)
N <- 1000       # Sample size
R <- 100        # Repetitions

# Define AR(2) models
ar2_models <- list(
  AR2_M1 = list(ar = c(0.1, 0.8), type = "AR2_M1"),
  AR2_M2 = list(ar = c(-0.8, 0.1), type = "AR2_M2"),
  AR2_M3 = list(ar = c(0.1, -0.8), type = "AR2_M3"),
  AR2_M4 = list(ar = c(-0.8, -0.1), type = "AR2_M4")
)

# ----------------------------------------------------------
# Helper: Jensen–Shannon Divergence
# ----------------------------------------------------------
Jensen_Shannon <- function(p, q) {
  m <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12)/(m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12)/(m + 1e-12))))
  return(js)
}

# ----------------------------------------------------------
# Extended Complexities for Rényi / Tsallis / Fisher
# ----------------------------------------------------------

# Jensen–Shannon type complexity for any entropy measure
GeneralizedComplexity <- function(prob, entropy_value) {
  Pe <- rep(1 / length(prob), length(prob))  # Uniform reference
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * entropy_value)
}

# Fisher–Information–based complexity (Fisher × Jensen–Shannon)
FisherBasedComplexity <- function(prob, fisher_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * fisher_value)
}

# ----------------------------------------------------------
# Main Loop Function
# ----------------------------------------------------------
generate_ar2_data <- function(ar_coef, model_name, n, r_iterations = R) {
  Output <- NULL
  for (r in 1:r_iterations) {
    ts_data <- arima.sim(model = list(ar = ar_coef), n = n)
    
    # Probability distribution from ordinal patterns
    ProbTS <- OPprob(ts_data, emb = D)
    
    # Entropy measures (normalized by default in StatOrdPattHxC)
    Hs <- HShannon(ProbTS)
    Hr <- HRenyi(ProbTS, beta = 1.5)
    Ht <- HTsallis(ProbTS, beta =  1.5)
    Hf <- HFisher(ProbTS)
    
    # Jensen–Shannon Statistical Complexity (Shannon baseline)
    C_Shannon <- StatComplexity(ProbTS)
    
    # Jensen–Shannon based complexities for Rényi and Tsallis entropies
    C_Renyi <- GeneralizedComplexity(ProbTS, Hr)
    C_Tsallis <- GeneralizedComplexity(ProbTS, Ht)
    
    # Fisher-based complexity
    C_Fisher <- FisherBasedComplexity(ProbTS, Hf)
    
    # Combine results
    Output <- rbind(Output, c(Hs, Hr, Ht, Hf, 
                              C_Shannon, C_Renyi, C_Tsallis, C_Fisher,
                              n, model_name))
  }
  
  # Assemble results
  Output <- as.data.frame(Output, stringsAsFactors = FALSE)
  names(Output) <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher",
                     "C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher",
                     "n", "Model")
  Output[, 1:8] <- lapply(Output[, 1:8], as.numeric)
  Output$n <- as.factor(Output$n)
  return(Output)
}

# ----------------------------------------------------------
# Run across all AR(2) models
# ----------------------------------------------------------
all_ar2_data <- data.frame()
for (i in seq_along(ar2_models)) {
  model_name <- names(ar2_models)[i]
  ar_coef <- ar2_models[[i]]$ar
  temp_data <- generate_ar2_data(ar_coef, model_name, N)
  all_ar2_data <- rbind(all_ar2_data, temp_data)
}

# ----------------------------------------------------------
# Save results
# ----------------------------------------------------------
write.csv(all_ar2_data, "AR2_All_Entropy_Complexity_Final.csv", row.names = FALSE)

cat("✅ Normalized Shannon, Rényi, Tsallis, Fisher entropies and respective complexities saved as 'AR2_All_Entropy_Complexity_Final.csv'\n")

##################################################################################################################

# ----------------------------------------------------------
# Entropy–Complexity Plane Visualization with Central Point
# ----------------------------------------------------------
library(ggplot2)
library(dplyr)
library(MASS)      # For kernel density estimation (kde2d)

# ----------------------------------------------------------
# 1. Load the CSV data
# ----------------------------------------------------------
data <- read.csv("AR2_All_Entropy_Complexity_Final.csv")

# Convert Model to factor for grouping
data$Model <- as.factor(data$Model)

# ----------------------------------------------------------
# Helper: Get axis ranges with padding
# ----------------------------------------------------------
p <- ggplot(data, aes(x = H_Shannon, y = C_Shannon, color = Model)) +
  geom_density_2d_filled(alpha = 0.3, show.legend = FALSE) +
  geom_point(alpha = 0.4, size = 1.5) +
  #geom_point(data = central_points, aes(x = H_star, y = C_star, color = Model),
  #           shape = 7, size = 5, stroke = 1.5) +  # star marker for (H, C)
  labs(title = "Entropy–Complexity Plane (H × C) for AR(2) Models",
       x = "Shannon Entropy (H)",
       y = "Statistical Complexity (C)") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))

print(p)

# Save plot
ggsave("HxC_Plane_AR2_Points.pdf", plot = p, width = 10, height = 7)




# ----------------------------------------------------------
# 2. Compute Central Point (H*, C*) for each model
# ----------------------------------------------------------
# Using MEDIAN (Chagas et al., 2022 approach)
central_points <- data %>%
  group_by(Model) %>%
  summarise(
    H_star = median(H_Shannon),
    C_star = median(C_Shannon),
    .groups = "drop"
  )

print("Central Points (H*, C*) for Each Model:")
print(central_points)

# ----------------------------------------------------------
# 3. H×C Plane: Density Contour Plot with Central Points
# ----------------------------------------------------------
p1 <- ggplot(data, aes(x = H_Shannon, y = C_Shannon, color = Model)) +
  geom_density_2d_filled(alpha = 0.4, show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = central_points, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5, stroke = 1.5) +  # Triangle marker for (H*, C*)
  labs(title = "Entropy–Complexity Plane (H × C) for AR(2) Models",
       x = "Shannon Entropy (H)",
       y = "Statistical Complexity (C)") +
  theme_minimal() +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14))

print(p1)

# Save plot
ggsave("HxC_Plane_AR2_with_Central_Points.pdf", plot = p1, width = 10, height = 7)

# ----------------------------------------------------------
# 4. Faceted H×C Plane: Separate Contour Plots per Model
# ----------------------------------------------------------
p2 <- ggplot(data, aes(x = H_Shannon, y = C_Shannon)) +
  geom_density_2d_filled(alpha = 0.6) +
  geom_point(alpha = 0.4, size = 1, color = "black") +
  geom_point(data = central_points, aes(x = H_star, y = C_star),
             color = "red", shape = 17, size = 4) +
  facet_wrap(~ Model, scales = "free") +
  labs(title = "H × C Plane by Model (with Central Points)",
       x = "Entropy (H)",
       y = "Complexity (C)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

print(p2)

# Save faceted plot
ggsave("HxC_Plane_AR2_Faceted.pdf", plot = p2, width = 12, height = 8)

# ----------------------------------------------------------
# 5. Alternative: Hexbin Density Plot
# ----------------------------------------------------------
library(hexbin)

p3 <- ggplot(data, aes(x = H_Shannon, y = C_Shannon)) +
  geom_hex(bins = 30, alpha = 0.7) +
  geom_point(data = central_points, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5, stroke = 1.5) +
  scale_fill_viridis_c() +
  facet_wrap(~ Model) +
  labs(title = "Hexbin Density: H × C Plane with Central Points",
       x = "Entropy (H)",
       y = "Complexity (C)") +
  theme_minimal()

print(p3)

# Save hexbin plot
ggsave("HxC_Plane_AR2_Hexbin.pdf", plot = p3, width = 12, height = 8)

# ----------------------------------------------------------
# 6. Overlay of All Models (Scatterplot with Centroids)
# ----------------------------------------------------------
p4 <- ggplot(data, aes(x = H_Shannon, y = C_Shannon, color = Model)) +
  geom_point(alpha = 0.3, size = 1.5) +
  stat_ellipse(level = 0.95, linetype = "dashed", size = 1) +  # 95% confidence ellipse
  geom_point(data = central_points, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "H × C Plane: All AR(2) Models with Central Points",
       x = "Normalized Shannon Entropy (H)",
       y = "Statistical Complexity (C)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p4)

# Save overlay plot
ggsave("HxC_Plane_AR2_Overlay.pdf", plot = p4, width = 10, height = 7)

# ----------------------------------------------------------
# 7. Summary Table of Central Points
# ----------------------------------------------------------
write.csv(central_points, "Central_Points_HxC.csv", row.names = FALSE)

cat("✅ All H×C plane plots saved:\n")
cat("   - HxC_Plane_AR2_with_Central_Points.pdf\n")
cat("   - HxC_Plane_AR2_Faceted.pdf\n")
cat("   - HxC_Plane_AR2_Hexbin.pdf\n")
cat("   - HxC_Plane_AR2_Overlay.pdf\n")
cat("✅ Central points (H*, C*) saved as 'Central_Points_HxC.csv'\n")

##########################################################################################################
# ----------------------------------------------------------
# Multi-Entropy Complexity Plane Visualization
# Includes: Shannon, Rényi, Tsallis, and Fisher
# ----------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining multiple plots

# ----------------------------------------------------------
# 1. Load CSV Data
# ----------------------------------------------------------
data <- read.csv("AR2_All_Entropy_Complexity_Final.csv")
data$Model <- as.factor(data$Model)

# ----------------------------------------------------------
# 2. Compute Central Points (H*, C*) for Each Entropy Type
# ----------------------------------------------------------

# Shannon
central_shannon <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Shannon),
            C_star = median(C_Shannon),
            Entropy_Type = "Shannon",
            .groups = "drop")

# Rényi
central_renyi <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Renyi),
            C_star = median(C_Renyi),
            Entropy_Type = "Renyi",
            .groups = "drop")

# Tsallis
central_tsallis <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Tsallis),
            C_star = median(C_Tsallis),
            Entropy_Type = "Tsallis",
            .groups = "drop")

# Fisher (optional - uses Fisher information as "entropy proxy")
central_fisher <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Fisher),
            C_star = median(C_Fisher),
            Entropy_Type = "Fisher",
            .groups = "drop")

# Combine all central points
central_all <- bind_rows(central_shannon, central_renyi, central_tsallis, central_fisher)

# Save central points
write.csv(central_all, "Central_Points_All_Entropies.csv", row.names = FALSE)

cat("✅ Central points for all entropy types saved as 'Central_Points_All_Entropies.csv'\n")

# ----------------------------------------------------------
# 3. Prepare Long-Format Data for Faceted Plotting
# ----------------------------------------------------------
data_long <- data %>%
  pivot_longer(cols = c(H_Shannon, H_Renyi, H_Tsallis, H_Fisher),
               names_to = "Entropy_Type", values_to = "Entropy") %>%
  pivot_longer(cols = c(C_Shannon, C_Renyi, C_Tsallis, C_Fisher),
               names_to = "Complexity_Type", values_to = "Complexity") %>%
  filter(
    (Entropy_Type == "H_Shannon" & Complexity_Type == "C_Shannon") |
      (Entropy_Type == "H_Renyi" & Complexity_Type == "C_Renyi") |
      (Entropy_Type == "H_Tsallis" & Complexity_Type == "C_Tsallis") |
      (Entropy_Type == "H_Fisher" & Complexity_Type == "C_Fisher")
  ) %>%
  mutate(Entropy_Type = gsub("H_", "", Entropy_Type))

# ----------------------------------------------------------
# 4. Multi-Panel H×C Plot: All Entropy Types
# ----------------------------------------------------------
p_multi <- ggplot(data_long, aes(x = Entropy, y = Complexity, color = Model)) +
  geom_point(alpha = 0.4, size = 1) +
  stat_ellipse(level = 0.95, linetype = "dashed", size = 0.8) +
  geom_point(data = central_all, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 4, stroke = 1.2) +
  facet_wrap(~ Entropy_Type, scales = "free", ncol = 2) +
  labs(title = "Entropy–Complexity Planes: Shannon, Rényi, Tsallis, Fisher",
       x = "Entropy (H)",
       y = "Complexity (C)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(strip.text = element_text(face = "bold", size = 12),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))

print(p_multi)

# Save multi-panel plot
ggsave("HxC_Multi_Entropy_Faceted.pdf", plot = p_multi, width = 12, height = 10)

# ----------------------------------------------------------
# 5. Individual Plots for Each Entropy Type
# ----------------------------------------------------------

# Shannon
p_shannon <- ggplot(data, aes(x = H_Shannon, y = C_Shannon, color = Model)) +
  geom_density_2d_filled(alpha = 0.4, show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = central_shannon, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Shannon Entropy–Complexity Plane", x = "H (Shannon)", y = "C (Shannon)") +
  theme_minimal(base_size = 12, base_family = "serif")

# Rényi
p_renyi <- ggplot(data, aes(x = H_Renyi, y = C_Renyi, color = Model)) +
  geom_density_2d_filled(alpha = 0.4, show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = central_renyi, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Rényi Entropy–Complexity Plane", x = "H (Rényi)", y = "C (Rényi)") +
  theme_minimal(base_size = 12, base_family = "serif")

# Tsallis
p_tsallis <- ggplot(data, aes(x = H_Tsallis, y = C_Tsallis, color = Model)) +
  geom_density_2d_filled(alpha = 0.4, show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = central_tsallis, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Tsallis Entropy–Complexity Plane", x = "H (Tsallis)", y = "C (Tsallis)") +
  theme_minimal(base_size = 12, base_family = "serif")

# Fisher
p_fisher <- ggplot(data, aes(x = H_Fisher, y = C_Fisher, color = Model)) +
  geom_density_2d_filled(alpha = 0.4, show.legend = FALSE) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = central_fisher, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Fisher Information–Complexity Plane", x = "Fisher Info", y = "C (Fisher)") +
  theme_minimal(base_size = 12, base_family = "serif")

# Combine using patchwork
combined_plot <- (p_shannon | p_renyi) / (p_tsallis | p_fisher) +
  plot_annotation(title = "Comparison of Entropy–Complexity Planes",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))

print(combined_plot)

# Save combined plot
ggsave("HxC_All_Entropies_Combined.pdf", plot = combined_plot, width = 14, height = 12)

# ----------------------------------------------------------
# 6. Overlay Plot: Compare Central Points Across Entropy Types
# ----------------------------------------------------------
p_overlay <- ggplot(central_all, aes(x = H_star, y = C_star, color = Model, shape = Entropy_Type)) +
  geom_point(size = 5, stroke = 1.5) +
  labs(title = "Central Points (H*, C*) Across Entropy Types",
       x = "Entropy (H*)",
       y = "Complexity (C*)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "right")

print(p_overlay)

# Save overlay of central points
ggsave("Central_Points_Overlay_All_Entropies.pdf", plot = p_overlay, width = 10, height = 7)

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------
cat("✅ All entropy–complexity plane plots saved:\n")
cat("   - HxC_Multi_Entropy_Faceted.pdf (faceted by entropy type)\n")
cat("   - HxC_All_Entropies_Combined.pdf (2x2 grid layout)\n")
cat("   - Central_Points_Overlay_All_Entropies.pdf (comparison of central points)\n")
cat("✅ Central points saved as 'Central_Points_All_Entropies.csv'\n")
##########################################################################################################

#PCA, Contour plot#####################################################################
# ----------------------------------------------------------
# Multi-Entropy H×C Plane with Contour Plots + PCA Analysis
# ----------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggfortify)  # For autoplot PCA

# ----------------------------------------------------------
# 1. Load CSV Data
# ----------------------------------------------------------
data <- read.csv("AR2_All_Entropy_Complexity_Final.csv")
data$Model <- as.factor(data$Model)

# ----------------------------------------------------------
# 2. Compute Central Points for Each Entropy Type
# ----------------------------------------------------------
central_shannon <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Shannon), C_star = median(C_Shannon),
            Entropy_Type = "Shannon", .groups = "drop")

central_renyi <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Renyi), C_star = median(C_Renyi),
            Entropy_Type = "Renyi", .groups = "drop")

central_tsallis <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Tsallis), C_star = median(C_Tsallis),
            Entropy_Type = "Tsallis", .groups = "drop")

central_fisher <- data %>%
  group_by(Model) %>%
  summarise(H_star = median(H_Fisher), C_star = median(C_Fisher),
            Entropy_Type = "Fisher", .groups = "drop")

central_all <- bind_rows(central_shannon, central_renyi, central_tsallis, central_fisher)

# ----------------------------------------------------------
# 3. Contour Plots (Filled Density) for Each Entropy Type
# ----------------------------------------------------------

# Shannon with filled contours
p_shannon_contour <- ggplot(data, aes(x = H_Shannon, y = C_Shannon)) +
  geom_density_2d_filled(alpha = 0.7) +
  geom_point(aes(color = Model), alpha = 0.5, size = 1.5) +
  geom_point(data = central_shannon, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Shannon H×C (Contour)", x = "H (Shannon)", y = "C (Shannon)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "none")

# Rényi with filled contours
p_renyi_contour <- ggplot(data, aes(x = H_Renyi, y = C_Renyi)) +
  geom_density_2d_filled(alpha = 0.7) +
  geom_point(aes(color = Model), alpha = 0.5, size = 1.5) +
  geom_point(data = central_renyi, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Rényi H×C (Contour)", x = "H (Rényi)", y = "C (Rényi)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "none")

# Tsallis with filled contours
p_tsallis_contour <- ggplot(data, aes(x = H_Tsallis, y = C_Tsallis)) +
  geom_density_2d_filled(alpha = 0.7) +
  geom_point(aes(color = Model), alpha = 0.5, size = 1.5) +
  geom_point(data = central_tsallis, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Tsallis H×C (Contour)", x = "H (Tsallis)", y = "C (Tsallis)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "none")

# Fisher with filled contours
p_fisher_contour <- ggplot(data, aes(x = H_Fisher, y = C_Fisher)) +
  geom_density_2d_filled(alpha = 0.7) +
  geom_point(aes(color = Model), alpha = 0.5, size = 1.5) +
  geom_point(data = central_fisher, aes(x = H_star, y = C_star, color = Model),
             shape = 17, size = 5) +
  labs(title = "Fisher H×C (Contour)", x = "Fisher Info", y = "C (Fisher)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "right")

# Combine contour plots
contour_combined <- (p_shannon_contour | p_renyi_contour) / 
  (p_tsallis_contour | p_fisher_contour) +
  plot_annotation(title = "Contour Density: H×C Planes with Central Points",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))

print(contour_combined)

# Save combined contour plot
ggsave("HxC_Contour_All_Entropies.pdf", plot = contour_combined, width = 14, height = 12)

# ----------------------------------------------------------
# 4. PCA Analysis on All Entropy-Complexity Features
# ----------------------------------------------------------

# Select feature columns for PCA
features <- dplyr::select(data, H_Shannon, H_Renyi, H_Tsallis, H_Fisher,
                          C_Shannon, C_Renyi, C_Tsallis, C_Fisher)


# Perform PCA (scale features for comparability)
pca_result <- prcomp(features, scale. = TRUE, center = TRUE)

# Summary of PCA
cat("\n========== PCA Summary ==========\n")
print(summary(pca_result))

# Variance explained by each principal component
variance_explained <- summary(pca_result)$importance[2, ]
cat("\nVariance Explained by Each PC:\n")
print(variance_explained)

# ----------------------------------------------------------
# 5. PCA Biplot: PC1 vs PC2
# ----------------------------------------------------------
p_pca <- autoplot(pca_result, data = data, colour = 'Model', 
                  loadings = TRUE, loadings.label = TRUE,
                  loadings.label.size = 4, loadings.colour = "black",
                  size = 2, alpha = 0.6) +
  labs(title = "PCA: Entropy-Complexity Feature Space",
       subtitle = paste0("PC1: ", round(variance_explained[1]*100, 1), "% | ",
                         "PC2: ", round(variance_explained[2]*100, 1), "%")) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

print(p_pca)

# Save PCA plot
ggsave("PCA_Entropy_Complexity_Features.pdf", plot = p_pca, width = 10, height = 8)

# ----------------------------------------------------------
# 6. PCA Loadings Table (Contribution of Each Feature)
# ----------------------------------------------------------
loadings <- as.data.frame(pca_result$rotation[, 1:4])  # First 4 PCs
loadings$Feature <- rownames(loadings)
loadings <- loadings[, c("Feature", "PC1", "PC2", "PC3", "PC4")]

cat("\n========== PCA Loadings (First 4 PCs) ==========\n")
print(loadings)

# Save loadings to CSV
write.csv(loadings, "PCA_Loadings.csv", row.names = FALSE)

# ----------------------------------------------------------
# 7. Scree Plot: Variance Explained by Each PC
# ----------------------------------------------------------
scree_data <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance = variance_explained * 100
)

p_scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "skyblue", alpha = 0.8) +
  geom_line(aes(group = 1), color = "red", size = 1) +
  geom_point(color = "red", size = 3) +
  labs(title = "Scree Plot: Variance Explained by Principal Components",
       x = "Principal Component",
       y = "Variance Explained (%)") +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_scree)

# Save scree plot
ggsave("PCA_Scree_Plot.pdf", plot = p_scree, width = 8, height = 6)

# ----------------------------------------------------------
# 8. PCA: PC1 vs PC3 (Alternative View)
# ----------------------------------------------------------
p_pca_13 <- autoplot(pca_result, data = data, colour = 'Model', 
                     x = 1, y = 3, size = 2, alpha = 0.6) +
  labs(title = "PCA: PC1 vs PC3") +
  theme_minimal(base_size = 12, base_family = "serif")

print(p_pca_13)

ggsave("PCA_PC1_vs_PC3.pdf", plot = p_pca_13, width = 10, height = 8)

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------
cat("\n✅ All visualizations and analyses saved:\n")
cat("   - HxC_Contour_All_Entropies.pdf (filled contour plots)\n")
cat("   - PCA_Entropy_Complexity_Features.pdf (biplot with loadings)\n")
cat("   - PCA_Scree_Plot.pdf (variance explained)\n")
cat("   - PCA_PC1_vs_PC3.pdf (alternative PCA view)\n")
cat("   - PCA_Loadings.csv (feature contributions to PCs)\n")
##########################################################################################################
# median point in the PCA#############################################################
##
# ----------------------------------------------------------
# Compute PCA as before
# ----------------------------------------------------------
features <- data %>%
  dplyr::select(H_Shannon, H_Renyi, H_Tsallis, H_Fisher,
                C_Shannon, C_Renyi, C_Tsallis, C_Fisher)

pca_result <- prcomp(features, scale. = TRUE, center = TRUE)
scores <- as.data.frame(pca_result$x)  # PCA scores for each observation
scores$Model <- data$Model  # Add model info

# ----------------------------------------------------------
# Compute Median Point per Model (PC1*, PC2*)
# ----------------------------------------------------------
pca_medians <- scores %>%
  group_by(Model) %>%
  summarise(PC1_star = median(PC1),
            PC2_star = median(PC2),
            .groups = 'drop')

print(pca_medians)

# ----------------------------------------------------------
# Enhanced PCA Plot with Medians
# ----------------------------------------------------------
library(ggplot2)

# Base PCA scatter plot
pca_plot <- ggplot(scores, aes(x = PC1, y = PC2, color = Model)) +
  geom_point(alpha = 0.6, size = 2) +
  stat_ellipse(level = 0.95, linetype = "dashed", size = 0.8) +
  geom_point(data = pca_medians, aes(x = PC1_star, y = PC2_star, color = Model),
             shape = 17, size = 5, stroke = 1.3) +  # triangle marker for median
  geom_text(data = pca_medians, aes(x = PC1_star, y = PC2_star, label = Model),
            vjust = -1.2, hjust = 0.5, size = 4, fontface = "bold") +
  labs(title = "PCA of Entropy–Complexity Features",
       subtitle = "Including cluster medians (PC1*, PC2*) per AR(2) model",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14))

print(pca_plot)

# Save plot
ggsave("PCA_Entropy_Complexity_with_MedianPoints.pdf",
       plot = pca_plot, width = 10, height = 8)

# ----------------------------------------------------------
# Save Median Data
# ----------------------------------------------------------
write.csv(pca_medians, "PCA_Central_Points.csv", row.names = FALSE)
cat("✅ PCA plot with median central points saved as 'PCA_Entropy_Complexity_with_MedianPoints.pdf'\n")
cat("✅ PCA medians table saved as 'PCA_Central_Points.csv'\n")

##########################################################################################################
# ----------------------------------------------------------
# Median Point and Stat Ellipse Plot for Each Entropy Type
# ----------------------------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)

# ----------------------------------------------------------
# 1. Load Data
# ----------------------------------------------------------
data <- read.csv("AR2_All_Entropy_Complexity_Final.csv")
data$Model <- as.factor(data$Model)

# ----------------------------------------------------------
# 2. Compute Median Points (H*, C*) per Model for all entropy types
# ----------------------------------------------------------
summarise_by_model <- function(df, H, C, entropy_type) {
  df %>%
    group_by(Model) %>%
    summarise(
      H_star = median({{ H }}, na.rm = TRUE),
      C_star = median({{ C }}, na.rm = TRUE),
      Entropy_Type = entropy_type,
      .groups = "drop"
    )
}

# Apply function to each entropy type
median_shannon  <- summarise_by_model(data, H_Shannon, C_Shannon, "Shannon")
median_renyi    <- summarise_by_model(data, H_Renyi, C_Renyi, "Renyi")
median_tsallis  <- summarise_by_model(data, H_Tsallis, C_Tsallis, "Tsallis")
median_fisher   <- summarise_by_model(data, H_Fisher, C_Fisher, "Fisher")

# Combine all medians
medians_all <- bind_rows(median_shannon, median_renyi, median_tsallis, median_fisher)

# Check results
print(medians_all)


# ----------------------------------------------------------
# 3. Function for H×C plane with Medians + Ellipses
# ----------------------------------------------------------
plot_entropy_complexity <- function(df, H_col, C_col, centroids, title, H_label, C_label) {
  ggplot(df, aes(x = {{ H_col }}, y = {{ C_col }}, color = Model)) +
    geom_point(alpha = 0.5, size = 1.5) +
    stat_ellipse(level = 0.95, linewidth = 0.8, linetype = "dashed") +
    geom_point(data = centroids, aes(x = H_star, y = C_star, color = Model),
               shape = 17, size = 4, stroke = 1.2) +
    geom_text(data = centroids, aes(x = H_star, y = C_star, label = Model),
              vjust = -1, fontface = "bold", size = 4) +
    labs(title = title, x = H_label, y = C_label) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "bottom")
}

# ----------------------------------------------------------
# 4. Plots for Each Entropy Type
# ----------------------------------------------------------
p_shannon <- plot_entropy_complexity(
  df = data,
  H_col = H_Shannon,
  C_col = C_Shannon,
  centroids = median_shannon,
  title = "Shannon Entropy–Complexity Plane",
  H_label = "H (Shannon)",
  C_label = "C (Shannon)"
)

p_renyi <- plot_entropy_complexity(
  df = data,
  H_col = H_Renyi,
  C_col = C_Renyi,
  centroids = median_renyi,
  title = "Rényi Entropy–Complexity Plane",
  H_label = "H (Rényi)",
  C_label = "C (Rényi)"
)

p_tsallis <- plot_entropy_complexity(
  df = data,
  H_col = H_Tsallis,
  C_col = C_Tsallis,
  centroids = median_tsallis,
  title = "Tsallis Entropy–Complexity Plane",
  H_label = "H (Tsallis)",
  C_label = "C (Tsallis)"
)

p_fisher <- plot_entropy_complexity(
  df = data,
  H_col = H_Fisher,
  C_col = C_Fisher,
  centroids = median_fisher,
  title = "Fisher Information–Complexity Plane",
  H_label = "H (Fisher)",
  C_label = "C (Fisher)"
)

# ----------------------------------------------------------
# 5. Combine All into One Visualization
# ----------------------------------------------------------
combined_plot <- (p_shannon | p_renyi) / (p_tsallis | p_fisher) +
  plot_annotation(title = "Entropy–Complexity Planes with Median Points and Stat Ellipses",
                  theme = theme(plot.title = element_text(face = "bold", size = 16)))

print(combined_plot)

# Save plot to file
ggsave("Entropy_Complexity_Medians_StatEllipses.pdf", plot = combined_plot,
       width = 14, height = 12)

# ----------------------------------------------------------
# 6. Export Median Points
# ----------------------------------------------------------
write.csv(medians_all, "Entropy_Complexity_Medians.csv", row.names = FALSE)
cat("✅ Saved combined entropy–complexity plot with median points and ellipses.\n")
cat("✅ Median coordinates saved as 'Entropy_Complexity_Medians.csv'\n")
##########################################################################################################

# EDA for Entropy–Complexity Data
# ===============================

library(dplyr)
library(ggplot2)
library(GGally)
library(cluster)
library(factoextra)
library(randomForest)
library(Rtsne)
library(data.table)

# ----------- 1. Load and Prepare Data ---------------
data <- fread("AR2_All_Entropy_Complexity_Final.csv")
data$Model <- as.factor(data$Model)

features <- data %>%
  dplyr::select(H_Shannon, H_Renyi, H_Tsallis, H_Fisher, 
                C_Shannon, C_Renyi, C_Tsallis, C_Fisher)

# ----------- 2. Correlation and Pairwise EDA --------
# Correlation matrix
corr_matrix <- cor(features, use = "complete.obs")
print("Correlation Matrix:")
print(round(corr_matrix, 2))

# Scatterplots and density for all features
ggpairs(features, aes(color = data$Model, alpha = 0.5)) +
  ggtitle("Pairwise Feature Relationships (Colored by Model)")
ggsave("pairs_entropy_complexity.pdf")

# ----------- 3. Clustering --------------------------
# Distance matrix and hierarchical clustering
dist_matrix <- dist(scale(features))
hc <- hclust(dist_matrix, method = "ward.D2")
plot(hc, hang=-1, labels=data$Model, main="Hierarchical Clustering on Features")

# Cut dendrogram into k=4 clusters (matching M1–M4)
clusters <- cutree(hc, k=4)
table(clusters, data$Model)

# Silhouette analysis
sil <- silhouette(clusters, dist_matrix)
fviz_silhouette(sil)
ggsave("silhouette_clusters.pdf")

# ----------- 4. Feature Importance (Random Forest) ---
rf <- randomForest(Model ~ ., 
                   data=data.frame(Model=data$Model, features), 
                   importance=TRUE, ntree=500)
print(rf)
varImpPlot(rf)
ggsave("varImpPlot_rf.pdf")

# ----------- 5. t-SNE Visualization -----------------
set.seed(2024)
tsne_out <- Rtsne(as.matrix(scale(features)), dims=2, perplexity=15, verbose=TRUE)
tsne_df <- data.frame(
  X = tsne_out$Y[,1],
  Y = tsne_out$Y[,2],
  Model = data$Model
)
ggplot(tsne_df, aes(x=X, y=Y, color=Model)) +
  geom_point(alpha = 0.6, size=2) +
  theme_minimal() +
  labs(title="t-SNE: Entropy–Complexity Features")
ggsave("tSNE_entropy_complexity.pdf")

# ----------- 6. Boxplots by Model -------------------
long <- pivot_longer(data, cols = c(H_Shannon, H_Renyi, H_Tsallis, H_Fisher,
                                    C_Shannon, C_Renyi, C_Tsallis, C_Fisher),
                     names_to = "Feature", values_to = "Value")
ggplot(long, aes(x=Model, y=Value, fill=Model)) +
  geom_boxplot(alpha=0.7) +
  facet_wrap(~Feature, ncol=4, scales="free") +
  theme_minimal() +
  labs(title="Distributions of Entropy and Complexity Features by Model")
ggsave("boxplots_entropy_complexity.pdf")

# ----------- 7. Optional: Feature ANOVA -------------
feature_names <- colnames(features)
for (f in feature_names) {
  cat("\nANOVA for", f, "by Model\n")
  print(summary(aov(features[[f]] ~ data$Model)))
}

