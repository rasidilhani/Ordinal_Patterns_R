library(readxl)      # Read Excel files
library(dplyr)       # Data manipulation
library(ggplot2)     # Create plots
library(reshape2)    # Reshape data (melt function)
library(ggcorrplot)  # Correlation plots
library(writexl)     # Write Excel files
library(here)          # Handle file paths)
library(cluster)      # For silhouette analysis
library(factoextra)    # For PCA visualization
library(GGally)

# ==============================================================================
# 1. SET PATHS
# ==============================================================================

df <- here("Data", "Hybrid model_data", "D4_Data", "Combined_All_Models_HC_Results_bySampleSize_D4.xlsx")
output_dir <- here("Hybrid_Analysis", "Clustering_Results", "Clustering_Plots", "n5000")

# Create output folder
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. LOAD AND PREPARE DATA
# ==============================================================================

# Read data
hc_data_all <- read_excel(df, sheet = 2)

# Select only the columns we need (Class + 9 features)
hc_data <- hc_data_all %>%
  select(Class, H_Shannon, C_Shannon, H_Renyi, C_Renyi, 
         H_Tsallis, C_Tsallis, H_Fisher, C_Fisher, Disequilibrium)

# Convert Class to factor (categorical variable)
hc_data$Class <- as.factor(hc_data$Class)

class_colors <- c(
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


class_shapes <- c(
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
# 3. CALCULATE SUMMARY STATISTICS
# ==============================================================================

# Calculate mean and SD for each feature by class
summary_stats <- hc_data %>%
  group_by(Class) %>%
  summarise(
    n = n(),
    H_Shannon_mean = mean(H_Shannon), H_Shannon_sd = sd(H_Shannon),
    C_Shannon_mean = mean(C_Shannon), C_Shannon_sd = sd(C_Shannon),
    H_Renyi_mean = mean(H_Renyi), H_Renyi_sd = sd(H_Renyi),
    C_Renyi_mean = mean(C_Renyi), C_Renyi_sd = sd(C_Renyi),
    H_Tsallis_mean = mean(H_Tsallis), H_Tsallis_sd = sd(H_Tsallis),
    C_Tsallis_mean = mean(C_Tsallis), C_Tsallis_sd = sd(C_Tsallis),
    H_Fisher_mean = mean(H_Fisher), H_Fisher_sd = sd(H_Fisher),
    C_Fisher_mean = mean(C_Fisher), C_Fisher_sd = sd(C_Fisher),
    Disequilibrium_mean = mean(Disequilibrium), Disequilibrium_sd = sd(Disequilibrium)
  )
print(summary_stats)
# Save summary to Excel
write_xlsx(summary_stats, file.path(output_dir, "Summary_Statistics.xlsx"))

# ==============================================================================
# 4. FIGURE 1: BOXPLOTS BY FEATUREs
# I checked the variability of the features for different class of ARMA models. Figure one shows the 
# boxplots of the features for each class. The results was quite clear. but the variable Fisher has 
# such a different scale that it makes it almost impossible to have an idea of the variability 
# of the other features.In order to check how the features vary according to the class, and since the 
# measures are all nonnegative (we know this from the summary statistics), I produce boxplots 
# in logarithmic scale. This way, we can better visualize the variability of the features, 
# especially for those with a wide range of values. (Figure 2).
# ==============================================================================

# Reshape data from wide to long format (needed for ggplot)
hc_data_long <- melt(hc_data, id.vars = "Class", 
                     variable.name = "Features", value.name = "Value")

# Create boxplot by taking into account the variability of the measures
p1 <- ggplot(hc_data_long, aes(x = Features, y = Value)) +
  geom_boxplot(aes(col=Features), notch = TRUE) +
  facet_wrap(~ Class, scales = "free_y") +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  labs(title = "Variability of the features", x = "Class", y = "Value")

print(p1)
ggsave(file.path(output_dir, "Figure_1_Boxplots_Variability.pdf"), p1, width = 12, height = 8)

# Create boxplot by taking into account the variability of the measures in logarithmic scale
p2 <- ggplot(hc_data_long, aes(x = Features, y = Value)) +
  geom_boxplot(aes(col=Features), notch = TRUE) +
  scale_y_log10() +
  facet_wrap(~ Class, scales = "free_y") +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") +
  labs(title = "Variability of the features by logarithmic scale", x = "Class", y = "Value")

print(p2)
ggsave(file.path(output_dir, "Figure_2_Boxplots_Variability by logarithmic scale.pdf"), p2, width = 12, height = 8)


# ==============================================================================
# 5. FIGURE 3: CORRELATION MATRIX
# ==============================================================================

# Calculate correlation (exclude Class column)
cor_matrix <- cor(hc_data[, -1])

# Create correlation plot
p3 <- ggcorrplot(cor_matrix, 
                 method = "circle", 
                 type = "lower", 
                 lab = TRUE, 
                 title = "Correlation Matrix") 
  #theme_bw(base_size = 12, base_family = "serif") 
  
print(p3)
ggsave(file.path(output_dir, "Figure_3_Correlation.pdf"), p3, width = 10, height = 10)

# Save correlation to Excel
write_xlsx(as.data.frame(cor_matrix), file.path(output_dir, "Correlation_Matrix.xlsx"))

# ==============================================================================
# 6. FIGURE 4: SCATTER PLOT (ALL POINTS, NO CLASS LABELS)
# ==============================================================================

# Use PCA to reduce 9 features to 2 dimensions for visualization
pca_result <- prcomp(hc_data[, -1], center = TRUE, scale. = TRUE)
pca_data <- as.data.frame(pca_result$x[, 1:2])  # Take first 2 components

# Calculate variance explained by each component
var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Create scatter plot (no class colors)
p4 <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 2.5, alpha = 0.6, color = "darkblue") +
  labs(title = "All Data Points",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw(base_size = 12, base_family = "serif") 

print(p4)
ggsave(file.path(output_dir, "Figure_4_Scatter_No_Class.pdf"), p4, width = 10, height = 8)

# ==============================================================================
# 7. FIGURE 5: DENSITY PLOTS BY CLASS
# ==============================================================================

p5 <- ggplot(hc_data_long, aes(x = Value, fill = Class)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Features, scales = "free", ncol = 3) +
  theme_bw(base_size = 12, base_family = "serif")  +
  labs(title = "Density Plots by Class", x = "Value", y = "Density")

print(p5)
ggsave(file.path(output_dir, "Figure_5_Density.pdf"), p5, width = 14, height = 10)

# ==============================================================================
# 8. K-means clustering
# STEP 1 - PREPARE DATA FOR CLUSTERING
# ==============================================================================
## Extract only features (remove Class column)
features <- hc_data[, -1]  # Remove first column (Class)

# Standardize features (important for k-means!)
# This makes all features have mean=0 and sd=1
features_scaled <- scale(features)

# ==============================================================================
# STEP 2 - FIND OPTIMAL NUMBER OF CLUSTERS
# ==============================================================================

# Test k from 2 to 10
k_values <- 2:15

set.seed(1234567890, kind = "Mersenne-Twister")

# Method 1: Within Sum of Squares (WSS) - Elbow Method
wss <- numeric(length(k_values))

for (i in 1:length(k_values)) {
  k <- k_values[i]
  kmeans_result <- kmeans(features_scaled, centers = k, nstart = 25)
  wss[i] <- kmeans_result$tot.withinss
  cat(sprintf("k=%d: WSS=%.2f\n", k, wss[i]))
}

# Method 2: Silhouette coefficients
silhouette_scores <- numeric(length(k_values))

for (i in 1:length(k_values)) {
  k <- k_values[i]
  kmeans_result <- kmeans(features_scaled, centers = k, nstart = 25)
  sil <- silhouette(kmeans_result$cluster, dist(features_scaled))
  silhouette_scores[i] <- mean(sil[, 3])
  cat(sprintf("k=%d: Silhouette=%.3f\n", k, silhouette_scores[i]))
}

# Find optimal k (highest silhouette)
optimal_k <- k_values[which.max(silhouette_scores)]
cat(sprintf("\nOptimal k = %d (highest silhouette)\n\n", optimal_k))


# ==============================================================================
# 9. FIGURE 6: ELBOW PLOT
# ==============================================================================

elbow_data <- data.frame(k = k_values, WSS = wss)

p6 <- ggplot(elbow_data, aes(x = k, y = WSS)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "blue", size = 3) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = k_values) +
  labs(title = "Elbow Method - Finding Optimal k",
       subtitle = paste("Optimal k =", optimal_k),
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares") +
  theme_bw(base_size = 12, base_family = "serif")

print(p6)
ggsave(file.path(output_dir, "Figure_6_Elbow_Plot.pdf"), p6, width = 10, height = 6)

# ==============================================================================
# 10. FIGURE 7: SILHOUETTE PLOT
# ==============================================================================

cat("### Creating Figure 7: Silhouette Scores ###\n")

silhouette_data <- data.frame(k = k_values, Silhouette = silhouette_scores)

p7 <- ggplot(silhouette_data, aes(x = k, y = Silhouette)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_point(color = "darkgreen", size = 3) +
  geom_point(data = silhouette_data[silhouette_data$k == optimal_k, ],
             aes(x = k, y = Silhouette), color = "red", size = 5) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = k_values) +
  labs(title = "Silhouette Score by k",
       subtitle = paste("Optimal k =", optimal_k),
       x = "Number of Clusters (k)",
       y = "Average Silhouette Score") +
  theme_bw(base_size = 12, base_family = "serif")

print(p7)
ggsave(file.path(output_dir, "Figure_7_Silhouette_Scores.pdf"), p7, width = 10, height = 6)

# ==============================================================================
# 11. RUN K-MEANS WITH OPTIMAL k
# ==============================================================================

#cat(sprintf("\n### STEP 3: RUN K-MEANS (k=%d) ###\n", optimal_k))

# Run k-means clustering
# nstart=25 means try 25 different random starting points
set.seed(1234567890, kind = "Mersenne-Twister")
kmeans_final <- kmeans(features_scaled, centers = optimal_k, nstart = 25)

cat("K-means complete!\n")
cat("Cluster sizes:\n")
print(table(kmeans_final$cluster))

# Add cluster assignments to data
hc_data$Cluster <- as.factor(kmeans_final$cluster)

# ==============================================================================
# 12. COMPARE CLUSTERS WITH TRUE CLASSES
# ==============================================================================

#cat("\n### STEP 4: COMPARE WITH TRUE CLASSES ###\n")

# Cross-tabulation: Cluster vs Class
comparison_table <- table(Cluster = hc_data$Cluster, Class = hc_data$Class)

cat("Cross-tabulation (Cluster vs Class):\n")
print(comparison_table)

# Save to Excel
write_xlsx(list(
  Cluster_vs_Class = as.data.frame.matrix(comparison_table),
  Cluster_Summary = data.frame(
    Cluster = 1:optimal_k,
    Size = as.numeric(table(kmeans_final$cluster)),
    Within_SS = kmeans_final$withinss
  ),
  Cluster_Assignments = hc_data[, c("Class", "Cluster")]
), file.path(output_dir, "K_Means_Results.xlsx"))

# ==============================================================================
# 13. FIGURE 8: PCA PLOT WITH CLUSTERS
# =====================================================pca_result <- prcomp(features_scaled, center = FALSE, scale. = FALSE)
pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Cluster = hc_data$Cluster,
  Class = hc_data$Class
)

var_explained <- summary(pca_result)$importance[2, 1:2] * 100

# Cluster centers in PCA space
centers_pca <- predict(pca_result, kmeans_final$centers)
centers_df <- data.frame(
  PC1 = centers_pca[, 1],
  PC2 = centers_pca[, 2],
  Cluster = factor(1:optimal_k)
)

p8 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.6) +
  # Add cluster centers as asterisks
  geom_point(data = centers_df, aes(x = PC1, y = PC2), 
             color = "black", size = 3, shape = 8, stroke = 2.5,
             inherit.aes = FALSE) +
  # Add labels for centers
  geom_text(data = centers_df, aes(x = PC1, y = PC2, label = Cluster),
            color = "white", size = 4, fontface = "bold",
            inherit.aes = FALSE) +
  scale_color_brewer(palette = "Set1") +
  labs(title = paste("K-Means Clustering (k =", optimal_k, ")"),
       subtitle = "Black asterisks (*) = cluster centers",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.position = "right")

print(p8)
ggsave(file.path(output_dir, "Figure_8_Clusters_PCA.pdf"), p8, width = 10, height = 8)

# ==============================================================================
# 14. FIGURE 9: PCA PLOT WITH TRUE CLASSES
# ==============================================================================

#cat("### Creating Figure 9: True Classes ###\n")

p9 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = class_colors) +
  labs(title = "True Classes (PCA Projection)",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.position = "right")

print(p9)
ggsave(file.path(output_dir, "Figure_9_True_Classes_PCA.pdf"), p9, width = 12, height = 8)

# ==============================================================================
# 15. FIGURE 10: CLUSTER SIZES
# ==============================================================================

#cat("### Creating Figure 10: Cluster Sizes ###\n")

cluster_sizes <- data.frame(
  Cluster = factor(1:optimal_k),
  Size = as.numeric(table(kmeans_final$cluster))
)

p10 <- ggplot(cluster_sizes, aes(x = Cluster, y = Size, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Size), vjust = -0.5, size = 5) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Cluster Sizes",
       x = "Cluster",
       y = "Number of Points") +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.position = "none")

print(p10)
ggsave(file.path(output_dir, "Figure_10_Cluster_Sizes.pdf"), p10, width = 10, height = 6)

# ==============================================================================
# 16. FIGURE 11: HEATMAP OF CLUSTER CENTERS
# ==============================================================================

#cat("### Creating Figure 11: Cluster Centers Heatmap ###\n")

# Get cluster centers
centers_original <- kmeans_final$centers
rownames(centers_original) <- paste("Cluster", 1:optimal_k)
colnames(centers_original) <- colnames(features)

par(family = "serif") 
# Display heatmap on screen
heatmap(centers_original,
        scale = "column",
        main = "Cluster Centers (Standardized Features)",
        col = colorRampPalette(c("#313695", "#4575B4", "#ABD9E9", 
                                 "#FFFFBF", "#FDAE61", "#F46D43", "#A50026"))(100),
        margins = c(10, 10),
        cexRow = 1.2,
        cexCol = 1.0)

# Create heatmap
pdf(file.path(output_dir, "Figure_11_Cluster_Centers_Heatmap.pdf"), width = 10, height = 6)

dev.off()

#------------------------------------------------------------------
# cross Cluster validity: Compare k-means clusters with true classes
#--------------------------------------------------------------------
# ==============================================================================
# 17. CROSS-TABULATE CLUSTERS VS TRUE CLASSES
# ==============================================================================

# Cross-tabulation: rows=clusters, columns=true classes
cluster_table <- table(Cluster = hc_data$Cluster, Class = hc_data$Class)
print("Cluster composition by model class:")
print(cluster_table)

# Show majority class per cluster
majority_class <- apply(cluster_table, 1, function(row) {
  colnames(cluster_table)[which.max(row)]
})
cat("\nMajority class per cluster:\n")
print(majority_class)

# Cluster purity (% majority class)
purity <- apply(cluster_table, 1, function(row) {
  max(row) / sum(row) * 100
})
cat("\nCluster purity (% majority class):\n")
print(round(purity, 1))

#-----------------------------------------------------
#-----------------------------------------------------------
#18. Shannon scatter plot
# ==============================================================================
# NEW FIGURE: SCATTER PLOT - RAW DATA (H_Shannon vs C_Shannon) WITH COLORS & SHAPES
# ==============================================================================

p0 <- ggplot(hc_data, aes(x = H_Shannon, y = C_Shannon, color = Class, shape = Class)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = class_colors) +
  scale_shape_manual(values = class_shapes) +
  labs(title = "Shannon Entropy vs Complexity",
       x = expression(italic(H)), 
       y = expression(italic(C))) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

print(p0)
ggsave(file.path(output_dir, "Figure0_HC_Scatter_Colors_Shapesn5000.pdf"), 
       p0, width = 12, height = 10)

#----------------------------------------------------------
#19. K-means clustering with k=15 

## Extract only features (remove Class column)
features <- hc_data[, -1]  # Remove first column (Class)

# Standardize features (important for k-means!)
# This makes all features have mean=0 and sd=1
features_scaled <- scale(features)

# Test k from 2 to 15
k_values <- 2:15

set.seed(1234567890, kind = "Mersenne-Twister")
kmeans_result <- kmeans(features_scaled, centers = k_values, nstart = 25)

sil <- silhouette(kmeans_result$cluster, dist(features_scaled))
plot(sil)
mean(sil[, 3])  # average silhouette width
sil_width <- c()

for (k in 2:15) {
  km <- kmeans(features_scaled, centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(features_scaled))
  sil_width[k] <- mean(ss[, 3])
}

plot(2:15, sil_width[2:15], type = "b")



