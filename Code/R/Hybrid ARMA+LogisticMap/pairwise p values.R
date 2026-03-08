library(readxl)      # Read Excel files
library(dplyr)       # Data manipulation
library(ggplot2)     # Create plots
library(reshape2)    # Reshape data (melt function)
library(ggcorrplot)  # Correlation plots
library(writexl)     # Write Excel files
library(here)        # Handle file paths
library(cluster)     # For silhouette analysis
library(factoextra)  # For PCA visualization
library(GGally)
library(dendextend)  # For dendrogram coloring
library(pvclust)     # For bootstrap p-values on clusters

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
  "ARMA(2,2)" = "#D55E00",
  "ARMA(1,1)" = "#0072B2",
  "AR(2)" = "#009E73",
  "AR(1)" = "#CC79A7",
  "MA(2)" = "#E69F00",
  "MA(1)" = "#56B4E9",
  "Logistic" = "#000000",
  "Hybrid_ARMA(2,2)" = "#8B4513",
  "Hybrid_ARMA(1,1)" = "#4B0082",
  "Hybrid_AR(2)" = "#696969",
  "Hybrid_AR(1)" = "#20B2AA",
  "Hybrid_MA(2)" = "#B22222",
  "Hybrid_MA(1)" = "#4169E1",
  "Logistic_r3_6" = "#2F4F4F",
  "Sine_Wave" = "#228B22",
  "Logistic_Sine_Combined" = "#8A2BE2"
)

# ==============================================================================
# 3. PREPARE DATA FOR CLUSTERING
# ==============================================================================

# Extract features (remove Class column)
features <- hc_data[, -1]

# Standardize features
features_scaled <- scale(features)

# Set seed for reproducibility
set.seed(123)

# ==============================================================================
# 4. PAIRWISE MULTIVARIATE TESTS BETWEEN CLASSES
# ==============================================================================

cat("\n### PAIRWISE COMPARISONS BETWEEN CLASSES ###\n")

# Get unique classes
classes <- levels(hc_data$Class)
n_classes <- length(classes)

# Initialize matrix for p-values
pairwise_pvalues <- matrix(NA, nrow = n_classes, ncol = n_classes)
rownames(pairwise_pvalues) <- classes
colnames(pairwise_pvalues) <- classes
diag(pairwise_pvalues) <- 1  # Diagonal is 1 (same class)

# Perform pairwise MANOVA (or Hotelling's T2)
library(ICSNP)  # For Hotelling's T2 test

cat("\nCalculating pairwise p-values (Hotelling's T2 test)...\n")

for (i in 1:(n_classes - 1)) {
  for (j in (i + 1):n_classes) {
    
    class1 <- classes[i]
    class2 <- classes[j]
    
    # Extract data for these two classes
    data1 <- features_scaled[hc_data$Class == class1, ]
    data2 <- features_scaled[hc_data$Class == class2, ]
    
    # Hotelling's T2 test
    test_result <- tryCatch({
      HotellingsT2(data1, data2)
    }, error = function(e) {
      list(p.value = NA)
    })
    
    # Store p-value
    pairwise_pvalues[i, j] <- test_result$p.value
    pairwise_pvalues[j, i] <- test_result$p.value
    
    cat(sprintf("%s vs %s: p = %.4f %s\n", 
                class1, class2, test_result$p.value,
                ifelse(test_result$p.value < 0.001, "***",
                       ifelse(test_result$p.value < 0.01, "**",
                              ifelse(test_result$p.value < 0.05, "*", "")))))
  }
}

# Save p-value matrix to Excel
write_xlsx(list(
  Pairwise_P_Values = as.data.frame(pairwise_pvalues)
), file.path(output_dir, "Pairwise_Class_PValues.xlsx"))

cat("\n✅ Pairwise p-values saved to Excel\n")

# ==============================================================================
# 5. VISUALIZE P-VALUE MATRIX
# ==============================================================================

cat("\n### Creating p-value heatmap ###\n")

# Convert to long format for ggplot
pvalue_long <- melt(pairwise_pvalues)
names(pvalue_long) <- c("Class1", "Class2", "P_Value")

# Create heatmap
p_pvalue <- ggplot(pvalue_long, aes(x = Class1, y = Class2, fill = -log10(P_Value))) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(P_Value < 0.001, "***",
                               ifelse(P_Value < 0.01, "**",
                                      ifelse(P_Value < 0.05, "*", "")))),
            size = 3, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 1.3, # -log10(0.05) = 1.3
                       name = "-log10(p)") +
  labs(title = "Pairwise Class Differences (Hotelling's T² Test)",
       subtitle = "* p<0.05, ** p<0.01, *** p<0.001",
       x = "Class", y = "Class") +
  theme_bw(base_size = 10, base_family = "serif") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_pvalue)
ggsave(file.path(output_dir, "Pairwise_PValue_Heatmap.pdf"), 
       p_pvalue, width = 12, height = 10)

# ==============================================================================
# 6. HIERARCHICAL CLUSTERING & DENDROGRAM
# ==============================================================================

cat("\n### HIERARCHICAL CLUSTERING ###\n")

# Calculate distance matrix
dist_matrix <- dist(features_scaled, method = "euclidean")

# Perform hierarchical clustering (different methods)
hc_complete <- hclust(dist_matrix, method = "complete")
hc_average <- hclust(dist_matrix, method = "average")
hc_ward <- hclust(dist_matrix, method = "ward.D2")

# Convert to dendrogram objects
dend_complete <- as.dendrogram(hc_complete)
dend_average <- as.dendrogram(hc_average)
dend_ward <- as.dendrogram(hc_ward)

# Color branches by true class
library(dendextend)

# Assign colors based on class
class_cols <- class_colors[as.character(hc_data$Class)]

# Method 1: Simple dendrogram with class colors
pdf(file.path(output_dir, "Dendrogram_Complete_Linkage.pdf"), 
    width = 14, height = 8)
par(family = "serif")
plot(hc_complete, 
     main = "Hierarchical Clustering Dendrogram (Complete Linkage)",
     xlab = "Observations",
     ylab = "Distance",
     sub = "",
     cex = 0.5)
# Color labels by class
labels_colors(dend_complete) <- class_cols[order.dendrogram(dend_complete)]
plot(dend_complete, main = "Dendrogram Colored by Class")
dev.off()

# Method 2: Ward's method (better clustering)
pdf(file.path(output_dir, "Dendrogram_Ward_Linkage.pdf"), 
    width = 14, height = 8)
par(family = "serif")
plot(hc_ward, 
     main = "Hierarchical Clustering Dendrogram (Ward's Method)",
     xlab = "Observations",
     ylab = "Distance",
     sub = "",
     cex = 0.5,
     hang = -1)
dev.off()

cat("✅ Dendrograms saved\n")

# ==============================================================================
# 7. DENDROGRAM WITH COLORED BRANCHES BY K-MEANS CLUSTERS
# ==============================================================================

# First, run k-means to get cluster assignments
cat("\n### Running K-means for dendrogram coloring ###\n")

# Find optimal k using silhouette
k_values <- 2:10
silhouette_scores <- numeric(length(k_values))

for (i in 1:length(k_values)) {
  k <- k_values[i]
  km_temp <- kmeans(features_scaled, centers = k, nstart = 100)
  sil <- silhouette(km_temp$cluster, dist_matrix)
  silhouette_scores[i] <- mean(sil[, 3])
}

optimal_k <- k_values[which.max(silhouette_scores)]
cat(sprintf("Optimal k = %d\n", optimal_k))

# Run k-means with optimal k
kmeans_final <- kmeans(features_scaled, centers = optimal_k, nstart = 100)

# Color dendrogram by k-means clusters
dend_colored <- color_branches(dend_ward, k = optimal_k)

pdf(file.path(output_dir, "Dendrogram_Ward_Colored_by_KMeans.pdf"), 
    width = 14, height = 8)
par(family = "serif")
plot(dend_colored, 
     main = paste0("Hierarchical Dendrogram (Ward) - Colored by ", 
                   optimal_k, " K-Means Clusters"),
     ylab = "Height")
dev.off()

cat("✅ K-means colored dendrogram saved\n")

# ==============================================================================
# 8. DENDROGRAM WITH CLASS LABELS
# ==============================================================================

# Create dendrogram with class labels instead of observation numbers
labels(dend_ward) <- as.character(hc_data$Class)[order.dendrogram(dend_ward)]
labels_colors(dend_ward) <- class_colors[hc_data$Class][order.dendrogram(dend_ward)]

pdf(file.path(output_dir, "Dendrogram_Ward_with_Class_Labels.pdf"), 
    width = 16, height = 10)
par(family = "serif", mar = c(10, 4, 4, 2))
plot(dend_ward, 
     main = "Hierarchical Clustering - Class Labels",
     ylab = "Height",
     cex = 0.4)
dev.off()

cat("✅ Class-labeled dendrogram saved\n")

# ==============================================================================
# 9. PVCLUST: BOOTSTRAP P-VALUES FOR HIERARCHICAL CLUSTERING
# ==============================================================================

cat("\n### Bootstrap P-values for Clusters (pvclust) ###\n")
cat("This may take a few minutes...\n")

# pvclust requires data as columns = variables, rows = observations
# We need to transpose
features_t <- t(features_scaled)

# Run pvclust (100 bootstrap replicates)
pvc <- pvclust(features_t, method.hclust = "ward.D2", 
               method.dist = "euclidean", nboot = 100)

# Plot with p-values
pdf(file.path(output_dir, "Dendrogram_with_Bootstrap_PValues.pdf"), 
    width = 14, height = 8)
par(family = "serif")
plot(pvc, 
     main = "Hierarchical Clustering with Bootstrap P-values",
     print.pv = "au",  # Show AU (Approximately Unbiased) p-values
     print.num = FALSE)
pvrect(pvc, alpha = 0.95)  # Draw rectangles around significant clusters
dev.off()

cat("✅ Bootstrap p-value dendrogram saved\n")

# Save pvclust results
pvclust_summary <- data.frame(
  Cluster_Height = pvc$hclust$height,
  AU_PValue = pvc$edges[, "au"],  # Approximately Unbiased p-value
  BP_PValue = pvc$edges[, "bp"]   # Bootstrap Probability
)

write_xlsx(list(
  PVClust_Results = pvclust_summary
), file.path(output_dir, "PVClust_Bootstrap_Results.xlsx"))

# ==============================================================================
# 10. COPHENETIC CORRELATION (Dendrogram Quality)
# ==============================================================================

cat("\n### Dendrogram Quality Metrics ###\n")

# Calculate cophenetic correlation for different methods
cophenetic_complete <- cor(dist_matrix, cophenetic(hc_complete))
cophenetic_average <- cor(dist_matrix, cophenetic(hc_average))
cophenetic_ward <- cor(dist_matrix, cophenetic(hc_ward))

cat(sprintf("Cophenetic correlation (Complete): %.3f\n", cophenetic_complete))
cat(sprintf("Cophenetic correlation (Average):  %.3f\n", cophenetic_average))
cat(sprintf("Cophenetic correlation (Ward):     %.3f\n", cophenetic_ward))
cat("\nHigher is better (>0.75 is good)\n")

# Save quality metrics
quality_metrics <- data.frame(
  Method = c("Complete", "Average", "Ward"),
  Cophenetic_Correlation = c(cophenetic_complete, cophenetic_average, cophenetic_ward)
)

write_xlsx(list(
  Dendrogram_Quality = quality_metrics
), file.path(output_dir, "Dendrogram_Quality_Metrics.xlsx"))

cat("\n✅ ALL ANALYSES COMPLETE! ✅\n")
cat(sprintf("Results saved to: %s\n", output_dir))