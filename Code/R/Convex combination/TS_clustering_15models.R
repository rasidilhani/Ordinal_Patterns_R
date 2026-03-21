#---------------------------------------------------------------

# ══════════════════════════════════════════════════════════════════════════════
#  Time Series Clustering using Ordinal-Pattern Features
#  Includes: Summary stats, Boxplots, Correlation, PCA, K-means
# ══════════════════════════════════════════════════════════════════════════════

library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(writexl)
library(here)
library(cluster)
library(factoextra)


# Global theme: serif font, size 11
theme_set(
  theme_bw(base_size = 11, base_family = "serif") +
    theme(
      plot.title = element_text(size = 11, family = "serif", face = "bold"),
      axis.title = element_text(size = 11, family = "serif"),
      axis.text  = element_text(size = 10, family = "serif"),
      legend.text = element_text(size = 9, family = "serif"),
      legend.title = element_text(size = 10, family = "serif", face = "bold")
    )
)
# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS AND PATHS
# ══════════════════════════════════════════════════════════════════════════════
D     <- 4
n_val <- 10000    

hc_data_path <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))
output_dir   <- here("Results", "Convex_combination", "Clustering",
                     paste0("n", n_val))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL COLORS AND SHAPES (for visualisation only)
# ══════════════════════════════════════════════════════════════════════════════
model_names <- c(
  "ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
  "ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)",
  "AR2+Logistic(w=0.1)", "AR2+Sine(w=0.8)",
  "MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)",
  "MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"
)

# Colors/shapes (unchanged)
model_colors <- c(
  "ARMA(2,2)" = "black",
  "AR(2)" = "tomato", "MA(2)" = "navy",
  "Logistic" = "dodgerblue", "Sine" = "forestgreen",
  setNames(colorRampPalette(c("#A8D8EA","#3B9AB2"))(3),
           c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)")),
  "AR2+Logistic(w=0.1)" = "#D4A017",
  "AR2+Sine(w=0.8)" = "#7EC8A4",
  setNames(colorRampPalette(c("#F5AAAA","#C0392B"))(2),
           c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  setNames(colorRampPalette(c("#D9B8E8","#8E44AD"))(3),
           c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"))
)

model_shapes <- c(
  "ARMA(2,2)" = 16, "AR(2)" = 17, "MA(2)" = 15,
  "Logistic" = 18, "Sine" = 8,
  setNames(rep(17,3), c("ARMA+Sine(w=0.1)","ARMA+Sine(w=0.2)","ARMA+Sine(w=0.3)")),
  "AR2+Logistic(w=0.1)" = 15, "AR2+Sine(w=0.8)" = 18,
  setNames(rep(1,2), c("MA2+Logistic(w=0.2)","MA2+Logistic(w=0.7)")),
  setNames(rep(2,3), c("MA2+Sine(w=0.4)","MA2+Sine(w=0.6)","MA2+Sine(w=0.8)"))
)

# ══════════════════════════════════════════════════════════════════════════════
# 14 FEATURES (FULL SET)
# ══════════════════════════════════════════════════════════════════════════════
features <- c(
  "H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher",
  "C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher",
  "Disequilibrium",
  "Var_H_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Var_H_Fisher", "Var_C_Shannon"
)

# ══════════════════════════════════════════════════════════════════════════════
#  LOAD & PREPARE DATA
# ══════════════════════════════════════════════════════════════════════════════
hc_data <- read_xlsx(hc_data_path) %>%
  select(Model, all_of(features)) %>%
  rename(Class = Model) %>%
  mutate(Class = factor(Class, levels = model_names))

# ══════════════════════════════════════════════════════════════════════════════
#  SUMMARY STATISTICS
# ══════════════════════════════════════════════════════════════════════════════
summary_stats <- hc_data %>%
  group_by(Class) %>%
  summarise(across(all_of(features), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

write_xlsx(summary_stats, file.path(output_dir, "Summary_Statistics.xlsx"))

# ══════════════════════════════════════════════════════════════════════════════
#  VARIABILITY CHECK — BOXPLOTS (regular + log scale)
# ══════════════════════════════════════════════════════════════════════════════

hc_long <- melt(hc_data, id.vars = "Class",
                variable.name = "Feature",
                value.name = "Value")

# --- Boxplot (normal scale)
p_box_normal <- ggplot(hc_long, aes(x = Feature, y = Value)) +
  geom_boxplot(aes(color = Feature), notch = TRUE) +
  facet_wrap(~ Class, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(title = "Feature Variability by Class",
       x = "Feature", y = "Value")
p_box_normal

ggsave(file.path(output_dir, "Boxplots_Normal.pdf"),
       p_box_normal, width = 12, height = 8)

# --- Boxplot (log10 scale)
p_box_log <- ggplot(hc_long, aes(x = Feature, y = Value)) +
  geom_boxplot(aes(color = Feature), notch = TRUE) +
  scale_y_log10() +
  facet_wrap(~ Class, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(title = "Feature Variability by Class (Log Scale)",
       x = "Feature", y = "log10(Value)")
p_box_log

ggsave(file.path(output_dir, "Boxplots_LogScale.pdf"),
       p_box_log, width = 12, height = 8)


# ══════════════════════════════════════════════════════════════════════════════
#  CORRELATION MATRIX
# ══════════════════════════════════════════════════════════════════════════════
cor_matrix <- cor(hc_data[, features], use="complete.obs")

p2 <- ggcorrplot(cor_matrix, method="circle", type="lower", lab=TRUE)
p2
ggsave(file.path(output_dir, "Correlation.pdf"), p2, width=10, height=10)
write_xlsx(as.data.frame(cor_matrix), file.path(output_dir, "Correlation.xlsx"))

# ══════════════════════════════════════════════════════════════════════════════
#  SHANNON HC PLANE (OPTIONAL VISUAL CHECK)
# ══════════════════════════════════════════════════════════════════════════════
p3 <- ggplot(hc_data, aes(H_Shannon, C_Shannon, color=Class, shape=Class)) +
  geom_point(size=2.5, alpha=0.7) +
  scale_color_manual(values=model_colors) +
  scale_shape_manual(values=model_shapes) +
  theme_bw() +
  labs(title="Shannon Entropy–Complexity Plane")
p3

ggsave(file.path(output_dir, "Shannon_HC.pdf"), p3, width=12, height=8)

# ══════════════════════════════════════════════════════════════════════════════
#  K-MEANS CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════
features_scaled <- scale(hc_data[, features])

k_values <- 2:15
wss_scores <- sil_scores <- numeric(length(k_values))

set.seed(123)

for (i in seq_along(k_values)) {
  k <- k_values[i]
  km <- kmeans(features_scaled, centers=k, nstart=25)
  wss_scores[i] <- km$tot.withinss
  sil_scores[i] <- mean(silhouette(km$cluster, dist(features_scaled))[,3])
}

optimal_k <- k_values[which.max(sil_scores)]

# --- Elbow plot
p4 <- ggplot(data.frame(k=k_values, WSS=wss_scores), aes(k,WSS)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = optimal_k, linetype="dashed", color="red") +
  theme_bw() +
  labs(title="Elbow Method", subtitle=paste("Optimal k =", optimal_k))
p4
ggsave(file.path(output_dir, "Elbow.pdf"), p4, width=10, height=6)

# --- Silhouette plot
p5 <- ggplot(data.frame(k=k_values, Sil=sil_scores), aes(k,Sil)) +
  geom_line(color="darkgreen") + geom_point(color="darkgreen") +
  geom_point(data=data.frame(k=optimal_k, Sil=max(sil_scores)),
             aes(k,Sil), color="red", size=4) +
  theme_bw() +
  labs(title="Silhouette Scores", subtitle=paste("Optimal k =", optimal_k))
p5
ggsave(file.path(output_dir, "Silhouette.pdf"), p5, width=10, height=6)

# ══════════════════════════════════════════════════════════════════════════════
#  FINAL K-MEANS
# ══════════════════════════════════════════════════════════════════════════════
km_final <- kmeans(features_scaled, centers=optimal_k, nstart=25)
hc_data$Cluster <- as.factor(km_final$cluster)

# ══════════════════════════════════════════════════════════════════════════════
#  CLUSTER SIZE BAR PLOT
# ══════════════════════════════════════════════════════════════════════════════

cluster_sizes <- as.data.frame(table(hc_data$Cluster))
colnames(cluster_sizes) <- c("Cluster", "Size")

n_clusters <- length(unique(hc_data$Cluster))

my_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_clusters)

p_cluster_size <- ggplot(cluster_sizes, aes(x = Cluster, y = Size, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = my_colors) +
  geom_text(aes(label = Size), vjust = -0.5, size = 4) +
  labs(title = "Cluster Sizes", x = "Cluster", y = "Count") +
  theme_bw()

p_cluster_size

ggsave(file.path(output_dir, "Cluster_Sizes.pdf"),
       p_cluster_size, width = 10, height = 6)


# ══════════════════════════════════════════════════════════════════════════════
#  CLUSTER HEATMAP (Feature Means)
# ══════════════════════════════════════════════════════════════════════════════

cluster_feature_means <- hc_data %>%
  group_by(Cluster) %>%
  summarise(across(all_of(features), mean))

# reshape to long format
heatmap_long <- cluster_feature_means %>%
  pivot_longer(cols = all_of(features),
               names_to = "Feature",
               values_to = "MeanValue")

p_heatmap <- ggplot(heatmap_long, aes(x = Feature, y = Cluster, fill = MeanValue)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C") +
  labs(title = "Cluster Heatmap of Feature Means",
       x = "Feature", y = "Cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_heatmap

ggsave(file.path(output_dir, "Cluster_Heatmap.pdf"),
       p_heatmap, width = 14, height = 8)


# ══════════════════════════════════════════════════════════════════════════════
#  CLUSTER HEATMAP (Z-SCALED FEATURE MEANS)
# ══════════════════════════════════════════════════════════════════════════════

library(tidyr)
library(dplyr)

# 1. Compute cluster means (raw)
cluster_feature_means <- hc_data %>%
  group_by(Cluster) %>%
  summarise(across(all_of(features), mean), .groups = "drop")

# 2. Z-scale each feature so heatmap is interpretable
cluster_feature_means_scaled <- cluster_feature_means %>%
  mutate(across(all_of(features),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

# 3. Convert to long format
heatmap_long <- cluster_feature_means_scaled %>%
  select(Cluster, ends_with("_z")) %>%
  pivot_longer(cols = -Cluster,
               names_to = "Feature",
               values_to = "ZValue")

# Clean feature names
heatmap_long$Feature <- gsub("_z", "", heatmap_long$Feature)

# 4. Plot heatmap
p_heatmap <- ggplot(heatmap_long,
                    aes(x = Feature, y = Cluster, fill = ZValue)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "C") +
  labs(title = "Cluster Heatmap of Z-Scaled Feature Means",
       x = "Feature",
       y = "Cluster",
       fill = "Z-Score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_heatmap

ggsave(file.path(output_dir, "Cluster_Heatmap_Zscaled.pdf"),
       p_heatmap, width = 14, height = 8)

# ══════════════════════════════════════════════════════════════════════════════
#  FINAL SILHOUETTE ANALYSIS (for chosen k)
# ══════════════════════════════════════════════════════════════════════════════

sil_final <- silhouette(km_final$cluster, dist(features_scaled))

# Silhouette plot
p_sil_final <- fviz_silhouette(sil_final) +
  ggtitle(paste("Silhouette Plot for k =", optimal_k))
p_sil_final

ggsave(file.path(output_dir, "Silhouette_Final_k.pdf"),
       p_sil_final, width = 10, height = 7)

# Save average silhouette width
avg_sil <- mean(sil_final[, 3])
writeLines(paste("Average silhouette width:", round(avg_sil, 4)),
           file.path(output_dir, "Silhouette_Final_k.txt"))

# ══════════════════════════════════════════════════════════════════════════════
#  PCA VISUALISATION
# ══════════════════════════════════════════════════════════════════════════════
pca_result <- prcomp(features_scaled)
pca_df <- data.frame(PC1=pca_result$x[,1],
                     PC2=pca_result$x[,2],
                     Cluster=hc_data$Cluster,
                     Class=hc_data$Class)

# PCA - clusters
p6 <- ggplot(pca_df, aes(PC1,PC2,color=Cluster)) +
  geom_point(alpha=0.6) +
  theme_bw() +
  labs(title=paste("K-means clustering (k =", optimal_k, ")"))
p6
ggsave(file.path(output_dir, "PCA_Clusters.pdf"), p6, width=10, height=8)

# PCA - true classes
p7 <- ggplot(pca_df, aes(PC1,PC2,color=Class)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=model_colors) +
  theme_bw() +
  labs(title="True Classes (PCA projection)")
p7
ggsave(file.path(output_dir, "PCA_Classes.pdf"), p7, width=12, height=8)

# ══════════════════════════════════════════════════════════════════════════════
#  STEP-BY-STEP CLUSTER CENTER MOVEMENT (PCA snapshots)
# ══════════════════════════════════════════════════════════════════════════════

set.seed(123)

# Step 0: initial centers (randomly chosen)
km_init <- kmeans(features_scaled, centers = optimal_k, nstart = 1, iter.max = 1)

# Step 1: after one iteration
km_step1 <- kmeans(features_scaled, centers = km_init$centers, nstart = 1, iter.max = 1)

# Step Final: your final result km_final already computed

# Project all three sets of centers into PCA space
cent_init  <- as.data.frame(predict(pca_result, km_init$centers))
cent_step1 <- as.data.frame(predict(pca_result, km_step1$centers))
cent_final <- as.data.frame(predict(pca_result, km_final$centers))

cent_init$Step  <- "Initial"
cent_step1$Step <- "Iteration 1"
cent_final$Step <- "Final"

cent_all <- rbind(cent_init, cent_step1, cent_final)
cent_all$Cluster <- rep(1:optimal_k, 3)

pA <- ggplot() +
  geom_point(data = pca_df, aes(PC1, PC2), color = "grey80", alpha = 0.4) +
  geom_point(data = cent_init, aes(PC1, PC2), color = "blue", size = 4, shape = 8) +
  geom_point(data = cent_final, aes(PC1, PC2), color = "red", size = 4, shape = 17) +
  labs(title = "Cluster Centres: Initial vs Final",
       subtitle = "Blue = Initial, Red = Final") +
  theme_bw()
pA
ggsave(file.path(output_dir, "PCA_Centres_Initial_Final.pdf"),
       pA, width = 10, height = 8)

pB <- ggplot() +
  geom_point(data = pca_df, aes(PC1, PC2), color = "grey80", alpha = 0.4) +
  geom_point(data = cent_all, aes(PC1, PC2, color = Step, shape = Step), size = 4) +
  labs(title = "Movement of Cluster Centres in PCA Space",
       subtitle = "Initial → Iteration 1 → Final") +
  scale_color_manual(values = c("Initial"="blue", "Iteration 1"="orange", "Final"="red")) +
  theme_bw()
pB

ggsave(file.path(output_dir, "PCA_Centre_Movement.pdf"),
       pB, width = 10, height = 8)


library(tidyr)
library(ggrepel)

# ----- reshape centres wider (initial, iteration1, final) -----
cent_wide <- cent_all %>%
  pivot_wider(
    names_from = Step,
    values_from = c(PC1, PC2)
  )

# ----- PCA Movement Plot including cluster points -----
pC <- ggplot() +
  
  # All points colored by their cluster
  geom_point(data = pca_df,
             aes(PC1, PC2, color = Cluster),
             alpha = 0.35, size = 1.8) +
  
  # Arrow from initial → final for each cluster
  geom_segment(
    data = cent_wide,
    aes(x = PC1_Initial, y = PC2_Initial,
        xend = PC1_Final, yend = PC2_Final,
        color = factor(Cluster)),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 1
  ) +
  
  # Highlight initial centres
  geom_point(
    data = cent_wide,
    aes(PC1_Initial, PC2_Initial, color = factor(Cluster)),
    size = 4, shape = 1, stroke = 1.2
  ) +
  
  # Highlight final centres
  geom_point(
    data = cent_wide,
    aes(PC1_Final, PC2_Final, color = factor(Cluster)),
    size = 4, shape = 16
  ) +
  
  # Label final centres
  geom_label_repel(
    data = cent_wide,
    aes(PC1_Final, PC2_Final, label = Cluster, fill = factor(Cluster)),
    size = 4, color = "white", show.legend = FALSE
  ) +
  
  labs(
    title = "Movement of Cluster Centres in PCA Space",
    subtitle = "Cluster points (faded), initial centres (circle), final centres (solid), arrows show movement",
    color = "Cluster"
  ) +
  theme_bw() +
  theme(legend.position = "right")
pC

# Save figure
ggsave(
  file.path(output_dir, "PCA_Centre_Movement_Arrows_WithPoints.pdf"),
  pC, width = 11, height = 8
)


# ══════════════════════════════════════════════════════════════════════════════
#  SAVE RESULTS
# ══════════════════════════════════════════════════════════════════════════════
write_xlsx(
  list(
    Cluster_Assignments = hc_data[, c("Class","Cluster")],
    Cluster_Centers     = as.data.frame(km_final$centers),
    Silhouette_Scores   = data.frame(k_values, sil_scores),
    WSS                 = data.frame(k_values, wss_scores)
  ),
  file.path(output_dir, "Clustering_Results.xlsx")
)

message("✓ Clustering completed. All outputs saved.")
# End of the code

#----------------------------------------------------------------------------------
#------------------------------------------------------------------------------------
# ══════════════════════════════════════════════════════════════════════════════
#  Clustering Analysis — Convex Combination Models
#  Includes: Summary stats, Boxplots, Correlation, PCA, K-Means
# ══════════════════════════════════════════════════════════════════════════════
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggcorrplot)
library(writexl)
library(here)
library(cluster)
library(factoextra)

# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS AND PATHS
# ══════════════════════════════════════════════════════════════════════════════
D     <- 4
n_val <- 10000    

hc_data_path <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))
output_dir   <- here("Results", "Convex_combination", "Clustering",
                     paste0("n", n_val))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL NAMES, COLORS, SHAPES
# ══════════════════════════════════════════════════════════════════════════════
model_names <- c(
  "ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
  "ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)",
  "AR2+Logistic(w=0.1)", "AR2+Sine(w=0.8)",
  "MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)",
  "MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"
)

model_colors <- c(
  "ARMA(2,2)" = "black",
  "AR(2)"     = "tomato",
  "MA(2)"     = "navy",
  "Logistic"  = "dodgerblue",
  "Sine"      = "forestgreen",
  setNames(colorRampPalette(c("#A8D8EA", "#3B9AB2"))(3),
           c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)")),
  "AR2+Logistic(w=0.1)" = "#D4A017",
  "AR2+Sine(w=0.8)"     = "#7EC8A4",
  setNames(colorRampPalette(c("#F5AAAA", "#C0392B"))(2),
           c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  setNames(colorRampPalette(c("#D9B8E8", "#8E44AD"))(3),
           c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"))
)

model_shapes <- c(
  "ARMA(2,2)" = 16, "AR(2)" = 17, "MA(2)" = 15,
  "Logistic"  = 18, "Sine"  = 8,
  setNames(rep(17, 3), c("ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)",
                         "ARMA+Sine(w=0.3)")),
  "AR2+Logistic(w=0.1)" = 15,
  "AR2+Sine(w=0.8)"     = 18,
  setNames(rep(1, 2), c("MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)")),
  setNames(rep(2, 3), c("MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)",
                        "MA2+Sine(w=0.8)"))
)

# ══════════════════════════════════════════════════════════════════════════════
#  FEATURES
# ══════════════════════════════════════════════════════════════════════════════
features <- c(
  "H_Shannon", "C_Shannon",
  "H_Renyi",   "C_Renyi",
  "H_Tsallis", "C_Tsallis",
  "H_Fisher",  "C_Fisher",
  "Disequilibrium"
)

# ══════════════════════════════════════════════════════════════════════════════
#  LOAD AND PREPARE DATA
# ══════════════════════════════════════════════════════════════════════════════
hc_data <- read_xlsx(hc_data_path, sheet = paste0("n", n_val)) %>%
  select(Model, all_of(features)) %>%
  rename(Class = Model) %>%
  mutate(Class = factor(Class, levels = model_names))

# ══════════════════════════════════════════════════════════════════════════════
#  SUMMARY STATISTICS
# ══════════════════════════════════════════════════════════════════════════════
summary_stats <- hc_data %>%
  group_by(Class) %>%
  summarise(
    n                  = n(),
    H_Shannon_mean     = mean(H_Shannon),     H_Shannon_sd     = sd(H_Shannon),
    C_Shannon_mean     = mean(C_Shannon),     C_Shannon_sd     = sd(C_Shannon),
    H_Renyi_mean       = mean(H_Renyi),       H_Renyi_sd       = sd(H_Renyi),
    C_Renyi_mean       = mean(C_Renyi),       C_Renyi_sd       = sd(C_Renyi),
    H_Tsallis_mean     = mean(H_Tsallis),     H_Tsallis_sd     = sd(H_Tsallis),
    C_Tsallis_mean     = mean(C_Tsallis),     C_Tsallis_sd     = sd(C_Tsallis),
    H_Fisher_mean      = mean(H_Fisher),      H_Fisher_sd      = sd(H_Fisher),
    C_Fisher_mean      = mean(C_Fisher),      C_Fisher_sd      = sd(C_Fisher),
    Disequilibrium_mean = mean(Disequilibrium), Disequilibrium_sd = sd(Disequilibrium)
  )

write_xlsx(summary_stats,
           file.path(output_dir, "Summary_Statistics.xlsx"))

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 1 — Boxplots by feature and class
# ══════════════════════════════════════════════════════════════════════════════
hc_long <- melt(hc_data, id.vars = "Class",
                variable.name = "Feature", value.name = "Value")

p1 <- ggplot(hc_long, aes(x = Feature, y = Value, color = Feature)) +
  geom_boxplot(notch = TRUE) +
  facet_wrap(~ Class, scales = "free_y") +
  labs(title = "Feature Variability by Class",
       x = "Feature", y = "Value") +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1),
        legend.position = "none")

ggsave(file.path(output_dir, "Figure_1_Boxplots.pdf"),
       p1, width = 12, height = 8, device = cairo_pdf)
print(p1)

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 2 — Boxplots (log scale)
# ══════════════════════════════════════════════════════════════════════════════
p2 <- p1 +
  scale_y_log10() +
  labs(title = "Feature Variability by Class (log scale)")

ggsave(file.path(output_dir, "Figure_2_Boxplots_Log.pdf"),
       p2, width = 12, height = 8, device = cairo_pdf)
print(p2)

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 3 — Correlation matrix
# ══════════════════════════════════════════════════════════════════════════════
cor_matrix <- cor(hc_data[, features], use = "complete.obs")

p3 <- ggcorrplot(cor_matrix,
                 method = "circle",
                 type   = "lower",
                 lab    = TRUE,
                 title  = "Correlation Matrix")

ggsave(file.path(output_dir, "Figure_3_Correlation.pdf"),
       p3, width = 10, height = 10, device = cairo_pdf)
print(p3)

write_xlsx(as.data.frame(cor_matrix),
           file.path(output_dir, "Correlation_Matrix.xlsx"))

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 4 — Shannon HC scatter plot (colored by class)
# ══════════════════════════════════════════════════════════════════════════════
p4 <- ggplot(hc_data, aes(x = H_Shannon, y = C_Shannon,
                          color = Class, shape = Class)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  labs(title = paste0("Shannon HC Plane  (n = ", n_val, ")"),
       x     = expression(italic(H)[Shannon]),
       y     = expression(italic(C)[Shannon])) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.text  = element_text(size = 8),
        legend.title = element_text(size = 9, face = "bold"))

ggsave(file.path(output_dir, "Figure_4_Shannon_HC.pdf"),
       p4, width = 12, height = 8, device = cairo_pdf)
print(p4)

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 5 — Density plots by class
# ══════════════════════════════════════════════════════════════════════════════
p5 <- ggplot(hc_long, aes(x = Value, fill = Class)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ Feature, scales = "free", ncol = 3) +
  labs(title = "Density Plots by Class", x = "Value", y = "Density") +
  theme_bw(base_size = 11, base_family = "serif")

ggsave(file.path(output_dir, "Figure_5_Density.pdf"),
       p5, width = 14, height = 10, device = cairo_pdf)
print(p5)

# ══════════════════════════════════════════════════════════════════════════════
#  K-MEANS CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════

# Scale features (mean = 0, sd = 1) — required before k-means
features_scaled <- scale(hc_data[, features])

set.seed(1234567890, kind = "Mersenne-Twister")
k_values <- 2:15

# Compute WSS (elbow method) and silhouette score for each k
wss_scores <- numeric(length(k_values))
sil_scores <- numeric(length(k_values))

for (i in seq_along(k_values)) {
  k  <- k_values[i]
  km <- kmeans(features_scaled, centers = k, nstart = 25)
  wss_scores[i] <- km$tot.withinss
  sil_scores[i] <- mean(silhouette(km$cluster, dist(features_scaled))[, 3])
}

# Optimal k = highest silhouette score
optimal_k <- k_values[which.max(sil_scores)]

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 6 — Elbow plot (WSS)
# ══════════════════════════════════════════════════════════════════════════════
elbow_df <- data.frame(k = k_values, WSS = wss_scores)

p6 <- ggplot(elbow_df, aes(x = k, y = WSS)) +
  geom_line(color = "blue", linewidth = 1.2) +
  geom_point(color = "blue", size = 3) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = k_values) +
  labs(title    = "Elbow Method",
       subtitle = paste("Optimal k =", optimal_k),
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares") +
  theme_bw(base_size = 12, base_family = "serif")

ggsave(file.path(output_dir, "Figure_6_Elbow.pdf"),
       p6, width = 10, height = 6, device = cairo_pdf)
print(p6)

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 7 — Silhouette scores
# ══════════════════════════════════════════════════════════════════════════════
sil_df <- data.frame(k = k_values, Silhouette = sil_scores)

p7 <- ggplot(sil_df, aes(x = k, y = Silhouette)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_point(color = "darkgreen", size = 3) +
  geom_point(data = filter(sil_df, k == optimal_k),
             aes(x = k, y = Silhouette),
             color = "red", size = 5) +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = k_values) +
  labs(title    = "Silhouette Scores",
       subtitle = paste("Optimal k =", optimal_k),
       x = "Number of Clusters (k)",
       y = "Average Silhouette Score") +
  theme_bw(base_size = 12, base_family = "serif")

ggsave(file.path(output_dir, "Figure_7_Silhouette.pdf"),
       p7, width = 10, height = 6, device = cairo_pdf)
print(p7)

# ══════════════════════════════════════════════════════════════════════════════
#  FIT FINAL K-MEANS WITH OPTIMAL k
# ══════════════════════════════════════════════════════════════════════════════
set.seed(1234567890, kind = "Mersenne-Twister")
km_final <- kmeans(features_scaled, centers = optimal_k, nstart = 25)

hc_data$Cluster <- as.factor(km_final$cluster)

# ══════════════════════════════════════════════════════════════════════════════
#  PCA (for 2D visualisation)
# ══════════════════════════════════════════════════════════════════════════════
pca_result  <- prcomp(features_scaled, center = FALSE, scale. = FALSE)
var_exp     <- summary(pca_result)$importance[2, 1:2] * 100

pca_df <- data.frame(
  PC1     = pca_result$x[, 1],
  PC2     = pca_result$x[, 2],
  Class   = hc_data$Class,
  Cluster = hc_data$Cluster
)

# Cluster centres in PCA space
centres_pca <- as.data.frame(predict(pca_result, km_final$centers))
centres_pca$Cluster <- factor(seq_len(optimal_k))

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 8 — PCA coloured by k-means cluster
# ══════════════════════════════════════════════════════════════════════════════
p8 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.6) +
  geom_point(data = centres_pca, aes(x = PC1, y = PC2),
             color = "black", size = 4, shape = 8, stroke = 2,
             inherit.aes = FALSE) +
  scale_color_brewer(palette = "Set1") +
  labs(title    = paste0("K-Means Clusters  (k = ", optimal_k, ")"),
       subtitle = "Black * = cluster centres",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_bw(base_size = 12, base_family = "serif")

ggsave(file.path(output_dir, "Figure_8_PCA_Clusters.pdf"),
       p8, width = 10, height = 8, device = cairo_pdf)
print(p8)

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 9 — PCA coloured by true class
# ══════════════════════════════════════════════════════════════════════════════
p9 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = model_colors) +
  labs(title = "True Classes  (PCA Projection)",
       x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.text = element_text(size = 8))

ggsave(file.path(output_dir, "Figure_9_PCA_TrueClasses.pdf"),
       p9, width = 12, height = 8, device = cairo_pdf)
print(p9)

# ══════════════════════════════════════════════════════════════════════════════
#  FIGURE 10 — Cluster sizes bar chart
# ══════════════════════════════════════════════════════════════════════════════
cluster_sizes <- data.frame(
  Cluster = factor(seq_len(optimal_k)),
  Size    = as.numeric(table(km_final$cluster))
)

p10 <- ggplot(cluster_sizes, aes(x = Cluster, y = Size, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Size), vjust = -0.5, size = 4.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Cluster Sizes", x = "Cluster", y = "Count") +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "Figure_10_Cluster_Sizes.pdf"),
       p10, width = 10, height = 6, device = cairo_pdf)
print(p10)

# ══════════════════════════════════════════════════════════════════════════════
#  CLUSTER VALIDITY — cross-tabulation vs true classes
# ══════════════════════════════════════════════════════════════════════════════
cluster_table <- table(Cluster = hc_data$Cluster, Class = hc_data$Class)

# Majority class per cluster
majority_class <- apply(cluster_table, 1, function(row)
  colnames(cluster_table)[which.max(row)])

# Cluster purity (% of majority class in each cluster)
purity <- apply(cluster_table, 1, function(row)
  round(max(row) / sum(row) * 100, 1))

# ══════════════════════════════════════════════════════════════════════════════
#  SAVE CLUSTERING RESULTS
# ══════════════════════════════════════════════════════════════════════════════
write_xlsx(
  list(
    Cluster_vs_Class  = as.data.frame.matrix(cluster_table),
    Majority_Purity   = data.frame(Cluster = names(majority_class),
                                   Majority_Class = majority_class,
                                   Purity_pct     = purity),
    Cluster_Summary   = data.frame(Cluster   = seq_len(optimal_k),
                                   Size      = as.numeric(table(km_final$cluster)),
                                   Within_SS = km_final$withinss),
    Assignments       = hc_data[, c("Class", "Cluster")]
  ),
  file.path(output_dir, "KMeans_Results.xlsx")
)

message("✓ All outputs saved to: ", output_dir)



#End of the code
