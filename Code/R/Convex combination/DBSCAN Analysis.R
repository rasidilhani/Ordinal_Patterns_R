library(readxl)
library(dplyr)
library(ggplot2)
library(dbscan)
library(factoextra)
library(writexl)
library(here)

# ─────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────
D        <- 4
n_val    <- 10000
eps_val  <- 0.8
minPts   <- 14

data_path  <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                   paste0("HC_Results_D", D, ".xlsx"))
output_dir <- here("Results", "Convex_combination", "Clustering", "DBSCAN",
                   paste0("n", n_val))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────
# FEATURES
# ─────────────────────────────────────────────────────────────
features <- c(
  "H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher",
  "C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher",
  "Disequilibrium",
  "Var_H_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Var_C_Shannon"
)

# ─────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────
hc_data <- read_xlsx(data_path) %>%
  select(Model, all_of(features)) %>%
  rename(Class = Model)

X <- scale(hc_data %>% select(all_of(features)))

# ─────────────────────────────────────────────────────────────
# kNN DISTANCE PLOT
# ─────────────────────────────────────────────────────────────
pdf(file.path(output_dir, paste0("kNN_D", D, ".pdf")), width = 7, height = 5)
par(family = "serif", bg = "white")
kNNdistplot(X, k = minPts)
abline(h = eps_val, col = "red", lty = 2)
dev.off()

# ─────────────────────────────────────────────────────────────
# DBSCAN
# ─────────────────────────────────────────────────────────────
db <- dbscan(X, eps = eps_val, minPts = minPts)

cat("Clusters:", length(unique(db$cluster)) - (0 %in% db$cluster), "\n")
cat("Noise points:", sum(db$cluster == 0), "\n")
print(table(hc_data$Class, db$cluster))

hc_data$Cluster <- db$cluster

# ─────────────────────────────────────────────────────────────
# PCA + ggplot CLUSTER PLOT
# ─────────────────────────────────────────────────────────────
pca <- prcomp(X)
var_ex <- round(summary(pca)$importance[2, 1:2] * 100, 1)

df  <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Cluster = factor(db$cluster),
  Class = hc_data$Class
)

p1 <- ggplot(df, aes(PC1, PC2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(
    title = paste0("DBSCAN (eps=", eps_val, ", minPts=", minPts, ")"),
    x = paste0("PC1 (", var_ex[1], "% variance)"),
    y = paste0("PC2 (", var_ex[2], "% variance)")
  ) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
print(p1)
ggsave(file.path(output_dir, paste0("DBSCAN_PCA_D", D, ".pdf")),
       p1, width = 8, height = 6)

# ─────────────────────────────────────────────────────────────
# fviz_cluster PLOT
# ─────────────────────────────────────────────────────────────
p2 <- fviz_cluster(
  list(data = X, cluster = db$cluster),
  geom = "point",
  ellipse = TRUE,
  main = paste0("fviz_cluster | eps=", eps_val)
) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
print(p2)
ggsave(file.path(output_dir, paste0("fviz_cluster_D", D, ".pdf")),
       p2, width = 8, height = 6)

# ─────────────────────────────────────────────────────────────
# SAVE CLUSTER ASSIGNMENTS
# ─────────────────────────────────────────────────────────────
write_xlsx(
  hc_data,
  file.path(output_dir, paste0("DBSCAN_clusters_D", D, "_n", n_val, ".xlsx"))
)

# ─────────────────────────────────────────────────────────────
# CLASS DISTRIBUTION PER CLUSTER
# ─────────────────────────────────────────────────────────────
p3 <- hc_data %>%
  count(Cluster, Class) %>%
  ggplot(aes(x = factor(Cluster), y = n, fill = Class)) +
  geom_col(position = "fill", color = "grey30", linewidth = 0.2) +
  labs(
    title = "Class Distribution per DBSCAN Cluster",
    x     = "Cluster",
    y     = "Proportion",
    fill  = "Class"
  ) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(nrow = 4))

print(p3)
ggsave(file.path(output_dir, paste0("DBSCAN_class_dist_D", D, ".pdf")),
       p3, width = 9, height = 7)

# ─────────────────────────────────────────────────────────────
# FEATURE BOXPLOTS PER CLUSTER
# ─────────────────────────────────────────────────────────────
p4 <- hc_data %>%
  tidyr::pivot_longer(all_of(features), names_to = "Feature", values_to = "Value") %>%
  ggplot(aes(x = factor(Cluster), y = Value, fill = factor(Cluster))) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 4) +
  labs(
    title = "Feature Profiles by Cluster",
    x     = "Cluster",
    y     = "Value"
  ) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p4)
ggsave(file.path(output_dir, paste0("DBSCAN_feature_boxplots_D", D, ".pdf")),
       p4, width = 14, height = 10)


