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
minPts   <- 13

data_path  <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                   paste0("HC_Results_D", D, ".xlsx"))
output_dir <- here("Results", "Convex_combination", "Clustering",
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


#########################################################################
library(readxl)
library(dplyr)
library(ggplot2)
library(dbscan)
library(writexl)
library(here)

# ── Parameters ────────────────────────────────────────────────────────────────
D          <- 4
n_val      <- 10000
eps_val    <- 0.8
minPts_val <- 14

# ── Paths ─────────────────────────────────────────────────────────────────────
data_path  <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                   paste0("HC_Results_D", D, ".xlsx"))
output_dir <- here("Results", "Convex_combination", "Clustering", "DBSCAN",
                   paste0("n", n_val))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ── Cluster colors (0 = noise in grey) ───────────────────────────────────────
cluster_palette <- c(
  "0" = "grey70",
  "1" = "#E41A1C",
  "2" = "#377EB8",
  "3" = "#4DAF4A",
  "4" = "#FF7F00",
  "5" = "#984EA3",
  "6" = "#A6761D",
  "7" = "#66C2A5",
  "8" = "#F781BF",
  "9" = "#999999"
)

# ── Features ──────────────────────────────────────────────────────────────────
features <- c(
  "H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher",
  "C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher",
  "Disequilibrium",
  "Var_H_Shannon", "Var_H_Renyi", "Var_H_Tsallis",
  "Var_C_Shannon"
)

# ── Common theme ──────────────────────────────────────────────────────────────
theme_clean <- theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(face = "bold")
  )

# ── Load and scale data ───────────────────────────────────────────────────────
hc_data <- read_xlsx(data_path) %>%
  select(Model, all_of(features)) %>%
  rename(Class = Model) %>%
  mutate(Class = factor(Class))

X <- scale(hc_data %>% select(all_of(features)))

# ── kNN distance plot ─────────────────────────────────────────────────────────
pdf(file.path(output_dir, paste0("kNN_D", D, ".pdf")), width = 7, height = 5)
par(family = "serif", bg = "white")
kNNdistplot(X, k = minPts_val)
abline(h = eps_val, col = "red", lty = 2, lwd = 1.5)
title(main = paste0("k-NN Distance Plot  (k = ", minPts_val, ",  D = ", D, ")"))
legend("topleft", legend = paste0("eps = ", eps_val),
       col = "red", lty = 2, lwd = 1.5, bty = "n")
dev.off()

# ── Run DBSCAN ────────────────────────────────────────────────────────────────
db <- dbscan(X, eps = eps_val, minPts = minPts_val)

cat("Clusters found :", length(unique(db$cluster)) - (0 %in% db$cluster), "\n")
cat("Noise points   :", sum(db$cluster == 0), "\n")

# ── Add cluster labels BEFORE building any plot data ─────────────────────────
hc_data <- hc_data %>%
  mutate(Cluster = factor(db$cluster))

# ── PCA (built AFTER cluster labels exist) ───────────────────────────────────
pca    <- prcomp(X, center = TRUE, scale. = FALSE)
var_ex <- round(summary(pca)$importance[2, 1:2] * 100, 1)

df_pca <- data.frame(
  PC1     = pca$x[, 1],
  PC2     = pca$x[, 2],
  Cluster = hc_data$Cluster,
  Class   = hc_data$Class
)

# ── Plot 1: Cluster sizes ─────────────────────────────────────────────────────
p_count <- hc_data %>%
  count(Cluster) %>%
  ggplot(aes(x = Cluster, y = n, fill = Cluster)) +
  geom_col(show.legend = FALSE, color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.4, family = "serif", size = 3.5) +
  scale_fill_manual(values = cluster_palette) +
  labs(title = "DBSCAN Cluster Sizes", x = "Cluster", y = "Count") +
  theme_clean

ggsave(file.path(output_dir, "DBSCAN_cluster_sizes.pdf"),
       p_count, width = 7, height = 5)

# ── Plot 2: PCA colored by cluster ───────────────────────────────────────────
p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_manual(values = cluster_palette) +
  labs(
    title = paste0("DBSCAN — PCA  (eps = ", eps_val, ",  minPts = ", minPts_val, ")"),
    x     = paste0("PC1 (", var_ex[1], "%)"),
    y     = paste0("PC2 (", var_ex[2], "%)"),
    color = "Cluster"
  ) +
  theme_clean +
  theme(legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(output_dir, "DBSCAN_PCA.pdf"),
       p_pca, width = 9, height = 7)

# ── Plot 3: Class distribution per cluster ────────────────────────────────────
p_dist <- hc_data %>%
  count(Cluster, Class) %>%
  ggplot(aes(x = Cluster, y = n, fill = Cluster)) +
  geom_col(position = "fill", color = "grey30", linewidth = 0.2) +
  scale_fill_manual(values = cluster_palette) +
  geom_text(aes(label = Class), position = position_fill(vjust = 0.5),
            size = 2.5, family = "serif") +
  labs(
    title = "Class Distribution per DBSCAN Cluster",
    x     = "Cluster", y = "Proportion", fill = "Cluster"
  ) +
  theme_clean +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "DBSCAN_class_distribution.pdf"),
       p_dist, width = 9, height = 7)

# ── Plot 4: Feature boxplots per cluster ──────────────────────────────────────
p_box <- hc_data %>%
  tidyr::pivot_longer(all_of(features), names_to = "Feature", values_to = "Value") %>%
  ggplot(aes(x = Cluster, y = Value, fill = Cluster)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = cluster_palette) +
  facet_wrap(~ Feature, scales = "free_y", ncol = 4) +
  labs(title = "Feature Profiles by Cluster", x = "Cluster", y = "Value") +
  theme_clean +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "DBSCAN_feature_boxplots.pdf"),
       p_box, width = 14, height = 10)

# ── Save Excel results ────────────────────────────────────────────────────────
write_xlsx(hc_data,
           file.path(output_dir, paste0("DBSCAN_clusters_D", D, ".xlsx")))

write_xlsx(hc_data %>% count(Cluster),
           file.path(output_dir, "DBSCAN_cluster_counts.xlsx"))

write_xlsx(
  as.data.frame.matrix(table(hc_data$Class, hc_data$Cluster)),
  file.path(output_dir, "DBSCAN_confusion_table.xlsx")
)

write_xlsx(
  hc_data %>% group_by(Cluster) %>% summarise(across(all_of(features), mean)),
  file.path(output_dir, "DBSCAN_cluster_feature_means.xlsx")
)

cat("All plots and Excel files saved to:", output_dir, "\n")