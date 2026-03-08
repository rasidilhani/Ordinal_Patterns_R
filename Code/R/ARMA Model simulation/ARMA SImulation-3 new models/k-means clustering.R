library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(StatOrdPattHxC)

# --- File Paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"

data("LinfLsup")

D <- 3

cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1","MA1_M1","MA1_M2","MA2_M1","ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4","MA1_M3","MA1_M4","MA2_M4","ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3","MA2_M2","MA2_M3","ARMA22_M2","ARMA22_M3")
)

sample_sizes <- c(1000, 5000)

model_colors <- c(
  "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
  "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
  "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
  "ARMA22_M3"="#fb9a99","AR1_M1"="#8dd3c7","AR1_M2"="red","ARMA11_M2"="#bebada",
  "MA1_M1"="#fb8072","MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69",
  "ARMA11_M4"="#fccde5","MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
)

model_shapes <- c(
  "ARMA11_M1"=21,"AR2_M1"=22,"MA1_M2"=23,"ARMA11_M3"=24,"AR2_M4"=25,"MA1_M4"=8,
  "ARMA22_M2"=15,"AR2_M3"=16,"MA2_M3"=17,"ARMA22_M1"=18,"ARMA22_M4"=19,"MA2_M2"=4,
  "ARMA22_M3"=3,"AR1_M1"=1,"AR1_M2"=2,"ARMA11_M2"=5,"MA1_M1"=6,"MA2_M1"=7,
  "AR1_M3"=9,"AR1_M4"=10,"ARMA11_M4"=11,"MA1_M3"=12,"MA2_M4"=13,"AR2_M2"=14
)

sheet_names <- excel_sheets(data_path)

# ----------------------------------------------------
# 1. Read all sheets and combine
# ----------------------------------------------------
# IMPORTANT: adjust column names below if needed.
# I assume columns: Model, n, H, C, Rep (or similar)
all_data <- purrr::map_dfr(
  sheet_names,
  ~ read_excel(data_path, sheet = .x) %>%
    mutate(Sheet = .x)
)

# Check the structure once (uncomment to inspect)
# str(all_data)
# head(all_data)

# ----------------------------------------------------
# 2. Filter only n = 5000
# ----------------------------------------------------
# If your sample-size column is not called "n", change it here.
df5000 <- all_data %>%
  filter(n == 5000)

# Again, check
# table(df5000$n)

# ----------------------------------------------------
# 3. Choose features for clustering
# ----------------------------------------------------
# Assuming the file has columns "H" and "C" (complexity-entropy plane).
# Add more features here if you want (e.g., extra entropies/complexities).
# ----------------------------------------------------
# 3. Choose features for clustering  (UPDATED)
# ----------------------------------------------------
feature_cols <- c("H_Shannon", "C_Shannon")

df5000 <- df5000 %>%
  select(Model, n, all_of(feature_cols), everything())

# Remove rows with NA in features (if any)
df5000 <- df5000 %>%
  drop_na(all_of(feature_cols))

# Matrix of features, scaled
X_all <- scale(df5000[, feature_cols])


# ----------------------------------------------------
# 4. K-means for ALL models together (n = 5000)
# ----------------------------------------------------
K_all <- length(unique(df5000$Model))  # one cluster per model type (you can change this)

set.seed(123)
km_all <- kmeans(X_all, centers = K_all, nstart = 50)

df5000$Cluster_all <- factor(km_all$cluster)

# Confusion table: how clusters match model types
cat("\n=== Confusion table: ALL MODELS (n=5000) ===\n")
print(table(df5000$Model, df5000$Cluster_all))

# PCA for 2D plotting
pc_all <- prcomp(X_all)
pc_df_all <- cbind(df5000, as.data.frame(pc_all$x[, 1:2]))

# True model types in PC space
p_true_all <- ggplot(pc_df_all, aes(x = PC1, y = PC2,
                                    color = Model,
                                    shape = Model)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = model_colors) +
  scale_shape_manual(values = model_shapes) +
  labs(title = "PC1–PC2: True model types (n = 5000)") +
  theme_minimal()

# K-means clusters (colour) vs true model types (shape)
p_clust_all <- ggplot(pc_df_all, aes(x = PC1, y = PC2,
                                     color = Cluster_all,
                                     shape = Model)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_shape_manual(values = model_shapes) +
  labs(title = "PC1–PC2: k-means clusters (colour) vs model types (shape)\nAll models, n = 5000") +
  theme_minimal()

print(p_true_all)
print(p_clust_all)

# ----------------------------------------------------
# 5. K-means separately for each CASE (n = 5000)
# ----------------------------------------------------
case_results <- list()

for (case_name in names(cases)) {
  cat("\n========================================\n")
  cat("Case:", case_name, "\n")
  cat("========================================\n")
  
  models_in_case <- cases[[case_name]]
  
  df_case <- df5000 %>%
    filter(Model %in% models_in_case)
  
  if (nrow(df_case) == 0) {
    cat("No rows for", case_name, "with n = 5000.\n")
    next
  }
  
  X_case <- scale(df_case[, feature_cols])
  
  K_case <- length(unique(df_case$Model))  # clusters = number of models in this case
  
  set.seed(123)
  km_case <- kmeans(X_case, centers = K_case, nstart = 50)
  
  df_case$Cluster_case <- factor(km_case$cluster)
  
  # Store for later if needed
  case_results[[case_name]] <- df_case
  
  # Confusion table for this case
  cat("\nConfusion table (Model vs Cluster) for", case_name, ":\n")
  print(table(df_case$Model, df_case$Cluster_case))
  
  # PCA
  pc_case <- prcomp(X_case)
  pc_df_case <- cbind(df_case, as.data.frame(pc_case$x[, 1:2]))
  
  # Plot: true model types
  p_true <- ggplot(pc_df_case, aes(x = PC1, y = PC2,
                                   color = Model,
                                   shape = Model)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    labs(title = paste0("PC1–PC2: True model types, ", case_name, " (n = 5000)")) +
    theme_minimal()
  
  # Plot: clusters vs model types
  p_clust <- ggplot(pc_df_case, aes(x = PC1, y = PC2,
                                    color = Cluster_case,
                                    shape = Model)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_shape_manual(values = model_shapes) +
    labs(title = paste0("PC1–PC2: k-means clusters vs model types, ",
                        case_name, " (n = 5000)")) +
    theme_minimal()
  
  print(p_true)
  print(p_clust)
}
# End of script