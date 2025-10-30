# --- Libraries ---
library(dplyr)
library(ggplot2)
library(GGally)
library(ggfortify)
library(tidyr)
library(patchwork)
library(readxl)

# --- Load Data ---
# Read Excel file
#setwd("C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Code/R/My_R_Work")
df <- read_excel("HC all results.xlsx")
if (!("Model" %in% names(df))) stop("Missing 'Model' column in CSV.")

# --- 1. Find Central Point (H*, C*) per Model ---
central_points <- df %>%
  group_by(Model) %>%
  summarise(
    H_star = median(H_Shannon, na.rm=TRUE),
    C_star = median(C_Shannon, na.rm=TRUE),
    i_central = which.min((H_Shannon - median(H_Shannon, na.rm=TRUE))^2 + (C_Shannon - median(C_Shannon, na.rm=TRUE))^2)[1],
    .groups = 'drop'
  )
write.csv(central_points, "Central_HC_points.csv", row.names=FALSE)

# --- 2. Density Contour Plot in H×C Plane + Central Points ---
p_density <- ggplot(df, aes(x = H_Shannon, y = C_Shannon, color=Model)) +
  geom_density_2d_filled(alpha=0.5) +
  geom_point(size = 1, alpha=0.5) +
  geom_point(data=central_points, aes(x=H_star,y=C_star,color=Model), shape=17, size=5) +
  labs(title="HC Entropy–Complexity Plane", x="Shannon Entropy (H)", y="Complexity (C)") +
  theme_minimal(base_size = 12, base_family = "serif")
ggsave("HC_density_contour.pdf", p_density, width=10, height=7)

# --- 3. Pairwise Feature Scatter Matrix ---
features <- df %>%
  dplyr::select(starts_with("H_"), starts_with("C_"), starts_with("Var_"))
ggpairs(df, columns = which(names(df) %in% colnames(features)), aes(color = Model, alpha = 0.6)) +
  labs(title = "Pairwise Feature Matrix by Model")
ggsave("EDA_feature_pairs.pdf", width=15, height=15)

# --- 4. Boxplots for Each Feature by Model ---
long_df <- df %>%
  pivot_longer(cols = starts_with("H_") | starts_with("C_") | starts_with("Var_"),
               names_to = "Feature", values_to = "Value")
p_box <- ggplot(long_df, aes(x=Model, y=Value, fill=Model)) +
  geom_boxplot(alpha=0.6) + facet_wrap(~Feature, scales="free", ncol=3) +
  theme_minimal(base_size = 12, base_family = "serif") + labs(title="Feature Distributions by Model")
ggsave("EDA_boxplots.pdf", p_box, width=16, height=12)

# --- 5. ANOVA & Kruskal-Wallis for Each Feature ---
feature_names <- colnames(features)
cat("ANOVA and Kruskal-Wallis tests (model discrimination):\n")
stat_results <- data.frame(Feature=character(), ANOVA_p=double(), KW_p=double())
for (f in feature_names) {
  aov_p  <- tryCatch(summary(aov(df[[f]] ~ df$Model))[[1]][["Pr(>F)"]][1], error=function(e) NA)
  kw_p   <- tryCatch(kruskal.test(df[[f]] ~ df$Model)$p.value, error=function(e) NA)
  stat_results <- rbind(stat_results, data.frame(Feature=f, ANOVA_p=aov_p, KW_p=kw_p))
}
write.csv(stat_results, "HC_ANOVA_KW_results.csv", row.names=FALSE)

# --- 6. PCA for All Features ---
pca_res <- prcomp(features, scale. = TRUE, center = TRUE)
summary_pca <- summary(pca_res)
autoplot(pca_res, data = df, colour = 'Model',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 4,
         loadings.colour = 'black', size=2, alpha=0.7) +
  labs(title="PCA: Entropy–Complexity Feature Space") + theme_minimal(base_size = 12, base_family = "serif")
ggsave("PCA_entropy_complexity.pdf", width=10, height=8)

# Scree plot
eig <- pca_res$sdev^2
p_scree <- ggplot(data.frame(PC=1:length(eig), Variance=100*eig/sum(eig)), aes(x=PC, y=Variance)) +
  geom_col(fill="skyblue") + geom_point(color="red", size=2) +
  geom_line(aes(group=1), color="red") +
  labs(title="PCA Scree Plot", x="Principal Component", y="Variance Explained (%)") +
  theme_minimal(base_size = 12, base_family = "serif")
ggsave("PCA_scree_plot.pdf", p_scree, width=8, height=5)

# Save PCA loadings
loadings <- as.data.frame(pca_res$rotation[,1:4])
loadings$Feature <- rownames(loadings)
write.csv(loadings, "PCA_loadings_first4PCs.csv", row.names=FALSE)

# --- 7. Describe and Report Discrimination/EDA ---
cat("--- PCA & EDA complete. Central points, ANOVA/KW stats, and plots saved.\n",
    "See Chagas et al. (2022), Brandmaier (2015), Boaretto et al. (2021) for methods.\n")

##############################################################################################
# --- LOAD LIBRARIES ---
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(GGally)
library(ggfortify)
library(readxl)

# --- 1. LOAD DATA ---
df <- read_excel("HC all results.xlsx")
df$Model <- as.factor(df$Model)
cat("Loaded data size:\n")
print(dim(df))
print(head(df))

# --- 2. COMPUTE MEDIAN (H*, C*) FOR EACH ENTROPY/COMPLEXITY + MODEL ---
summarise_by_median <- function(df, H, C, entropy_type) {
  df %>%
    group_by(Model) %>%
    summarise(
      H_star = median({{ H }}, na.rm=TRUE),
      C_star = median({{ C }}, na.rm=TRUE),
      Entropy_Type = entropy_type,
      n = dplyr::n(),
      .groups = 'drop'
    )
}
med_shannon <- summarise_by_median(df, H_Shannon, C_Shannon, "Shannon")
med_renyi   <- summarise_by_median(df, H_Renyi,   C_Renyi,   "Renyi")
med_tsallis <- summarise_by_median(df, H_Tsallis, C_Tsallis, "Tsallis")
med_fisher  <- summarise_by_median(df, H_Fisher,  C_Fisher,  "Fisher")
medians_all <- bind_rows(med_shannon, med_renyi, med_tsallis, med_fisher)
cat("\nMedian (H*, C*) points by model and entropy type:\n")
print(medians_all)

# --- 3. DENSITY + ELLIPSE PLOTS FOR EACH --- ENTROPY/COMPLEXITY PAIR ---
plot_entropy_complexity <- function(data, H_col, C_col, centroids, title, H_label, C_label) {
  p <- ggplot(data, aes(x = {{ H_col }}, y = {{ C_col }}, color = Model)) +
    geom_point(alpha = 0.4, size = 1) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    geom_density_2d_filled(alpha=0.3, show.legend=FALSE) +
    geom_point(data = centroids, aes(x = H_star, y = C_star, color = Model),
               shape = 8, size = 3, stroke = 1.4) +
    geom_text(data = centroids, aes(x = H_star, y = C_star, label = Model),
              fontface = "bold", vjust = -1, size = 4) +
    labs(title = title, x = H_label, y = C_label) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(
    plot.title = element_text(size=14, face="bold"),
    legend.position = "bottom")
  print(p)
  return(p)
}

p_shannon <- plot_entropy_complexity(df, H_Shannon, C_Shannon, med_shannon,
                                     "Shannon Entropy–Complexity", "H (Shannon)", "C (Shannon)")
p_renyi   <- plot_entropy_complexity(df, H_Renyi, C_Renyi, med_renyi,
                                     "Rényi Entropy–Complexity", "H (Rényi)", "C (Rényi)")
p_tsallis <- plot_entropy_complexity(df, H_Tsallis, C_Tsallis, med_tsallis,
                                     "Tsallis Entropy–Complexity", "H (Tsallis)", "C (Tsallis)")
p_fisher  <- plot_entropy_complexity(df, H_Fisher, C_Fisher, med_fisher,
                                     "Fisher Information–Complexity", "H (Fisher)", "C (Fisher)")

# For viewing combined grid:
print((p_shannon | p_renyi) / (p_tsallis | p_fisher))

# --- 4. ADVANCED EDA: PAIRWISE GRAPH, BOXPLOTS, VIOLINS, VARIANCE TESTS ---
feat_cols <- grep("^(H|C|Var)_", colnames(df), value=TRUE)
cat("\nFeature columns included in EDA:\n")
print(feat_cols)

pairmat <- ggpairs(df, columns = which(names(df) %in% feat_cols), 
                   aes(color = Model, alpha = 0.5)) +
  labs(title = "Pairwise Matrix: Entropy, Complexity, and Variances")
print(pairmat)

long_df <- pivot_longer(df, cols = feat_cols, names_to = "Feature", values_to = "Value")
p_box <- ggplot(long_df, aes(x=Model, y=Value, fill=Model)) +
  geom_boxplot(alpha=0.7) +
  facet_wrap(~Feature, ncol=4, scales="free") +
  theme_minimal(base_size = 12, base_family = "serif") +
  labs(title="Feature Distributions by Model")
print(p_box)

# --- 5. STATISTICAL TESTS (ANOVA, KRUSKAL-WALLIS) ---
stat_results <- data.frame(Feature=character(), ANOVA=double(), Kruskal=double())
for (f in feat_cols) {
  aov_p <- tryCatch(summary(aov(df[[f]] ~ df$Model))[[1]][["Pr(>F)"]][1], error=function(e) NA)
  kw_p  <- tryCatch(kruskal.test(df[[f]] ~ df$Model)$p.value, error=function(e) NA)
  stat_results <- rbind(stat_results, data.frame(Feature=f, ANOVA=aov_p, Kruskal=kw_p))
}
cat("\nStatistical discrimination p-values:\n")
print(stat_results)

# --- 6. PCA: FULL MULTIVARIATE SEPARATION ---
pca_model <- prcomp(df[, feat_cols], scale. = TRUE)
cat("\nPCA summary:\n")
print(summary(pca_model))
cat("\nPCA loadings (first 4 PCs):\n")
print(pca_model$rotation[,1:4])

pca_plot <- autoplot(pca_model, data=df, colour="Model", loadings=TRUE, loadings.label=TRUE, loadings.colour="black", size=2, alpha=0.6) +
  labs(title="PCA: Multifeature Separation") + theme_minimal(base_size = 12, base_family = "serif")
print(pca_plot)

pca_var <- pca_model$sdev^2
scree <- data.frame(PC=seq_along(pca_var), Var=100*pca_var/sum(pca_var))
p_scree <- ggplot(scree, aes(x=PC, y=Var)) + geom_col(fill="skyblue") + geom_point(color="red", size=2) +
  labs(title="PCA Scree Plot",x="Principal Component",y="Variance Explained (%)") + theme_minimal(base_size = 12, base_family = "serif")
print(p_scree)

# --- 7. STRUCTURAL: MODEL SEPARABILITY/FEATURE RANKING (Optional) ---
cat("\nRandom Forest variable importance (Model discrimination):\n")
rf <- randomForest(Model~., data=data.frame(Model=df$Model, df[, feat_cols]), importance=TRUE, ntree=500)
print(importance(rf))
varImpPlot(rf)

# --- (Now you can SAVE any object, after visual inspection above) ---
# Example save:
ggsave("All_HxC_planes.pdf", (p_shannon | p_renyi) / (p_tsallis | p_fisher), width=14, height=12)
write.csv(stat_results, "HC_feature_stat_tests.csv", row.names=FALSE)
write.csv(pca_model$rotation, "PCA_feature_loadings.csv", row.names=TRUE)
write.csv(medians_all, "HC_median_points_all.csv", row.names=FALSE)
