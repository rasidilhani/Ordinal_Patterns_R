# --- Libraries ---
library(dplyr)
library(ggplot2)
library(GGally)
library(ggfortify)
library(tidyr)
library(hexbin)
library(patchwork)
library(readxl)
library(reshape2)


# Load your Excel file
df <- readxl::read_excel("ARMA Model/Summary results.xlsx")

# List all entropy/complexity columns
entropy_cols <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher")
complexity_cols <- c("C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher")

# Ensure all relevant columns are numeric
for (col in c(entropy_cols, complexity_cols)) {
  df[[col]] <- as.numeric(df[[col]])
}

# Remove rows with ANY NA or bad values in selected columns
df <- df %>%
  filter(if_all(all_of(c(entropy_cols, complexity_cols)), ~!is.na(.) & is.finite(.)))

# how many rows removed
cat("Number of valid rows after cleaning:", nrow(df), "\n")

# Model and class parsing, assuming Model like 'AR2_M1' or 'MA2_M3'
df$Group <- sub("_M[1-4]$", "", df$Model)
df$Class <- sub("^.*_(M[1-4])", "\\1", df$Model)

plot_entropy_complexity <- function(df, h_col, c_col, main_title, h_lab, c_lab) {
  p <- ggplot(df, aes_string(x = h_col, y = c_col, color = "Class")) +
    geom_density_2d_filled(alpha = 0.5, contour_var = "density", show.legend = FALSE) +
    geom_point(size = 1.2, alpha = 0.5) +
    scale_color_manual(values = c("M1" = "green", "M2" = "orange",
                                  "M3" = "blue", "M4" = "red")) +
    labs(title = main_title, x = h_lab, y = c_lab, color = "Class") +
    theme_minimal(base_size = 12, base_family = "serif")
  print(p)
  p
}

# Example: plot for Shannon (all models shown together, colored by M class)
p_shannon <- plot_entropy_complexity(df, "H_Shannon", "C_Shannon",
                                     "Shannon Entropy–Complexity by Class", "H (Shannon)", "C (Shannon)")
ggsave("Shannon_HxC_by_class.pdf", p_shannon, width = 7.5, height = 6)

# If you want one plot per AR/MA/ARMA, use
for (mod in unique(df$Group)) {
  df_mod <- df[df$Group == mod, ]
  p <- plot_entropy_complexity(df_mod, "H_Shannon", "C_Shannon",
                               paste0(mod, ": Shannon Entropy–Complexity"), "H (Shannon)", "C (Shannon)")
  ggsave(paste0("Shannon_HxC_by_class_", mod, ".pdf"), p, width = 7.5, height = 6)
}

# Repeat for Renyi, Tsallis, Fisher as needed:
p_renyi <- plot_entropy_complexity(df, "H_Renyi", "C_Renyi",
                                   "Renyi Entropy–Complexity by Class", "H (Renyi)", "C (Renyi)")
ggsave("Renyi_HxC_by_class.pdf", p_renyi, width = 7.5, height = 6)

p_tsallis <- plot_entropy_complexity(df, "H_Tsallis", "C_Tsallis",
                                     "Tsallis Entropy–Complexity by Class", "H (Tsallis)", "C (Tsallis)")
ggsave("Tsallis_HxC_by_class.pdf", p_tsallis, width = 7.5, height = 6)

p_fisher <- plot_entropy_complexity(df, "H_Fisher", "C_Fisher",
                                    "Fisher Information–Complexity by Class", "H (Fisher)", "C (Fisher)")
ggsave("Fisher_HxC_by_class.pdf", p_fisher, width = 7.5, height = 6)


##############################################################################################
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load data
df <- readxl::read_excel("ARMA Model/Summary results.xlsx")

# Ensure numeric and clean data
entropy_cols <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher")
complexity_cols <- c("C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher")
df[entropy_cols] <- lapply(df[entropy_cols], as.numeric)
df[complexity_cols] <- lapply(df[complexity_cols], as.numeric)
df <- df %>%
  filter(if_all(all_of(c(entropy_cols, complexity_cols)), ~!is.na(.) & is.finite(.)))

# Add 'Group' (model class) and M class label columns
df$Group <- sub("_M[1-4]$", "", df$Model)
df$Class <- sub(".*_(M[1-4])$", "\\1", df$Model)

# Fixed colors for classes
class_colors <- c("M1" = "green", "M2" = "orange", "M3" = "blue", "M4" = "red")

# Function to reshape for a single group, include class as color
get_long_for_group <- function(df, group_label) {
  df_group <- filter(df, Group == group_label)
  data.frame(
    Class = rep(df_group$Class, 4),
    Type = rep(c("Shannon", "Renyi", "Tsallis", "Fisher"), each = nrow(df_group)),
    H = c(df_group$H_Shannon, df_group$H_Renyi, df_group$H_Tsallis, df_group$H_Fisher),
    C = c(df_group$C_Shannon, df_group$C_Renyi, df_group$C_Tsallis, df_group$C_Fisher)
  )
}

# Produce and save 6 separate graphs, one for each model group
for (group_label in unique(df$Group)) {
  df_long <- get_long_for_group(df, group_label)
  df_long$Type <- factor(df_long$Type, levels = c("Shannon", "Renyi", "Tsallis", "Fisher"))
  df_long$Class <- factor(df_long$Class, levels = c("M1", "M2", "M3", "M4"))
  p <- ggplot(df_long, aes(x = H, y = C, color = Class, fill = Class)) +
    geom_density_2d_filled(alpha = 0.4, show.legend = FALSE, contour_var = "density") +
    geom_point(size = 1.1, alpha = 0.40) +
    facet_wrap(~ Type, scales = "free", nrow = 2) +
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values = class_colors) +
    labs(title = paste(group_label, ": HxC Density by Class (Shannon/Renyi/Tsallis/Fisher)"),
         x = "Entropy H", y = "Complexity C") +
    theme_minimal(base_size = 12, base_family = "serif") 
  print(p)
  ggsave(paste0("HC_Density_", group_label, "_byClass.pdf"), p, width = 13, height = 7)
}

##############################################################################################
#Hexbin plots for all models together
library(readxl)
library(ggplot2)
library(dplyr)

df <- readxl::read_excel("ARMA Model/Summary results.xlsx")
entropy_cols <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher")
complexity_cols <- c("C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher")
df[entropy_cols] <- lapply(df[entropy_cols], as.numeric)
df[complexity_cols] <- lapply(df[complexity_cols], as.numeric)
df <- df %>%
  filter(if_all(all_of(c(entropy_cols, complexity_cols)), ~!is.na(.) & is.finite(.)))

df$Group <- sub("_M[1-4]$", "", df$Model)
df$Class <- sub(".*_(M[1-4])$", "\\1", df$Model)

get_long_for_group <- function(df, group_label) {
  df_group <- filter(df, Group == group_label)
  data.frame(
    Class = rep(df_group$Class, 4),
    Type = rep(c("Shannon", "Renyi", "Tsallis", "Fisher"), each = nrow(df_group)),
    H = c(df_group$H_Shannon, df_group$H_Renyi, df_group$H_Tsallis, df_group$H_Fisher),
    C = c(df_group$C_Shannon, df_group$C_Renyi, df_group$C_Tsallis, df_group$C_Fisher)
  )
}

for (group_label in unique(df$Group)) {
  df_long <- get_long_for_group(df, group_label)
  df_long$Type <- factor(df_long$Type, levels = c("Shannon", "Renyi", "Tsallis", "Fisher"))
  df_long$Class <- factor(df_long$Class, levels = c("M1", "M2", "M3", "M4"))
  p_hex <- ggplot(df_long, aes(x = H, y = C, color = Class, fill = Class)) +
    stat_binhex(bins = 30, aes(fill = after_stat(count))) +
    facet_wrap(~ Type, scales = "free", nrow = 2) +
    scale_fill_viridis_c(option = "C", guide = "colourbar") +
    labs(title = paste(group_label, ": Hexbin Counts by Class (All Entropies)"),
         x = "Entropy H", y = "Complexity C", fill = "Count") +
    theme_minimal(base_size = 12, base_family = "serif")
  print(p_hex)
  ggsave(paste0("HC_Hexbin_", group_label, ".pdf"), p_hex, width = 13, height = 7)
}

##############################################################################################
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

df <- readxl::read_excel("ARMA Model/Summary results.xlsx")
df$Model <- as.factor(df$Model)

# Fix: convert all entropy/complexity columns to numeric and drop NA/malformed rows
entropy_cols <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher")
complexity_cols <- c("C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher")
for (col in c(entropy_cols, complexity_cols)) {
  df[[col]] <- as.numeric(df[[col]])
}
df <- df %>% filter(if_all(all_of(c(entropy_cols, complexity_cols)), ~is.finite(.) & !is.na(.)))

cat("Loaded data size:\n")
print(dim(df))
print(head(df))


# --- 1. DEFINE MEDIAN COMPUTATION BY MODEL+ENTROPY ---
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

# --- 2. PLOTTING FUNCTION: HIGHLIGHT (H*, C*) POINT ---
plot_entropy_complexity <- function(data, H_col, C_col, centroids, title, H_label, C_label) {
  p <- ggplot(data, aes(x = {{ H_col }}, y = {{ C_col }}, color = Model)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_density_2d_filled(alpha=0.3, show.legend=FALSE) +
    stat_ellipse(level = 0.95, linetype = "dashed") +
    geom_point(data = centroids, aes(x = H_star, y = C_star, color = Model),
               shape = 8, size = 4, stroke = 1.5, fill = "yellow") +
    geom_text(data = centroids, aes(x = H_star, y = C_star, label = Model),
              fontface = "bold", vjust = -1, size = 4) +
    labs(title = title, x = H_label, y = C_label) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(
      plot.title = element_text(size=12, face="bold"),
      legend.position = "bottom")
  print(p)
  return(p)
}

# --- 3. PLOTS FOR ALL ENTROPY-TYPE PAIRS ---
p_shannon <- plot_entropy_complexity(df, H_Shannon, C_Shannon, med_shannon,
                                     "Shannon Entropy–Complexity", "H (Shannon)", "C (Shannon)")
p_renyi   <- plot_entropy_complexity(df, H_Renyi, C_Renyi, med_renyi,
                                     "Rényi Entropy–Complexity", "H (Rényi)", "C (Rényi)")
p_tsallis <- plot_entropy_complexity(df, H_Tsallis, C_Tsallis, med_tsallis,
                                     "Tsallis Entropy–Complexity", "H (Tsallis)", "C (Tsallis)")
p_fisher  <- plot_entropy_complexity(df, H_Fisher, C_Fisher, med_fisher,
                                     "Fisher Information–Complexity", "H (Fisher)", "C (Fisher)")

# View all together (patchwork)
(p_shannon | p_renyi) / (p_tsallis | p_fisher)
ggsave("HC_All_Entropy_Types_with_Medians.pdf", 
       (p_shannon | p_renyi) / (p_tsallis | p_fisher), 
       width = 14, height = 10)
##############################################################################################

#Emblaematic points for each model

library(dplyr)

get_CH_emblematic <- function(df, h_col, c_col, model_col = "Model") {
  # For each model, compute PC1/PC2 for (H, C)
  res <- df %>%
    group_by_at(model_col) %>%
    group_modify(~{
      dat <- select(.x, all_of(h_col), all_of(c_col))
      pcs <- prcomp(dat, scale. = FALSE)
      pcscores <- pcs$x
      N <- nrow(dat)
      ordering <- order(pcscores[,1])
      median_index <- ordering[ceiling((N + 1)/2)]
      pc_median <- pcscores[median_index,]
      # Project back to H-C plane (here, just pick the original row)
      em <- dat[median_index, , drop=FALSE]
      out <- tibble(
        H_star = em[[h_col]],
        C_star = em[[c_col]],
        Emblematic_row = median_index
      )
      return(out)
    }) %>% ungroup()
  return(res)
}

# Example for Shannon:
emblematic_shannon <- get_CH_emblematic(df, "H_Shannon", "C_Shannon")
print(emblematic_shannon)
# You can repeat for other entropy types as needed
emblematic_renyi <- get_CH_emblematic(df, "H_Renyi", "C_Renyi")
emblematic_tsallis <- get_CH_emblematic(df, "H_Tsallis", "C_Tsallis")
emblematic_fisher <- get_CH_emblematic(df, "H_Fisher", "C_Fisher")
print(emblematic_renyi)
print(emblematic_tsallis)
print(emblematic_fisher)
write.csv(emblematic_shannon, "Emblematic_Points_Shannon.csv", row.names = FALSE)
write.csv(emblematic_renyi, "Emblematic_Points_Renyi.csv", row.names = FALSE)
write.csv(emblematic_tsallis, "Emblematic_Points_Tsallis.csv", row.names = FALSE)
write.csv(emblematic_fisher, "Emblematic_Points_Fisher.csv", row.names = FALSE)
##############################################################################################
# Now you can use these emblematic points in your plots as needed
##############################################################################################
library(readxl)
library(ggplot2)
library(dplyr)

# Read emblematic results
emblem <- readxl::read_excel("ARMA Model/Emblematic results.xlsx")

# Ensure correct type for plotting
emblem$Entropy_Type <- factor(emblem$Entropy_Type, levels = c("Shannon", "Renyi", "Tsallis", "Fisher"))

# Plot all emblematic points colored by entropy type
p <- ggplot(emblem, aes(x = H_star, y = C_star, color = Entropy_Type, shape = Entropy_Type)) +
  geom_point(size = 3, stroke = 1.4) +
  geom_text(aes(label = Model), nudge_y = 0.01, size = 3.5, fontface = "bold", show.legend = FALSE) +
  labs(
    title = expression(paste("Median Points (", H^"*", ", ", C^"*", ") for Each Model and Entropy")),
    x = expression(H^"*"),
    y = expression(C^"*")
  ) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "bottom")
print(p)
ggsave("Median_HC_by_entropy.pdf", p, width = 9, height = 7)

# Optional: facet by entropy type if you want per-panel comparison
p_facet <- ggplot(emblem, aes(x = H_star, y = C_star, color = Model)) +
  geom_point(size = 3, stroke = 0.8) +
  geom_text(aes(label = Model), nudge_y = 0.01, size = 3, fontface = "bold", show.legend = FALSE) +
  facet_wrap(~Entropy_Type) +
  labs(
    title = expression(paste("Faceted Median Points per Entropy Type")),
    x = expression(H^"*"),
    y = expression(C^"*")
  ) +
  theme_minimal(base_size = 12, base_family = "serif")
print(p_facet)
ggsave("Median_HC_facet_entropy.pdf", p_facet, width = 12, height = 8)
##############################################################################################