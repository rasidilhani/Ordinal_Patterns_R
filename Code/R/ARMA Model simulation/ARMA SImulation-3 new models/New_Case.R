# =====================================================
# 📊 General HC Shannon Analysis for Any Case
# =====================================================

library(readxl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(viridis)
library(here)

# ---- 1. Define Case ----
# Just change this line for Case1, Case2, ..., Case6
case_name <- "New_Case1"

# ---- 2. Define file paths ----
data_file <- here("GitHub", "Ordinal_Patterns_R", "Data",
                  "New ARMA Time series results", case_name,
                  paste0("HC_Shannon_Results_", case_name, ".xlsx"))

output_dir <- here("GitHub", "Ordinal_Patterns_R", "Data",
                   "New ARMA Time series results", case_name)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

plots_dir  <- file.path(output_dir, "plots_HC")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

sheets_to_process <- c("HC_results_5000", "HC_results_10000")

# ---- 3. PCA helper function ----
get_pca_center <- function(df_model) {
  dat <- as.matrix(df_model[, c("H_Shannon", "C_Shannon")])
  pca <- prcomp(dat, center = TRUE, scale. = FALSE)
  u <- pca$x[, 1]
  N <- length(u)
  r <- order(u)
  median_idx <- r[ceiling((N + 1) / 2)]
  emblem_idx <- df_model$.orig_row[median_idx]
  list(
    H_star = df_model$H_Shannon[median_idx],
    C_star = df_model$C_Shannon[median_idx],
    Central_point_Row = emblem_idx
  )
}

# ---- 4. Define consistent color/shape mappings ----
model_colors <- c(
  "ARMA11_M1" = "#1b9e77",
  "AR2_M1"    = "#d95f02",
  "MA1_M2"    = "#7570b3",
  "ARMA11_M3" = "#e7298a",
  "AR2_M4"    = "#66a61e",
  "MA1_M4"    = "#e6ab02",
  "ARMA22_M2" = "#a6761d",
  "AR2_M3"    = "#666666",
  "MA2_M3"    = "#1f78b4",
  "ARMA22_M1" = "#b15928",
  "ARMA22_M4" = "#6a3d9a",
  "MA2_M2"    = "#33a02c",
  "ARMA22_M3" = "#fb9a99"
)

model_shapes <- c(
  "ARMA11_M1" = 21, "AR2_M1" = 22, "MA1_M2" = 23,
  "ARMA11_M3" = 24, "AR2_M4" = 25, "MA1_M4" = 8,
  "ARMA22_M2" = 15, "AR2_M3" = 16, "MA2_M3" = 17,
  "ARMA22_M1" = 18, "ARMA22_M4" = 19, "MA2_M2" = 4, "ARMA22_M3" = 3
)

# ---- 5. Process both sheets (5000 & 10000) ----
central_points_all <- list()

for (sheet_name in sheets_to_process) {
  message("Processing: ", sheet_name)
  
  df <- read_excel(data_file, sheet = sheet_name) %>%
    mutate(.orig_row = row_number())
  
  models <- unique(df$Model)
  colors <- model_colors[models]
  shapes <- model_shapes[models]
  
  # ---------------- Plot 1: Entropy–Complexity Scatter ----------------
  p1 <- ggplot(df, aes(x = H_Shannon, y = C_Shannon, color = Model, shape = Model)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(
      title = paste("Shannon Entropy–Complexity Plane (", sheet_name, ")", sep = ""),
      x = expression(italic(H)), y = expression(italic(C))
    ) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(file.path(plots_dir, paste0("HC_Scatter_", sheet_name, ".pdf")),
         p1, width = 8, height = 6)
  print(p1)
  
  # ---------------- Compute PCA Central Points ----------------
  centers <- data.frame(Model = character(),
                        Entropy_Type = character(),
                        H_star = numeric(),
                        C_star = numeric(),
                        Central_point_Row = integer(),
                        stringsAsFactors = FALSE)
  
  for (m in models) {
    sub <- df %>% filter(Model == m, !is.na(H_Shannon), !is.na(C_Shannon))
    if (nrow(sub) < 3) next
    cpt <- get_pca_center(sub)
    centers <- rbind(centers, data.frame(Model = m,
                                         Entropy_Type = sheet_name,
                                         H_star = cpt$H_star,
                                         C_star = cpt$C_star,
                                         Central_point_Row = cpt$Central_point_Row))
  }
  
  central_points_all[[sheet_name]] <- centers
  
  # ---------------- Plot 2: Central Points (H*, C*) ----------------
  p2 <- ggplot(centers, aes(x = H_star, y = C_star, color = Model, shape = Model, label = Model)) +
    geom_point(size = 3, stroke = 1.2, fill = "white") +
    geom_text(vjust = -0.8, size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(
      title = paste("Central Points (H*, C*) – ", sheet_name),
      x = expression(italic(H^"*")), y = expression(italic(C^"*"))
    ) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  ggsave(file.path(plots_dir, paste0("HC_CentralPoints_", sheet_name, ".pdf")),
         p2, width = 8, height = 6)
  print(p2)
  
  # ---------------- Plot 3: Heatmap ----------------
  p3 <- ggplot(df, aes(x = H_Shannon, y = C_Shannon)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, bins = 60, alpha = 0.6) +
    scale_fill_viridis_c(option = "plasma") +
    geom_point(aes(color = Model, shape = Model), size = 2) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    labs(
      title = paste("Heatmap of Shannon Entropy–Complexity Plane (", sheet_name, ")", sep = ""),
      x = expression(italic(H)), y = expression(italic(C))
    ) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(legend.position = "right", legend.title = element_blank())
  
  ggsave(file.path(plots_dir, paste0("HC_Heatmap_", sheet_name, ".pdf")),
         p3, width = 8, height = 6)
  print(p3)
}

# ---- 6. Save Central Points ----
wb <- createWorkbook()
for (sheet_name in names(central_points_all)) {
  addWorksheet(wb, paste0("Central_", sheet_name))
  writeData(wb, paste0("Central_", sheet_name), central_points_all[[sheet_name]])
}
out_excel <- file.path(output_dir, paste0("HC_Shannon_Central_Points_", case_name, ".xlsx"))
saveWorkbook(wb, out_excel, overwrite = TRUE)

write.csv(bind_rows(central_points_all, .id = "Sheet"),
          file.path(output_dir, paste0("HC_Shannon_Central_Points_", case_name, ".csv")),
          row.names = FALSE)

cat("\n✅ All plots and results saved for", case_name, "in:\n", output_dir, "\n")
#################################################################################################################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 