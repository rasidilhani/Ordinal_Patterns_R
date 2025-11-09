# =============================
# 📊 Entropy–Complexity Contour & Scatter Plots for Time Series Models
# =============================
library(here)
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(viridis)
library(gridExtra)

# ---- 1) Load Excel sheet ----
xlsx_path <- here("GitHub", "Ordinal_Patterns_R", "Data", "ARMA Time series results", "Summary results.xlsx")


sheet_index <- 4
df_raw <- read_excel(xlsx_path, sheet = sheet_index)
cat("✅ Loaded sheet", sheet_index, "with columns:\n")
print(colnames(df_raw))

# ---- 2) Detect main columns automatically ----
colnames_lower <- tolower(colnames(df_raw))

find_col <- function(patterns) {
  idx <- which(sapply(patterns, function(p) grepl(p, colnames_lower)))
  if (length(idx) == 0) return(NA_character_)
  return(colnames(df_raw)[idx[1]])
}

mapping <- list(
  Shannon_H = find_col(c("h_shannon", "shannon")),
  Shannon_C = find_col(c("c_shannon", "complexity_shannon")),
  Renyi_H   = find_col(c("h_renyi", "renyi")),
  Renyi_C   = find_col(c("c_renyi")),
  Tsallis_H = find_col(c("h_tsallis", "tsallis")),
  Tsallis_C = find_col(c("c_tsallis")),
  Fisher_H  = find_col(c("h_fisher", "fisher")),
  Fisher_C  = find_col(c("c_fisher"))
)

cat("\nDetected column mapping:\n")
print(mapping)

# ---- 3) Detect model column ----
model_col <- find_col(c("model", "series", "type"))
if (is.na(model_col)) stop("❌ Could not detect a model column (expected 'Model', 'Type', etc.).")

# Extract model names (e.g., AR1, AR2, MA1, MA2, ARMA11, ARMA22)
df_raw$Model <- toupper(df_raw[[model_col]])
df_raw$Model <- gsub("\\s+", "", df_raw$Model)

# ---- 4) Build tidy combined dataset ----
entropy_pairs <- list(
  Shannon = c(mapping$Shannon_H, mapping$Shannon_C),
  Renyi   = c(mapping$Renyi_H, mapping$Renyi_C),
  Tsallis = c(mapping$Tsallis_H, mapping$Tsallis_C),
  Fisher  = c(mapping$Fisher_H, mapping$Fisher_C)
)

all_entropy_data <- list()

for (entropy_name in names(entropy_pairs)) {
  pair <- entropy_pairs[[entropy_name]]
  if (any(is.na(pair))) next
  
  temp_df <- df_raw %>%
    select(all_of(c("Model", pair))) %>%
    rename(
      Entropy = !!pair[1],
      Complexity = !!pair[2]
    ) %>%
    mutate(
      Entropy = as.numeric(Entropy),
      Complexity = as.numeric(Complexity),
      Entropy_Type = entropy_name
    )
  
  all_entropy_data[[entropy_name]] <- temp_df
}

combined_df <- bind_rows(all_entropy_data)

# ---- 5) Save combined results ----
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Case 3 results"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

output_csv <- file.path(plot_dir, "Entropy_Complexity_AllModels_case3.csv")
write.csv(combined_df, output_csv, row.names = FALSE)
cat("\n✅ Combined entropy–complexity data saved as:", output_csv, "\n")

# ---- 6) Summary statistics ----
summary_stats <- combined_df %>%
  group_by(Entropy_Type, Model) %>%
  summarise(
    N = n(),
    Entropy_mean = mean(Entropy, na.rm = TRUE),
    Entropy_sd = sd(Entropy, na.rm = TRUE),
    Complexity_mean = mean(Complexity, na.rm = TRUE),
    Complexity_sd = sd(Complexity, na.rm = TRUE),
    .groups = "drop"
  )

summary_csv <- file.path(plot_dir, "Entropy_Complexity_Model_Summary_case3.csv")
write.csv(summary_stats, summary_csv, row.names = FALSE)
cat("✅ Summary statistics saved as:", summary_csv, "\n")

# ---- 7) Ensure output folder exists ----
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Case 3 results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 8) Clean data and ensure correct model types ----
combined_df <- combined_df %>%
  mutate(
    Entropy = as.numeric(Entropy),
    Complexity = as.numeric(Complexity),
    Model = toupper(trimws(Model))
  ) %>%
  filter(!is.na(Entropy), !is.na(Complexity))

# Try to extract clear model names like AR1, AR2, MA1, MA2, ARMA11, ARMA22
combined_df$Model <- gsub("[^A-Z0-9]", "", combined_df$Model)
combined_df$Model <- case_when(
  grepl("ARMA11", combined_df$Model) ~ "ARMA11",
  grepl("ARMA22", combined_df$Model) ~ "ARMA22",
  grepl("AR2", combined_df$Model) ~ "AR2",
  grepl("AR1", combined_df$Model) ~ "AR1",
  grepl("MA2", combined_df$Model) ~ "MA2",
  grepl("MA1", combined_df$Model) ~ "MA1",
  TRUE ~ combined_df$Model
)

# ---- 9) Define consistent shapes and colors ----
model_shapes <- c(
  "AR1" = 16, "AR2" = 17,
  "MA1" = 15, "MA2" = 0,
  "ARMA11" = 18, "ARMA22" = 1
)

model_colors <- c(
  "AR1" = "blue", "AR2" = "red",
  "MA1" = "brown", "MA2" = "yellow",
  "ARMA11" = "orange", "ARMA22" = "black"
)

# ---- 10) Plot for each entropy type ----
plots_list <- list()

for (etype in unique(combined_df$Entropy_Type)) {
  temp <- combined_df %>% filter(Entropy_Type == etype)
  
  cat("\nPlotting:", etype, "  ->  rows:", nrow(temp), "\n")
  
  p <- ggplot(temp, aes(x = Entropy, y = Complexity, color = Model, shape = Model)) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 1) +
    geom_density_2d(color = "grey60", alpha = 0.4, size = 0.4) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    scale_shape_manual(values = model_shapes, drop = FALSE) +
    labs(
      title = paste(etype),
      x = expression(H),
      y = expression(C)
    ) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(face = "bold")
    )
  
  plots_list[[etype]] <- p
  ggsave(filename = file.path(output_dir, paste0(etype, "_Plane_case3.pdf")),
         plot = p, width = 7, height = 5)
  
  cat("✅ Saved:", paste0(etype, "_Plane.pdf"), "\n")
}

# ---- 11) Create one combined 4-in-1 plot with a single shared legend ----

# Helper function to extract legend
get_legend <- function(a_plot) {
  tmp <- ggplot_gtable(ggplot_build(a_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

# Extract shared legend from first plot
shared_legend <- get_legend(plots_list[[1]])

# Remove legends from all individual plots
plots_no_legend <- lapply(plots_list, function(p) p + theme(legend.position = "none"))

# Arrange 4 plots (2×2 grid) + single legend at bottom
combined_plot <- grid.arrange(
  arrangeGrob(grobs = plots_no_legend, ncol = 2),
  shared_legend,
  ncol = 1,
  heights = c(10, 1),
  top = "Entropy–Complexity Planes"
)

# Show preview
print(combined_plot)

# Save combined plot
ggsave(filename = file.path(output_dir, "All_Entropy_Comparison_case3.pdf"),
       plot = combined_plot, width = 12, height = 9)

cat("\n✅ Saved: All_Entropy_Comparison.pdf with single legend\n")


# =============================
# =============================
# 📊 Emblematic (Median) H* x C* Points for Time Series Models
# --- Load libraries ---
library(here)
library(readxl)
library(openxlsx)
library(writexl)
library(dplyr)
library(ggplot2)
library(gridExtra)

# --- 2. Define file paths ---
data_path <- here("GitHub", "Ordinal_Patterns_R", "Data", "ARMA Time series results", "Summary results.xlsx")

#plot_dir  <- here("Plots", "ARMA plots")
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Case 3 results/Emblematic results"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- 3. Load data ---
sheet_index <- 4
df <- read_excel(data_path, sheet = sheet_index)
cat("✅ Loaded sheet", sheet_index, "with columns:\n")
print(colnames(df))

# --- 4. Function to compute emblematic (median) H*, C* points ---
get_CH_emblematic <- function(df, h_col, c_col, model_col = "Model", entropy_type = "Shannon") {
  df <- df %>%
    mutate(
      !!h_col := as.numeric(.data[[h_col]]),
      !!c_col := as.numeric(.data[[c_col]])
    ) %>%
    group_by(.data[[model_col]]) %>%
    group_modify(~{
      dat <- select(.x, all_of(h_col), all_of(c_col))
      pcs <- prcomp(dat, scale. = FALSE)
      pcscores <- pcs$x[, 1]
      N <- nrow(dat)
      median_index <- order(pcscores)[ceiling((N + 1) / 2)]
      em <- dat[median_index, , drop = FALSE]
      tibble(
        Model = .x[[model_col]][1],
        Entropy_Type = entropy_type,
        H_star = em[[h_col]],
        C_star = em[[c_col]],
        Emblematic_Row = median_index
      )
    }) %>%
    ungroup()
  return(df)
}

# --- 5. Compute emblematic rows for each entropy ---
em_shannon <- get_CH_emblematic(df, "H_Shannon", "C_Shannon", entropy_type = "Shannon")
em_renyi   <- get_CH_emblematic(df, "H_Renyi", "C_Renyi", entropy_type = "Renyi")
em_tsallis <- get_CH_emblematic(df, "H_Tsallis", "C_Tsallis", entropy_type = "Tsallis")
em_fisher  <- get_CH_emblematic(df, "H_Fisher", "C_Fisher", entropy_type = "Fisher")

# --- 6. Combine results ---
emblematic_all <- bind_rows(em_shannon, em_renyi, em_tsallis, em_fisher)
cat("\n🔹 Emblematic (median) results with raw indices:\n")
print(emblematic_all)

# --- 7. Save emblematic results ---
# Add as new sheet to Excel
wb <- loadWorkbook(data_path)
if ("Emblematic_Rows" %in% names(wb)) removeWorksheet(wb, "Emblematic_Rows")
addWorksheet(wb, "Emblematic_Rows")
writeData(wb, "Emblematic_Rows", emblematic_all)
saveWorkbook(wb, data_path, overwrite = TRUE)

# Also save CSV
output_csv <- file.path(plot_dir, "Emblematic_All_Raw_case3.csv")
write.csv(emblematic_all, output_csv, row.names = FALSE)
cat("\n✅ Saved emblematic points to Excel and CSV.\n")

# --- 8. Define colors & shapes for all 6 models ---
model_colors <-  c(
  "AR2_M2" = "blue",
  "AR2_M3" = "red",
  "MA2_M2" = "purple",
  "MA2_M3" = "#b15928",
  "ARMA22_M2" = "green",
  "ARMA22_M3" = "black"
)

model_shapes <- c(
  "AR2_M2" = 16, "AR2_M3" = 1,
  "MA2_M2" = 17, "MA2_M3" = 2,
  "ARMA22_M2" = 18, "ARMA22_M3" = 5
)

# --- 9. Individual plots per entropy type ---
plots_list <- list()
unique_entropies <- unique(emblematic_all$Entropy_Type)

for (etype in unique_entropies) {
  df_plot <- filter(emblematic_all, Entropy_Type == etype)
  
  p_sep <- ggplot(df_plot, aes(x = H_star, y = C_star, color = Model, shape = Model)) +
    geom_point(size = 1.5, stroke = 1.2) +
    geom_text(aes(label = Model), vjust = -0.6, size = 2, show.legend = FALSE) +
    scale_shape_manual(values = model_shapes, drop = FALSE) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    labs(
      title = paste( etype),
      x = expression(H^"*"), y = expression(C^"*")
    ) +
    theme_minimal(base_size = 12, base_family = "serif") +
    theme(legend.position = "bottom")
  
  cat(paste0("\n👀 Preview for ", etype, " entropy plot...\n"))
  print(p_sep)
  
  ggsave(
    filename = file.path(plot_dir, paste0("Median_HC_", etype, "_Entropy_case3.pdf")),
    plot = p_sep, width = 8, height = 6
  )
  
  plots_list[[etype]] <- p_sep
}

# --- 10. Create 4-in-1 plot with shared legend ---
get_legend <- function(myplot) {
  tmp <- ggplot_gtable(ggplot_build(myplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

shared_legend <- get_legend(plots_list[[1]])
plots_no_leg <- lapply(plots_list, function(p) p + theme(legend.position = "none"))

p_grid <- grid.arrange(
  arrangeGrob(grobs = plots_no_leg, ncol = 2),
  shared_legend,
  ncol = 1,
  heights = c(10, 1),
  top = "Entropy–Complexity Median Comparison:Case 3"
)

cat("\n👀 Showing 4-in-1 plot preview...\n")
print(p_grid)

ggsave(
  filename = file.path(plot_dir, "Median_Comparison_4in1_case3.pdf"),
  plot = p_grid, width = 12, height = 10
)

cat("\n✅ All plots and emblematic results saved in:\n", plot_dir, "\n")
# =============================

####################################################################################
#Time series plot and related emblematic point in HC plane
####################################################################################
# ================================================================
# 📊 Entropy–Complexity Plane with Time-Series Insets
# ================================================================

library(readxl)
library(ggplot2)
library(dplyr)
library(here)

# --- Paths ---
data_path <- here("GitHub", "Ordinal_Patterns_R", "Data", "ARMA Time series results", "All cases time series.xlsx")
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Case 3 results/Time series and median points"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- Read data ---
hc_med <- read_excel(data_path, sheet = 5)
ts_data <- read_excel(data_path, sheet = 6)

cat("✅ Loaded HC median points & time series sheets.\n")
print(colnames(ts_data))  # check structure

# --- Define models and colors ---
models_to_plot <- c("AR2_M2", "AR2_M3", 
                    "MA2_M2", "MA2_M3", 
                    "ARMA22_M2", "ARMA22_M3") 

manual_colors <- c(
  "AR2_M2" = "blue",
  "AR2_M3" = "red",
  "MA2_M2" = "purple",
  "MA2_M3" = "#b15928",
  "ARMA22_M2" = "green",
  "ARMA22_M3" = "black"
)

# --- Determine correct column names dynamically ---
# Try to find suitable column names automatically
col_lower <- tolower(colnames(ts_data))
time_col <- colnames(ts_data)[grepl("time", col_lower)][1]
value_col <- colnames(ts_data)[grepl("value|series|data|y", col_lower)][1]

if (is.na(time_col) | is.na(value_col)) {
  stop("❌ Could not automatically detect 'Time' or 'Value' column. Please rename columns in Excel (e.g., 'Time', 'Value').")
}

# --- Plot time series per model ---
for (mdl in models_to_plot) {
  df_model <- ts_data %>% filter(Model == mdl)
  
  if (nrow(df_model) == 0) {
    warning(paste("Skipping model", mdl, "- no data found"))
    next
  }
  
  p <- ggplot(df_model, aes(x = .data[[time_col]], y = .data[[value_col]])) +
    geom_line(color = manual_colors[mdl], size = 0.9) +
    labs(title = paste("Time Series -", mdl), x = "Time", y = "Value") +
    theme_minimal(base_size = 12, base_family = "serif")
  
  print(p)  # 👀 preview before saving
  ggsave(filename = file.path(plot_dir, paste0("TS_plot_", mdl, ".pdf")),
         plot = p, width = 8, height = 4)
}

# --- Plot Entropy–Complexity Median Points ---
p_hc <- ggplot(hc_med, aes(x = H_star, y = C_star, color = Model, label = Model)) + 
  geom_point(size = 4) + 
  geom_text(nudge_y = 0.02, fontface = "bold", show.legend = FALSE) + 
  scale_color_manual(values = manual_colors) + 
  labs(title = "Entropy–Complexity Plane (Median Points)", 
       x = expression(H^"*"), y = expression(C^"*")) + 
  theme_minimal(base_size = 12, base_family = "serif") + 
  theme(legend.position = "bottom")

print(p_hc)
ggsave(filename = file.path(plot_dir, "HC_plane_median_points_case3.pdf"), 
       plot = p_hc, width = 8, height = 5)


# ================================================================
# =========================================================
# ===========================================================
# 📊 Combined Entropy–Complexity Plane + Time Series Insets
# ===========================================================
# ===========================================================
# 📊 Entropy–Complexity Plane (Center) + Time Series Insets (Top/Bottom)
# ===========================================================
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(here)

# --- 1️⃣ Paths ---
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Plots/ARMA plots/Case 3 results/Time series and median points"
data_path <- here("GitHub", "Ordinal_PatternS_R", "Data", "ARMA Time series results", "All cases time series.xlsx")

if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- 2️⃣ Load Data ---
hc_med <- read_excel(data_path, sheet = 5)
ts_data <- read_excel(data_path, sheet = 6)

cat("✅ Loaded Excel sheets successfully.\n")
print(colnames(ts_data))  # 🧐 check structure

# --- 3️⃣ Automatically detect time/value columns ---
col_lower <- tolower(colnames(ts_data))
time_col <- colnames(ts_data)[grepl("time", col_lower)][1]
value_col <- colnames(ts_data)[grepl("value|series|data|y", col_lower)][1]

if (is.na(time_col) | is.na(value_col)) {
  stop("❌ Could not detect 'Time' or 'Value' column in your time series sheet. Please rename them clearly (e.g., 'Time', 'Value').")
}

# --- 4️⃣ Define models and colors ---
models_to_plot <- c("AR2_M2", "AR2_M3", 
                    "MA2_M2", "MA2_M3", 
                    "ARMA22_M2", "ARMA22_M3")

manual_colors <- c(
  "AR2_M2" = "blue",
  "AR2_M3" = "red",
  "MA2_M2" = "purple",
  "MA2_M3" = "#b15928",
  "ARMA22_M2" = "green",
  "ARMA22_M3" = "black"
)

# --- 5️⃣ Create small time-series plots ---
ts_plots <- lapply(models_to_plot, function(mdl) {
  df <- ts_data %>% filter(Model == mdl)
  if (nrow(df) == 0) {
    warning(paste("⚠️ Skipping", mdl, "- no data found"))
    return(ggplot() + theme_void() + ggtitle(paste(mdl, "(No Data)")))
  }
  ggplot(df, aes(x = .data[[time_col]], y = .data[[value_col]])) +
    geom_line(color = manual_colors[mdl], size = 0.6) +
    labs(title = mdl) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 9, hjust = 0.5)
    )
})

# --- 6️⃣ Central H–C plane (main plot, no legend) ---
p_center <- ggplot(hc_med, aes(x = H_star, y = C_star, color = Model, label = Model)) +
  geom_point(size = 5) +
  geom_text(nudge_y = 0.02, fontface = "bold", size = 4, show.legend = FALSE) +
  scale_color_manual(values = manual_colors) +
  labs(title = "H–C Plane (Median Points)",
       x = expression(H^"*"), y = expression(C^"*")) +
  theme_minimal(base_size = 14, base_family = "serif") +
  theme(legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"))

# --- 7️⃣ Arrange layout: 3 top, 1 middle (main), 3 bottom ---
# Top 3 → models 1–3
# Bottom 3 → models 4–6
top_row    <- arrangeGrob(grobs = ts_plots[1:3], ncol = 3)
bottom_row <- arrangeGrob(grobs = ts_plots[4:6], ncol = 3)

combined <- grid.arrange(
  top_row,
  p_center,
  bottom_row,
  ncol = 1,
  heights = c(2, 3.5, 2),
  top = textGrob("Median Points in the H–C Plane and Associated Time Series (Case 3)",
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# --- 8️⃣ Save output ---
output_file <- file.path(plot_dir, "Combined_HC_TS_center_case3.pdf")
ggsave(filename = output_file, plot = combined, width = 12, height = 12)

cat("\n✅ Saved combined figure to:", output_file, "\n")

# =============================