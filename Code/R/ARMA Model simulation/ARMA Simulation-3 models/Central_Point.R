# 📊 Central (Median) H* × C* Points for Shannon Entropy (3 Cases)
# ================================================================
# --- 1. Load libraries ---
library(here)
library(readxl)
library(openxlsx)
library(writexl)
library(dplyr)
library(ggplot2)
library(stats)

# --- 2. Define file paths ---
data_path <- data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/Summary results.xlsx"

# --- 3. Define output directory ---
plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Emblematic results"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# --- 4. Function: Compute central (median) Shannon H*, C* point with SemiLength ---
get_CH_emblematic <- function(df, h_col, c_col, var_col, model_col = "Model", entropy_type = "Shannon") {
  df <- df %>%
    mutate(
      !!h_col := as.numeric(.data[[h_col]]),
      !!c_col := as.numeric(.data[[c_col]]),
      !!var_col := as.numeric(.data[[var_col]])
    ) %>%
    group_by(.data[[model_col]]) %>%
    group_modify(~{
      dat <- select(.x, all_of(h_col), all_of(c_col))
      pcs <- prcomp(dat, scale. = FALSE)
      pcscores <- pcs$x[, 1]
      N <- nrow(dat)
      median_index <- order(pcscores)[ceiling((N + 1) / 2)]
      em <- dat[median_index, , drop = FALSE]
      var_val <- .x[[var_col]][median_index]
      
      # --- Compute SemiLength using given formula ---
      SemiLength <- sqrt(var_val) / sqrt(1000 - 3) * qnorm(1 - 0.05/2)
      
      tibble(
        Model = .x[[model_col]][1],
        Entropy_Type = entropy_type,
        H_Shannon = em[[h_col]],
        C_Shannon = em[[c_col]],
        Var_Shannon = var_val,
        SemiLength = SemiLength,
        Emblematic_Row = median_index
      )
    }) %>%
    ungroup()
  
  return(df)
}

# --- 5. Loop through 3 cases (Sheets 2–4) ---
case_sheets <- list("Case 1" = 2, "Case 2" = 3, "Case 3" = 4)

wb <- loadWorkbook(data_path)

for (case_name in names(case_sheets)) {
  sheet_index <- case_sheets[[case_name]]
  
  # --- Load data ---
  df <- read_excel(data_path, sheet = sheet_index)
  cat("\n✅ Loaded", case_name, "(Sheet", sheet_index, ") with columns:\n")
  print(colnames(df))
  
  # --- Compute emblematic Shannon results ---
  em_shannon <- get_CH_emblematic(
    df,
    h_col = "H_Shannon",
    c_col = "C_Shannon",
    var_col = "Var_Shannon",
    entropy_type = "Shannon"
  )
  
  # --- Add central point (median H*, C*) ---
  central_point <- em_shannon %>%
    summarise(Central_H = median(H_Shannon, na.rm = TRUE),
              Central_C = median(C_Shannon, na.rm = TRUE))
  
  em_shannon <- em_shannon %>%
    mutate(Central_H = central_point$Central_H,
           Central_C = central_point$Central_C)
  
  # --- Write to Excel as a separate sheet ---
  sheet_name <- paste0("Emblematic_", case_name)
  if (sheet_name %in% names(wb)) removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, em_shannon)
  
  # --- Also save as CSV ---
  output_csv <- file.path(plot_dir, paste0("Emblematic_", case_name, ".csv"))
  write.csv(em_shannon, output_csv, row.names = FALSE)
  
  cat("✅ Saved emblematic Shannon results for", case_name, "to Excel and CSV.\n")
}

# --- 6. Save workbook ---
saveWorkbook(wb, data_path, overwrite = TRUE)
cat("\n🎯 All three emblematic result sheets successfully saved in 'Summary results.xlsx'.\n")


##########################################################################################################

###############################################################################################################
############################################################################################################

# =======================================================
# 📊 HC Central Points + CI + Heatmaps (n = 1000 only)
# =======================================================

# --- Load libraries ---
library(ggplot2)
library(dplyr)
library(readr)
library(viridis)
library(StatOrdPattHxC)

# --- Define consistent color and shape palettes ---
model_colors <- c(
  "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
  "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
  "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
  "ARMA22_M3"="#fb9a99",
  "AR1_M1"="#8dd3c7","AR1_M2"="#ffffb3","ARMA11_M2"="#bebada","MA1_M1"="#fb8072",
  "MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69","AMA11_M4"="#fccde5",
  "MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
)

model_shapes <- c(
  "ARMA11_M1"=21,"AR2_M1"=22,"MA1_M2"=23,"ARMA11_M3"=24,"AR2_M4"=25,"MA1_M4"=8,
  "ARMA22_M2"=15,"AR2_M3"=16,"MA2_M3"=17,"ARMA22_M1"=18,"ARMA22_M4"=19,"MA2_M2"=4,
  "ARMA22_M3"=3,
  "AR1_M1"=1,"AR1_M2"=2,"ARMA11_M2"=5,"MA1_M1"=6,"MA2_M1"=7,"AR1_M3"=9,
  "AR1_M4"=10,"AMA11_M4"=11,"MA1_M3"=12,"MA2_M4"=13,"AR2_M2"=14
)

# --- Load LinfLsup boundaries ---
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# --- Define file directories ---
base_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results"
case_names <- c("Case 1", "Case 2", "Case 3")

# --- Loop through cases ---
for (case_name in case_names) {
  
  message("Processing ", case_name, " ...")
  
  # File and output directories
  data_path <- file.path(base_dir, case_name, paste0("Emblematic_", case_name, ".csv"))
  plot_dir  <- file.path(base_dir, case_name, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # --- Read emblematic data ---
  df <- read_csv(data_path, show_col_types = FALSE)
  
  # Ensure column names match expectations
  df <- df %>%
    rename(H_data = H_star, C_data = C_star) %>%
    mutate(
      Model = factor(Model),
      n = 1000   # fixed sample size
    )
  
  # Define focus region for boundaries
  h_range <- range(df$H_data, na.rm = TRUE)
  c_range <- range(df$C_data, na.rm = TRUE)
  Linf_focus <- LinfLsup_subset %>%
    filter(H >= min(h_range) - 0.05, H <= max(h_range) + 0.05,
           C >= min(c_range) - 0.05, C <= max(c_range) + 0.05)
  
  # ===================================================
  # 📍 Scatter Plot with Confidence Intervals
  # ===================================================
  p_scatter <- ggplot() +
    geom_line(data=subset(Linf_focus, Side=="Lower"), aes(H, C),
              linetype="dashed", color="gray40", linewidth=0.6) +
    geom_line(data=subset(Linf_focus, Side=="Upper"), aes(H, C),
              linetype="dashed", color="gray40", linewidth=0.6) +
    geom_errorbarh(
      data=df, aes(y=C_data, xmin=H_data - SemiLength, xmax=H_data + SemiLength, color=Model),
      height=0.002, linewidth=0.7
    ) +
    geom_point(
      data=df, aes(H_data, C_data, fill=Model, shape=Model),
      size=4, color="black"
    ) +
    geom_text(
      data=df, aes(H_data, C_data, label=Model),
      vjust=-1, fontface="bold", size=3
    ) +
    scale_fill_manual(values=model_colors) +
    scale_color_manual(values=model_colors) +
    scale_shape_manual(values=model_shapes) +
    labs(
      title=paste0("Entropy–Complexity Plane (", case_name, ", n=1000)"),
      x=expression(italic(H)), y=expression(italic(C))
    ) +
    theme_minimal(base_size=13, base_family="serif") +
    theme(legend.position="bottom", legend.title=element_blank())
  
  ggsave(file.path(plot_dir, paste0("HC_Scatter_CI_", case_name, "_n1000.pdf")),
         p_scatter, width=8, height=6)
  print(p_scatter)
  
  # ===================================================
  # 🔥 Heatmap of Central Points
  # ===================================================
  p_heat <- ggplot() +
    stat_density_2d(
      data=df, aes(x=H_data, y=C_data, fill=after_stat(level)),
      geom="polygon", bins=60, alpha=0.6
    ) +
    geom_point(data=df, aes(H_data, C_data, shape=Model, color=Model), size=3) +
    scale_fill_viridis_c(option="plasma") +
    scale_color_manual(values=model_colors) +
    scale_shape_manual(values=model_shapes) +
    labs(
      title=paste("Heatmap of Central Points (", case_name, ", n=1000)", sep=""),
      x=expression(italic(H)), y=expression(italic(C))
    ) +
    theme_minimal(base_size=13, base_family="serif") +
    theme(legend.position="bottom", legend.title=element_blank())
  
  ggsave(file.path(plot_dir, paste0("HC_Heatmap_", case_name, "_n1000.pdf")),
         p_heat, width=8, height=6)
  print(p_heat)
}

message("✅ All plots generated successfully for n = 1000.")
###############################################################################################################

