########################################################################################
# In this script, we simulate time series data from various ARMA models, compute their
# entropy and complexity measures, and save the results in Excel files.
#######################################################################
# --- Required Packages ---
library(StatOrdPattHxC)
library(dplyr)
library(writexl)

# --- Parameters ---
set.seed(1234567890, kind = "Mersenne-Twister")
D <- 3                     # Embedding dimension
#N <- c(500, 1000)          # Sample sizes
N <- c(5000, 10000)
R <- 100                   # Number of replications

# --- Model Definitions ---
ar1_models <- list(
  AR1_M1 = list(ar = c(0.8), type = "AR1_M1"),
  AR1_M2 = list(ar = c(0.1), type = "AR1_M2"),
  AR1_M3 = list(ar = c(-0.8), type = "AR1_M3"),
  AR1_M4 = list(ar = c(-0.1), type = "AR1_M4")
)
ar2_models <- list(
  AR2_M1 = list(ar = c(0.1, 0.8), type = "AR2_M1"),
  AR2_M2 = list(ar = c(-0.8, 0.1), type = "AR2_M2"),
  AR2_M3 = list(ar = c(0.1, -0.8), type = "AR2_M3"),
  AR2_M4 = list(ar = c(-0.8, -0.1), type = "AR2_M4")
)
ma1_models <- list(
  MA1_M1 = list(ma = c(0.8), type = "MA1_M1"),
  MA1_M2 = list(ma = c(0.1), type = "MA1_M2"),
  MA1_M3 = list(ma = c(-0.8), type = "MA1_M3"),
  MA1_M4 = list(ma = c(-0.1), type = "MA1_M4")
)
ma2_models <- list(
  MA2_M1 = list(ma = c(0.1, 0.8), type = "MA2_M1"),
  MA2_M2 = list(ma = c(-0.8, 0.1), type = "MA2_M2"),
  MA2_M3 = list(ma = c(0.1, -0.8), type = "MA2_M3"),
  MA2_M4 = list(ma = c(-0.8, -0.1), type = "MA2_M4")
)
arma11_models <- list(
  ARMA11_M1 = list(ar = c(0.8), ma = c(0.8), type = "ARMA11_M1"),
  ARMA11_M2 = list(ar = c(0.1), ma = c(0.1), type = "ARMA11_M2"),
  ARMA11_M3 = list(ar = c(-0.8), ma = c(-0.8), type = "ARMA11_M3"),
  ARMA11_M4 = list(ar = c(-0.1), ma = c(-0.1), type = "ARMA11_M4")
)
arma22_models <- list(
  ARMA22_M1 = list(ar = c(0.1, 0.8), ma = c(0.1, 0.8), type = "ARMA22_M1"),
  ARMA22_M2 = list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), type = "ARMA22_M2"),
  ARMA22_M3 = list(ar = c(0.1, -0.8), ma = c(0.1, -0.8), type = "ARMA22_M3"),
  ARMA22_M4 = list(ar = c(-0.8, -0.1), ma = c(-0.8, -0.1), type = "ARMA22_M4")
)

# --- Combine all models ---
all_models <- c(ar1_models, ar2_models, ma1_models, ma2_models, arma11_models, arma22_models)

# --- Function to simulate data and compute entropy/complexity ---
generate_model_data <- function(model, n, D, R) {
  results <- list()
  ts_store <- list()
  
  for (r in 1:R) {
    ts_data <- arima.sim(model = model, n = n)
    ProbTS <- OPprob(ts_data, emb = D)
    
    # Entropy and complexity
    Hs <- HShannon(ProbTS)
    Cs <- StatComplexity(ProbTS)
    
    # Variances (keep negative)
    Var_HD <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(ProbTS, n - 2))
    a_ratio <- Var_HD / Var_HI
    Var_CI <- suppressWarnings(varC(ProbTS, n - 2))
    Var_CD <- a_ratio * Var_CI
    
    # Semi-lengths: NA if variance <= 0
    SemiLengthH <- ifelse(!is.finite(Var_HD) | Var_HD <= 0, NA,
                          sqrt(Var_HD) / sqrt(n - 3) * qnorm(1 - 0.05 / 2))
    SemiLengthC <- ifelse(!is.finite(Var_CD) | Var_CD <= 0, NA,
                          sqrt(Var_CD) / sqrt(n - 3) * qnorm(1 - 0.05 / 2))
    
    results[[r]] <- data.frame(
      Model = model$type,
      n = n,
      Rep = r,
      H_Shannon = Hs,
      C_Shannon = Cs,
      Var_H = Var_HD,
      Var_C = Var_CD,
      SemiLength_H = SemiLengthH,
      SemiLength_C = SemiLengthC
    )
    
    ts_store[[r]] <- data.frame(
      Model = model$type,
      n = n,
      Rep = r,
      Value = as.numeric(ts_data)
    )
  }
  
  list(
    summary = bind_rows(results),
    timeseries = bind_rows(ts_store)
  )
}

# --- Run simulations for all models and sample sizes ---
all_summary <- list()
all_timeseries <- list()

for (n_val in N) {
  message("Simulating for n = ", n_val)
  for (model in all_models) {
    res <- generate_model_data(model, n = n_val, D = D, R = R)
    # Save each model and n combination as separate sheet
    sheet_name <- paste0(model$type, "_n", n_val)
    all_summary[[sheet_name]] <- res$summary
    all_timeseries[[sheet_name]] <- res$timeseries
  }
}

# --- Output paths ---
summary_path <- "Entropy_Complexity_Results_all_Models.xlsx"
ts_path <- "TimeSeries_Data_all_Models.xlsx"

# --- Save Excel files ---
write_xlsx(all_summary, path = summary_path)
write_xlsx(all_timeseries, path = ts_path)

cat("✅ All simulations complete!\n")
cat("Entropy–Complexity for all model results saved at:", summary_path, "\n")
cat("Time series data for all models saved at:", ts_path, "\n")

########################################################################################
# In this script, we compute the emblematic (central) Shannon entropy and complexity points
# for various ARMA models from precomputed results stored in an Excel file. The results
# are saved back to the same Excel file in new sheets.
########################################################################################
# --- Load Required Packages ---
library(readxl)
library(openxlsx)
library(dplyr)
library(stats)

# --- Define file paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
output_path <- data_path  # same file

# --- Load sheet names ---
all_sheets <- excel_sheets(data_path)

# --- Helper function: calculate emblematic point for each model ---
get_emblematic_point <- function(df, model_name, n_value) {
  # ensure numeric
  df <- df %>%
    mutate(across(c(H_Shannon, C_Shannon, Var_H, Var_C), as.numeric))
  
  # PCA-based central point
  dat <- df[, c("H_Shannon", "C_Shannon")]
  pcs <- prcomp(dat, scale. = FALSE)
  pcscores <- pcs$x[, 1]
  N <- nrow(dat)
  median_index <- order(pcscores)[ceiling((N + 1) / 2)]
  
  em <- dat[median_index, , drop = FALSE]
  varH <- df$Var_H[median_index]
  varC <- df$Var_C[median_index]
  
  # --- SemiLength: keep NA if variance < 0 ---
  SemiLength_H <- ifelse(varH >= 0, sqrt(varH) / sqrt(n_value - 3) * qnorm(1 - 0.05/2), NA)
  SemiLength_C <- ifelse(varC >= 0, sqrt(varC) / sqrt(n_value - 3) * qnorm(1 - 0.05/2), NA)
  
  tibble(
    Model = model_name,
    n = n_value,
    Emblematic_Row = median_index,
    H_Star = em$H_Shannon,
    C_Star = em$C_Shannon,
    Var_H = varH,
    Var_C = varC,
    SemiLength_H = SemiLength_H,
    SemiLength_C = SemiLength_C
  )
}

# --- Storage lists for n=5000 and n=10000 ---
results_1000 <- list()
results_5000 <- list()

# --- Loop through all sheets ---
for (sheet in all_sheets) {
  df <- read_excel(data_path, sheet = sheet)
  
  if (!all(c("H_Shannon", "C_Shannon", "Var_H", "Var_C") %in% names(df))) {
    message("⚠️ Skipping ", sheet, " — missing columns.")
    next
  }
  
  # Extract model name and n value
  parts <- strsplit(sheet, "_")[[1]]
  model_name <- paste(parts[1], parts[2], sep = "_")
  n_value <- as.numeric(gsub("n", "", parts[length(parts)]))
  
  message("✅ Processing ", sheet, " (n = ", n_value, ") ...")
  
  # Compute emblematic point
  em <- get_emblematic_point(df, model_name, n_value)
  
  # Store results
  if (n_value == 1000) results_1000[[model_name]] <- em
  if (n_value == 5000) results_5000[[model_name]] <- em
}

# --- Combine results ---
df_1000 <- bind_rows(results_1000)
df_5000 <- bind_rows(results_5000)

if (nrow(df_1000) == 0 & nrow(df_5000) == 0) {
  cat("\n⚠️ No results were generated. Check sheet naming or n values.\n")
} else {
  # --- Save to same Excel file ---
  wb <- loadWorkbook(output_path)
  
  if ("CentralPoints_n1000" %in% names(wb)) removeWorksheet(wb, "CentralPoints_n1000")
  if ("CentralPoints_n5000" %in% names(wb)) removeWorksheet(wb, "CentralPoints_n5000")
  
  addWorksheet(wb, "CentralPoints_n1000")
  addWorksheet(wb, "CentralPoints_n5000")
  
  writeData(wb, "CentralPoints_n1000", df_1000)
  writeData(wb, "CentralPoints_n5000", df_5000)
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  
  cat("\n🎯 Central point results saved to the same Excel file:\n", output_path, "\n")
}

##########################################################################################
# In this script, we generate scatter plots of entropy vs. complexity for various ARMA models.
# The data is read from an Excel file, and the plots are saved in a structured directory
# based on cases and sample sizes.
####################################################################################
# --- Required Packages ---
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(StatOrdPattHxC)

# --- File Paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"

# --- Load LinfLsup boundaries ---
data("LinfLsup")

# --- Embedding Dimension ---
D <- 3  # You can adjust this if needed

# --- Cases ---
cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1","MA1_M1","MA1_M2","MA2_M1","ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4","MA1_M3","MA1_M4","MA2_M4","ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3","MA2_M2","MA2_M3","ARMA22_M2","ARMA22_M3")
)

sample_sizes <- c(1000, 5000)

# --- Define Colors and Shapes ---
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

# --- Main Loop ---
sheet_names <- excel_sheets(data_path)

for(case_name in names(cases)) {
  case_models <- cases[[case_name]]
  case_plot_dir <- file.path(base_plot_dir, case_name, "Plots")
  dir.create(case_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(n_val in sample_sizes) {
    
    # Read all sheets for this case and sample size
    all_data <- lapply(case_models, function(model_name) {
      sheet_n <- paste0(model_name, "_n", n_val)
      if(sheet_n %in% sheet_names) {
        df <- read_excel(data_path, sheet = sheet_n)
        df$Model <- model_name
        df$n <- n_val
        return(df)
      } else {
        warning(paste("Skipping", sheet_n, "- sheet not found"))
        return(NULL)
      }
    }) %>% bind_rows()
    
    if(nrow(all_data) == 0) next
    
    # --- Focus Linf/Lsup Boundaries for the Data Range ---
    LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))
    h_range <- range(all_data$H_Shannon, na.rm = TRUE)
    c_range <- range(all_data$C_Shannon, na.rm = TRUE)
    
    Linf_focus <- LinfLsup_subset %>%
      filter(H >= min(h_range) - 0.05, H <= max(h_range) + 0.05,
             C >= min(c_range) - 0.05, C <= max(c_range) + 0.05)
    
    # --- Scatter Plot ---
    p_scatter <- ggplot() +
      geom_line(data=subset(Linf_focus, Side=="Lower"), aes(H, C),
                linetype="dashed", color="gray40", linewidth=0.6) +
      geom_line(data=subset(Linf_focus, Side=="Upper"), aes(H, C),
                linetype="dashed", color="gray40", linewidth=0.6) +
      geom_point(data=all_data, aes(H_Shannon, C_Shannon, color=Model, shape=Model), size=3) +
      scale_color_manual(values=model_colors) +
      scale_shape_manual(values=model_shapes) +
      labs(title=paste0("Entropy–Complexity Scatter with Boundaries (", case_name, ", n=", n_val, ")"),
           x=expression(italic(H)), y=expression(italic(C))) +
      theme_minimal(base_family="serif", base_size=13) +
      theme(
        legend.position="bottom",
        legend.title=element_blank(),
        panel.grid.minor=element_blank()
      )
    
    # --- Save Plot ---
    ggsave(file.path(case_plot_dir, paste0("Scatter_", case_name, "_n", n_val, ".pdf")),
           p_scatter, width=8, height=6)
    
    message(paste("✅ Scatter plot generated for", case_name, "n=", n_val))
  }
}

message("🎉 All scatter plots completed with Linf–Lsup boundaries!")

######################################################################################
# In this script, we generate scatter plots with confidence intervals 
# for the central points of various ARMA models. The central points are read from an
# Excel file, and the plots are saved in a structured directory based on cases.
######################################################################################
# --- Required Libraries ---
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(StatOrdPattHxC)   

# --- File Paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"

# --- Define Cases ---
cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1","MA1_M1","MA1_M2","MA2_M1","ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4","MA1_M3","MA1_M4","MA2_M4","ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3","MA2_M2","MA2_M3","ARMA22_M2","ARMA22_M3")
)

# --- Color and Shape Palettes ---
model_colors <- c(
  "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
  "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
  "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
  "ARMA22_M3"="#fb9a99",
  "AR1_M1"="#8dd3c7","AR1_M2"="red","ARMA11_M2"="#bebada","MA1_M1"="#fb8072",
  "MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69","ARMA11_M4"="#fccde5",
  "MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
)

model_shapes <- c(
  "ARMA11_M1"=21,"AR2_M1"=22,"MA1_M2"=23,"ARMA11_M3"=24,"AR2_M4"=25,"MA1_M4"=8,
  "ARMA22_M2"=15,"AR2_M3"=16,"MA2_M3"=17,"ARMA22_M1"=18,"ARMA22_M4"=19,"MA2_M2"=4,
  "ARMA22_M3"=3,
  "AR1_M1"=1,"AR1_M2"=2,"ARMA11_M2"=5,"MA1_M1"=6,"MA2_M1"=7,"AR1_M3"=9,
  "AR1_M4"=10,"ARMA11_M4"=11,"MA1_M3"=12,"MA2_M4"=13,"AR2_M2"=14
)

# --- Load LinfLsup boundaries and choose embedding dimension ---
data("LinfLsup")
D <- 3   # same as in your second script
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# --- Function to Read and Clean Sheet ---
read_central <- function(sheet_name){
  df <- read_excel(data_path, sheet = sheet_name)
  names(df) <- trimws(names(df))                 # remove leading/trailing spaces
  df <- df %>% rename_all(~gsub("\\s+","",.))    # remove internal spaces
  return(df)
}

# --- Loop over Sheets and Cases ---
sheets <- c("CentralPoints_n1000","CentralPoints_n5000")

for(sheet_name in sheets){
  
  # infer n from sheet name
  n_val <- ifelse(grepl("1000", sheet_name), 1000, 5000)
  df <- read_central(sheet_name)
  
  for(case_name in names(cases)){
    case_models <- cases[[case_name]]
    
    # Create Plot Folder
    plot_dir <- file.path(base_plot_dir, case_name, "Plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Filter Case Models
    df_case <- df %>% filter(Model %in% case_models)
    if(nrow(df_case) == 0) next
    
    # --- Subset Linf/Lsup to the H*, C* range for this case ---
    h_range <- range(df_case$H_Star, na.rm = TRUE)
    c_range <- range(df_case$C_Star, na.rm = TRUE)
    
    Linf_focus <- LinfLsup_subset %>%
      filter(H >= min(h_range) - 0.05, H <= max(h_range) + 0.05,
             C >= min(c_range) - 0.05, C <= max(c_range) + 0.05)
    
    # --- Central Points with CI: keep only rows with valid CI lengths ---
    df_ci <- df_case %>%
      filter(!is.na(SemiLength_H) & SemiLength_H > 0 &
               !is.na(SemiLength_C) & SemiLength_C > 0)
    
    # --- Plot with Linf/Lsup boundaries + central points + CI ---
    p_central <- ggplot() +
      # boundaries (same style as in your scatter code)
      geom_line(data = subset(Linf_focus, Side == "Lower"),
                aes(H, C), linetype = "dashed",
                color = "gray40", linewidth = 0.6) +
      geom_line(data = subset(Linf_focus, Side == "Upper"),
                aes(H, C), linetype = "dashed",
                color = "gray40", linewidth = 0.6) +
      # central points
      geom_point(data = df_case,
                 aes(H_Star, C_Star, color = Model, shape = Model),
                 size = 4) +
      # --- NEW: Label each model next to its point ---
      geom_text_repel(
        data = df_case,
        aes(H_Star, C_Star, label = Model, color = Model),
        size = 3.5,
        box.padding = 0.3,
        point.padding = 0.25,
        segment.color = "grey40",
        max.overlaps = Inf
      ) +
      # horizontal CI
      geom_errorbarh(data = df_ci,
                     aes(y = C_Star,
                         xmin = H_Star - SemiLength_H,
                         xmax = H_Star + SemiLength_H,
                         color = Model),
                     height = 0.002, linewidth = 0.8) +
      # vertical CI
      geom_errorbar(data = df_ci,
                    aes(x = H_Star,
                        ymin = C_Star - SemiLength_C,
                        ymax = C_Star + SemiLength_C,
                        color = Model),
                    width = 0.002, linewidth = 0.8) +
      scale_color_manual(values = model_colors) +
      scale_shape_manual(values = model_shapes) +
      labs(title = paste0("Central Points with CI and Boundaries (",
                          case_name, ", n=", n_val, ")"),
           x = expression(italic(H)*"*"),
           y = expression(italic(C)*"*")) +
      theme_minimal(base_family = "serif", base_size = 13) +
      theme(
        legend.position = "bottom",
        legend.title    = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # --- Save Plot ---
    ggsave(file.path(plot_dir,
                     paste0("CentralPoints_",
                            case_name, "_n", n_val, ".pdf")),
           p_central, width = 8, height = 6)
  }
}

message("🎉 All central points plots completed with Linf–Lsup boundaries!")


###############################################################################
#In this plotting script, we generate time series plots for the emblematic points of 
#various ARMA models. We read the central points from an Excel file, extract the emblematic
#row for each model, and plot the corresponding time series data. The plots are saved in
#a structured directory based on cases.  
###############################################################
# 📚 Required Packages
###############################################################
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(StatOrdPattHxC)
library(stringr)
library(ggrepel)


#ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/TimeSeries_by_Case.xlsx"
#base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results"
ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/TimeSeries_by_Case_n1000_n5000.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"



cases <- list(
  Case1 = "Case1",
  Case2 = "Case2",
  Case3 = "Case3"
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



###############################################################
# ⭐ FINAL FUNCTION — ASCENDING ORDER + 3 TOP / 3 RIGHT / 3 BOTTOM
###############################################################
plot_case_timeseries <- function(case_name, n_val) {
  
  message("------------------------------------------------")
  message("Processing ", case_name, " (n = ", n_val, ")")
  
  # Output folder
  plot_dir <- file.path(base_plot_dir, case_name, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  # -------------------------------------------------------------
  # 1. Load Data
  # -------------------------------------------------------------
  df <- read_excel(ts_case_path, sheet = case_name) %>%
    filter(n == n_val) %>%
    mutate(Model = as.character(Model))
  
  if (nrow(df) == 0) stop("❌ No data found for this case and n.")
  
  # Clean TimeSeries
  clean_ts <- function(x) {
    x <- gsub("[c()]", "", x)
    x <- gsub("[[:space:]]+", "", x)
    as.numeric(strsplit(x, ",")[[1]])
  }
  
  ts_list_raw <- lapply(df$TimeSeries, clean_ts)
  names(ts_list_raw) <- df$Model
  
  # -------------------------------------------------------------
  # 2. Sort by H_Star ASCENDING (lowest → highest)
  # -------------------------------------------------------------
  df_sorted <- df %>% arrange(H_Star)
  ts_sorted <- ts_list_raw[df_sorted$Model]
  
  n_ts <- length(ts_sorted)
  message("Detected ", n_ts, " time series.")
  
  # -------------------------------------------------------------
  # 3. Pick TOP, MIDDLE, BOTTOM groups
  # -------------------------------------------------------------
  if (n_ts == 9) {
    top_idx    <- 1:3     # lowest H
    right_idx  <- 4:6     # mid H
    bottom_idx <- 7:9     # highest H
  } else if (n_ts == 6) {
    top_idx    <- 1:3
    right_idx  <- NULL
    bottom_idx <- 4:6
  } else {
    stop("❌ Expected 6 or 9 time series.")
  }
  
  ts_top    <- ts_sorted[top_idx]
  if (!is.null(right_idx)) ts_right <- ts_sorted[right_idx]
  ts_bottom <- ts_sorted[bottom_idx]
  
  # -------------------------------------------------------------
  # 4. Central H–C Scatter Plot (with labels + CI)
  # -------------------------------------------------------------
  data("LinfLsup")
  Linf_focus <- subset(LinfLsup, Dimension == "3") %>%
    filter(H >= min(df$H_Star)-0.05,
           H <= max(df$H_Star)+0.05,
           C >= min(df$C_Star)-0.05,
           C <= max(df$C_Star)+0.05)
  
  p_center <- ggplot() +
    geom_line(data=subset(Linf_focus, Side=="Lower"), aes(H, C),
              color="gray40", linetype="dashed") +
    geom_line(data=subset(Linf_focus, Side=="Upper"), aes(H, C),
              color="gray40", linetype="dashed") +
    
    geom_point(data=df_sorted, aes(H_Star, C_Star, color=Model, shape=Model), size=3) +
    
    geom_errorbarh(data=df_sorted,
                   aes(y=C_Star,
                       xmin=H_Star - SemiLength_H,
                       xmax=H_Star + SemiLength_H,
                       color=Model),
                   height=0.002, linewidth=0.6) +
    
    geom_errorbar(data=df_sorted,
                  aes(x=H_Star,
                      ymin=C_Star - SemiLength_C,
                      ymax=C_Star + SemiLength_C,
                      color=Model),
                  width=0.002, linewidth=0.6) +
    
    # labels next to points
    geom_text_repel(
      data=df_sorted,
      aes(H_Star, C_Star, label=Model, color=Model),
      family="serif", size=4, box.padding=0.4
    ) +
    
    scale_color_manual(values=model_colors) +
    scale_shape_manual(values=model_shapes) +
    
    labs(title=paste0("Central H–C Scatter — ", case_name, " (n=",n_val,")"),
         x=expression(italic(H)), y=expression(italic(C))) +
    
    theme_minimal(base_family="serif", base_size=14) +
    theme(legend.position="none")
  
  
  # -------------------------------------------------------------
  # 5. Time-Series plots — ASCENDING ORDER split into groups
  # -------------------------------------------------------------
  make_ts_plot <- function(ts_data, name) {
    df_ts <- data.frame(x=seq_along(ts_data), y=ts_data)
    
    ggplot(df_ts, aes(x=x, y=y)) +
      geom_line(color=model_colors[name], linewidth=0.5) +
      #scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000)) +
      scale_x_continuous(limits=c(min(df_ts$x), max(df_ts$x))) +
      scale_y_continuous(limits=c(min(df_ts$y), max(df_ts$y))) +
      labs(title=name) +
      theme_minimal(base_family="serif", base_size=8) +
      theme(
        axis.title.x = element_blank(),    
        axis.title.y = element_blank(),    
        plot.margin = margin(2,4,2,2)
      )
  }
  
  ts_top_plots    <- mapply(make_ts_plot, ts_top,    names(ts_top),    SIMPLIFY=FALSE)
  ts_bottom_plots <- mapply(make_ts_plot, ts_bottom, names(ts_bottom), SIMPLIFY=FALSE)
  if (!is.null(right_idx)) {
    ts_right_plots <- mapply(make_ts_plot, ts_right, names(ts_right), SIMPLIFY=FALSE)
  }
  
  # -------------------------------------------------------------
  # 6. Arrange layout 3–3–3
  # -------------------------------------------------------------
  if (n_ts == 9) {
    combined <- grid.arrange(
      arrangeGrob(grobs=ts_top_plots,    ncol=3),      # top
      arrangeGrob(
        p_center,
        arrangeGrob(grobs=ts_right_plots, ncol=1),     # right side
        ncol=2, widths=c(3,1.2)
      ),
      arrangeGrob(grobs=ts_bottom_plots, ncol=3),      # bottom
      ncol=1,
      heights=c(0.9, 5, 0.9)
    )
  } else {
    combined <- grid.arrange(            # Case 3: 6 TS only
      arrangeGrob(grobs=ts_top_plots, ncol=3),
      p_center,
      arrangeGrob(grobs=ts_bottom_plots, ncol=3),
      ncol=1,
      heights=c(0.9, 5, 0.9)
    )
  }
  
  # -------------------------------------------------------------
  # 7. Save PDF
  # -------------------------------------------------------------
  out_file <- file.path(plot_dir,
                        paste0("CentralPoints_TimeSeries_",
                               case_name, "_n", n_val, ".pdf"))
  
  ggsave(out_file, combined, width=14, height=12)
  
  message("✅ Saved: ", out_file)
}



###############################################################
# 🔁 RUN ALL CASES × SAMPLE SIZES
###############################################################
for (case_name in names(cases)) {
  for (n_val in sample_sizes) {
    plot_case_timeseries(case_name, n_val)
  }
}

message("🎉 All plots successfully generated!")

####################################################################################
# In this script, we generate plots of three selected time series points
#from various ARMA models based on their central H* values. The time series data
#is read from an Excel file, and the plots are saved in a structured directory
#based on cases and sample sizes.
#####################################################################################
###############################################################
# 📚 Required libraries
###############################################################
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(stringr)
library(ggrepel)
library(StatOrdPattHxC)

###############################################################
# 📁 File paths
###############################################################
#ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/TimeSeries_by_Case_n1000_n5000.xlsx"
#base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results"
ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/TimeSeries_by_Case_n1000_n5000.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"

cases <- list(
  Case1 = "Case1",
  Case2 = "Case2",
  Case3 = "Case3"
)

sample_sizes <- c(1000, 5000)

###############################################################
# 🎨 Colors & Shapes
###############################################################
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

###############################################################
# 🔧 Helper functions
###############################################################

clean_ts <- function(x) {
  x <- gsub("[c()]", "", x)
  x <- gsub("[[:space:]]+", "", x)
  as.numeric(strsplit(x, ",")[[1]])
}

###############################################################
# ⭐ 1) Plot 3 selected points + TS panels
###############################################################

plot_three_points_with_ts <- function(case_name, n_val) {
  
  message(">>> Creating 3-point plot for ", case_name, " n=", n_val)
  
  plot_dir <- file.path(base_plot_dir, case_name, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  df <- read_excel(ts_case_path, sheet=case_name) %>%
    filter(n == n_val) %>%
    mutate(Model = as.character(Model))
  
  if (nrow(df) < 3) {
    warning("Not enough models for ", case_name, " n=", n_val)
    return()
  }
  
  df_sorted <- df %>% arrange(H_Star)
  idx_top <- 1
  idx_mid <- ceiling(nrow(df_sorted)/2)
  idx_bottom <- nrow(df_sorted)
  
  df_sel <- df_sorted[c(idx_top, idx_mid, idx_bottom), ]
  df_sel$Position <- factor(c("Top","Middle","Bottom"))
  
  ts_all <- lapply(df_sel$TimeSeries, clean_ts)
  names(ts_all) <- df_sel$Model
  
  # --- Cropped feasible region ---
  data("LinfLsup")
  Linf_focus <- subset(LinfLsup, Dimension=="3") %>%
    filter(
      H >= min(df_sel$H_Star) - 0.05,
      H <= max(df_sel$H_Star) + 0.05,
      C >= min(df_sel$C_Star) - 0.05,
      C <= max(df_sel$C_Star) + 0.05
    )
  
  # --- Central H–C plot ---
  p_center <- ggplot() +
    geom_line(data=subset(Linf_focus, Side=="Lower"), aes(H,C), color="gray40") +
    geom_line(data=subset(Linf_focus, Side=="Upper"), aes(H,C), color="gray40") +
    
    geom_point(data=df_sel, aes(H_Star, C_Star, color=Model, shape=Model), size=3) +
    
    geom_errorbarh(data=df_sel,
                   aes(y=C_Star, xmin=H_Star-SemiLength_H, xmax=H_Star+SemiLength_H, color=Model),
                   height=0.002, linewidth=0.6) +
    
    geom_errorbar(data=df_sel,
                  aes(x=H_Star, ymin=C_Star-SemiLength_C, ymax=C_Star+SemiLength_C, color=Model),
                  width=0.002, linewidth=0.6) +
    
    geom_text_repel(
      data=df_sel, aes(H_Star, C_Star, label=Model, color=Model),
      family="serif", size=4
    ) +
    
    scale_color_manual(values=model_colors) +
    scale_shape_manual(values=model_shapes) +
    labs(title=paste0("Selected Points – ", case_name, " n=", n_val),
         x=expression(italic(H)), y=expression(italic(C))) +
    theme_minimal(base_family="serif", base_size=14) +
    theme(legend.position="none")
  
  # --- TS panels ---
  make_ts_plot <- function(ts_data, model_name, pos_label) {
    df_ts <- data.frame(x=seq_along(ts_data), y=ts_data)
    ggplot(df_ts, aes(x,y)) +
      geom_line(color=model_colors[model_name], linewidth=0.5) +
      scale_x_continuous(breaks=c(0,1000,2000,3000,4000,5000)) +
      scale_y_continuous(limits=c(min(df_ts$y), max(df_ts$y))) +
      labs(title=paste0(pos_label,": ",model_name)) +
      theme_minimal(base_family="serif", base_size=8) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
  }
  
  ts_plots <- mapply(make_ts_plot, ts_all, names(ts_all), df_sel$Position, SIMPLIFY=FALSE)
  ts_panel <- arrangeGrob(grobs=ts_plots, ncol=1)
  
  combined <- grid.arrange(p_center, ts_panel, ncol=2, widths=c(3,1.7))
  
  out_file <- file.path(plot_dir, paste0("new_centralpoint_timeseries_",case_name,"_n",n_val,".pdf"))
  ggsave(out_file, combined, width=14, height=10)
  
  message("✔ Saved: ", out_file)
}



for (case_name in names(cases)) {
  for (n_val in sample_sizes) {
    plot_three_points_with_ts(case_name, n_val)
  }
}

message("🎉 All 3-point plots completed successfully!")

###############################################################
# ⭐ NEW FUNCTION: Combined scatter plot (9 points)
###############################################################

plot_combined_threepoint_scatter <- function() {
  
  message(">>> Creating combined scatter plot (9 points total)")
  
  # Storage for selected points from all cases
  combined_points <- list()
  
  # Loop through Case1, Case2, Case3
  for (case_name in names(cases)) {
    
    df <- read_excel(ts_case_path, sheet=case_name) %>%
      mutate(Model = as.character(Model))
    
    # Sort ascending by H
    df_sorted <- df %>% arrange(H_Star)
    
    # Select top/mid/bottom
    idx_top <- 1
    idx_mid <- ceiling(nrow(df_sorted)/2)
    idx_bottom <- nrow(df_sorted)
    
    sel <- df_sorted[c(idx_top, idx_mid, idx_bottom), ]
    
    # Assign coefficient class based on CASE
    coef_class <- dplyr::case_when(
      case_name == "Case1" ~ "Positive",
      case_name == "Case2" ~ "Negative",
      case_name == "Case3" ~ "Mixed",
      TRUE ~ "Other"
    )
    
    sel$CoefClass <- coef_class
    sel$Case <- case_name
    
    combined_points[[case_name]] <- sel
  }
  
  # Combine 3 sets of selected points (9 rows total)
  df_all_sel <- bind_rows(combined_points)
  
  # --- CROP FEASIBLE REGION BASED ON THE 9 POINTS ---
  data("LinfLsup")
  Linf_focus <- subset(LinfLsup, Dimension=="3") %>%
    filter(
      H >= min(df_all_sel$H_Star) - 0.05,
      H <= max(df_all_sel$H_Star) + 0.05,
      C >= min(df_all_sel$C_Star) - 0.05,
      C <= max(df_all_sel$C_Star) + 0.05
    )
  
  # --- COMBINED SCATTER PLOT ---
  p <- ggplot() +
    geom_line(data=subset(Linf_focus,Side=="Lower"), aes(H,C), color="gray40") +
    geom_line(data=subset(Linf_focus,Side=="Upper"), aes(H,C), color="gray40") +
    
    geom_point(data=df_all_sel,
               aes(H_Star, C_Star, color=CoefClass, shape=Model),
               size=4, alpha=0.9) +
    
    geom_errorbarh(data=df_all_sel,
                   aes(y=C_Star,
                       xmin=H_Star-SemiLength_H,
                       xmax=H_Star+SemiLength_H,
                       color=CoefClass),
                   height=0.002, linewidth=0.7) +
    
    geom_errorbar(data=df_all_sel,
                  aes(x=H_Star,
                      ymin=C_Star-SemiLength_C,
                      ymax=C_Star+SemiLength_C,
                      color=CoefClass),
                  width=0.002, linewidth=0.7) +
    
    geom_text_repel(data=df_all_sel,
                    aes(H_Star, C_Star, label=Model, color=CoefClass),
                    family="serif", size=4, box.padding=0.4) +
    
    # Coefficient class color scheme
    scale_color_manual(values=c(
      "Positive"="darkblue",
      "Negative"="darkgreen",
      "Mixed"="orange"
    )) +
    
    scale_shape_manual(values=model_shapes) +
    
    labs(title="Combined H–C Scatter Plot (3 points from Case1–3)",
         #subtitle="Positive vs Negative vs Mixed ARMA coefficients",
         x=expression(italic(H)),
         y=expression(italic(C)),
         color="Coefficient Sign",
         shape="Model") +
    
    theme_minimal(base_family="serif", base_size=14)
  
  # Save plot
  out_file <- file.path(base_plot_dir, "combined_new_scatterplot.pdf")
  ggsave(out_file, p, width=10, height=7)
  
  message("✔ Saved combined scatter plot: ", out_file)
}
plot_combined_threepoint_scatter()

#####################################################################################
# Selected poins scatter plots with Linf–Lsup boundaries
################################################################################
###############################################################
# 📦 Libraries
###############################################################
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(StatOrdPattHxC)

###############################################################
# 📁 Paths
###############################################################
#full_data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/Entropy_Complexity_Results_all_Models.xlsx"
#output_dir     <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results"

full_data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"
###############################################################
# ⭐ Selected Models per Case
###############################################################
selected_models <- list(
  Case1 = c("ARMA11_M1", "AR1_M1", "AR1_M2"),
  Case2 = c("ARMA11_M3", "MA1_M3", "MA1_M4"),
  Case3 = c("ARMA22_M2", "ARMA22_M3", "MA2_M2")
)

sample_sizes <- c(1000, 5000)

###############################################################
# 🎨 Colors & Shapes
###############################################################
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

###############################################################
# ⭐ FUNCTION — Plot Selected Points for Each Case × n
###############################################################
plot_selected_scatter <- function(case_name, n_val) {
  
  message(">>> Processing ", case_name, " for n = ", n_val)
  
  models_this_case <- selected_models[[case_name]]
  df_list <- list()
  
  # --- Load data from each model sheet ---
  for (model_name in models_this_case) {
    
    sheet_to_read <- paste0(model_name, "_n", n_val)
    message("     Reading sheet: ", sheet_to_read)
    
    df_model <- read_excel(full_data_path, sheet = sheet_to_read)
    df_model$Model <- model_name   # ensure model column
    
    df_list[[model_name]] <- df_model
  }
  
  df <- bind_rows(df_list)
  
  # Rename to simpler H, C
  df <- df %>%
    rename(
      H = H_Shannon,
      C = C_Shannon
    )
  
  # --- Feasible region crop ---
  data("LinfLsup")
  Linf_focus <- subset(LinfLsup, Dimension=="3") %>%
    filter(
      H >= min(df$H) - 0.05,
      H <= max(df$H) + 0.05,
      C >= min(df$C) - 0.05,
      C <= max(df$C) + 0.05
    )
  
  # --- Scatter plot (legend ON, no labels inside plot) ---
  p <- ggplot() +
    geom_line(data=subset(Linf_focus, Side=="Lower"), aes(H,C), color="gray40") +
    geom_line(data=subset(Linf_focus, Side=="Upper"), aes(H,C), color="gray40") +
    
    geom_point(
      data=df,
      aes(H, C, color=Model, shape=Model),
      size=4
    ) +
    
    scale_color_manual(values=model_colors) +
    scale_shape_manual(values=model_shapes) +
    
    labs(
      title = paste0("Selected Models — ", case_name, " (n=", n_val, ")"),
      x = expression(italic(H)),
      y = expression(italic(C)),
      color = "Model",
      shape = "Model"
    ) +
    
    theme_minimal(base_family="serif", base_size=14)
  
  # Save output
  outfile <- file.path(output_dir,
                       paste0("Scatter_Selected_", case_name, "_n", n_val, ".pdf"))
  
  ggsave(outfile, p, width = 8, height = 6)
  message("✔ Saved: ", outfile)
}

###############################################################
# ⭐ RUN FOR ALL CASES × SAMPLE SIZES
###############################################################
for (case_name in names(selected_models)) {
  for (n_val in sample_sizes) {
    plot_selected_scatter(case_name, n_val)
  }
}

message("\n🎉 All 6 scatter plots created successfully (legend ON, no repeated labels)!\n")

#----------------------------------------------------------------------------------
#-------------------------------------------------------------------------
# In this script we generate a density of the Entropy-Complexity plane
# for a selected cases (Case1-3) and sample size (n = 5000 and 10000).
# The density visualizes the density of points corresponding to different
# time series models (AR, MA, ARMA) in the Entropy-Complexity space.
#--------------------------------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
base_output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/ARMA Plots n1000_n5000"

#--------------------------------------------------------
# Case definitions (models per case)
#--------------------------------------------------------
cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1",
            "MA1_M1","MA1_M2","MA2_M1",
            "ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4",
            "MA1_M3","MA1_M4","MA2_M4",
            "ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3",
            "MA2_M2","MA2_M3",
            "ARMA22_M2","ARMA22_M3")
)

#--------------------------------------------------------
# Master color palette per model
#--------------------------------------------------------
model_colors <- c(
  "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
  "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
  "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
  "ARMA22_M3"="#fb9a99","AR1_M1"="#8dd3c7","AR1_M2"="#ff0000","ARMA11_M2"="#bebada",
  "MA1_M1"="#fb8072","MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69",
  "ARMA11_M4"="#fccde5","MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
)

#--------------------------------------------------------
# Sample sizes
#--------------------------------------------------------
sample_sizes <- c(1000, 5000)

#--------------------------------------------------------
# Feasible region (Linf / Lsup) for D = 3
#--------------------------------------------------------
data("LinfLsup")
D <- 3
Linf_all <- subset(LinfLsup, Side == "Lower" & Dimension == as.character(D))
Lsup_all <- subset(LinfLsup, Side == "Upper" & Dimension == as.character(D))

#--------------------------------------------------------
# Main loops: over cases and sample sizes
#--------------------------------------------------------
for (case_name in names(cases)) {
  
  models_case <- cases[[case_name]]
  
  for (n_val in sample_sizes) {
    
    message("===== Processing ", case_name, ", n = ", n_val, " =====")
    
    #--------------------------------------------------------
    # Read Excel sheets for this case & n
    #--------------------------------------------------------
    df_list <- lapply(models_case, function(m) {
      sheet_name <- paste0(m, "_n", n_val)
      message("  Reading sheet: ", sheet_name)
      read_excel(data_path, sheet = sheet_name) |>
        mutate(Model = m, n = n_val)
    })
    
    df_case <- bind_rows(df_list) |>
      dplyr::rename(
        H = H_Shannon,
        C = C_Shannon
      ) |>
      dplyr::filter(is.finite(H), is.finite(C))
    
    if (nrow(df_case) == 0) {
      warning("No data for ", case_name, " n = ", n_val, ", skipping.")
      next
    }
    
    #--------------------------------------------------------
    # Classify Families (AR / MA / ARMA) and fine Types (each model)
    #--------------------------------------------------------
    df_case <- df_case |>
      mutate(
        Family = case_when(
          str_starts(Model, "ARMA") ~ "ARMA",
          str_starts(Model, "AR")   ~ "AR",
          str_starts(Model, "MA")   ~ "MA",
          TRUE                      ~ "Other"
        ),
        # Fine-grained type: each model in this case is its own "type"
        Type = factor(Model, levels = models_case)
      )
    
    df_case$Family <- factor(df_case$Family,
                             levels = c("AR","MA","ARMA","Other"))
    
    #--------------------------------------------------------
    # Feasible region (crop to H-range of this case)
    #--------------------------------------------------------
    H_min <- min(df_case$H)
    H_max <- max(df_case$H)
    C_min <- min(df_case$C)
    C_max <- max(df_case$C)
    
    Linf_crop <- Linf_all |>
      dplyr::filter(H >= H_min, H <= H_max)
    Lsup_crop <- Lsup_all |>
      dplyr::filter(H >= H_min, H <= H_max)
    
    #--------------------------------------------------------
    # Dynamic bandwidth for density smoothing
    #--------------------------------------------------------
    Hx <- diff(range(df_case$H))
    Cx <- diff(range(df_case$C))
    h_vec <- c(Hx / 5, Cx / 5)   # adaptive smoothing
    
    #--------------------------------------------------------
    # Colors per fine model Type — use your predefined palette
    #--------------------------------------------------------
    type_cols <- model_colors[models_case]
    
    #--------------------------------------------------------
    # Plot
    #--------------------------------------------------------
    p <- ggplot(df_case, aes(H, C)) +
      
      # Feasible region boundaries
      geom_line(data = Linf_crop, aes(H, C),
                colour = "black", linewidth = 0.8) +
      geom_line(data = Lsup_crop, aes(H, C),
                colour = "black", linewidth = 0.8) +
      
      # Background points
      geom_point(alpha = 0.07, size = 0.8, colour = "black") +
      
      # Density polygons per model Type
      stat_density_2d(
        aes(fill = Type, alpha = after_stat(level)),
        geom = "polygon",
        contour = TRUE,
        bins = 10,
        colour = NA,
        h = h_vec
      ) +
      
      # Contour lines per model Type
      stat_density_2d(
        aes(color = Type),
        contour = TRUE,
        bins = 6,
        linewidth = 0.25,
        h = h_vec
      ) +
      
      scale_fill_manual(values = type_cols, name = "Model (Type)") +
      scale_color_manual(values = type_cols, guide = "none") +
      scale_alpha(range = c(0.10, 0.55), guide = "none") +
      
      labs(
        title = paste("Density plot by Model Type —", case_name, "(n =", n_val, ")"),
        x = expression(italic(H)),
        y = expression(italic(C))
      ) +
      
      coord_cartesian(
        xlim = c(H_min, H_max),
        ylim = c(C_min - 0.02, C_max + 0.02)
      ) +
      
      theme_minimal(base_size = 16) +
      theme(
        panel.grid      = element_blank(),
        legend.position = "right"
      )
    
    print(p)
    
    #--------------------------------------------------------
    # Save plot
    #--------------------------------------------------------
    case_output_dir <- file.path(base_output_dir, case_name, "Plots")
    dir.create(case_output_dir, recursive = TRUE, showWarnings = FALSE)
    
    output_plot_path <- file.path(
      case_output_dir,
      paste0("HC_Density_ModelTypes_", case_name, "_n", n_val, ".pdf")
    )
    
    ggsave(output_plot_path, p, width = 10, height = 8)
    message("  Saved: ", output_plot_path, "\n")
  }
}

message("🎉 All density plots completed for all cases, sample sizes, and model types!")

#-------------------------------------------------------------------------------------------
#end of script
#----------------------------------------------------------------------------------------
#This script generates 2D histograms of the Entropy-Complexity plane
#for selected cases (Case1-3) and sample sizes (n = 1000 and 5000).
#The histograms visualize the distribution of points corresponding to different
#time series models (AR, MA, ARMA) in the Entropy-Complexity space.
#_--------------------------------------------------------------------------------
    library(readxl)
      library(dplyr)
      library(ggplot2)
      library(stringr)
      library(StatOrdPattHxC)
      
      #--------------------------------------------------------
      # CHOOSE PLOT OPTION HERE
      #   1 = 2D histogram (all points) + density contours by model type
      #   2 = 2D histogram faceted by model type (one panel per model)
      #   3 = BIVARIATE NORMAL ELLIPSES per model type
      #--------------------------------------------------------
      plot_option <- 3   # <--- SET THIS TO 3 for bivariate normal model view
      
      #--------------------------------------------------------
      # Paths
      #--------------------------------------------------------
      data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results n1000_n5000/HC_Results_all_Models_n1000_n5000.xlsx"
      base_output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Plots/ARMA plots/ARMA Plots n1000_n5000"
      
      #--------------------------------------------------------
      # Case definitions (models per case)
      #--------------------------------------------------------
      cases <- list(
        Case1 = c("AR1_M1","AR1_M2","AR2_M1",
                  "MA1_M1","MA1_M2","MA2_M1",
                  "ARMA11_M1","ARMA11_M2","ARMA22_M1"),
        Case2 = c("AR1_M3","AR1_M4","AR2_M4",
                  "MA1_M3","MA1_M4","MA2_M4",
                  "ARMA11_M3","ARMA11_M4","ARMA22_M4"),
        Case3 = c("AR2_M2","AR2_M3",
                  "MA2_M2","MA2_M3",
                  "ARMA22_M2","ARMA22_M3")
      )
      
      #--------------------------------------------------------
      # Master color palette per model (for contours/ellipses)
      #--------------------------------------------------------
      model_colors <- c(
        "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
        "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
        "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
        "ARMA22_M3"="#fb9a99","AR1_M1"="#8dd3c7","AR1_M2"="#ff0000","ARMA11_M2"="#bebada",
        "MA1_M1"="#fb8072","MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69",
        "ARMA11_M4"="#fccde5","MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
      )
      
      #--------------------------------------------------------
      # Sample sizes
      #--------------------------------------------------------
      sample_sizes <- c(1000, 5000)
      
      #--------------------------------------------------------
      # Feasible region (Linf / Lsup) for D = 3
      #--------------------------------------------------------
      data("LinfLsup")
      D <- 3
      Linf_all <- subset(LinfLsup, Side == "Lower" & Dimension == as.character(D))
      Lsup_all <- subset(LinfLsup, Side == "Upper" & Dimension == as.character(D))
      
      #--------------------------------------------------------
      # Main loops: over cases and sample sizes
      #--------------------------------------------------------
      for (case_name in names(cases)) {
        
        models_case <- cases[[case_name]]
        
        for (n_val in sample_sizes) {
          
          message("===== Processing ", case_name, ", n = ", n_val, " =====")
          
          #--------------------------------------------------------
          # Read Excel sheets for this case & n
          #--------------------------------------------------------
          df_list <- lapply(models_case, function(m) {
            sheet_name <- paste0(m, "_n", n_val)
            message("  Reading sheet: ", sheet_name)
            read_excel(data_path, sheet = sheet_name) |>
              mutate(Model = m, n = n_val)
          })
          
          df_case <- bind_rows(df_list) |>
            dplyr::rename(
              H = H_Shannon,
              C = C_Shannon
            ) |>
            dplyr::filter(is.finite(H), is.finite(C))
          
          if (nrow(df_case) == 0) {
            warning("No data for ", case_name, " n = ", n_val, ", skipping.")
            next
          }
          
          #--------------------------------------------------------
          # Classify Families (AR / MA / ARMA) and fine Types (each model)
          #--------------------------------------------------------
          df_case <- df_case |>
            mutate(
              Family = case_when(
                str_starts(Model, "ARMA") ~ "ARMA",
                str_starts(Model, "AR")   ~ "AR",
                str_starts(Model, "MA")   ~ "MA",
                TRUE                      ~ "Other"
              ),
              # Fine-grained type: each model in this case is its own "type"
              Type = factor(Model, levels = models_case)
            )
          
          df_case$Family <- factor(df_case$Family,
                                   levels = c("AR","MA","ARMA","Other"))
          
          #--------------------------------------------------------
          # Feasible region (crop to H-range of this case)
          #--------------------------------------------------------
          H_min <- min(df_case$H)
          H_max <- max(df_case$H)
          C_min <- min(df_case$C)
          C_max <- max(df_case$C)
          
          Linf_crop <- Linf_all |>
            dplyr::filter(H >= H_min, H <= H_max)
          Lsup_crop <- Lsup_all |>
            dplyr::filter(H >= H_min, H <= H_max)
          
          #--------------------------------------------------------
          # Dynamic bandwidth for density smoothing (used in Option 1)
          #--------------------------------------------------------
          Hx <- diff(range(df_case$H))
          Cx <- diff(range(df_case$C))
          h_vec <- c(Hx / 5, Cx / 5)   # adaptive smoothing
          
          #--------------------------------------------------------
          # Colors per fine model Type
          #--------------------------------------------------------
          type_cols <- model_colors[models_case]
          # lighter, semi-transparent versions for ellipse fill in option 3
          fill_cols <- grDevices::adjustcolor(type_cols, alpha.f = 0.25)
          
          #--------------------------------------------------------
          # PLOT: three options
          #--------------------------------------------------------
          if (plot_option == 1) {
            #======================================================
            # OPTION 1:
            # 2D histogram of ALL points + density contours by Type
            #======================================================
            
            p <- ggplot(df_case, aes(H, C)) +
              
              # 2D histogram of ALL points
              stat_bin_2d(
                bins = 20,
                aes(fill = after_stat(count))
              ) +
              
              # colour scale for counts
              scale_fill_viridis_c(
                option = "magma",
                trans  = "sqrt",
                name   = "Count"
              ) +
              
              # Feasible region boundaries
              geom_line(
                data        = Linf_crop,
                aes(H, C),
                inherit.aes = FALSE,
                colour      = "black",
                linewidth   = 0.8
              ) +
              geom_line(
                data        = Lsup_crop,
                aes(H, C),
                inherit.aes = FALSE,
                colour      = "black",
                linewidth   = 0.8
              ) +
              
              # Density contours per model Type (smooth overlay)
              stat_density_2d(
                aes(color = Type),
                contour   = TRUE,
                bins      = 15,
                linewidth = 0.35,
                h         = h_vec
              ) +
              
              scale_color_manual(values = type_cols, name = "Model (Type)") +
              
              labs(
                title = paste("2D histogram + density —", case_name, "(n =", n_val, ")"),
                x     = expression(italic(H)),
                y     = expression(italic(C))
              ) +
              
              coord_cartesian(
                xlim = c(H_min, H_max),
                ylim = c(C_min - 0.02, C_max + 0.02)
              ) +
              
              theme_minimal(base_size = 16) +
              theme(
                panel.grid      = element_blank(),
                legend.position = "right"
              )
            
            plot_name_stub <- "HC_2DHistDensity_ModelTypes"
            
          } else if (plot_option == 2) {
            #======================================================
            # OPTION 2:
            # 2D histogram FACETED by model Type (one panel/model)
            #======================================================
            
            p <- ggplot(df_case, aes(H, C)) +
              
              # 2D histogram within each facet (per Type)
              stat_bin_2d(
                bins = 40,
                aes(fill = after_stat(count))
              ) +
              
              scale_fill_viridis_c(
                option = "magma",
                trans  = "sqrt",
                name   = "Count"
              ) +
              
              # Feasible region boundaries (same in all facets)
              geom_line(
                data        = Linf_crop,
                aes(H, C),
                inherit.aes = FALSE,
                colour      = "black",
                linewidth   = 0.6
              ) +
              geom_line(
                data        = Lsup_crop,
                aes(H, C),
                inherit.aes = FALSE,
                colour      = "black",
                linewidth   = 0.6
              ) +
              
              facet_wrap(~ Type, ncol = 3) +
              
              labs(
                title = paste("2D histograms by Model Type —", case_name, "(n =", n_val, ")"),
                x     = expression(italic(H)),
                y     = expression(italic(C))
              ) +
              
              coord_cartesian(
                xlim = c(H_min, H_max),
                ylim = c(C_min - 0.02, C_max + 0.02)
              ) +
              
              theme_minimal(base_size = 14) +
              theme(
                panel.grid = element_blank()
              )
            
            plot_name_stub <- "HC_2DHistFacets_ModelTypes"
            
          } else if (plot_option == 3) {
            #======================================================
            # OPTION 3 (NEW):
            # BIVARIATE NORMAL ELLIPSES per model Type
            #
            # stat_ellipse(type = "norm") fits a bivariate normal
            # to each group (Type) and draws a constant-density
            # ellipse (here: 95% level).
            #======================================================
            
            p <- ggplot(df_case, aes(H, C)) +
              
              # light scatter (to see the raw clouds)
              geom_point(
                aes(color = Type),
                alpha = 0.15,
                size  = 0.7
              ) +
              
              # 95% bivariate normal ellipse per model Type
              stat_ellipse(
                aes(color = Type, fill = Type),
                type      = "norm",     # assumes bivariate normal within each Type
                level     = 0.95,       # 95% ellipse
                linewidth = 0.9,
                alpha     = 0.35
              ) +
              
              # Feasible region boundaries
              geom_line(
                data        = Linf_crop,
                aes(H, C),
                inherit.aes = FALSE,
                colour      = "black",
                linewidth   = 0.8
              ) +
              geom_line(
                data        = Lsup_crop,
                aes(H, C),
                inherit.aes = FALSE,
                colour      = "black",
                linewidth   = 0.8
              ) +
              
              scale_color_manual(values = type_cols, name = "Model (Type)") +
              scale_fill_manual(values  = fill_cols, guide = "none") +
              
              labs(
                title = paste("Bivariate normal ellipses by Model Type —", case_name, "(n =", n_val, ")"),
                x     = expression(italic(H)),
                y     = expression(italic(C))
              ) +
              
              coord_cartesian(
                xlim = c(H_min, H_max),
                ylim = c(C_min - 0.02, C_max + 0.02)
              ) +
              
              theme_minimal(base_size = 16) +
              theme(
                panel.grid      = element_blank(),
                legend.position = "right"
              )
            
            plot_name_stub <- "HC_BivNormEllipses_ModelTypes"
            
          } else {
            stop("Unknown plot_option. Use 1, 2, or 3.")
          }
          
          #--------------------------------------------------------
          # Print plot to device
          #--------------------------------------------------------
          print(p)
          
          #--------------------------------------------------------
          # Save plot
          #--------------------------------------------------------
          case_output_dir <- file.path(base_output_dir, case_name, "Plots")
          dir.create(case_output_dir, recursive = TRUE, showWarnings = FALSE)
          
          output_plot_path <- file.path(
            case_output_dir,
            paste0(plot_name_stub, "_", case_name, "_n", n_val, ".pdf")
          )
          
          ggsave(output_plot_path, p, width = 10, height = 8)
          message("  Saved: ", output_plot_path, "\n")
        }
      }
      
      message("🎉 All plots completed for all cases, sample sizes, and model types!")
      #-------------------------------------------------------------------------------------------