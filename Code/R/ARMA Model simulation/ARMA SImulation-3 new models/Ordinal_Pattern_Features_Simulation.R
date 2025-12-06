#-------------------------------------------------------------------------
# Extended Entropy–Complexity Simulation Script
# Shannon, Rényi, Tsallis, Fisher + variances + semi-lengths
#-----------------------------------------------------------------------

#-----------------------------
# Libraries
#-----------------------------
library(StatOrdPattHxC)
library(dplyr)
library(writexl)

#-----------------------------
# Global Parameters
#-----------------------------
set.seed(1234567890, kind = "Mersenne-Twister")

D    <- 3                 # Embedding dimension
N    <- c(5000, 10000)    # Sample sizes
R    <- 100               # Number of replications
BETA <- 1.5               # Parameter for Rényi / Tsallis

# Make J available globally (used internally by some functions)
J <- factorial(D)
assign("J", J, envir = .GlobalEnv)

#======================================================================
# Helper Functions: Jensen–Shannon, Generalized Complexities
#======================================================================

# Jensen–Shannon divergence between two discrete pmfs p, q
Jensen_Shannon <- function(p, q) {
  m  <- 0.5 * (p + q)
  js <- 0.5 * sum(ifelse(p == 0, 0, p * log((p + 1e-12)/(m + 1e-12)))) +
    0.5 * sum(ifelse(q == 0, 0, q * log((q + 1e-12)/(m + 1e-12))))
  return(js)
}

# Generalized statistical complexity based on a given entropy value
# (for Rényi or Tsallis)
GeneralizedComplexity <- function(prob, entropy_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * entropy_value)
}

# Fisher-based statistical complexity
FisherBasedComplexity <- function(prob, fisher_value) {
  Pe <- rep(1 / length(prob), length(prob))
  JS <- Jensen_Shannon(prob, Pe)
  return(JS * fisher_value)
}

#======================================================================
# Model Definitions: AR, MA, ARMA
#======================================================================

ar1_models <- list(
  AR1_M1 = list(ar = c(0.8),  type = "AR1_M1"),
  AR1_M2 = list(ar = c(0.1),  type = "AR1_M2"),
  AR1_M3 = list(ar = c(-0.8), type = "AR1_M3"),
  AR1_M4 = list(ar = c(-0.1), type = "AR1_M4")
)

ar2_models <- list(
  AR2_M1 = list(ar = c(0.1, 0.8),   type = "AR2_M1"),
  AR2_M2 = list(ar = c(-0.8, 0.1),  type = "AR2_M2"),
  AR2_M3 = list(ar = c(0.1, -0.8),  type = "AR2_M3"),
  AR2_M4 = list(ar = c(-0.8, -0.1), type = "AR2_M4")
)

ma1_models <- list(
  MA1_M1 = list(ma = c(0.8),  type = "MA1_M1"),
  MA1_M2 = list(ma = c(0.1),  type = "MA1_M2"),
  MA1_M3 = list(ma = c(-0.8), type = "MA1_M3"),
  MA1_M4 = list(ma = c(-0.1), type = "MA1_M4")
)

ma2_models <- list(
  MA2_M1 = list(ma = c(0.1, 0.8),   type = "MA2_M1"),
  MA2_M2 = list(ma = c(-0.8, 0.1),  type = "MA2_M2"),
  MA2_M3 = list(ma = c(0.1, -0.8),  type = "MA2_M3"),
  MA2_M4 = list(ma = c(-0.8, -0.1), type = "MA2_M4")
)

arma11_models <- list(
  ARMA11_M1 = list(ar = c(0.8),  ma = c(0.8),  type = "ARMA11_M1"),
  ARMA11_M2 = list(ar = c(0.1),  ma = c(0.1),  type = "ARMA11_M2"),
  ARMA11_M3 = list(ar = c(-0.8), ma = c(-0.8), type = "ARMA11_M3"),
  ARMA11_M4 = list(ar = c(-0.1), ma = c(-0.1), type = "ARMA11_M4")
)

arma22_models <- list(
  ARMA22_M1 = list(ar = c(0.1, 0.8),   ma = c(0.1, 0.8),   type = "ARMA22_M1"),
  ARMA22_M2 = list(ar = c(-0.8, 0.1),  ma = c(-0.8, 0.1),  type = "ARMA22_M2"),
  ARMA22_M3 = list(ar = c(0.1, -0.8),  ma = c(0.1, -0.8),  type = "ARMA22_M3"),
  ARMA22_M4 = list(ar = c(-0.8, -0.1), ma = c(-0.8, -0.1), type = "ARMA22_M4")
)

# Combine all models into one list
all_models <- c(
  ar1_models, ar2_models,
  ma1_models, ma2_models,
  arma11_models, arma22_models
)

#======================================================================
# Core Simulation Function
#======================================================================

generate_model_data <- function(model, n, D, R, beta = BETA) {
  results   <- list()
  ts_store  <- list()
  z_alpha   <- qnorm(1 - 0.05 / 2)   # z for 95% CI
  
  for (r in 1:R) {
    #-----------------------------
    # Simulate time series
    #-----------------------------
    ts_data <- arima.sim(model = model, n = n)
    ProbTS  <- OPprob(ts_data, emb = D)
    
    #-----------------------------
    # Entropies
    #-----------------------------
    Hs <- HShannon(ProbTS)
    Hr <- HRenyi(ProbTS, beta = beta)
    Ht <- HTsallis(ProbTS, beta = beta)
    Hf <- HFisher(ProbTS)
    
    #-----------------------------
    # Complexities
    #-----------------------------
    C_Shannon <- StatComplexity(ProbTS)          # Shannon-based
    C_Renyi   <- GeneralizedComplexity(ProbTS, Hr)
    C_Tsallis <- GeneralizedComplexity(ProbTS, Ht)
    C_Fisher  <- FisherBasedComplexity(ProbTS, Hf)
    
    #-----------------------------
    # Variances of entropies via sigma2q
    # ent codes: "S", "R", "T", "F"
    #-----------------------------
    Var_HS <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "S"))
    Var_HR <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "R", beta = beta))
    Var_HT <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "T", beta = beta))
    Var_HF <- suppressWarnings(sigma2q(ts_data, emb = D, ent = "F"))
    
    #-----------------------------
    # Shannon complexity variance (distance-based)
    # (as in your previous system)
    #-----------------------------
    Var_HI <- suppressWarnings(asymptoticVarHShannonMultinomial(ProbTS, n - 2))
    a_ratio <- Var_HS / Var_HI
    Var_CI  <- suppressWarnings(varC(ProbTS, n - 2))
    Var_CS  <- a_ratio * Var_CI   # distance-based variance for C_Shannon
    
    #-----------------------------
    # Semi-lengths for entropies & Shannon complexity
    #-----------------------------
    Semi_HS <- ifelse(!is.finite(Var_HS) | Var_HS <= 0, NA,
                      sqrt(Var_HS) / sqrt(n - 3) * z_alpha)
    Semi_HR <- ifelse(!is.finite(Var_HR) | Var_HR <= 0, NA,
                      sqrt(Var_HR) / sqrt(n - 3) * z_alpha)
    Semi_HT <- ifelse(!is.finite(Var_HT) | Var_HT <= 0, NA,
                      sqrt(Var_HT) / sqrt(n - 3) * z_alpha)
    Semi_HF <- ifelse(!is.finite(Var_HF) | Var_HF <= 0, NA,
                      sqrt(Var_HF) / sqrt(n - 3) * z_alpha)
    
    Semi_CS <- ifelse(!is.finite(Var_CS) | Var_CS <= 0, NA,
                      sqrt(Var_CS) / sqrt(n - 3) * z_alpha)
    
    #-----------------------------
    # Store summary row
    #-----------------------------
    results[[r]] <- data.frame(
      Model   = model$type,
      n       = n,
      Rep     = r,
      
      # Entropies
      H_Shannon = Hs,
      H_Renyi   = Hr,
      H_Tsallis = Ht,
      H_Fisher  = Hf,
      
      # Complexities
      C_Shannon = C_Shannon,
      C_Renyi   = C_Renyi,
      C_Tsallis = C_Tsallis,
      C_Fisher  = C_Fisher,
      
      # Variances of entropies
      Var_H_Shannon = Var_HS,
      Var_H_Renyi   = Var_HR,
      Var_H_Tsallis = Var_HT,
      Var_H_Fisher  = Var_HF,
      
      # Variance of Shannon complexity
      Var_C_Shannon = Var_CS,
      
      # Semi-lengths for entropies
      SemiLength_H_Shannon = Semi_HS,
      SemiLength_H_Renyi   = Semi_HR,
      SemiLength_H_Tsallis = Semi_HT,
      SemiLength_H_Fisher  = Semi_HF,
      
      # Semi-length for Shannon complexity
      SemiLength_C_Shannon = Semi_CS
    )
    
    #-----------------------------
    # Store raw time series
    #-----------------------------
    ts_store[[r]] <- data.frame(
      Model = model$type,
      n     = n,
      Rep   = r,
      Value = as.numeric(ts_data)
    )
  }
  
  list(
    summary    = bind_rows(results),
    timeseries = bind_rows(ts_store)
  )
}

#======================================================================
# Run Simulations for All Models and Sample Sizes
#======================================================================

all_summary    <- list()
all_timeseries <- list()

for (n_val in N) {
  message("Simulating for n = ", n_val)
  
  for (model in all_models) {
    message("  Model: ", model$type)
    res <- generate_model_data(model, n = n_val, D = D, R = R, beta = BETA)
    
    sheet_name <- paste0(model$type, "_n", n_val)
    all_summary[[sheet_name]]    <- res$summary
    all_timeseries[[sheet_name]] <- res$timeseries
  }
}

#======================================================================
# Save Results to Excel
#======================================================================

summary_path <- "HC_Results_all_Models_extended_entropies.xlsx"
ts_path      <- "TimeSeries_Data_all_Models_extended.xlsx"

write_xlsx(all_summary,    path = summary_path)
write_xlsx(all_timeseries, path = ts_path)

cat("✅ All simulations complete!\n")
cat("📁 Extended entropy–complexity results saved at: ", summary_path, "\n")
cat("📁 Time series data for all models saved at: ", ts_path, "\n")
#-------------------------------------------------------------------------------
## End of Extended Entropy–Complexity Simulation Script
#-------------------------------------------------------------------------------
#======================================================================
# Central Point Extraction for Extended Entropies/Complexities
#======================================================================
# --- Load Required Packages ---
library(readxl)
library(openxlsx)
library(dplyr)
library(stats)

# --- Define file paths ---
data_path   <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/HC_Results_all_Models_extended_entropies.xlsx"
output_path <- data_path  # same file

# --- Load sheet names ---
all_sheets <- excel_sheets(data_path)

# --- Helper: PCA-based emblematic point, extended but backward-compatible ---
get_emblematic_point <- function(df, model_name, n_value) {
  
  # We ONLY require H_Shannon and C_Shannon to exist
  if (!all(c("H_Shannon", "C_Shannon") %in% names(df))) {
    stop("Sheet for model ", model_name, " is missing H_Shannon or C_Shannon.")
  }
  
  # Harmonize Shannon variance column names if present
  # (old files: Var_H, Var_C; new files: Var_H_Shannon, Var_C_Shannon)
  if ("Var_H_Shannon" %in% names(df) && !"Var_H" %in% names(df)) {
    df$Var_H <- df$Var_H_Shannon
  }
  if ("Var_C_Shannon" %in% names(df) && !"Var_C" %in% names(df)) {
    df$Var_C <- df$Var_C_Shannon
  }
  
  # Ensure numeric for relevant columns (if they exist)
  numeric_cols <- intersect(
    c("H_Shannon", "C_Shannon",
      "H_Renyi",   "C_Renyi",
      "H_Tsallis", "C_Tsallis",
      "H_Fisher",  "C_Fisher",
      "Var_H", "Var_C",
      "Var_H_Renyi", "Var_H_Tsallis", "Var_H_Fisher"),
    names(df)
  )
  
  df <- df %>%
    mutate(across(all_of(numeric_cols), as.numeric))
  
  # --- PCA-based central point (Shannon H–C) ---
  dat <- df[, c("H_Shannon", "C_Shannon")]
  pcs <- prcomp(dat, scale. = FALSE)
  pcscores <- pcs$x[, 1]
  N <- nrow(dat)
  median_index <- order(pcscores)[ceiling((N + 1) / 2)]
  
  row_em <- df[median_index, , drop = FALSE]
  
  # --- Extract variances (may be NA if not present) ---
  varHS <- if ("Var_H" %in% names(row_em)) row_em$Var_H else NA_real_
  varCS <- if ("Var_C" %in% names(row_em)) row_em$Var_C else NA_real_
  
  varHR <- if ("Var_H_Renyi"   %in% names(row_em)) row_em$Var_H_Renyi   else NA_real_
  varHT <- if ("Var_H_Tsallis" %in% names(row_em)) row_em$Var_H_Tsallis else NA_real_
  varHF <- if ("Var_H_Fisher"  %in% names(row_em)) row_em$Var_H_Fisher  else NA_real_
  
  # --- Compute semi-lengths from variances (95% CI) ---
  z_alpha <- qnorm(1 - 0.05 / 2)
  semi_fun <- function(v) {
    ifelse(is.finite(v) & v > 0,
           sqrt(v) / sqrt(n_value - 3) * z_alpha,
           NA_real_)
  }
  
  Semi_HS <- semi_fun(varHS)
  Semi_CS <- semi_fun(varCS)
  Semi_HR <- semi_fun(varHR)
  Semi_HT <- semi_fun(varHT)
  Semi_HF <- semi_fun(varHF)
  
  # --- Safe extract helper: if a column is missing, return NA ---
  get_col <- function(df_row, nm) {
    if (nm %in% names(df_row)) df_row[[nm]] else NA_real_
  }
  
  # --- Return emblematic row with all available measures ---
  tibble(
    Model = model_name,
    n     = n_value,
    Emblematic_Row = median_index,
    
    # Entropies (Shannon always present; others if available)
    H_Shannon = get_col(row_em, "H_Shannon"),
    H_Renyi   = get_col(row_em, "H_Renyi"),
    H_Tsallis = get_col(row_em, "H_Tsallis"),
    H_Fisher  = get_col(row_em, "H_Fisher"),
    
    # Complexities (Shannon always; others if available)
    C_Shannon = get_col(row_em, "C_Shannon"),
    C_Renyi   = get_col(row_em, "C_Renyi"),
    C_Tsallis = get_col(row_em, "C_Tsallis"),
    C_Fisher  = get_col(row_em, "C_Fisher"),
    
    # Variances
    Var_H_Shannon = varHS,
    Var_H_Renyi   = varHR,
    Var_H_Tsallis = varHT,
    Var_H_Fisher  = varHF,
    Var_C_Shannon = varCS,
    
    # Semi-lengths (95% CI)
    SemiLength_H_Shannon = Semi_HS,
    SemiLength_H_Renyi   = Semi_HR,
    SemiLength_H_Tsallis = Semi_HT,
    SemiLength_H_Fisher  = Semi_HF,
    SemiLength_C_Shannon = Semi_CS
  )
}

# --- Storage lists for n = 5000 and n = 10000 ---
results_5000  <- list()
results_10000 <- list()

# --- Loop through all sheets ---
for (sheet in all_sheets) {
  df <- read_excel(data_path, sheet = sheet)
  
  # Require ONLY H_Shannon & C_Shannon to exist
  if (!all(c("H_Shannon", "C_Shannon") %in% names(df))) {
    message("⚠️ Skipping ", sheet, " — missing H_Shannon or C_Shannon.")
    next
  }
  
  # Extract model name and n value from sheet name, e.g. "AR1_M1_n5000"
  parts <- strsplit(sheet, "_")[[1]]
  if (length(parts) < 3) {
    message("⚠️ Skipping ", sheet, " — unexpected sheet naming.")
    next
  }
  model_name <- paste(parts[1], parts[2], sep = "_")
  n_value    <- as.numeric(gsub("n", "", parts[length(parts)]))
  
  if (is.na(n_value)) {
    message("⚠️ Skipping ", sheet, " — unable to parse n value.")
    next
  }
  
  message("✅ Processing ", sheet, " (n = ", n_value, ") ...")
  
  # Compute emblematic point (now robust to missing extended columns)
  em <- get_emblematic_point(df, model_name, n_value)
  
  # Store results by n
  if (n_value == 5000)  results_5000[[model_name]]  <- em
  if (n_value == 10000) results_10000[[model_name]] <- em
}

# --- Combine results into data frames ---
df_5000  <- bind_rows(results_5000)
df_10000 <- bind_rows(results_10000)

if (nrow(df_5000) == 0 & nrow(df_10000) == 0) {
  cat("\n⚠️ No results were generated. Check sheet naming or n values.\n")
} else {
  # --- Save to same Excel file ---
  wb <- loadWorkbook(output_path)
  
  if ("CentralPoints_n5000"  %in% names(wb)) removeWorksheet(wb, "CentralPoints_n5000")
  if ("CentralPoints_n10000" %in% names(wb)) removeWorksheet(wb, "CentralPoints_n10000")
  
  addWorksheet(wb, "CentralPoints_n5000")
  addWorksheet(wb, "CentralPoints_n10000")
  
  writeData(wb, "CentralPoints_n5000",  df_5000)
  writeData(wb, "CentralPoints_n10000", df_10000)
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  
  cat("\n🎯 Central point results saved to:\n", output_path, "\n")
}
#======================================================================
# End of Central Point Extraction for Extended Entropies/Complexities
#======================================================================
#This script generates scatter plots of entropy vs. complexity for various ARMA models,
#using extended entropy measures (Shannon, Rényi, Tsallis, Fisher).
#======================================================================
# --- Required Packages ---
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(StatOrdPattHxC)

# --- File Paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/HC_Results_all_Models_extended_entropies.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

# --- Load LinfLsup boundaries (used ONLY for Shannon) ---
data("LinfLsup")

# --- Embedding Dimension ---
D <- 3  # You can adjust this if needed

# --- Cases ---
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

sample_sizes <- c(5000, 10000)

# --- Define Colors and Shapes ---
model_colors <- c(
  "ARMA11_M1"="#1b9e77","AR2_M1"="#d95f02","MA1_M2"="#7570b3","ARMA11_M3"="#e7298a",
  "AR2_M4"="#66a61e","MA1_M4"="#e6ab02","ARMA22_M2"="#a6761d","AR2_M3"="#666666",
  "MA2_M3"="#1f78b4","ARMA22_M1"="#b15928","ARMA22_M4"="#6a3d9a","MA2_M2"="#33a02c",
  "ARMA22_M3"="#fb9a99","AR1_M1"="#8dd3c7","AR1_M2"="#ff0000","ARMA11_M2"="#bebada",
  "MA1_M1"="#fb8072","MA2_M1"="#80b1d3","AR1_M3"="#fdb462","AR1_M4"="#b3de69",
  "ARMA11_M4"="#fccde5","MA1_M3"="#d9d9d9","MA2_M4"="#bc80bd","AR2_M2"="#ccebc5"
)

model_shapes <- c(
  "ARMA11_M1"=21,"AR2_M1"=22,"MA1_M2"=23,"ARMA11_M3"=24,"AR2_M4"=25,"MA1_M4"=8,
  "ARMA22_M2"=15,"AR2_M3"=16,"MA2_M3"=17,"ARMA22_M1"=18,"ARMA22_M4"=19,"MA2_M2"=4,
  "ARMA22_M3"=3,"AR1_M1"=1,"AR1_M2"=2,"ARMA11_M2"=5,"MA1_M1"=6,"MA2_M1"=7,
  "AR1_M3"=9,"AR1_M4"=10,"ARMA11_M4"=11,"MA1_M3"=12,"MA2_M4"=13,"AR2_M2"=14
)

# --- Entropy configurations ---
# col_H / col_C: column names in the data
# add_bounds: whether to add Linf–Lsup boundaries
entropy_configs <- list(
  Shannon = list(
    col_H = "H_Shannon",
    col_C = "C_Shannon",
    add_bounds = TRUE,   # <- ONLY Shannon has boundaries
    xlab = expression(italic(H)[S]),
    ylab = expression(italic(C)[S])
  ),
  Renyi = list(
    col_H = "H_Renyi",
    col_C = "C_Renyi",
    add_bounds = FALSE,  # <- NO boundaries
    xlab = expression(italic(H)[R]),
    ylab = expression(italic(C)[R])
  ),
  Tsallis = list(
    col_H = "H_Tsallis",
    col_C = "C_Tsallis",
    add_bounds = FALSE,  # <- NO boundaries
    xlab = expression(italic(H)[T]),
    ylab = expression(italic(C)[T])
  ),
  Fisher = list(
    col_H = "H_Fisher",
    col_C = "C_Fisher",
    add_bounds = FALSE,  # <- NO boundaries
    xlab = expression(italic(H)[F]),
    ylab = expression(italic(C)[F])
  )
)

# --- Main Loop ---
sheet_names <- excel_sheets(data_path)

for(case_name in names(cases)) {
  case_models   <- cases[[case_name]]
  case_plot_dir <- file.path(base_plot_dir, case_name, "Plots")
  dir.create(case_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(n_val in sample_sizes) {
    
    # Read all sheets for this case and sample size
    all_data <- lapply(case_models, function(model_name) {
      sheet_n <- paste0(model_name, "_n", n_val)
      if (sheet_n %in% sheet_names) {
        df       <- read_excel(data_path, sheet = sheet_n)
        df$Model <- model_name
        df$n     <- n_val
        return(df)
      } else {
        warning(paste("Skipping", sheet_n, "- sheet not found"))
        return(NULL)
      }
    }) %>% bind_rows()
    
    if (nrow(all_data) == 0) next
    
    # Subset Linf/Lsup once per D
    LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))
    
    #----------------------------------------------------------
    # Loop over entropy types: Shannon, Renyi, Tsallis, Fisher
    #----------------------------------------------------------
    for(ent_name in names(entropy_configs)) {
      cfg   <- entropy_configs[[ent_name]]
      H_col <- cfg$col_H
      C_col <- cfg$col_C
      
      # If this entropy is not present in the data, skip
      if (!all(c(H_col, C_col) %in% names(all_data))) {
        message("⚠️ Skipping ", ent_name, " for ", case_name, " n=", n_val,
                " — missing columns (", H_col, ", ", C_col, ").")
        next
      }
      
      # Remove rows with NA in these columns
      df_ent <- all_data %>%
        filter(is.finite(.data[[H_col]]), is.finite(.data[[C_col]]))
      
      if (nrow(df_ent) == 0) {
        message("⚠️ No finite data for ", ent_name, " (", case_name, ", n=", n_val, ").")
        next
      }
      
      # --- Optionally focus Linf/Lsup boundaries (Shannon only) ---
      if (cfg$add_bounds) {
        h_range <- range(df_ent[[H_col]], na.rm = TRUE)
        c_range <- range(df_ent[[C_col]], na.rm = TRUE)
        
        Linf_focus <- LinfLsup_subset %>%
          filter(H >= min(h_range) - 0.05, H <= max(h_range) + 0.05,
                 C >= min(c_range) - 0.05, C <= max(c_range) + 0.05)
      } else {
        Linf_focus <- NULL
      }
      
      # --- Scatter Plot for this entropy type ---
      p_scatter <- ggplot() +
        # Boundaries ONLY when add_bounds = TRUE (so only Shannon)
        {
          if (!is.null(Linf_focus)) {
            list(
              geom_line(
                data = subset(Linf_focus, Side == "Lower"),
                aes(H, C),
                linetype = "dashed", color = "gray40", linewidth = 0.6
              ),
              geom_line(
                data = subset(Linf_focus, Side == "Upper"),
                aes(H, C),
                linetype = "dashed", color = "gray40", linewidth = 0.6
              )
            )
          } else {
            NULL
          }
        } +
        geom_point(
          data = df_ent,
          aes(x = .data[[H_col]], y = .data[[C_col]],
              color = Model, shape = Model),
          size = 3
        ) +
        scale_color_manual(values = model_colors) +
        scale_shape_manual(values = model_shapes) +
        labs(
          title = paste0(ent_name, " Entropy–Complexity Scatter",
                         if (cfg$add_bounds) " with Boundaries " else " ",
                         "(", case_name, ", n=", n_val, ")"),
          x = cfg$xlab,
          y = cfg$ylab
        ) +
        theme_minimal(base_family = "serif", base_size = 13) +
        theme(
          legend.position = "bottom",
          legend.title    = element_blank(),
          panel.grid.minor = element_blank()
        )
      
      # --- Save Plot ---
      file_name <- paste0("Scatter_", ent_name, "_", case_name, "_n", n_val, ".pdf")
      ggsave(file.path(case_plot_dir, file_name),
             p_scatter, width = 8, height = 6)
      
      message("✅ ", ent_name, " scatter plot generated for ", case_name, " n=", n_val)
    } # end entropy loop
  }   # end n loop
}     # end case loop

message("🎉 All entropy–complexity scatter plots completed (LinfLsup only for Shannon)!")

#======================================================================
# End of Entropy–Complexity Scatter Plot Generation Script
#======================================================================
#This script generates plots of central points with confidence intervals
#for various ARMA models using extended entropy measures (Shannon, Rényi, Tsallis
#and Fisher).
#======================================================================
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(StatOrdPattHxC)
library(ggrepel)

# --- File Paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/HC_Results_all_Models_extended_entropies.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

# --- Define Cases ---
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

# --- Load LinfLsup boundaries (USED ONLY for Shannon) ---
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# --- Read and harmonize central points sheet ---
read_central <- function(sheet_name){
  df <- read_excel(data_path, sheet = sheet_name)
  names(df) <- trimws(names(df))
  df <- df %>% rename_all(~gsub("\\s+","",.))
  
  # Backwards compatibility for older Shannon names
  if ("H_Star" %in% names(df) && !"H_Shannon" %in% names(df)) df$H_Shannon <- df$H_Star
  if ("C_Star" %in% names(df) && !"C_Shannon" %in% names(df)) df$C_Shannon <- df$C_Star
  if ("SemiLength_H" %in% names(df) && !"SemiLength_H_Shannon" %in% names(df)) df$SemiLength_H_Shannon <- df$SemiLength_H
  if ("SemiLength_C" %in% names(df) && !"SemiLength_C_Shannon" %in% names(df)) df$SemiLength_C_Shannon <- df$SemiLength_C
  
  return(df)
}

# --- Entropy configurations ---
entropy_configs <- list(
  Shannon = list(
    H_col      = "H_Shannon",
    C_col      = "C_Shannon",
    semiH_col  = "SemiLength_H_Shannon",
    semiC_col  = "SemiLength_C_Shannon",
    add_bounds = TRUE,   # <-- ONLY SHANNON gets Linf-Lsup
    xlab       = expression(italic(H)[S]^"*"),
    ylab       = expression(italic(C)[S]^"*")
  ),
  Renyi = list(
    H_col      = "H_Renyi",
    C_col      = "C_Renyi",
    semiH_col  = "SemiLength_H_Renyi",
    semiC_col  = NULL,
    add_bounds = FALSE,
    xlab       = expression(italic(H)[R]^"*"),
    ylab       = expression(italic(C)[R]^"*")
  ),
  Tsallis = list(
    H_col      = "H_Tsallis",
    C_col      = "C_Tsallis",
    semiH_col  = "SemiLength_H_Tsallis",
    semiC_col  = NULL,
    add_bounds = FALSE,
    xlab       = expression(italic(H)[T]^"*"),
    ylab       = expression(italic(C)[T]^"*")
  ),
  Fisher = list(
    H_col      = "H_Fisher",
    C_col      = "C_Fisher",
    semiH_col  = "SemiLength_H_Fisher",
    semiC_col  = NULL,
    add_bounds = FALSE,  # <-- NEVER add boundaries for Fisher
    xlab       = expression(italic(H)[F]^"*"),
    ylab       = expression(italic(C)[F]^"*")
  )
)

# --- Loop over Sheets and Cases ---
sheets <- c("CentralPoints_n5000","CentralPoints_n10000")

for(sheet_name in sheets){
  
  n_val <- ifelse(grepl("5000", sheet_name), 5000, 10000)
  df <- read_central(sheet_name)
  
  for(case_name in names(cases)){
    
    case_models <- cases[[case_name]]
    df_case <- df %>% filter(Model %in% case_models)
    if(nrow(df_case) == 0) next
    
    plot_dir <- file.path(base_plot_dir, case_name, "Plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Loop through entropy types
    for(ent_name in names(entropy_configs)){
      
      cfg <- entropy_configs[[ent_name]]
      
      # Required columns
      needed <- c("Model", cfg$H_col, cfg$C_col, cfg$semiH_col)
      if (!all(needed %in% names(df_case))) {
        message("⚠️ Missing columns for ", ent_name, " in ", case_name, ", skipping.")
        next
      }
      
      df_ent <- df_case %>%
        mutate(
          H_val  = .data[[cfg$H_col]],
          C_val  = .data[[cfg$C_col]],
          Semi_H = .data[[cfg$semiH_col]],
          Semi_C = if (!is.null(cfg$semiC_col) && cfg$semiC_col %in% names(df_case))
            .data[[cfg$semiC_col]] else NA_real_
        ) %>% filter(is.finite(H_val), is.finite(C_val))
      
      # CI subset
      df_ci <- df_ent %>% filter(!is.na(Semi_H) & Semi_H > 0)
      
      # Shannon only: Linf–Lsup boundaries
      if (cfg$add_bounds) {
        h_range <- range(df_ent$H_val)
        c_range <- range(df_ent$C_val)
        Linf_focus <- LinfLsup_subset %>%
          filter(H >= h_range[1]-0.05, H <= h_range[2]+0.05,
                 C >= c_range[1]-0.05, C <= c_range[2]+0.05)
      } else {
        Linf_focus <- NULL
      }
      
      # Build plot
      p <- ggplot() +
        {
          if (!is.null(Linf_focus)) {
            list(
              geom_line(data=subset(Linf_focus, Side=="Lower"),
                        aes(H,C), linetype="dashed", color="gray40"),
              geom_line(data=subset(Linf_focus, Side=="Upper"),
                        aes(H,C), linetype="dashed", color="gray40")
            )
          }
        } +
        geom_point(data=df_ent, aes(H_val, C_val, color=Model, shape=Model), size=4) +
        geom_text_repel(data=df_ent,
                        aes(H_val, C_val, label=Model, color=Model),
                        size=3.5, max.overlaps=Inf) +
        geom_errorbarh(data=df_ci,
                       aes(y=C_val, xmin=H_val-Semi_H, xmax=H_val+Semi_H, color=Model),
                       height=0.001) +
        {
          if (!is.null(cfg$semiC_col)) {
            geom_errorbar(data=df_ci %>% filter(!is.na(Semi_C)&Semi_C>0),
                          aes(x=H_val, ymin=C_val-Semi_C, ymax=C_val+Semi_C, color=Model),
                          width=0.001)
          }
        } +
        scale_color_manual(values=model_colors) +
        scale_shape_manual(values=model_shapes) +
        labs(
          title=paste0(ent_name," Central Points (",case_name,", n=",n_val,")"),
          x=cfg$xlab, y=cfg$ylab
        ) +
        theme_minimal(base_size=13, base_family="serif") +
        theme(legend.position="bottom",
              legend.title=element_blank())
      
      # Save
      out_file <- paste0("CentralPoints_",ent_name,"_",case_name,"_n",n_val,".pdf")
      ggsave(file.path(plot_dir, out_file), p, width=8, height=6)
      
      message("✅ Saved ", ent_name," plot for ",case_name," n=",n_val)
    }
  }
}

message("🎉 All central point plots completed — LinfLsup ONLY for Shannon!")


#======================================================================
# End of Central Points Plot Generation Script
#======================================================================
#This script processes time series data for different ARMA model cases,
#extracts emblematic points, and generates plots with confidence intervals.
#======================================================================
# --- Required Packages ---
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(StatOrdPattHxC)
library(stringr)
library(ggrepel)

# --- Paths ---
ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/TimeSeries_by_Case_extended.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

# --- Cases (by sheet name) ---
cases <- list(
  Case1 = "Case1",
  Case2 = "Case2",
  Case3 = "Case3"
)

sample_sizes <- c(5000, 10000)

# --- Colors / Shapes (per model) ---
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

# --- LinfLsup: ONLY for Shannon H–C plane ---
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

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
  # 1. Load Data (by case sheet) and filter n
  # -------------------------------------------------------------
  df <- read_excel(ts_case_path, sheet = case_name) %>%
    filter(n == n_val) %>%
    mutate(Model = as.character(Model))
  
  if (nrow(df) == 0) {
    stop("❌ No data found for this case and n.")
  }
  
  # -------------------------------------------------------------
  # 1a. Harmonize names: old (H_Star) vs new (H_Shannon)
  # -------------------------------------------------------------
  if ("H_Star" %in% names(df) && !"H_Shannon" %in% names(df)) {
    df$H_Shannon <- df$H_Star
  }
  if ("C_Star" %in% names(df) && !"C_Shannon" %in% names(df)) {
    df$C_Shannon <- df$C_Star
  }
  if ("SemiLength_H" %in% names(df) && !"SemiLength_H_Shannon" %in% names(df)) {
    df$SemiLength_H_Shannon <- df$SemiLength_H
  }
  if ("SemiLength_C" %in% names(df) && !"SemiLength_C_Shannon" %in% names(df)) {
    df$SemiLength_C_Shannon <- df$SemiLength_C
  }
  
  required <- c("Model", "H_Shannon", "C_Shannon",
                "SemiLength_H_Shannon", "SemiLength_C_Shannon", "TimeSeries")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("❌ Missing required columns in sheet ", case_name, ": ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Clean TimeSeries (stored as character "c( ... )")
  clean_ts <- function(x) {
    x <- gsub("[c()]", "", x)
    x <- gsub("[[:space:]]+", "", x)
    as.numeric(strsplit(x, ",")[[1]])
  }
  
  ts_list_raw <- lapply(df$TimeSeries, clean_ts)
  names(ts_list_raw) <- df$Model
  
  # -------------------------------------------------------------
  # 2. Sort by Shannon H_Shannon ASCENDING (lowest → highest)
  # -------------------------------------------------------------
  df_sorted <- df %>% arrange(H_Shannon)
  ts_sorted <- ts_list_raw[df_sorted$Model]
  
  n_ts <- length(ts_sorted)
  message("Detected ", n_ts, " time series.")
  
  # -------------------------------------------------------------
  # 3. Pick TOP, MIDDLE, BOTTOM groups based on Shannon H
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
  # 4. Central Shannon H–C Scatter Plot (with LinfLsup + CI)
  # -------------------------------------------------------------
  # Focus Linf–Lsup region **for Shannon only**
  h_range <- range(df_sorted$H_Shannon, na.rm = TRUE)
  c_range <- range(df_sorted$C_Shannon, na.rm = TRUE)
  
  Linf_focus <- LinfLsup_subset %>%
    filter(H >= h_range[1] - 0.05,
           H <= h_range[2] + 0.05,
           C >= c_range[1] - 0.05,
           C <= c_range[2] + 0.05)
  
  p_center <- ggplot() +
    # Linf–Lsup boundaries (Shannon EC plane)
    geom_line(data = subset(Linf_focus, Side == "Lower"),
              aes(H, C), color = "gray40", linetype = "dashed") +
    geom_line(data = subset(Linf_focus, Side == "Upper"),
              aes(H, C), color = "gray40", linetype = "dashed") +
    
    # Central Shannon points
    geom_point(data = df_sorted,
               aes(H_Shannon, C_Shannon, color = Model, shape = Model),
               size = 3) +
    
    # CIs in H and C (Shannon only)
    geom_errorbarh(data = df_sorted,
                   aes(y = C_Shannon,
                       xmin = H_Shannon - SemiLength_H_Shannon,
                       xmax = H_Shannon + SemiLength_H_Shannon,
                       color = Model),
                   height = 0.002, linewidth = 0.6) +
    
    geom_errorbar(data = df_sorted,
                  aes(x = H_Shannon,
                      ymin = C_Shannon - SemiLength_C_Shannon,
                      ymax = C_Shannon + SemiLength_C_Shannon,
                      color = Model),
                  width = 0.002, linewidth = 0.6) +
    
    # labels next to points
    geom_text_repel(
      data = df_sorted,
      aes(H_Shannon, C_Shannon, label = Model, color = Model),
      family = "serif", size = 4, box.padding = 0.4
    ) +
    
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    
    labs(title = paste0("Central Shannon H–C Scatter — ", case_name, " (n=", n_val, ")"),
         x = expression(italic(H)[S]),
         y = expression(italic(C)[S])) +
    
    theme_minimal(base_family = "serif", base_size = 14) +
    theme(legend.position = "none")
  
  # -------------------------------------------------------------
  # 5. Time-Series plots — ASCENDING ORDER split into groups
  # -------------------------------------------------------------
  make_ts_plot <- function(ts_data, name) {
    df_ts <- data.frame(x = seq_along(ts_data), y = ts_data)
    
    ggplot(df_ts, aes(x = x, y = y)) +
      geom_line(color = model_colors[name], linewidth = 0.5) +
      scale_x_continuous(limits = c(min(df_ts$x), max(df_ts$x))) +
      scale_y_continuous(limits = c(min(df_ts$y), max(df_ts$y))) +
      labs(title = name) +
      theme_minimal(base_family = "serif", base_size = 8) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin  = margin(2, 4, 2, 2)
      )
  }
  
  ts_top_plots    <- mapply(make_ts_plot, ts_top,    names(ts_top),    SIMPLIFY = FALSE)
  ts_bottom_plots <- mapply(make_ts_plot, ts_bottom, names(ts_bottom), SIMPLIFY = FALSE)
  if (!is.null(right_idx)) {
    ts_right_plots <- mapply(make_ts_plot, ts_right, names(ts_right), SIMPLIFY = FALSE)
  }
  
  # -------------------------------------------------------------
  # 6. Arrange layout 3–3–3 (or 3–center–3)
  # -------------------------------------------------------------
  if (n_ts == 9) {
    combined <- grid.arrange(
      arrangeGrob(grobs = ts_top_plots,    ncol = 3),      # top
      arrangeGrob(
        p_center,
        arrangeGrob(grobs = ts_right_plots, ncol = 1),     # right side
        ncol = 2, widths = c(3, 1.2)
      ),
      arrangeGrob(grobs = ts_bottom_plots, ncol = 3),      # bottom
      ncol = 1,
      heights = c(0.9, 5, 0.9)
    )
  } else {
    combined <- grid.arrange(            # Case3: 6 TS only
      arrangeGrob(grobs = ts_top_plots, ncol = 3),
      p_center,
      arrangeGrob(grobs = ts_bottom_plots, ncol = 3),
      ncol = 1,
      heights = c(0.9, 5, 0.9)
    )
  }
  
  # -------------------------------------------------------------
  # 7. Save PDF
  # -------------------------------------------------------------
  out_file <- file.path(
    plot_dir,
    paste0("CentralPoints_TimeSeries_Shannon_", case_name, "_n", n_val, ".pdf")
  )
  
  ggsave(out_file, combined, width = 14, height = 12)
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

message("🎉 All plots successfully generated (Shannon H–C with LinfLsup only)!")
#======================================================================
# End of Time Series Plot Generation Script
#======================================================================
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
full_data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/HC_Results_all_Models_extended_entropies.xlsx"
output_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

###############################################################
# ⭐ Selected Models per Case
###############################################################
selected_models <- list(
  Case1 = c("ARMA11_M1", "AR1_M1", "AR1_M2"),
  Case2 = c("ARMA11_M3", "MA1_M3", "MA1_M4"),
  Case3 = c("ARMA22_M2", "ARMA22_M3", "MA2_M2")
)

sample_sizes <- c(5000, 10000)

###############################################################
# 🎨 Colors & Shapes (unchanged)
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
# 🔢 Entropies to plot (column names in each sheet)
###############################################################
# Assumes these columns exist in each model sheet:
#   H_Shannon, C_Shannon
#   H_Renyi,   C_Renyi
#   H_Tsallis, C_Tsallis
#   H_Fisher,  C_Fisher
entropy_configs <- list(
  Shannon = list(
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    add_bounds = TRUE   # use Linf–Lsup feasible region
  ),
  Renyi = list(
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    add_bounds = FALSE
  ),
  Tsallis = list(
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    add_bounds = FALSE
  ),
  Fisher = list(
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    add_bounds = FALSE   # no feasible region
  )
)

###############################################################
# 🔺 Feasible region (for Shannon)
###############################################################
data("LinfLsup")
LinfLsup_D3 <- subset(LinfLsup, Dimension == "3")

###############################################################
# ⭐ FUNCTION — Plot Selected Points for Each Case × n × entropy
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
  
  # --- Loop over entropies ---
  for (ent_name in names(entropy_configs)) {
    
    cfg <- entropy_configs[[ent_name]]
    H_col <- cfg$H_col
    C_col <- cfg$C_col
    
    # Skip if columns not present
    if (!(H_col %in% names(df)) || !(C_col %in% names(df))) {
      message("  ⚠️ Columns ", H_col, " or ", C_col,
              " not found for ", ent_name, " — skipping.")
      next
    }
    
    # Create a simple H, C version for plotting
    df_ent <- df %>%
      mutate(
        H = .data[[H_col]],
        C = .data[[C_col]]
      ) %>%
      filter(is.finite(H), is.finite(C))
    
    if (nrow(df_ent) == 0) {
      message("  ⚠️ No finite values for ", ent_name,
              " in ", case_name, " (n=", n_val, "), skipping.")
      next
    }
    
    # --- Feasible region crop (only for Shannon) ---
    if (isTRUE(cfg$add_bounds)) {
      Linf_focus <- LinfLsup_D3 %>%
        filter(
          H >= min(df_ent$H) - 0.05,
          H <= max(df_ent$H) + 0.05,
          C >= min(df_ent$C) - 0.05,
          C <= max(df_ent$C) + 0.05
        )
    } else {
      Linf_focus <- NULL
    }
    
    # --- Scatter plot (legend ON, no labels inside plot) ---
    p <- ggplot() +
      {
        if (!is.null(Linf_focus)) {
          list(
            geom_line(data = subset(Linf_focus, Side=="Lower"),
                      aes(H, C), color="gray40"),
            geom_line(data = subset(Linf_focus, Side=="Upper"),
                      aes(H, C), color="gray40")
          )
        }
      } +
      geom_point(
        data = df_ent,
        aes(H, C, color=Model, shape=Model),
        size=4
      ) +
      scale_color_manual(values=model_colors) +
      scale_shape_manual(values=model_shapes) +
      labs(
        title = paste0("Selected Models — ", case_name,
                       " (", ent_name, ", n=", n_val, ")"),
        x = expression(italic(H)),
        y = expression(italic(C)),
        color = "Model",
        shape = "Model"
      ) +
      theme_minimal(base_family="serif", base_size=14)
    
    # Save output
    outfile <- file.path(
      output_dir,
      paste0("Scatter_Selected_", case_name, "_", ent_name, "_n", n_val, ".pdf")
    )
    
    ggsave(outfile, p, width = 8, height = 6)
    message("  ✔ Saved: ", outfile)
  }
}

###############################################################
# ⭐ RUN FOR ALL CASES × SAMPLE SIZES
###############################################################
for (case_name in names(selected_models)) {
  for (n_val in sample_sizes) {
    plot_selected_scatter(case_name, n_val)
  }
}

message("\n🎉 All scatter plots for selected cases and all entropies created successfully!\n")
#----------------------------------------------------------------------
# End of Selected Models Scatter Plot Generation Script
#======================================================================
#This script processes time series data for different ARMA model cases,
#extracts emblematic points, and generates plots with confidence intervals.
#======================================================================
# --- Required Packages ---

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
ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/TimeSeries_by_Case_extended.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

cases <- list(
  Case1 = "Case1",
  Case2 = "Case2",
  Case3 = "Case3"
)

sample_sizes <- c(5000, 10000)

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
# 🔧 Helper: clean time series
###############################################################
clean_ts <- function(x) {
  x <- gsub("[c()]", "", x)
  x <- gsub("[[:space:]]+", "", x)
  as.numeric(strsplit(x, ",")[[1]])
}

###############################################################
# 🔧 Entropy configurations
###############################################################
# For each entropy, define column names and whether to add LinfLsup
entropy_configs <- list(
  Shannon = list(
    H_col     = "H_Shannon",
    C_col     = "C_Shannon",
    semiH_col = "SemiLength_H_Shannon",
    semiC_col = "SemiLength_C_Shannon",
    add_bounds = TRUE,
    xlab = expression(italic(H)[S]),
    ylab = expression(italic(C)[S])
  ),
  Renyi = list(
    H_col     = "H_Renyi",
    C_col     = "C_Renyi",
    semiH_col = "SemiLength_H_Renyi",
    semiC_col = NULL,  # no C variance for Renyi
    add_bounds = FALSE,
    xlab = expression(italic(H)[R]),
    ylab = expression(italic(C)[R])
  ),
  Tsallis = list(
    H_col     = "H_Tsallis",
    C_col     = "C_Tsallis",
    semiH_col = "SemiLength_H_Tsallis",
    semiC_col = NULL,
    add_bounds = FALSE,
    xlab = expression(italic(H)[T]),
    ylab = expression(italic(C)[T])
  ),
  Fisher = list(
    H_col     = "H_Fisher",
    C_col     = "C_Fisher",
    semiH_col = "SemiLength_H_Fisher",
    semiC_col = NULL,
    add_bounds = FALSE,  # NEVER LinfLsup for Fisher
    xlab = expression(italic(H)[F]),
    ylab = expression(italic(C)[F])
  )
)

# LinfLsup for Shannon only
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

###############################################################
# 🔧 Helper: harmonise Shannon column names (H_Star → H_Shannon, etc.)
###############################################################
harmonise_shannon_names <- function(df) {
  if ("H_Star" %in% names(df) && !"H_Shannon" %in% names(df)) {
    df$H_Shannon <- df$H_Star
  }
  if ("C_Star" %in% names(df) && !"C_Shannon" %in% names(df)) {
    df$C_Shannon <- df$C_Star
  }
  if ("SemiLength_H" %in% names(df) && !"SemiLength_H_Shannon" %in% names(df)) {
    df$SemiLength_H_Shannon <- df$SemiLength_H
  }
  if ("SemiLength_C" %in% names(df) && !"SemiLength_C_Shannon" %in% names(df)) {
    df$SemiLength_C_Shannon <- df$SemiLength_C
  }
  df
}

###############################################################
# ⭐ 1) 3 selected points + TS panels, for a given ENTROPY
###############################################################
plot_three_points_with_ts_entropy <- function(case_name, n_val, ent_name) {
  
  cfg <- entropy_configs[[ent_name]]
  
  message(">>> [", ent_name, "] 3-point plot for ", case_name, " n=", n_val)
  
  plot_dir <- file.path(base_plot_dir, case_name, "Plots")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  df <- read_excel(ts_case_path, sheet = case_name) %>%
    filter(n == n_val) %>%
    mutate(Model = as.character(Model))
  
  if (nrow(df) < 3) {
    warning("Not enough models for ", case_name, " n=", n_val)
    return()
  }
  
  # Harmonise Shannon names if needed
  df <- harmonise_shannon_names(df)
  
  # Check entropy-specific columns exist
  needed <- c("Model", "TimeSeries", cfg$H_col, cfg$C_col, cfg$semiH_col)
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    message("⚠️ Skipping ", ent_name, " for ", case_name, " n=", n_val,
            " — missing columns: ", paste(missing_cols, collapse = ", "))
    return()
  }
  
  # Build entropy-specific working columns
  df <- df %>%
    mutate(
      H_val    = .data[[cfg$H_col]],
      C_val    = .data[[cfg$C_col]],
      Semi_H   = .data[[cfg$semiH_col]],
      Semi_C   = if (!is.null(cfg$semiC_col) && cfg$semiC_col %in% names(df))
        .data[[cfg$semiC_col]] else NA_real_
    ) %>%
    filter(is.finite(H_val), is.finite(C_val))
  
  if (nrow(df) < 3) {
    message("⚠️ Not enough finite points for ", ent_name,
            " in ", case_name, " n=", n_val)
    return()
  }
  
  # Sort ascending by entropy H
  df_sorted <- df %>% arrange(H_val)
  
  idx_top    <- 1
  idx_mid    <- ceiling(nrow(df_sorted) / 2)
  idx_bottom <- nrow(df_sorted)
  
  df_sel <- df_sorted[c(idx_top, idx_mid, idx_bottom), ]
  df_sel$Position <- factor(c("Top","Middle","Bottom"), levels = c("Top","Middle","Bottom"))
  
  # Time series for the 3 selected models
  ts_all <- lapply(df_sel$TimeSeries, clean_ts)
  names(ts_all) <- df_sel$Model
  
  # --- LinfLsup (ONLY when cfg$add_bounds == TRUE, i.e. Shannon) ---
  if (cfg$add_bounds) {
    h_range <- range(df_sel$H_val, na.rm = TRUE)
    c_range <- range(df_sel$C_val, na.rm = TRUE)
    
    Linf_focus <- LinfLsup_subset %>%
      filter(
        H >= h_range[1] - 0.05,
        H <= h_range[2] + 0.05,
        C >= c_range[1] - 0.05,
        C <= c_range[2] + 0.05
      )
  } else {
    Linf_focus <- NULL
  }
  
  # --- Central H–C plot for this entropy ---
  p_center <- ggplot() +
    {
      if (!is.null(Linf_focus)) {
        list(
          geom_line(data = subset(Linf_focus, Side == "Lower"),
                    aes(H, C), color = "gray40", linetype = "dashed"),
          geom_line(data = subset(Linf_focus, Side == "Upper"),
                    aes(H, C), color = "gray40", linetype = "dashed")
        )
      } else {
        NULL
      }
    } +
    geom_point(data = df_sel,
               aes(H_val, C_val, color = Model, shape = Model),
               size = 3) +
    geom_errorbarh(data = df_sel,
                   aes(y = C_val,
                       xmin = H_val - Semi_H,
                       xmax = H_val + Semi_H,
                       color = Model),
                   height = 0.002, linewidth = 0.6) +
    {
      if (!all(is.na(df_sel$Semi_C))) {
        geom_errorbar(
          data = df_sel %>% filter(!is.na(Semi_C) & Semi_C > 0),
          aes(x = H_val,
              ymin = C_val - Semi_C,
              ymax = C_val + Semi_C,
              color = Model),
          width = 0.002, linewidth = 0.6
        )
      } else NULL
    } +
    geom_text_repel(
      data = df_sel,
      aes(H_val, C_val, label = Model, color = Model),
      family = "serif", size = 4
    ) +
    scale_color_manual(values = model_colors) +
    scale_shape_manual(values = model_shapes) +
    labs(
      title = paste0(ent_name, " – Selected Points (", case_name, ", n=", n_val, ")"),
      x = cfg$xlab,
      y = cfg$ylab
    ) +
    theme_minimal(base_family = "serif", base_size = 14) +
    theme(legend.position = "none")
  
  # --- TS panels for the 3 selected models ---
  make_ts_plot <- function(ts_data, model_name, pos_label) {
    df_ts <- data.frame(x = seq_along(ts_data), y = ts_data)
    ggplot(df_ts, aes(x, y)) +
      geom_line(color = model_colors[model_name], linewidth = 0.5) +
      scale_x_continuous(limits = c(min(df_ts$x), max(df_ts$x))) +
      scale_y_continuous(limits = c(min(df_ts$y), max(df_ts$y))) +
      labs(title = paste0(pos_label, ": ", model_name)) +
      theme_minimal(base_family = "serif", base_size = 8) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
      )
  }
  
  ts_plots <- mapply(make_ts_plot, ts_all, names(ts_all), df_sel$Position, SIMPLIFY = FALSE)
  ts_panel <- arrangeGrob(grobs = ts_plots, ncol = 1)
  
  combined <- grid.arrange(p_center, ts_panel, ncol = 2, widths = c(3, 1.7))
  
  out_file <- file.path(
    plot_dir,
    paste0("centralpoint_timeseries_", ent_name, "_", case_name, "_n", n_val, ".pdf")
  )
  ggsave(out_file, combined, width = 14, height = 10)
  
  message("✔ Saved: ", out_file)
}

###############################################################
# ⭐ 2) Combined scatter plot (3 points × 3 cases = 9) per ENTROPY
###############################################################
plot_combined_threepoint_scatter_entropy <- function(ent_name) {
  
  cfg <- entropy_configs[[ent_name]]
  message(">>> [", ent_name, "] Creating combined 9-point scatter")
  
  combined_points <- list()
  
  for (case_name in names(cases)) {
    df <- read_excel(ts_case_path, sheet = case_name) %>%
      mutate(Model = as.character(Model))
    
    df <- harmonise_shannon_names(df)
    
    needed <- c("Model", cfg$H_col, cfg$C_col, cfg$semiH_col)
    missing <- setdiff(needed, names(df))
    if (length(missing) > 0) {
      message("⚠️ Skipping ", case_name, " for ", ent_name,
              " — missing: ", paste(missing, collapse = ", "))
      next
    }
    
    df <- df %>%
      mutate(
        H_val  = .data[[cfg$H_col]],
        C_val  = .data[[cfg$C_col]],
        Semi_H = .data[[cfg$semiH_col]],
        Semi_C = if (!is.null(cfg$semiC_col) && cfg$semiC_col %in% names(df))
          .data[[cfg$semiC_col]] else NA_real_
      ) %>%
      filter(is.finite(H_val), is.finite(C_val))
    
    if (nrow(df) < 3) next
    
    df_sorted <- df %>% arrange(H_val)
    
    idx_top    <- 1
    idx_mid    <- ceiling(nrow(df_sorted) / 2)
    idx_bottom <- nrow(df_sorted)
    
    sel <- df_sorted[c(idx_top, idx_mid, idx_bottom), ]
    
    sel$CoefClass <- dplyr::case_when(
      case_name == "Case1" ~ "Positive",
      case_name == "Case2" ~ "Negative",
      case_name == "Case3" ~ "Mixed",
      TRUE ~ "Other"
    )
    sel$Case <- case_name
    
    combined_points[[case_name]] <- sel
  }
  
  df_all_sel <- bind_rows(combined_points)
  if (nrow(df_all_sel) == 0) {
    message("⚠️ No points for combined scatter (", ent_name, ").")
    return()
  }
  
  # LinfLsup only for Shannon combined scatter
  if (cfg$add_bounds) {
    h_range <- range(df_all_sel$H_val, na.rm = TRUE)
    c_range <- range(df_all_sel$C_val, na.rm = TRUE)
    
    Linf_focus <- LinfLsup_subset %>%
      filter(
        H >= h_range[1] - 0.05,
        H <= h_range[2] + 0.05,
        C >= c_range[1] - 0.05,
        C <= c_range[2] + 0.05
      )
  } else {
    Linf_focus <- NULL
  }
  
  p <- ggplot() +
    {
      if (!is.null(Linf_focus)) {
        list(
          geom_line(data = subset(Linf_focus, Side == "Lower"),
                    aes(H, C), color = "gray40"),
          geom_line(data = subset(Linf_focus, Side == "Upper"),
                    aes(H, C), color = "gray40")
        )
      } else NULL
    } +
    geom_point(
      data = df_all_sel,
      aes(H_val, C_val, color = CoefClass, shape = Model),
      size = 4, alpha = 0.9
    ) +
    geom_errorbarh(
      data = df_all_sel,
      aes(y = C_val,
          xmin = H_val - Semi_H,
          xmax = H_val + Semi_H,
          color = CoefClass),
      height = 0.002, linewidth = 0.7
    ) +
    {
      if (!all(is.na(df_all_sel$Semi_C))) {
        geom_errorbar(
          data = df_all_sel %>% filter(!is.na(Semi_C) & Semi_C > 0),
          aes(x = H_val,
              ymin = C_val - Semi_C,
              ymax = C_val + Semi_C,
              color = CoefClass),
          width = 0.002, linewidth = 0.7
        )
      } else NULL
    } +
    geom_text_repel(
      data = df_all_sel,
      aes(H_val, C_val, label = Model, color = CoefClass),
      family = "serif", size = 4, box.padding = 0.4
    ) +
    scale_color_manual(values = c(
      "Positive" = "darkblue",
      "Negative" = "darkgreen",
      "Mixed"    = "orange",
      "Other"    = "grey40"
    )) +
    scale_shape_manual(values = model_shapes) +
    labs(
      title = paste0(ent_name, " – Combined H–C Scatter (3 points per Case)"),
      x = cfg$xlab,
      y = cfg$ylab,
      color = "Coefficient Sign",
      shape = "Model"
    ) +
    theme_minimal(base_family = "serif", base_size = 14)
  
  out_file <- file.path(base_plot_dir,
                        paste0("combined_scatter_", ent_name, ".pdf"))
  ggsave(out_file, p, width = 10, height = 7)
  message("✔ Saved combined scatter for ", ent_name, ": ", out_file)
}

###############################################################
# 🔁 RUN: 3-point TS layouts for ALL entropies / cases / n
###############################################################
for (ent_name in names(entropy_configs)) {
  for (case_name in names(cases)) {
    for (n_val in sample_sizes) {
      plot_three_points_with_ts_entropy(case_name, n_val, ent_name)
    }
  }
}

###############################################################
# 🔁 RUN: Combined 9-point scatter for ALL entropies
###############################################################
for (ent_name in names(entropy_configs)) {
  plot_combined_threepoint_scatter_entropy(ent_name)
}

message("🎉 All entropy-extended TS + scatter plots generated!")
###############################################################
#======================================================================
# Start of Ordinal Patterns Features Plot Generation Script
#======================================================================
# This script reads ordinal pattern features from Excel files,
# It gives heatmaps and scatter plots for different entropies and models.
#======================================================================
# --- Required Packages ---
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(StatOrdPattHxC)

#--------------------------------------------------------
# Paths
#--------------------------------------------------------
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Data/Ordinal_patterns_Features n5000_n10000/HC_Results_all_Models_extended_entropies.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

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
sample_sizes <- c(5000, 10000)

#--------------------------------------------------------
# Feasible region (Linf / Lsup) for D = 3
#   👉 Will be used ONLY for Shannon heatmaps
#--------------------------------------------------------
data("LinfLsup")
D <- 3
Linf_all <- subset(LinfLsup, Side == "Lower" & Dimension == as.character(D))
Lsup_all <- subset(LinfLsup, Side == "Upper" & Dimension == as.character(D))

#--------------------------------------------------------
# Entropy configurations
#   H_col / C_col: column names for that entropy
#   add_bounds: whether to add Linf–Lsup (only Shannon)
#--------------------------------------------------------
entropy_configs <- list(
  Shannon = list(
    H_col = "H_Shannon",
    C_col = "C_Shannon",
    add_bounds = TRUE,
    xlab = expression(italic(H)[S]),
    ylab = expression(italic(C)[S])
  ),
  Renyi = list(
    H_col = "H_Renyi",
    C_col = "C_Renyi",
    add_bounds = FALSE,
    xlab = expression(italic(H)[R]),
    ylab = expression(italic(C)[R])
  ),
  Tsallis = list(
    H_col = "H_Tsallis",
    C_col = "C_Tsallis",
    add_bounds = FALSE,
    xlab = expression(italic(H)[T]),
    ylab = expression(italic(C)[T])
  ),
  Fisher = list(
    H_col = "H_Fisher",
    C_col = "C_Fisher",
    add_bounds = FALSE,  # ⛔ No LinfLsup for Fisher
    xlab = expression(italic(H)[F]),
    ylab = expression(italic(C)[F])
  )
)

#--------------------------------------------------------
# Main loops: cases × sample sizes × entropies
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
      df <- read_excel(data_path, sheet = sheet_name)
      df$Model <- m
      df$n     <- n_val
      df
    })
    
    df_all <- bind_rows(df_list)
    
    # If completely empty, skip
    if (nrow(df_all) == 0) {
      warning("No data for ", case_name, " n = ", n_val, ", skipping all entropies.")
      next
    }
    
    #--------------------------------------------------------
    # Entropy loop: Shannon, Renyi, Tsallis, Fisher
    #--------------------------------------------------------
    for (ent_name in names(entropy_configs)) {
      
      cfg <- entropy_configs[[ent_name]]
      H_col <- cfg$H_col
      C_col <- cfg$C_col
      
      message("  --- Entropy: ", ent_name, " ---")
      
      # Check required columns exist for this entropy
      if (!all(c(H_col, C_col) %in% names(df_all))) {
        message("  ⚠️ Skipping ", ent_name, " — missing columns (", H_col, ", ", C_col, ").")
        next
      }
      
      #--------------------------------------------------------
      # Build entropy-specific H, C; classify Family & Type
      #--------------------------------------------------------
      df_case <- df_all |>
        mutate(
          H = .data[[H_col]],
          C = .data[[C_col]]
        ) |>
        filter(is.finite(H), is.finite(C))
      
      if (nrow(df_case) == 0) {
        message("  ⚠️ No finite H,C for ", ent_name, " in ", case_name, " n=", n_val)
        next
      }
      
      # Classify AR/MA/ARMA family and fine Type (per model)
      df_case <- df_case |>
        mutate(
          Family = case_when(
            str_starts(Model, "ARMA") ~ "ARMA",
            str_starts(Model, "AR")   ~ "AR",
            str_starts(Model, "MA")   ~ "MA",
            TRUE                      ~ "Other"
          ),
          Type = factor(Model, levels = models_case)
        )
      
      df_case$Family <- factor(df_case$Family,
                               levels = c("AR","MA","ARMA","Other"))
      
      #--------------------------------------------------------
      # Feasible region cropping (ONLY if cfg$add_bounds == TRUE)
      #--------------------------------------------------------
      H_min <- min(df_case$H)
      H_max <- max(df_case$H)
      C_min <- min(df_case$C)
      C_max <- max(df_case$C)
      
      if (cfg$add_bounds) {
        Linf_crop <- Linf_all |>
          dplyr::filter(H >= H_min, H <= H_max)
        Lsup_crop <- Lsup_all |>
          dplyr::filter(H >= H_min, H <= H_max)
      } else {
        Linf_crop <- NULL
        Lsup_crop <- NULL
      }
      
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
        
        # Feasible region boundaries (ONLY for Shannon)
        {
          if (!is.null(Linf_crop) && !is.null(Lsup_crop)) {
            list(
              geom_line(data = Linf_crop, aes(H, C),
                        colour = "black", linewidth = 0.8),
              geom_line(data = Lsup_crop, aes(H, C),
                        colour = "black", linewidth = 0.8)
            )
          } else NULL
        } +
        
        # Background points
        geom_point(alpha = 0.07, size = 0.8, colour = "black") +
        
        # Density polygons per model Type
        stat_density_2d(
          aes(fill = Type, alpha = after_stat(level)),
          geom = "polygon",
          contour = TRUE,
          bins = 12,
          colour = NA,
          h = h_vec
        ) +
        
        # Contour lines per model Type
        stat_density_2d(
          aes(color = Type),
          contour = TRUE,
          bins = 8,
          linewidth = 0.25,
          h = h_vec
        ) +
        
        scale_fill_manual(values = type_cols, name = "Model (Type)") +
        scale_color_manual(values = type_cols, guide = "none") +
        scale_alpha(range = c(0.10, 0.55), guide = "none") +
        
        labs(
          title = paste("Heatmap by Model Type —", ent_name, case_name, "(n =", n_val, ")"),
          x = cfg$xlab,
          y = cfg$ylab
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
      case_output_dir <- file.path(base_plot_dir, case_name, "Plots")
      dir.create(case_output_dir, recursive = TRUE, showWarnings = FALSE)
      
      output_plot_path <- file.path(
        case_output_dir,
        paste0("HC_Heatmap_", ent_name, "_ModelTypes_",
               case_name, "_n", n_val, ".pdf")
      )
      
      ggsave(output_plot_path, p, width = 10, height = 8)
      message("  Saved: ", output_plot_path, "\n")
    } # end entropy loop
  }   # end n loop
}     # end case loop

message("🎉 All heatmaps completed for all cases, sample sizes, and entropies!")
#---------------------------------------------------------------------------------
# End of HC Heatmap Generation Script
#---------------------------------------------------------------------------------
#This script reads ordinal pattern features from Excel files,
#and generates scatter plots with error bars for different entropies and models. All the model
#configurations, color schemes, and plot settings are defined within the script.
#---------------------------------------------------------------------------------
library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)
library(StatOrdPattHxC)
library(ggrepel)

# --- File Paths ---
data_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/Ordinal_patterns_Features n5000_n10000/HC_Results_all_Models_extended_entropies.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

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

# --- Load LinfLsup boundaries (USED ONLY for Shannon) ---
data("LinfLsup")
D <- 3
LinfLsup_subset <- subset(LinfLsup, Dimension == as.character(D))

# --- Read and harmonize central points sheet ---
read_central <- function(sheet_name){
  df <- read_excel(data_path, sheet = sheet_name)
  names(df) <- trimws(names(df))
  df <- df %>% rename_all(~gsub("\\s+","",.))
  
  # Backwards compatibility for older Shannon names
  if ("H_Star" %in% names(df) && !"H_Shannon" %in% names(df)) df$H_Shannon <- df$H_Star
  if ("C_Star" %in% names(df) && !"C_Shannon" %in% names(df)) df$C_Shannon <- df$C_Star
  if ("SemiLength_H" %in% names(df) && !"SemiLength_H_Shannon" %in% names(df)) df$SemiLength_H_Shannon <- df$SemiLength_H
  if ("SemiLength_C" %in% names(df) && !"SemiLength_C_Shannon" %in% names(df)) df$SemiLength_C_Shannon <- df$SemiLength_C
  
  return(df)
}

# --- Entropy configurations ---
entropy_configs <- list(
  Shannon = list(
    H_col      = "H_Shannon",
    C_col      = "C_Shannon",
    semiH_col  = "SemiLength_H_Shannon",
    semiC_col  = "SemiLength_C_Shannon",
    add_bounds = TRUE,   # <-- ONLY SHANNON gets Linf-Lsup
    xlab       = expression(italic(H)[S]^"*"),
    ylab       = expression(italic(C)[S]^"*")
  ),
  Renyi = list(
    H_col      = "H_Renyi",
    C_col      = "C_Renyi",
    semiH_col  = "SemiLength_H_Renyi",
    semiC_col  = NULL,
    add_bounds = FALSE,
    xlab       = expression(italic(H)[R]^"*"),
    ylab       = expression(italic(C)[R]^"*")
  ),
  Tsallis = list(
    H_col      = "H_Tsallis",
    C_col      = "C_Tsallis",
    semiH_col  = "SemiLength_H_Tsallis",
    semiC_col  = NULL,
    add_bounds = FALSE,
    xlab       = expression(italic(H)[T]^"*"),
    ylab       = expression(italic(C)[T]^"*")
  ),
  Fisher = list(
    H_col      = "H_Fisher",
    C_col      = "C_Fisher",
    semiH_col  = "SemiLength_H_Fisher",
    semiC_col  = NULL,
    add_bounds = FALSE,  # <-- NEVER add boundaries for Fisher
    xlab       = expression(italic(H)[F]^"*"),
    ylab       = expression(italic(C)[F]^"*")
  )
)

# --- Sheets: central points for each n ---
sheets <- c("CentralPoints_n5000","CentralPoints_n10000")

for(sheet_name in sheets){
  
  n_val <- ifelse(grepl("5000", sheet_name), 5000, 10000)
  message(">>> Processing sheet: ", sheet_name, " (n = ", n_val, ")")
  
  df <- read_central(sheet_name)
  
  # Use ALL models present in the sheet
  df_all <- df %>% filter(!is.na(Model))
  if (nrow(df_all) == 0) {
    message("⚠️ No models found in sheet ", sheet_name, ", skipping.")
    next
  }
  
  # Output directory for "AllModels"
  plot_dir <- file.path(base_plot_dir, "AllModels", "Plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Loop through entropy types
  for(ent_name in names(entropy_configs)){
    
    cfg <- entropy_configs[[ent_name]]
    
    # Required columns
    needed <- c("Model", cfg$H_col, cfg$C_col, cfg$semiH_col)
    if (!all(needed %in% names(df_all))) {
      message("⚠️ Missing columns for ", ent_name, " in sheet ", sheet_name, ", skipping.")
      next
    }
    
    df_ent <- df_all %>%
      mutate(
        H_val  = .data[[cfg$H_col]],
        C_val  = .data[[cfg$C_col]],
        Semi_H = .data[[cfg$semiH_col]],
        Semi_C = if (!is.null(cfg$semiC_col) && cfg$semiC_col %in% names(df_all))
          .data[[cfg$semiC_col]] else NA_real_
      ) %>%
      filter(is.finite(H_val), is.finite(C_val))
    
    if (nrow(df_ent) == 0) {
      message("⚠️ No finite H/C for ", ent_name, " in sheet ", sheet_name, ", skipping.")
      next
    }
    
    # Rows with H-CI
    df_ci <- df_ent %>% filter(!is.na(Semi_H) & Semi_H > 0)
    
    # Shannon only: Linf–Lsup boundaries
    if (cfg$add_bounds) {
      h_range <- range(df_ent$H_val)
      c_range <- range(df_ent$C_val)
      Linf_focus <- LinfLsup_subset %>%
        filter(H >= h_range[1]-0.05, H <= h_range[2]+0.05,
               C >= c_range[1]-0.05, C <= c_range[2]+0.05)
    } else {
      Linf_focus <- NULL
    }
    
    # --- Build extended scatter plot for ALL models ---
    p <- ggplot() +
      {
        if (!is.null(Linf_focus)) {
          list(
            geom_line(data=subset(Linf_focus, Side=="Lower"),
                      aes(H,C), linetype="dashed", color="gray40"),
            geom_line(data=subset(Linf_focus, Side=="Upper"),
                      aes(H,C), linetype="dashed", color="gray40")
          )
        }
      } +
      geom_point(
        data=df_ent,
        aes(H_val, C_val, color=Model, shape=Model),
        size=4
      ) +
      # Optional labels – comment out if it gets too cluttered
      geom_text_repel(
        data=df_ent,
        aes(H_val, C_val, label=Model, color=Model),
        size=3.5, max.overlaps = Inf
      ) +
      geom_errorbarh(
        data=df_ci,
        aes(y=C_val, xmin=H_val-Semi_H, xmax=H_val+Semi_H, color=Model),
        height=0.001
      ) +
      {
        if (!is.null(cfg$semiC_col)) {
          geom_errorbar(
            data=df_ci %>% filter(!is.na(Semi_C) & Semi_C > 0),
            aes(x=H_val, ymin=C_val-Semi_C, ymax=C_val+Semi_C, color=Model),
            width=0.001
          )
        }
      } +
      scale_color_manual(values=model_colors) +
      scale_shape_manual(values=model_shapes) +
      labs(
        title = paste0(ent_name," Central Points (All Models, n=", n_val, ")"),
        x = cfg$xlab,
        y = cfg$ylab
      ) +
      theme_minimal(base_size=13, base_family="serif") +
      theme(
        legend.position = "bottom",
        legend.title    = element_blank()
      )
    
    # --- Save ---
    out_file <- paste0("CentralPoints_", ent_name, "_AllModels_n", n_val, ".pdf")
    ggsave(file.path(plot_dir, out_file), p, width=8, height=6)
    
    message("✅ Saved ", ent_name, " plot for ALL MODELS, n = ", n_val)
  }
}

message("🎉 All extended scatter plots completed — LinfLsup ONLY for Shannon, ALL models included!")
#--------------------------------------------------------
# End of Ordinal Patterns Features Extended Scatter Plot Script

# This script reads time series data and ordinal pattern features from Excel files,
# and generates combined scatter plots with error bars and time series rings
# for different entropies and models. All the model configurations, color schemes,
# and plot settings are defined within the script.
#--------------------------------------------------------
###############################################################
# 📦 Libraries
###############################################################
library(readxl)
library(StatOrdPattHxC)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(grid)

###############################################################
# 📁 Paths and settings
###############################################################
ts_case_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Data/Ordinal_patterns_Features n5000_n10000/TimeSeries_by_Case_extended.xlsx"
base_plot_dir <- "C:/Users/UserA1/Documents/GitHub/Ordinal_PatternS_R/Plots/ARMA plots/Ordinal_patterns_Features Plots n5000_n10000"

cases        <- c("Case1", "Case2", "Case3")
sample_sizes <- c(5000, 10000)

if (!dir.exists(base_plot_dir)) dir.create(base_plot_dir, recursive = TRUE)

###############################################################
# 🎨 Colors & Shapes (your mapping)
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
# 🔧 Helpers
###############################################################

clean_ts <- function(x) {
  x <- gsub("[c()]", "", x)
  x <- gsub("[[:space:]]+", "", x)
  as.numeric(strsplit(x, ",")[[1]])
}

harmonise_shannon_names <- function(df) {
  nms <- names(df)
  if ("H_Star" %in% nms && !"H_Shannon" %in% nms) {
    df$H_Shannon <- df$H_Star
  }
  if ("C_Star" %in% nms && !"C_Shannon" %in% nms) {
    df$C_Shannon <- df$C_Star
  }
  if ("SemiLength_H" %in% nms && !"SemiLength_H_Shannon" %in% nms) {
    df$SemiLength_H_Shannon <- df$SemiLength_H
  }
  if ("SemiLength_C" %in% nms && !"SemiLength_C_Shannon" %in% nms) {
    df$SemiLength_C_Shannon <- df$SemiLength_C
  }
  df
}

sanitize_names <- function(df, max_len = 200) {
  nms <- names(df)
  too_long <- nchar(nms, type = "bytes") > max_len
  nms[too_long] <- substr(nms[too_long], 1, max_len)
  nms <- make.names(nms, unique = TRUE)
  names(df) <- nms
  df
}

###############################################################
# 🔧 Linf–Lsup for D = 3 (Shannon)
###############################################################
data("LinfLsup")
LinfLsup_D3 <- LinfLsup[LinfLsup$Dimension == "3", ]

###############################################################
# ⭐ MAIN FUNCTION
###############################################################
plot_combined_selected_with_ts <- function(n_val) {
  
  message(">>> Creating layout like example for n = ", n_val)
  combined_list <- list()
  
  # 1) select 3 (min/mid/max H) per case
  for (case_name in cases) {
    df <- read_excel(ts_case_path, sheet = case_name)
    df <- sanitize_names(df)
    df <- harmonise_shannon_names(df)
    
    if (!("n" %in% names(df))) stop("Column 'n' missing in ", case_name)
    
    df <- df[df$n == n_val, , drop = FALSE]
    if (nrow(df) < 3) next
    
    needed <- c("Model","TimeSeries",
                "H_Shannon","C_Shannon",
                "SemiLength_H_Shannon","SemiLength_C_Shannon")
    if (length(setdiff(needed, names(df))) > 0) next
    
    df$Model  <- as.character(df$Model)
    df$H_val  <- df$H_Shannon
    df$C_val  <- df$C_Shannon
    df$Semi_H <- df$SemiLength_H_Shannon
    df$Semi_C <- df$SemiLength_C_Shannon
    
    df <- df[is.finite(df$H_val) & is.finite(df$C_val), , drop = FALSE]
    if (nrow(df) < 3) next
    
    o <- order(df$H_val)
    df_sorted <- df[o, ]
    
    idx_top    <- 1
    idx_mid    <- ceiling(nrow(df_sorted) / 2)
    idx_bottom <- nrow(df_sorted)
    
    sel <- df_sorted[c(idx_top, idx_mid, idx_bottom), ]
    sel$Case <- case_name
    combined_list[[case_name]] <- sel
  }
  
  if (!length(combined_list)) {
    warning("No points for n = ", n_val)
    return(invisible(NULL))
  }
  
  df_all <- do.call(rbind, combined_list)
  rownames(df_all) <- NULL
  
  ## global ascending H
  df_all <- df_all[order(df_all$H_val), ]
  rownames(df_all) <- NULL
  
  df_all$Model <- factor(df_all$Model, levels = names(model_colors))
  
  ts_list <- lapply(df_all$TimeSeries, clean_ts)
  names(ts_list) <- as.character(df_all$Model)
  
  ###########################################################
  # Build TS ggplots in H-ascending order
  ###########################################################
  make_ts_plot <- function(ts_data, model_name, case_name) {
    df_ts <- data.frame(x = seq_along(ts_data), y = ts_data)
    ggplot(df_ts, aes(x, y)) +
      geom_line(color = model_colors[model_name], linewidth = 0.4) +
      labs(title = model_name, x = NULL, y = NULL) +
      theme_minimal(base_family = "serif", base_size = 8) +
      theme(
        plot.title   = element_text(hjust = 0.5, size = 8),
        axis.text.x  = element_text(size = 6),
        axis.text.y  = element_text(size = 6),
        panel.grid.major = element_line(linewidth = 0.1),
        panel.grid.minor = element_blank()
      )
  }
  
  ts_plots <- mapply(
    FUN        = make_ts_plot,
    ts_data    = ts_list,
    model_name = as.character(df_all$Model),
    case_name  = df_all$Case,
    SIMPLIFY   = FALSE
  )
  
  # pad up to 9 if needed
  if (length(ts_plots) < 9) {
    for (k in (length(ts_plots)+1):9) {
      ts_plots[[k]] <- ggplot() + theme_void()
    }
  }
  
  # assign: top = 1–3, right = 4–6, bottom = 7–9
  top_row    <- gridExtra::arrangeGrob(ts_plots[[1]], ts_plots[[2]], ts_plots[[3]],
                                       ncol = 3)
  right_col  <- gridExtra::arrangeGrob(ts_plots[[4]], ts_plots[[5]], ts_plots[[6]],
                                       ncol = 1)
  bottom_row <- gridExtra::arrangeGrob(ts_plots[[7]], ts_plots[[8]], ts_plots[[9]],
                                       ncol = 3)
  
  ###########################################################
  # Central H–C scatter
  ###########################################################
  h_range <- range(df_all$H_val, na.rm = TRUE)
  c_range <- range(df_all$C_val, na.rm = TRUE)
  
  Linf_focus <- LinfLsup_D3[
    LinfLsup_D3$H >= h_range[1] - 0.05 &
      LinfLsup_D3$H <= h_range[2] + 0.05 &
      LinfLsup_D3$C >= c_range[1] - 0.05 &
      LinfLsup_D3$C <= c_range[2] + 0.05,
  ]
  
  p_center <- ggplot() +
    geom_line(
      data = subset(Linf_focus, Side == "Lower"),
      aes(H, C),
      color = "grey40", linetype = "dashed", linewidth = 0.4
    ) +
    geom_line(
      data = subset(Linf_focus, Side == "Upper"),
      aes(H, C),
      color = "grey40", linetype = "dashed", linewidth = 0.4
    ) +
    geom_point(
      data = df_all,
      aes(H_val, C_val, color = Model, shape = Model),
      size = 3
    ) +
    geom_errorbarh(
      data = df_all,
      aes(y = C_val,
          xmin = H_val - Semi_H,
          xmax = H_val + Semi_H,
          color = Model),
      height = 0.002
    ) +
    geom_errorbar(
      data = df_all,
      aes(x = H_val,
          ymin = C_val - Semi_C,
          ymax = C_val + Semi_C,
          color = Model),
      width = 0.002
    ) +
    geom_text_repel(
      data = df_all,
      aes(H_val, C_val, label = Model, color = Model),
      size = 3,
      family = "serif",
      box.padding = 0.3,
      max.overlaps = Inf
    ) +
    scale_color_manual(values = model_colors, drop = FALSE) +
    scale_shape_manual(values = model_shapes, drop = FALSE) +
    labs(
      title = paste0("Central Shannon H–C Scatter — n = ", n_val),
      x = expression(italic(H)[S]),
      y = expression(italic(C)[S])
    ) +
    theme_minimal(base_family = "serif", base_size = 12) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 14),
      legend.position = "none"
    )
  
  # centre + right column (like your example)
  mid_row <- gridExtra::arrangeGrob(
    p_center, right_col,
    ncol = 2,
    widths = c(3.2, 1.3)   # adjust these to change relative sizes
  )
  
  # final layout: top / middle / bottom
  full_plot <- gridExtra::arrangeGrob(
    top_row,
    mid_row,
    bottom_row,
    ncol    = 1,
    heights = c(1.2, 3.5, 1.2)  # central band taller; tweak if needed
  )
  
  ###########################################################
  # Save
  ###########################################################
  out_file <- file.path(
    base_plot_dir,
    paste0("Shannon_layout_like_example_n", n_val, ".pdf")
  )
  ggsave(out_file, full_plot, width = 12, height = 9)
  message("✔ Saved layout for n = ", n_val, ": ", out_file)
}

###############################################################
# 🚀 RUN
###############################################################
for (n_val in sample_sizes) {
  plot_combined_selected_with_ts(n_val)
}

message("🎉 All figures created in “top / centre+right / bottom” style!")
