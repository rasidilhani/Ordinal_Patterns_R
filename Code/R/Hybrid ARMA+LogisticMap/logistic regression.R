#---------------------------------------------------------------------------------------------
# This script generates a feature analysis of the models. It includes correlation analysis and 
#prepares data for multinomial logistic regression. It uses the same features as the previous 
#script (H and C for Shannon, Rényi, Tsallis, Fisher, and Disequilibrium).
#----------------------------------------------------------------------------------------------
library(nnet)
library(caret)
library(corrplot)
library(readxl)
library(dplyr)
library(tidyr)

# ==============================================================================
# PARAMETERS AND PATHS
# ==============================================================================
D <- 4  # Embedding dimension

# Paths
hc_data_path <- here("Data", "Hybrid model_data", "D4_Data", "Combined_All_Models_HC_Results_D4.xlsx")
output_dir <- here("Hybrid_Analysis", "Multinomial_Logistic_Regression")
# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Sample sizes to process
sample_sizes <- c(1000, 5000)

# ==============================================================================
# DEFINE MODEL NAMES
# ==============================================================================
model_names <- c(
  "ARMA(2,2)", "ARMA(1,1)", "AR(2)", "AR(1)", "MA(2)", "MA(1)",
  "Logistic",
  "Hybrid_ARMA(2,2)", "Hybrid_ARMA(1,1)", "Hybrid_AR(2)", 
  "Hybrid_AR(1)", "Hybrid_MA(2)", "Hybrid_MA(1)", "Logistic_r3_6","Sine_Wave","Logistic_Sine_Combined"
)

# ==============================================================================
# LOAD DATA
# ==============================================================================
cat("Loading data...\n")
hc_data_all <- read_excel(hc_data_path, sheet = 1)

# Select relevant columns
hc_data <- hc_data_all %>%
  select(Model, n, rep,
         H_Shannon, C_Shannon,
         H_Renyi, C_Renyi,
         H_Tsallis, C_Tsallis,
         H_Fisher, C_Fisher,
         Disequilibrium)

# Rename for consistency
#names(hc_data)[names(hc_data) == disequil_col] <- "Disequilibrium"

# Ensure Model is factor with correct order
hc_data$Model <- factor(hc_data$Model, levels = model_names)
# ==============================================================================
# DEFINE FEATURES
# ==============================================================================
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "Disequilibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

print(features)

# ==============================================================================
# ANALYSIS FOR EACH SAMPLE SIZE
# ==============================================================================

for (n_val in sample_sizes) {
  
  cat(sprintf("\n\n================================================================================\n"))
  cat(sprintf("ANALYSIS FOR SAMPLE SIZE n = %d\n", n_val))
  cat(sprintf("================================================================================\n\n"))
  
  # Create subdirectory for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # PREPARE DATA
  # --------------------------------------------------------------------------
  
  analysis_data <- hc_data %>%
    filter(n == n_val) %>%
    select(Model, all_of(features)) %>%
    na.omit()  # Remove any rows with missing values
  
  cat(sprintf("Data after filtering and removing NA: %d rows\n", nrow(analysis_data)))
  cat(sprintf("Models in dataset: %d\n", length(unique(analysis_data$Model))))
  
  # Check class distribution
  cat("\nClass distribution:\n")
  print(table(analysis_data$Model))
  
  # --------------------------------------------------------------------------
  # CORRELATION ANALYSIS
  # --------------------------------------------------------------------------
  
  cat("\n### CORRELATION ANALYSIS ###\n")
  
  # Calculate correlation matrix
  cor_matrix <- cor(analysis_data[, features], use = "complete.obs")
  
  cat("\nCorrelation Matrix:\n")
  print(round(cor_matrix, 3))
  
  # Save correlation matrix as CSV
  write.csv(round(cor_matrix, 3),
            file = file.path(n_dir, paste0("Correlation_Matrix_n", n_val, ".csv")),
            row.names = TRUE)
  
  # Save correlation plot as PDF
  pdf(file.path(n_dir, paste0("Correlation_Matrix_n", n_val, ".pdf")),
      width = 10, height = 10)
  corrplot(cor_matrix, 
           method = "color", 
           addCoef.col = "black",
           tl.col = "black",
           tl.srt = 45,
           number.cex = 0.7,
           title = paste0("Feature Correlation Matrix (n = ", n_val, ")"),
           mar = c(0, 0, 2, 0))
  dev.off()
  
  cat(sprintf("\n✅ Correlation matrix saved as PDF\n"))
  
  # Identify highly correlated pairs (|r| > 0.8)
  high_cor <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
  if (nrow(high_cor) > 0) {
    cat("\nHighly correlated feature pairs (|r| > 0.8):\n")
    for (i in 1:nrow(high_cor)) {
      row_idx <- high_cor[i, 1]
      col_idx <- high_cor[i, 2]
      if (row_idx < col_idx) {  # Avoid duplicates
        cat(sprintf("  %s <-> %s: %.3f\n", 
                    rownames(cor_matrix)[row_idx],
                    colnames(cor_matrix)[col_idx],
                    cor_matrix[row_idx, col_idx]))
      }
    }
  } else {
    cat("\nNo highly correlated feature pairs (|r| > 0.8) found.\n")
  }
  
  # --------------------------------------------------------------------------
  # FIT MULTINOMIAL LOGISTIC REGRESSION MODEL
  # --------------------------------------------------------------------------
  
  cat("\n### MULTINOMIAL LOGISTIC REGRESSION MODEL ###\n")
  
  # Fit the model
  cat("\nFitting multinomial logistic regression model...\n")
  
  set.seed(42)  # For reproducibility
  
  # Choose the baseline for logistic regression
  analysis_data$Model <- relevel(analysis_data$Model, ref = "ARMA(2,2)")
  
  # Now fit the model
  multinom_model <- multinom(
    Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
      C_Shannon + C_Renyi + C_Tsallis + C_Fisher + Disequilibrium, 
    data = analysis_data,
    maxit = 500,
    trace = FALSE
  )
  
  cat("\nBaseline class:", levels(analysis_data$Model)[1], "\n")
  # --------------------------------------------------------------------------
  # MODEL SUMMARY AND STATISTICS
  # --------------------------------------------------------------------------
  
  # Capture model summary
  sink(file.path(n_dir, paste0("Model_Summary_n", n_val, ".txt")))
  cat("================================================================================\n")
  cat(sprintf("MULTINOMIAL LOGISTIC REGRESSION MODEL SUMMARY (n = %d)\n", n_val))
  cat("================================================================================\n\n")
  print(summary(multinom_model))
  cat("\n\n")
  sink()
  
  # Get AIC
  model_aic <- AIC(multinom_model)
  cat(sprintf("\nModel AIC: %.2f\n", model_aic))
  
  # --------------------------------------------------------------------------
  # MODEL PREDICTIONS AND ACCURACY
  # --------------------------------------------------------------------------
  
  # Make predictions on the same data (training accuracy)
  predictions <- predict(multinom_model, newdata = analysis_data, type = "class")
  
  # Create confusion matrix
  conf_matrix <- confusionMatrix(predictions, analysis_data$Model)
  
  # Extract overall accuracy
  overall_accuracy <- conf_matrix$overall['Accuracy']
  cat(sprintf("Overall Model Accuracy: %.4f (%.2f%%)\n", 
              overall_accuracy, overall_accuracy * 100))
  
  # Per-class accuracy (sensitivity)
  per_class_accuracy <- conf_matrix$byClass[, 'Sensitivity']
  
  cat("\nPer-Class Accuracy (Sensitivity):\n")
  print(round(per_class_accuracy, 4))
  
  # --------------------------------------------------------------------------
  # SAVE CONFUSION MATRIX
  # --------------------------------------------------------------------------
  
  # Save confusion matrix
  sink(file.path(n_dir, paste0("Confusion_Matrix_n", n_val, ".txt")))
  cat("================================================================================\n")
  cat(sprintf("CONFUSION MATRIX (n = %d)\n", n_val))
  cat("================================================================================\n\n")
  print(conf_matrix)
  sink()
  
  # Save confusion matrix table as CSV
  conf_table <- as.data.frame.matrix(table(Predicted = predictions, Actual = analysis_data$Model))
  write.csv(conf_table,
            file = file.path(n_dir, paste0("Confusion_Matrix_Table_n", n_val, ".csv")),
            row.names = TRUE)
  
  # --------------------------------------------------------------------------
  # SAVE MODEL COEFFICIENTS
  # --------------------------------------------------------------------------
  
  # Extract coefficients
  coefficients <- summary(multinom_model)$coefficients
  std_errors <- summary(multinom_model)$standard.errors
  
  # Calculate z-values and p-values
  z_values <- coefficients / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # Save coefficients
  write.csv(coefficients,
            file = file.path(n_dir, paste0("Model_Coefficients_n", n_val, ".csv")),
            row.names = TRUE)
  
  write.csv(std_errors,
            file = file.path(n_dir, paste0("Model_StdErrors_n", n_val, ".csv")),
            row.names = TRUE)
  
  write.csv(p_values,
            file = file.path(n_dir, paste0("Model_PValues_n", n_val, ".csv")),
            row.names = TRUE)
  
  # --------------------------------------------------------------------------
  # SAVE PERFORMANCE METRICS
  # --------------------------------------------------------------------------
  
  # Create performance summary
  performance_summary <- data.frame(
    Metric = c("Sample Size", "Number of Observations", "Number of Classes", 
               "AIC", "Overall Accuracy", "Kappa"),
    Value = c(n_val, 
              nrow(analysis_data),
              length(unique(analysis_data$Model)),
              round(model_aic, 2),
              round(overall_accuracy, 4),
              round(conf_matrix$overall['Kappa'], 4))
  )
  
  write.csv(performance_summary,
            file = file.path(n_dir, paste0("Performance_Summary_n", n_val, ".csv")),
            row.names = FALSE)
  
  cat("\n### PERFORMANCE SUMMARY ###\n")
  print(performance_summary)
  
  # --------------------------------------------------------------------------
  # SAVE PER-CLASS METRICS
  # --------------------------------------------------------------------------
  
  # Extract per-class metrics
  per_class_metrics <- data.frame(
    Model = rownames(conf_matrix$byClass),
    Sensitivity = conf_matrix$byClass[, 'Sensitivity'],
    Specificity = conf_matrix$byClass[, 'Specificity'],
    Pos_Pred_Value = conf_matrix$byClass[, 'Pos Pred Value'],
    Neg_Pred_Value = conf_matrix$byClass[, 'Neg Pred Value'],
    Precision = conf_matrix$byClass[, 'Precision'],
    Recall = conf_matrix$byClass[, 'Recall'],
    F1 = conf_matrix$byClass[, 'F1'],
    Balanced_Accuracy = conf_matrix$byClass[, 'Balanced Accuracy']
  )
  
  write.csv(per_class_metrics,
            file = file.path(n_dir, paste0("Per_Class_Metrics_n", n_val, ".csv")),
            row.names = FALSE)
  
  cat("\n✅ All results saved for n =", n_val, "\n")
  
  # --------------------------------------------------------------------------
  # SAVE R MODEL OBJECT
  # --------------------------------------------------------------------------
  
  saveRDS(multinom_model,
          file = file.path(n_dir, paste0("Multinom_Model_n", n_val, ".rds")))
  
  cat(sprintf("\n✅ Model object saved as RDS file\n"))
  
}
