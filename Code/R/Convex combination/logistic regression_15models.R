# ══════════════════════════════════════════════════════════════════════════════
#  Multinomial Logistic Regression — Convex Combination Models
# ══════════════════════════════════════════════════════════════════════════════
library(nnet)
library(caret)
library(corrplot)
library(readxl)
library(dplyr)
library(here)

# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
D            <- 4
sample_sizes <- c(1000, 10000)

# ══════════════════════════════════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════════════════════════════════
hc_data_path <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))
output_dir   <- here("Results", "Convex_combination", "Multinomial_Logistic")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL NAMES
# ══════════════════════════════════════════════════════════════════════════════
model_names <- c(
  "ARMA(2,2)", "AR(2)", "MA(2)", "Logistic", "Sine",
  "ARMA+Sine(w=0.1)", "ARMA+Sine(w=0.2)", "ARMA+Sine(w=0.3)",
  "AR2+Logistic(w=0.1)", "AR2+Sine(w=0.8)",
  "MA2+Logistic(w=0.2)", "MA2+Logistic(w=0.7)",
  "MA2+Sine(w=0.4)", "MA2+Sine(w=0.6)", "MA2+Sine(w=0.8)"
)

# ══════════════════════════════════════════════════════════════════════════════
#  FEATURES
# ══════════════════════════════════════════════════════════════════════════════
features <- c(
  "H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher",
  "C_Shannon", "C_Renyi", "C_Tsallis", "C_Fisher",
  "Disequilibrium"
)

# ══════════════════════════════════════════════════════════════════════════════
#  LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
hc_data <- bind_rows(
  read_xlsx(hc_data_path, sheet = "n1000"),
  read_xlsx(hc_data_path, sheet = "n10000")
) %>%
  select(Model, n, Rep, all_of(features)) %>%
  mutate(Model = factor(Model, levels = model_names))

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN LOOP
# ══════════════════════════════════════════════════════════════════════════════
for (n_val in sample_sizes) {
  
  # Output folder for this sample size
  n_dir <- file.path(output_dir, paste0("n", n_val))
  dir.create(n_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ── Prepare analysis data ─────────────────────────────────────────────────
  analysis_data <- hc_data %>%
    filter(n == n_val) %>%
    select(Model, all_of(features)) %>%
    na.omit() %>%
    mutate(Model = relevel(Model, ref = "ARMA(2,2)"))  # baseline class
  
  # ── Correlation matrix ────────────────────────────────────────────────────
  cor_matrix <- cor(analysis_data[, features], use = "complete.obs")
  
  write.csv(round(cor_matrix, 3),
            file.path(n_dir, paste0("Correlation_Matrix_n", n_val, ".csv")),
            row.names = TRUE)
  
  pdf(file.path(n_dir, paste0("Correlation_Matrix_n", n_val, ".pdf")),
      width = 10, height = 10)
  corrplot(cor_matrix,
           method      = "color",
           addCoef.col = "black",
           tl.col      = "black",
           tl.srt      = 45,
           number.cex  = 0.7,
           title       = paste0("Correlation Matrix  (n = ", n_val, ")"),
           mar         = c(0, 0, 2, 0))
  dev.off()
  
  # ── Fit multinomial logistic regression ───────────────────────────────────
  set.seed(42)
  
  model_fit <- multinom(
    Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher +
      C_Shannon + C_Renyi + C_Tsallis + C_Fisher + Disequilibrium,
    data  = analysis_data,
    maxit = 500,
    trace = FALSE
  )
  
  # ── Model summary → text file ─────────────────────────────────────────────
  sink(file.path(n_dir, paste0("Model_Summary_n", n_val, ".txt")))
  print(summary(model_fit))
  sink()
  
  # ── Predictions and confusion matrix ─────────────────────────────────────
  predictions <- predict(model_fit, newdata = analysis_data, type = "class")
  conf_matrix <- confusionMatrix(predictions, analysis_data$Model)
  
  sink(file.path(n_dir, paste0("Confusion_Matrix_n", n_val, ".txt")))
  print(conf_matrix)
  sink()
  
  write.csv(
    as.data.frame.matrix(table(Predicted = predictions,
                               Actual    = analysis_data$Model)),
    file.path(n_dir, paste0("Confusion_Matrix_Table_n", n_val, ".csv")),
    row.names = TRUE
  )
  
  # ── Coefficients, standard errors, p-values ───────────────────────────────
  coef_mat <- summary(model_fit)$coefficients
  se_mat   <- summary(model_fit)$standard.errors
  p_mat    <- 2 * (1 - pnorm(abs(coef_mat / se_mat)))
  
  write.csv(coef_mat, file.path(n_dir, paste0("Coefficients_n", n_val, ".csv")),
            row.names = TRUE)
  write.csv(se_mat,   file.path(n_dir, paste0("StdErrors_n",    n_val, ".csv")),
            row.names = TRUE)
  write.csv(p_mat,    file.path(n_dir, paste0("PValues_n",      n_val, ".csv")),
            row.names = TRUE)
  
  # ── Performance summary ───────────────────────────────────────────────────
  perf_df <- data.frame(
    Metric = c("Sample Size", "Observations", "Classes",
               "AIC", "Overall Accuracy", "Kappa"),
    Value  = c(n_val,
               nrow(analysis_data),
               length(unique(analysis_data$Model)),
               round(AIC(model_fit), 2),
               round(conf_matrix$overall["Accuracy"], 4),
               round(conf_matrix$overall["Kappa"],    4))
  )
  write.csv(perf_df,
            file.path(n_dir, paste0("Performance_Summary_n", n_val, ".csv")),
            row.names = FALSE)
  
  # ── Per-class metrics ─────────────────────────────────────────────────────
  per_class_df <- data.frame(
    Model             = rownames(conf_matrix$byClass),
    Sensitivity       = conf_matrix$byClass[, "Sensitivity"],
    Specificity       = conf_matrix$byClass[, "Specificity"],
    Precision         = conf_matrix$byClass[, "Precision"],
    Recall            = conf_matrix$byClass[, "Recall"],
    F1                = conf_matrix$byClass[, "F1"],
    Balanced_Accuracy = conf_matrix$byClass[, "Balanced Accuracy"]
  )
  write.csv(per_class_df,
            file.path(n_dir, paste0("Per_Class_Metrics_n", n_val, ".csv")),
            row.names = FALSE)
  
  # ── Save model object ─────────────────────────────────────────────────────
  saveRDS(model_fit,
          file.path(n_dir, paste0("Multinom_Model_n", n_val, ".rds")))
  
  message("✓ Done  n = ", n_val, "  →  ", n_dir)
}

message("\n✅ Analysis complete!")