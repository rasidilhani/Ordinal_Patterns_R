# ══════════════════════════════════════════════════════════════════════════════
#  Multinomial Logistic Regression — Detailed Statistical Output
# ══════════════════════════════════════════════════════════════════════════════
#install.packages("VGAM")
# ── Libraries ─────────────────────────────────────────────
library(VGAM)
library(readxl)
library(dplyr)
library(here)
library(caret)

# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
D <- 4
n_val <- 10000

# ══════════════════════════════════════════════════════════════════════════════
#  PATHS
# ══════════════════════════════════════════════════════════════════════════════
hc_data_path <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))

output_dir <- here("Results", "Convex_combination", "Multinom_Detailed")
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

# ── Filter n = 10000 ──────────────────────────────────────
analysis_data <- hc_data %>%
  filter(n == n_val) %>%
  select(Model, all_of(features)) %>%
  na.omit()

# ══════════════════════════════════════════════════════════════════════════════
#  CASE PROCESSING SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
case_summary <- data.frame(
  Total = nrow(hc_data),
  Used = nrow(analysis_data),
  Missing = nrow(hc_data) - nrow(analysis_data)
)

write.csv(case_summary,
          file.path(output_dir, "Case_Processing_Summary.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  FIT MODEL (VGAM)
# ══════════════════════════════════════════════════════════════════════════════
model_fit <- vglm(
  Model ~ .,
  family = multinomial(refLevel = "ARMA(2,2)"),
  data = analysis_data
)

# ══════════════════════════════════════════════════════════════════════════════
#  MODEL FITTING INFORMATION
# ══════════════════════════════════════════════════════════════════════════════
loglik_full <- logLik(model_fit)

null_model <- vglm(Model ~ 1,
                   family = multinomial(refLevel = "ARMA(2,2)"),
                   data = analysis_data)

loglik_null <- logLik(null_model)

model_fit_info <- data.frame(
  LogLik_Full = as.numeric(loglik_full),
  LogLik_Null = as.numeric(loglik_null),
  AIC = AIC(model_fit)
)

write.csv(model_fit_info,
          file.path(output_dir, "Model_Fit_Info.csv"),
          row.names = FALSE)

# ── Likelihood Ratio Test (FIXED) ─────────────────────────
lr_stat <- 2 * (loglik_full - loglik_null)

df <- length(coef(model_fit)) - length(coef(null_model))

p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)

lr_results <- data.frame(
  LR_Statistic = as.numeric(lr_stat),
  DF = df,
  P_Value = p_value
)

write.csv(lr_results,
          file.path(output_dir, "Likelihood_Ratio_Test.csv"),
          row.names = FALSE)
# ══════════════════════════════════════════════════════════════════════════════
#  COMPUTE ALL RESULTS (NO SAVING YET)
# ══════════════════════════════════════════════════════════════════════════════

# Log-likelihoods
loglik_full <- logLik(model_fit)

null_model <- vglm(Model ~ 1,
                   family = multinomial(refLevel = "ARMA(2,2)"),
                   data = analysis_data)

loglik_null <- logLik(null_model)

# Likelihood Ratio Test
lr_stat <- 2 * (loglik_full - loglik_null)
df <- length(coef(model_fit)) - length(coef(null_model))
p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)

# Pseudo R²
pseudo_r2 <- 1 - (loglik_full / loglik_null)

# Goodness-of-fit
deviance_val <- deviance(model_fit)
df_residual  <- df.residual(model_fit)

# ── Predictions (VGAM FIX) ────────────────────────────────
pred_probs <- predict(model_fit, newdata = analysis_data, type = "response")
pred_class <- colnames(pred_probs)[max.col(pred_probs)]

# Confusion Matrix
conf_matrix <- confusionMatrix(as.factor(pred_class), analysis_data$Model)

# ── Case processing summary ──────────────────────────────
total_obs <- nrow(hc_data)
used_obs  <- nrow(analysis_data)
missing_obs <- total_obs - used_obs

# ══════════════════════════════════════════════════════════════════════════════
#  FINAL COMBINED TABLE (ONLY ONE OUTPUT)
# ══════════════════════════════════════════════════════════════════════════════

final_table <- data.frame(
  Sample_Size   = n_val,
  Total_Obs     = total_obs,
  Used_Obs      = used_obs,
  Missing_Obs   = missing_obs,
  Classes       = length(unique(analysis_data$Model)),
  
  LogLik_Full   = round(as.numeric(loglik_full), 3),
  LogLik_Null   = round(as.numeric(loglik_null), 3),
  AIC           = round(AIC(model_fit), 3),
  
  LR_Statistic  = round(as.numeric(lr_stat), 3),
  DF            = df,
  P_Value       = signif(p_value, 4),
  
  McFadden_R2   = round(as.numeric(pseudo_r2), 4),
  
  Deviance      = round(deviance_val, 3),
  DF_Residual   = df_residual,
  
  Accuracy      = round(conf_matrix$overall["Accuracy"], 4),
  Kappa         = round(conf_matrix$overall["Kappa"], 4)
)

# Save ONLY this table
write.csv(final_table,
          file.path(output_dir, paste0("Final_Model_Summary_n", n_val, ".csv")),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
#  OPTIONAL (ONLY if you want extra detail)
# ══════════════════════════════════════════════════════════════════════════════

# Confusion matrix (optional)
write.csv(
  as.data.frame.matrix(table(Predicted = pred_class,
                             Actual = analysis_data$Model)),
  file.path(output_dir, paste0("Confusion_Matrix_n", n_val, ".csv"))
)

# Coefficients (optional)
coef_table <- coef(summary(model_fit))
write.csv(coef_table,
          file.path(output_dir, paste0("Estimates_n", n_val, ".csv")))

message("\n✅ Clean analysis complete — results saved in ONE table.")
#-------------------------------------------------------------------------------

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
#sample_sizes <- c(1000, 10000)
sample_sizes <- 10000

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