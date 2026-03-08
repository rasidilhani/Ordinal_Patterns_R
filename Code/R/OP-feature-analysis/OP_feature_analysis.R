install.packages("nnet")      # For multinomial logistic regression
install.packages("dplyr")     
install.packages("readxl")    
install.packages("here")      

# Load packages
library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)
library(corrplot)

#here()

# read data from directory
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

#dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

str(df)
head(df)

# Check unique models
unique(df$Model)

# Select the features to use for the model
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

# Check correlations (multicollinearity check)
cor_matrix <- cor(df[, features], use = "complete.obs")
print("Correlation Matrix:")
print(round(cor_matrix, 3))
corrplot(cor_matrix, method = "color", addCoef.col = "black")  # Visualize


# Create analysis dataset with Model and features
analysis_data <- df %>%
  select(Model, all_of(features))

# Check for missing values
colSums(is.na(analysis_data))

# Remove rows with missing values if any
analysis_data <- na.omit(analysis_data)

# Convert Model to factor
analysis_data$Model <- as.factor(analysis_data$Model)

# Check levels
#levels(analysis_data$Model)

# Set a reference level (I considered AR1_M1 as baseline)
analysis_data$Model <- relevel(analysis_data$Model, ref = "AR1_M1")

# Verify
#levels(analysis_data$Model)

# Fit the model
multinom_model <- multinom(Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
                             C_Shannon + H_Disequlibrium + C_Renyi + 
                             C_Tsallis + C_Fisher, 
                           data = analysis_data)

# View summary
Model_Summary <- summary(multinom_model)
print(Model_Summary)

# Get z-values (coefficient/standard error)
z_values <- summary(multinom_model)$coefficients / summary(multinom_model)$standard.errors

# Get p-values (2-tailed test)
p_values <- (1 - pnorm(abs(z_values))) * 2

# View p-values
print(round(p_values, 4))

# Identify significant predictors (p < 0.05)
significant <- p_values < 0.05
print(significant)

# Predictions
predictions <- predict(multinom_model, type = "class")

# Confusion matrix
table(Predicted = predictions, Actual = analysis_data$Model)

# Accuracy
accuracy <- mean(predictions == analysis_data$Model)
cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")



####-----------------------------------------------------------------------------------
# Save the results into excel
#-----------------------------------------------------------------------------------------
# Install packages (run once)
# install.packages("nnet")      
# install.packages("dplyr")     
# install.packages("readxl")    
# install.packages("here")
# install.packages("writexl")

# Load packages
library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)

# Read data from directory
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

str(df)
head(df)

# Check unique models
unique(df$Model)

# Select the features to use for the model
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

# Create analysis dataset with Model and features
analysis_data <- df %>%
  select(Model, all_of(features))

# Check for missing values
colSums(is.na(analysis_data))

# Remove rows with missing values if any
analysis_data <- na.omit(analysis_data)

# Convert Model to factor
analysis_data$Model <- as.factor(analysis_data$Model)

# Set a reference level (I considered AR1_M1 as baseline)
analysis_data$Model <- relevel(analysis_data$Model, ref = "AR1_M1")

# Fit the model
multinom_model <- multinom(Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
                             C_Shannon + H_Disequlibrium + C_Renyi + 
                             C_Tsallis + C_Fisher, 
                           data = analysis_data)

# View summary
Model_Summary <- summary(multinom_model)
print(Model_Summary)

# Get z-values (coefficient/standard error)
z_values <- summary(multinom_model)$coefficients / summary(multinom_model)$standard.errors

# Get p-values (2-tailed test)
p_values <- (1 - pnorm(abs(z_values))) * 2

# View p-values
print(round(p_values, 4))

# Identify significant predictors (p < 0.05)
significant <- p_values < 0.05
print(significant)

# Predictions
predictions <- predict(multinom_model, type = "class")

# Confusion matrix
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data$Model)
print(confusion_matrix)

# Accuracy
accuracy <- mean(predictions == analysis_data$Model)
cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")

# ====== SAVE RESULTS TO EXCEL ======

# Create Results folder
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

# 1. Coefficients
coefficients_df <- as.data.frame(coef(multinom_model))
coefficients_df <- cbind(Model = rownames(coefficients_df), coefficients_df)
rownames(coefficients_df) <- NULL

# 2. Standard Errors
std_errors_df <- as.data.frame(Model_Summary$standard.errors)
std_errors_df <- cbind(Model = rownames(std_errors_df), std_errors_df)
rownames(std_errors_df) <- NULL

# 3. Z-values
z_values_df <- as.data.frame(z_values)
z_values_df <- cbind(Model = rownames(z_values_df), z_values_df)
rownames(z_values_df) <- NULL

# 4. P-values
p_values_df <- as.data.frame(round(p_values, 4))
p_values_df <- cbind(Model = rownames(p_values_df), p_values_df)
rownames(p_values_df) <- NULL

# 5. Significant indicators
significant_df <- as.data.frame(significant)
significant_df <- cbind(Model = rownames(significant_df), significant_df)
rownames(significant_df) <- NULL

# 6. Confusion Matrix
confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted_Model = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

# 7. Accuracy Summary
model_accuracy <- analysis_data %>%
  mutate(Predicted = predictions,
         Correct = Model == Predicted) %>%
  group_by(Model) %>%
  summarise(
    Total = n(),
    Correct = sum(Correct),
    Accuracy_Percent = round(mean(Correct) * 100, 2)
  )

overall_accuracy <- data.frame(
  Model = "OVERALL",
  Total = nrow(analysis_data),
  Correct = sum(predictions == analysis_data$Model),
  Accuracy_Percent = round(accuracy * 100, 2)
)

accuracy_summary <- bind_rows(model_accuracy, overall_accuracy)

# 8. Model Information
model_info <- data.frame(
  Information = c("Reference Level", "Number of Models", "Number of Observations", 
                  "Number of Features", "Overall Accuracy (%)", "AIC"),
  Value = c(levels(analysis_data$Model)[1], 
            length(levels(analysis_data$Model)), 
            nrow(analysis_data),
            length(features),
            round(accuracy * 100, 2),
            round(multinom_model$AIC, 2))
)

# Save all results in one Excel file with multiple sheets
all_results <- list(
  "Model_Info" = model_info,
  "Accuracy_Summary" = accuracy_summary,
  "Confusion_Matrix" = confusion_df,
  "Coefficients" = coefficients_df,
  "Standard_Errors" = std_errors_df,
  "Z_Values" = z_values_df,
  "P_Values" = p_values_df,
  "Significant_p005" = significant_df
)

write_xlsx(all_results, 
           path = here("Results", "OP_feature_analysis", "Multinomial_Results.xlsx"))

cat("\n========================================\n")
cat("Results saved to:\n")
cat(here("Results", "OP_feature_analysis", "Multinomial_Results.xlsx"), "\n")
cat("========================================\n")

#------------------------------------------------------------------------------------------------------------------
#Based on the case study results, we can try different modeling strategies to improve classification performance. 
#----------------------------------------------------------------------------------------------
library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)
library(corrplot)
library(tidyr)

df <- read_excel(
  here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"),
  sheet = "D3n1000"
)

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

features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

for (case_name in names(cases)) {
  cat("\n====================\n", case_name, "\n====================\n")
  
  models_in_case <- cases[[case_name]]
  
  analysis_data <- df %>%
    filter(Model %in% models_in_case) %>%
    select(Model, all_of(features)) %>%
    na.omit()
  
  if (nrow(analysis_data) < length(models_in_case)) {
    warning("Not enough observations for ", case_name, " – skipping.")
    next
  }
  
  analysis_data$Model <- factor(analysis_data$Model)
  
  baseline <- dplyr::case_when(
    case_name == "Case1" ~ "AR1_M1",
    case_name == "Case2" ~ "AR1_M3",
    case_name == "Case3" ~ "AR2_M2",
    TRUE ~ levels(analysis_data$Model)[1]
  )
  
  if (!baseline %in% levels(analysis_data$Model)) {
    baseline <- levels(analysis_data$Model)[1]
  }
  
  analysis_data$Model <- relevel(analysis_data$Model, ref = baseline)
  
  formula_case <- as.formula(paste("Model ~", paste(features, collapse = " + ")))
  multinom_model <- multinom(formula_case, data = analysis_data)
  
  model_summary <- summary(multinom_model)
  coefs <- model_summary$coefficients
  ses   <- model_summary$standard.errors
  z_vals <- coefs / ses
  p_vals <- (1 - pnorm(abs(z_vals))) * 2
  
  preds <- predict(multinom_model, type = "class")
  cm    <- table(Predicted = preds, Actual = analysis_data$Model)
  acc   <- mean(preds == analysis_data$Model)
  
  cat("Accuracy for", case_name, ":", round(acc * 100, 2), "%\n")
  
  ## ---------- Create all export sheets for this case ----------
  
  # 1. Model Info
  model_info <- data.frame(
    Case = case_name,
    Baseline = baseline,
    N_observations = nrow(analysis_data),
    N_models = length(unique(analysis_data$Model)),
    Accuracy = round(acc, 4),
    row.names = NULL
  )
  
  # 2. Coefficients (wide format)
  coef_df <- as.data.frame(coefs)
  coef_df$Comparison <- rownames(coefs)
  
  # 3. Standard Errors (wide format)
  se_df <- as.data.frame(ses)
  se_df$Comparison <- rownames(ses)
  
  # 4. Z-values (wide format)
  z_df <- as.data.frame(z_vals)
  z_df$Comparison <- rownames(z_vals)
  
  # 5. P-values (wide format)
  p_df <- as.data.frame(p_vals)
  p_df$Comparison <- rownames(p_vals)
  
  # 6. Significance (TRUE/FALSE)
  sig_df <- as.data.frame(p_vals < 0.05)
  sig_df$Comparison <- rownames(sig_df)
  colnames(sig_df)[1:(ncol(sig_df)-1)] <- paste0(colnames(sig_df)[1:(ncol(sig_df)-1)], "_significant")
  
  # 7. Confusion Matrix
  cm_df <- as.data.frame(cm)
  names(cm_df) <- c("Predicted", "Actual", "Freq")
  cm_df$Case     <- case_name
  cm_df$Baseline <- baseline
  
  # 8. Accuracy Summary (with per-class accuracy)
  class_acc <- prop.table(table(preds, analysis_data$Model), 2)
  acc_long <- as.data.frame(class_acc)
  names(acc_long) <- c("Predicted", "Actual", "Accuracy")
  acc_long$Case <- case_name
  
  ## ---------- Named list of ALL sheets for this case ----------
  all_case_sheets <- list(
    "Model_Info"              = model_info,
    "Coefficients"            = coef_df,
    "Standard_Errors"         = se_df,
    "Z_Values"                = z_df,
    "P_Values"                = p_df,
    "Significant"             = sig_df,
    "Confusion_Matrix"        = cm_df,
    "Class_Accuracy"          = acc_long
  )
  
  ## ---------- Save ONE Excel file per case ----------
  excel_file <- here("Results", "OP_feature_analysis", 
                     paste0("Multinomial_Results_", case_name, ".xlsx"))
  
  write_xlsx(all_case_sheets, path = excel_file)
  
  cat("✅ Saved", case_name, "results to:", excel_file, "\n")
}

cat("\n🎉 All case analyses complete! Check Results/OP_feature_analysis/ for 3 Excel files.\n")

#------------------------------------------------------------------------------------------------------------------
# Case study with interaction terms
#------------------------------------------------------------------------------------------------------------------
library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)
library(corrplot)
library(tidyr)

df <- read_excel(
  here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"),
  sheet = "D3n1000"
)

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

features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

for (case_name in names(cases)) {
  cat("\n====================\n", case_name, "\n====================\n")
  
  models_in_case <- cases[[case_name]]
  
  analysis_data <- df %>%
    filter(Model %in% models_in_case) %>%
    select(Model, all_of(features)) %>%
    na.omit()
  
  if (nrow(analysis_data) < length(models_in_case)) {
    warning("Not enough observations for ", case_name, " – skipping.")
    next
  }
  
  analysis_data$Model <- factor(analysis_data$Model)
  
  baseline <- dplyr::case_when(
    case_name == "Case1" ~ "AR1_M1",
    case_name == "Case2" ~ "ARMA11_M3",
    case_name == "Case3" ~ "AR2_M2",
    TRUE ~ levels(analysis_data$Model)[1]
  )
  
  if (!baseline %in% levels(analysis_data$Model)) {
    baseline <- levels(analysis_data$Model)[1]
  }
  
  analysis_data$Model <- relevel(analysis_data$Model, ref = baseline)
  
  # ✅ CHANGED: Model formula with interaction terms
  multinom_model <- multinom(Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
                               C_Shannon + H_Disequlibrium + C_Renyi + 
                               C_Tsallis + C_Fisher +
                               H_Fisher:C_Fisher +      # Interaction
                               H_Shannon:C_Shannon +    # Interaction
                               H_Tsallis:C_Tsallis +    # Interaction
                               H_Renyi:C_Renyi,         # Interaction
                             data = analysis_data)
  
  model_summary <- summary(multinom_model)
  coefs <- model_summary$coefficients
  ses   <- model_summary$standard.errors
  z_vals <- coefs / ses
  p_vals <- (1 - pnorm(abs(z_vals))) * 2
  
  preds <- predict(multinom_model, type = "class")
  cm    <- table(Predicted = preds, Actual = analysis_data$Model)
  acc   <- mean(preds == analysis_data$Model)
  
  cat("Accuracy for", case_name, ":", round(acc * 100, 2), "%\n")
  
  ## ---------- Create all export sheets for this case ----------
  
  # 1. Model Info (updated filename to indicate interactions)
  model_info <- data.frame(
    Case = case_name,
    Baseline = baseline,
    N_observations = nrow(analysis_data),
    N_models = length(unique(analysis_data$Model)),
    Accuracy = round(acc, 4),
    Model_Type = "Multinomial_with_Interactions",
    row.names = NULL
  )
  
  # 2. Coefficients (wide format) - now includes interaction terms
  coef_df <- as.data.frame(coefs)
  coef_df$Comparison <- rownames(coefs)
  
  # 3. Standard Errors
  se_df <- as.data.frame(ses)
  se_df$Comparison <- rownames(ses)
  
  # 4. Z-values
  z_df <- as.data.frame(z_vals)
  z_df$Comparison <- rownames(z_vals)
  
  # 5. P-values
  p_df <- as.data.frame(p_vals)
  p_df$Comparison <- rownames(p_vals)
  
  # 6. Significance (TRUE/FALSE)
  sig_df <- as.data.frame(p_vals < 0.05)
  sig_df$Comparison <- rownames(sig_df)
  colnames(sig_df)[1:(ncol(sig_df)-1)] <- paste0(colnames(sig_df)[1:(ncol(sig_df)-1)], "_significant")
  
  # 7. Confusion Matrix
  cm_df <- as.data.frame(cm)
  names(cm_df) <- c("Predicted", "Actual", "Freq")
  cm_df$Case     <- case_name
  cm_df$Baseline <- baseline
  
  # 8. Class Accuracy
  class_acc <- prop.table(table(preds, analysis_data$Model), 2)
  acc_long <- as.data.frame(class_acc)
  names(acc_long) <- c("Predicted", "Actual", "Accuracy")
  acc_long$Case <- case_name
  
  ## ---------- Named list of ALL sheets ----------
  all_case_sheets <- list(
    "Model_Info"              = model_info,
    "Coefficients"            = coef_df,
    "Standard_Errors"         = se_df,
    "Z_Values"                = z_df,
    "P_Values"                = p_df,
    "Significant"             = sig_df,
    "Confusion_Matrix"        = cm_df,
    "Class_Accuracy"          = acc_long
  )
  
  ## ---------- Save ONE Excel file per case ----------
  excel_file <- here("Results", "OP_feature_analysis", 
                     paste0("Multinomial_Interactions_", case_name, ".xlsx"))
  
  write_xlsx(all_case_sheets, path = excel_file)
  
  cat("✅ Saved INTERACTION model results for", case_name, "to:", excel_file, "\n")
}

cat("\n🎉 All interaction model analyses complete!\n")
cat("Check Results/OP_feature_analysis/ for 3 new Excel files with '_Interactions_' suffix.\n")

#------------------------------------------------------------------------------------------------------------------
# Modified model to exclude H_Fisher and C_Fisher based on significance results
#-------------------------------------------------------------------------------------------------
# # ===== STRATEGY 1: FEATURE SELECTION =====
# Keep only features significant for at least 50% of models

library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)

# Read data
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

# Original features
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

analysis_data <- df %>%
  select(Model, all_of(features)) %>%
  na.omit()

analysis_data$Model <- as.factor(analysis_data$Model)
analysis_data$Model <- relevel(analysis_data$Model, ref = "AR1_M1")

# Fit original model to identify important features
multinom_model_original <- multinom(Model ~ ., data = analysis_data)
z_values_original <- summary(multinom_model_original)$coefficients / 
  summary(multinom_model_original)$standard.errors
p_values_original <- (1 - pnorm(abs(z_values_original))) * 2

# Count significance for each feature
feature_importance <- colSums(p_values_original < 0.05, na.rm = TRUE)
print(feature_importance)

# Fit new model with selected features Total for out of 23
important_features <- names(feature_importance[feature_importance >= 12])
cat("Selected features:", paste(important_features, collapse = ", "), "\n")

# Create new dataset with selected features only
analysis_data_selected <- analysis_data %>%
  select(Model, all_of(important_features))

# Fit new model with selected features
multinom_model <- multinom(Model ~ ., data = analysis_data_selected)

# Get results
Model_Summary <- summary(multinom_model)
z_values <- Model_Summary$coefficients / Model_Summary$standard.errors
p_values <- (1 - pnorm(abs(z_values))) * 2
significant <- p_values < 0.05

predictions <- predict(multinom_model, type = "class")
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data_selected$Model)
accuracy <- mean(predictions == analysis_data_selected$Model)

cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")

# ===== SAVE RESULTS =====
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

coefficients_df <- as.data.frame(coef(multinom_model))
coefficients_df <- cbind(Model = rownames(coefficients_df), coefficients_df)
rownames(coefficients_df) <- NULL

std_errors_df <- as.data.frame(Model_Summary$standard.errors)
std_errors_df <- cbind(Model = rownames(std_errors_df), std_errors_df)
rownames(std_errors_df) <- NULL

z_values_df <- as.data.frame(z_values)
z_values_df <- cbind(Model = rownames(z_values_df), z_values_df)
rownames(z_values_df) <- NULL

p_values_df <- as.data.frame(round(p_values),4)
p_values_df <- cbind(Model = rownames(p_values_df), p_values_df)
rownames(p_values_df) <- NULL

significant_df <- as.data.frame(significant)
significant_df <- cbind(Model = rownames(significant_df), significant_df)
rownames(significant_df) <- NULL

confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted_Model = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

model_accuracy <- analysis_data_selected %>%
  mutate(Predicted = predictions, Correct = Model == Predicted) %>%
  group_by(Model) %>%
  summarise(Total = n(), Correct = sum(Correct), Accuracy_Percent = round(mean(Correct) * 100, 2))

overall_accuracy <- data.frame(Model = "OVERALL", Total = nrow(analysis_data_selected),
                               Correct = sum(predictions == analysis_data_selected$Model),
                               Accuracy_Percent = round(accuracy * 100, 2))
accuracy_summary <- bind_rows(model_accuracy, overall_accuracy)

model_info <- data.frame(
  Information = c("Strategy", "Reference Level", "Number of Models", "Number of Observations", 
                  "Number of Features", "Selected Features", "Overall Accuracy (%)", "AIC"),
  Value = c("Feature Selection (50% threshold)",
            levels(analysis_data_selected$Model)[1], 
            length(levels(analysis_data_selected$Model)), 
            nrow(analysis_data_selected),
            length(important_features),
            paste(important_features, collapse = ", "),
            round(accuracy * 100, 2),
            round(multinom_model$AIC, 2))
)

all_results <- list(
  "Model_Info" = model_info,
  "Accuracy_Summary" = accuracy_summary,
  "Confusion_Matrix" = confusion_df,
  "Coefficients" = coefficients_df,
  "Standard_Errors" = std_errors_df,
  "Z_Values" = z_values_df,
  "P_Values" = p_values_df,
  "Significant_p005" = significant_df
)

write_xlsx(all_results, 
           path = here("Results", "OP_feature_analysis", "Strategy1_Feature_Selection.xlsx"))

cat("\nResults saved to: Strategy1_Feature_Selection.xlsx\n")
#------------------------------------------------------------------------------------------------------------------

# ===== STRATEGY 2: INTERACTION TERMS =====
# Add interactions between important features

library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)

# Read data
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

analysis_data <- df %>%
  select(Model, all_of(features)) %>%
  na.omit()

analysis_data$Model <- as.factor(analysis_data$Model)
analysis_data$Model <- relevel(analysis_data$Model, ref = "AR1_M1")

# Fit model with interaction terms
multinom_model <- multinom(Model ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
                             C_Shannon + H_Disequlibrium + C_Renyi + 
                             C_Tsallis + C_Fisher +
                             H_Fisher:C_Fisher +      # Interaction
                             H_Shannon:C_Shannon +    # Interaction
                             H_Tsallis:C_Tsallis +   # Interaction
                             H_Renyi:C_Renyi,        # Interaction
                           data = analysis_data)

# Get results
Model_Summary <- summary(multinom_model)
z_values <- Model_Summary$coefficients / Model_Summary$standard.errors
p_values <- (1 - pnorm(abs(z_values))) * 2
significant <- p_values < 0.05

predictions <- predict(multinom_model, type = "class")
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data$Model)
accuracy <- mean(predictions == analysis_data$Model)

cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")

# ===== SAVE RESULTS =====
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

coefficients_df <- as.data.frame(coef(multinom_model))
coefficients_df <- cbind(Model = rownames(coefficients_df), coefficients_df)
rownames(coefficients_df) <- NULL

std_errors_df <- as.data.frame(round(Model_Summary$standard.errors),4)
std_errors_df <- cbind(Model = rownames(std_errors_df), std_errors_df)
rownames(std_errors_df) <- NULL

z_values_df <- as.data.frame(z_values)
z_values_df <- cbind(Model = rownames(z_values_df), z_values_df)
rownames(z_values_df) <- NULL

p_values_df <- as.data.frame(round(p_values, 4))
p_values_df <- cbind(Model = rownames(p_values_df), p_values_df)
rownames(p_values_df) <- NULL

significant_df <- as.data.frame(significant)
significant_df <- cbind(Model = rownames(significant_df), significant_df)
rownames(significant_df) <- NULL

confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted_Model = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

model_accuracy <- analysis_data %>%
  mutate(Predicted = predictions, Correct = Model == Predicted) %>%
  group_by(Model) %>%
  summarise(Total = n(), Correct = sum(Correct), Accuracy_Percent = round(mean(Correct) * 100, 2))

overall_accuracy <- data.frame(Model = "OVERALL", Total = nrow(analysis_data),
                               Correct = sum(predictions == analysis_data$Model),
                               Accuracy_Percent = round(accuracy * 100, 2))
accuracy_summary <- bind_rows(model_accuracy, overall_accuracy)

model_info <- data.frame(
  Information = c("Strategy", "Reference Level", "Number of Models", "Number of Observations", 
                  "Number of Features", "Interaction Terms", "Overall Accuracy (%)", "AIC"),
  Value = c("Interaction Terms",
            levels(analysis_data$Model)[1], 
            length(levels(analysis_data$Model)), 
            nrow(analysis_data),
            "9 main + 3 interactions",
            "H_Fisher:C_Fisher, H_Shannon:C_Shannon, H_Tsallis:C_Tsallis",
            round(accuracy * 100, 2),
            round(multinom_model$AIC, 2))
)

all_results <- list(
  "Model_Info" = model_info,
  "Accuracy_Summary" = accuracy_summary,
  "Confusion_Matrix" = confusion_df,
  "Coefficients" = coefficients_df,
  "Standard_Errors" = std_errors_df,
  "Z_Values" = z_values_df,
  "P_Values" = p_values_df,
  "Significant_p005" = significant_df
)

write_xlsx(all_results, 
           path = here("Results", "OP_feature_analysis", "Strategy2_Interaction_Terms.xlsx"))

cat("\nResults saved to: Strategy2_Interaction_Terms.xlsx\n")
#------------------------------------------------------------------------------------------------------------------

# ===== STRATEGY 3: FEATURE ENGINEERING =====
# Create composite features from existing ones

library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)

# Read data
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

analysis_data <- df %>%
  select(Model, all_of(features)) %>%
  na.omit()

# Create engineered features
analysis_data_engineered <- analysis_data %>%
  mutate(
    Entropy_Avg = (H_Shannon + H_Renyi + H_Tsallis) / 3,
    Complexity_Avg = (C_Shannon + C_Renyi + C_Tsallis) / 3,
    Fisher_Ratio = H_Fisher / (C_Fisher + 0.001)  # Avoid division by zero
    #Total_Info = H_Shannon + C_Shannon,
    #Entropy_Complexity_Diff = Entropy_Avg - Complexity_Avg
  )

analysis_data_engineered$Model <- as.factor(analysis_data_engineered$Model)
analysis_data_engineered$Model <- relevel(analysis_data_engineered$Model, ref = "AR1_M1")

# Fit model with engineered features
#multinom_model <- multinom(Model ~ H_Fisher + C_Fisher + Fisher_Ratio + 
#                             Entropy_Avg + Complexity_Avg + 
#                             H_Disequlibrium + Total_Info + 
 #                            Entropy_Complexity_Diff,
 #                          data = analysis_data_engineered)
multinom_model <- multinom(Model ~ H_Fisher + C_Fisher + Fisher_Ratio + 
                             Entropy_Avg + Complexity_Avg + 
                             H_Disequlibrium,
                           data = analysis_data_engineered)

# Get results
Model_Summary <- summary(multinom_model)
z_values <- Model_Summary$coefficients / Model_Summary$standard.errors
p_values <- (1 - pnorm(abs(z_values))) * 2
significant <- p_values < 0.05

predictions <- predict(multinom_model, type = "class")
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data_engineered$Model)
accuracy <- mean(predictions == analysis_data_engineered$Model)

cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")

# ===== SAVE RESULTS =====
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

coefficients_df <- as.data.frame(coef(multinom_model))
coefficients_df <- cbind(Model = rownames(coefficients_df), coefficients_df)
rownames(coefficients_df) <- NULL

std_errors_df <- as.data.frame(round(Model_Summary$standard.errors),4)
std_errors_df <- cbind(Model = rownames(std_errors_df), std_errors_df)
rownames(std_errors_df) <- NULL

z_values_df <- as.data.frame(z_values)
z_values_df <- cbind(Model = rownames(z_values_df), z_values_df)
rownames(z_values_df) <- NULL

p_values_df <- round(as.data.frame(p_values),4)
p_values_df <- cbind(Model = rownames(p_values_df), p_values_df)
rownames(p_values_df) <- NULL

significant_df <- as.data.frame(significant)
significant_df <- cbind(Model = rownames(significant_df), significant_df)
rownames(significant_df) <- NULL

confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted_Model = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

model_accuracy <- analysis_data_engineered %>%
  mutate(Predicted = predictions, Correct = Model == Predicted) %>%
  group_by(Model) %>%
  summarise(Total = n(), Correct = sum(Correct), Accuracy_Percent = round(mean(Correct) * 100, 2))

overall_accuracy <- data.frame(Model = "OVERALL", Total = nrow(analysis_data_engineered),
                               Correct = sum(predictions == analysis_data_engineered$Model),
                               Accuracy_Percent = round(accuracy * 100, 2))
accuracy_summary <- bind_rows(model_accuracy, overall_accuracy)

model_info <- data.frame(
  Information = c("Strategy", "Reference Level", "Number of Models", "Number of Observations", 
                  "Number of Features", "Engineered Features", "Overall Accuracy (%)", "AIC"),
  Value = c("Feature Engineering",
            levels(analysis_data_engineered$Model)[1], 
            length(levels(analysis_data_engineered$Model)), 
            nrow(analysis_data_engineered),
            "8 features",
            "Entropy_Avg, Complexity_Avg, Fisher_Ratio, Total_Info, Entropy_Complexity_Diff",
            round(accuracy * 100, 2),
            round(multinom_model$AIC, 2))
)

all_results <- list(
  "Model_Info" = model_info,
  "Accuracy_Summary" = accuracy_summary,
  "Confusion_Matrix" = confusion_df,
  "Coefficients" = coefficients_df,
  "Standard_Errors" = std_errors_df,
  "Z_Values" = z_values_df,
  "P_Values" = p_values_df,
  "Significant_p005" = significant_df
)

write_xlsx(all_results, 
           path = here("Results", "OP_feature_analysis", "Strategy3_Feature_Engineering.xlsx"))

cat("\nResults saved to: Strategy3_Feature_Engineering.xlsx\n")

#------------------------------------------------------------------------------------------------------------------
# Model reduction by combining similar models
#------------------------------------------------------------------------------------------------------------------
# Load packages
library(nnet)
library(dplyr)
library(readxl)
library(here)
library(writexl)
library(corrplot)

# Read data from directory
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

#str(df)
#head(df)

# Check original unique models
#cat("Original unique models:\n")
print(unique(df$Model))

# Create a function to extract base model type
extract_model_type <- function(model_name) {
  # Extract the model type (e.g., "AR1_M1" -> "AR(1)")
  if (grepl("^AR1_", model_name)) return("AR(1)")
  if (grepl("^AR2_", model_name)) return("AR(2)")
  if (grepl("^MA1_", model_name)) return("MA(1)")
  if (grepl("^MA2_", model_name)) return("MA(2)")
  if (grepl("^ARMA11_", model_name)) return("ARMA(1,1)")
  if (grepl("^ARMA22_", model_name)) return("ARMA(2,2)")
  return(NA)  # For any unmatched models
}

# Create new simplified Model column
df$Model_Type <- sapply(df$Model, extract_model_type)

# Check if there are any NA values (unmatched models)
if (any(is.na(df$Model_Type))) {
  cat("\nWarning: Some models could not be classified:\n")
  print(unique(df$Model[is.na(df$Model_Type)]))
}

# Remove rows with unclassified models
df <- df %>% filter(!is.na(Model_Type))

# Check new unique model types
cat("\nNew model types:\n")
print(unique(df$Model_Type))
cat("\nNumber of observations per model type:\n")
print(table(df$Model_Type))

# Select the features to use for the model
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

# Check correlations (multicollinearity check)
cor_matrix <- cor(df[, features], use = "complete.obs")
print("\nCorrelation Matrix:")
print(round(cor_matrix, 3))
corrplot(cor_matrix, method = "color", addCoef.col = "black")  # Visualize

# Create analysis dataset with Model_Type and features
analysis_data <- df %>%
  select(Model_Type, all_of(features))

# Check for missing values
cat("\nMissing values per column:\n")
print(colSums(is.na(analysis_data)))

# Remove rows with missing values if any
analysis_data <- na.omit(analysis_data)

# Convert Model_Type to factor
analysis_data$Model_Type <- as.factor(analysis_data$Model_Type)

# Check levels
cat("\nModel levels:\n")
print(levels(analysis_data$Model_Type))

# Set a reference level (AR(1) as baseline)
analysis_data$Model_Type <- relevel(analysis_data$Model_Type, ref = "AR(1)")

# Verify
cat("\nModel levels after releveling:\n")
print(levels(analysis_data$Model_Type))

# Fit the multinomial logistic regression model
multinom_model <- multinom(Model_Type ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
                             C_Shannon + H_Disequlibrium + C_Renyi + 
                             C_Tsallis + C_Fisher, 
                           data = analysis_data)

# View summary
Model_Summary <- summary(multinom_model)
print("\nModel Summary:")
print(Model_Summary)

# Get z-values (coefficient/standard error)
z_values <- summary(multinom_model)$coefficients / summary(multinom_model)$standard.errors

# Get p-values (2-tailed test)
p_values <- (1 - pnorm(abs(z_values))) * 2

# View p-values
cat("\nP-values:\n")
print(round(p_values, 4))

# Identify significant predictors (p < 0.05)
significant <- p_values < 0.05
cat("\nSignificant predictors (p < 0.05):\n")
print(significant)

# Predictions
predictions <- predict(multinom_model, type = "class")

# Confusion matrix
cat("\nConfusion Matrix:\n")
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data$Model_Type)
print(confusion_matrix)

# Accuracy
accuracy <- mean(predictions == analysis_data$Model_Type)
cat("\nOverall Accuracy:", round(accuracy * 100, 2), "%\n")

# ============================================================
# PREPARE RESULTS FOR EXCEL EXPORT
# ============================================================

# 1. Coefficients
coefficients_df <- as.data.frame(summary(multinom_model)$coefficients)
coefficients_df <- cbind(Model = rownames(coefficients_df), coefficients_df)
rownames(coefficients_df) <- NULL

# 2. Standard Errors
std_errors_df <- as.data.frame(summary(multinom_model)$standard.errors)
std_errors_df <- cbind(Model = rownames(std_errors_df), std_errors_df)
rownames(std_errors_df) <- NULL

# 3. Z-values
z_values_df <- as.data.frame(z_values)
z_values_df <- cbind(Model = rownames(z_values_df), z_values_df)
rownames(z_values_df) <- NULL

# 4. P-values
p_values_df <- as.data.frame(p_values)
p_values_df <- cbind(Model = rownames(p_values_df), p_values_df)
rownames(p_values_df) <- NULL

# 5. Significant predictors
significant_df <- as.data.frame(significant)
significant_df <- cbind(Model = rownames(significant_df), significant_df)
rownames(significant_df) <- NULL

# 6. Confusion Matrix
confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

# 7. Model Performance Summary
performance_summary <- data.frame(
  Metric = c("Overall Accuracy (%)", "Number of Observations", "Number of Model Types"),
  Value = c(round(accuracy * 100, 2), nrow(analysis_data), nlevels(analysis_data$Model_Type))
)

# 8. Correlation Matrix
cor_matrix_df <- as.data.frame(cor_matrix)
cor_matrix_df <- cbind(Feature = rownames(cor_matrix_df), cor_matrix_df)
rownames(cor_matrix_df) <- NULL

# Create a list of all results
results_list <- list(
  "Performance_Summary" = performance_summary,
  "Coefficients" = coefficients_df,
  "Standard_Errors" = std_errors_df,
  "Z_Values" = z_values_df,
  "P_Values" = p_values_df,
  "Significant_Predictors" = significant_df,
  "Confusion_Matrix" = confusion_df,
  "Correlation_Matrix" = cor_matrix_df
)

# Create Results directory if it doesn't exist
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

# Save to Excel
output_file <- here("Results", "OP_feature_analysis", "Strategy4_Grouped_Models.xlsx")
write_xlsx(results_list, output_file)

cat("\n============================================================\n")
cat("Results saved to:", output_file, "\n")
cat("============================================================\n")


#------------------------------------------------------------------------------------------------------------------
# This script discuss about GLM that use to classificy models----------

library(readxl)
library(here)
library(dplyr)
library(pROC)
library(ordinalGMifs)
library(ggplot2)

df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

# create BINARY response: 1=ARMA models, 0=others
df$Is_ARMA <- ifelse(grepl("ARMA", df$Model), 1, 0)
df$Is_ARMA <- factor(df$Is_ARMA)  # Make it a factor for glm

# Check data
table(df$Is_ARMA, df$Model)
str(df$Is_ARMA)

# Fit BINARY logistic regression
glm_model <- glm(Is_ARMA ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + 
                   C_Shannon + H_Disequlibrium + C_Renyi + 
                   C_Tsallis + C_Fisher, 
                 data = df, family = binomial(link = "logit"))

summary(glm_model)

# Check if ARMA perfectly separated
table(df$Is_ARMA)
prop.table(table(df$Is_ARMA))  # Balance?

# Visualize key predictors
ggplot(df, aes(x = H_Renyi, fill = Is_ARMA)) + geom_density(alpha=0.5)
ggplot(df, aes(x = C_Tsallis, fill = Is_ARMA)) + geom_density(alpha=0.5)

# Use only significant predictors (remove Fisher)
glm_simple <- glm(Is_ARMA ~ H_Renyi + H_Tsallis + C_Shannon + H_Disequlibrium + C_Tsallis, 
                  data = df, family = binomial())

# Odds ratios
exp(cbind(OR = coef(glm_simple), confint(glm_simple)))

# Predictions & classification
preds <- predict(glm_simple, type="response")
df$Pred_ARMA <- ifelse(preds > 0.5, 1, 0)
table(Actual = df$Is_ARMA, Predicted = df$Pred_ARMA)

# AUC (gold standard)
roc_obj <- roc(df$Is_ARMA, preds)
plot(roc_obj); auc(roc_obj)  # AUC > 0.9 = excellent
#------------------------------------------------------------------------------------------------------------------

#-------- MULTI-GLM FOR 6-WAY MODEL CLASSIFICATION ---------
#-------------------------------------------------------------------------------------------------
# Load libraries in CORRECT order
library(dplyr)
library(readxl) 
library(here)
library(broom)
library(writexl)

# RELOAD YOUR DATA FRESH (CRITICAL!)
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

model_families <- c("AR1", "AR2", "MA1", "MA2", "ARMA11", "ARMA22")

predictors <- c(
  "H_Shannon","H_Renyi","H_Tsallis","H_Fisher",
  "C_Shannon","H_Disequlibrium","C_Renyi",
  "C_Tsallis","C_Fisher"
)

df[predictors] <- scale(df[predictors])


print(head(df_test))

# Now create model family columns 
df$Is_AR1    <- as.numeric(grepl("AR1_", df$Model))
df$Is_AR2    <- as.numeric(grepl("AR2_", df$Model))
df$Is_MA1    <- as.numeric(grepl("MA1_", df$Model))
df$Is_MA2    <- as.numeric(grepl("MA2_", df$Model))
df$Is_ARMA11 <- as.numeric(grepl("ARMA11_", df$Model))
df$Is_ARMA22 <- as.numeric(grepl("ARMA22_", df$Model))

print("Model family columns created!")
print(table(df$Is_AR1))


# Fit 6 separate GLMs
models <- list()
for(family in model_families) {
  formula <- as.formula(paste0("Is_", family, " ~ H_Shannon + H_Renyi + H_Tsallis + H_Fisher + ",
                               "C_Shannon + H_Disequlibrium + C_Renyi + C_Tsallis + C_Fisher"))
  
  models[[family]] <- glm(formula, data = df, family = binomial())
  
  cat("\n=== GLM for", family, "===\n")
  print(summary(models[[family]]))
}

# Extract key results for comparison
results_summary <- bind_rows(
  lapply(model_families, function(family) {
    tidy(models[[family]], conf.int = FALSE) %>%
      filter(term != "(Intercept)") %>%
      mutate(Model_Family = family) %>%
      select(Model_Family, term, estimate, std.error, statistic, p.value)
  })
)

print("📊 KEY PREDICTORS BY MODEL FAMILY:")
print(results_summary %>%
        filter(p.value < 0.05) %>%
        arrange(Model_Family, p.value))

# Predict probabilities for 6 models
df$Pred_AR1    <- predict(models$AR1, newdata = df, type = "response")
df$Pred_AR2    <- predict(models$AR2, newdata = df, type = "response")
df$Pred_MA1    <- predict(models$MA1, newdata = df, type = "response")
df$Pred_MA2    <- predict(models$MA2, newdata = df, type = "response")
df$Pred_ARMA11 <- predict(models$ARMA11, newdata = df, type = "response")
df$Pred_ARMA22 <- predict(models$ARMA22, newdata = df, type = "response")

# classification based on highest predicted probability
pred_cols <- paste0("Pred_", model_families)
prob_mat <- as.matrix(df[, pred_cols])
df$Predicted_Model <- model_families[max.col(prob_mat)]

# Check accuracy by true model
df$true_family <- case_when(
  grepl("AR1_", df$Model) ~ "AR1",
  grepl("AR2_", df$Model) ~ "AR2",
  grepl("MA1_", df$Model) ~ "MA1",
  grepl("MA2_", df$Model) ~ "MA2",
  grepl("ARMA11_", df$Model) ~ "ARMA11",
  grepl("ARMA22_", df$Model) ~ "ARMA22"
)


table(True = df$true_family, Predicted = df$Predicted_Model)
accuracy <- mean(df$true_family == df$Predicted_Model, na.rm = TRUE)
cat("Overall 6-way classification accuracy:", round(accuracy * 100, 2), "%\n")

cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1","MA1_M1","MA1_M2","MA2_M1","ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4","MA1_M3","MA1_M4","MA2_M4","ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3","MA2_M2","MA2_M3","ARMA22_M2","ARMA22_M3")
)

case_results <- list()
for(case_name in names(cases)) {
  case_df <- df %>% filter(Model %in% cases[[case_name]])
  
  # Predict within case only
  pred_cols_case <- paste0("Pred_", intersect(model_families, 
                                              unique(gsub("_.+", "", case_df$Model))))
  case_fams <- intersect(model_families,
                         unique(gsub("_.+", "", case_df$Model)))
  
  prob_case <- as.matrix(case_df[, paste0("Pred_", case_fams)])
  case_df$Predicted_Case <- case_fams[max.col(prob_case)]
  
  
  acc_case <- mean(gsub("_.+", "", case_df$Model) == case_df$Predicted_Case, na.rm = TRUE)
  case_results[[case_name]] <- acc_case
  
  cat(case_name, "accuracy:", round(acc_case * 100, 1), "%\n")
}

# Save all model summaries to Excel
model_summaries <- lapply(model_families, function(family) {
  tidy(models[[family]]) %>%
    mutate(Model_Family = family)
})
names(model_summaries) <- paste0("GLM_", model_families)

write_xlsx(model_summaries, 
           here("Results", "OP_feature_analysis", "GLM_6way_Model_Discrimination.xlsx"))


#------------------------------------------------------------------------------------------------------------------

# This script discuss about ordinal GMIFS that use to classificy models----------

#-------------------------------------------------
#Stereotype

library(ordinalgmifs)
library(readxl)
library(here)

df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

# Example: Create ordinal response (1=simple AR/MA1, 2=moderate ARMA11/MA2, 3=complex AR2/ARMA22)
df$Model_Ordinal <- factor(df$Model_Name, 
                           levels = c("AR1_M1", "MA1_M1", "ARMA11_M1", "AR2_M1", "ARMA22_M1", ), 
                           ordered = TRUE)

# Prepare predictors (entropy features) and response
X <- as.matrix(df[, c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
                      "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")])  # Your features
y <- df$Model_Ordinal

# Fit cumulative logit model (default)
fit <- ordinal.gmifs(y ~ 1, x = X, probability.model = "Cumulative", link = "logit")

# Results
summary(fit)
plot(fit)
coef(fit)  # Selected features
predict(fit, newx = X, type = "class")  # Predicted classes


