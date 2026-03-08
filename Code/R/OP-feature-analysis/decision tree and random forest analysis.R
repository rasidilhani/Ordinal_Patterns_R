
# ===== DECISION TREE CLASSIFICATION =====

# Install packages if needed
# install.packages("rpart")
# install.packages("rpart.plot")

library(rpart)
library(rpart.plot)
library(dplyr)
library(readxl)
library(here)
library(writexl)

# Read data
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

# Select features
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

analysis_data <- df %>%
  select(Model, all_of(features)) %>%
  na.omit()

# Convert Model to factor
analysis_data$Model <- as.factor(analysis_data$Model)

cat("Number of observations:", nrow(analysis_data), "\n")
cat("Number of models:", length(levels(analysis_data$Model)), "\n")
cat("Number of features:", length(features), "\n\n")

# ===== FIT DECISION TREE MODEL =====
set.seed(123)  # For reproducibility

# Fit the tree
dt_model <- rpart(Model ~ ., 
                  data = analysis_data,
                  method = "class",
                  control = rpart.control(
                    minsplit = 20,      # Minimum observations for split
                    minbucket = 7,      # Minimum observations in terminal node
                    cp = 0.001,         # Complexity parameter
                    maxdepth = 30       # Maximum depth
                  ))

# Print tree summary
print(dt_model)

# Print complexity parameter table
cat("\n=== Complexity Parameter Table ===\n")
printcp(dt_model)

# ===== FIND OPTIMAL CP FOR PRUNING =====
# Get the CP with minimum cross-validation error
optimal_cp <- dt_model$cptable[which.min(dt_model$cptable[,"xerror"]), "CP"]
cat("\nOptimal CP:", optimal_cp, "\n")

# Prune the tree
dt_model_pruned <- prune(dt_model, cp = optimal_cp)

cat("\n=== Using Pruned Tree for Predictions ===\n")

# ===== PREDICTIONS =====
predictions <- predict(dt_model_pruned, type = "class")

# ===== CONFUSION MATRIX =====
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data$Model)
print(confusion_matrix)

# ===== ACCURACY =====
accuracy <- mean(predictions == analysis_data$Model)
cat("\nOverall Accuracy:", round(accuracy * 100, 2), "%\n\n")

# ===== VARIABLE IMPORTANCE =====
var_importance <- dt_model_pruned$variable.importance

if (!is.null(var_importance)) {
  importance_df <- data.frame(
    Feature = names(var_importance),
    Importance = var_importance,
    Relative_Importance = round(var_importance / sum(var_importance) * 100, 2)
  )
  importance_df <- importance_df[order(-importance_df$Importance), ]
  rownames(importance_df) <- NULL
  
  cat("=== VARIABLE IMPORTANCE ===\n")
  print(importance_df)
  cat("\n")
} else {
  importance_df <- data.frame(
    Feature = character(),
    Importance = numeric(),
    Relative_Importance = numeric()
  )
  cat("No variable importance available (tree too simple)\n")
}

# ===== CLASS-SPECIFIC ACCURACY =====
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

# ===== TREE STRUCTURE INFORMATION =====
tree_info <- data.frame(
  Metric = c("Number of Terminal Nodes", "Tree Depth", "Number of Splits"),
  Value = c(
    sum(dt_model_pruned$frame$var == "<leaf>"),
    max(rpart:::tree.depth(as.numeric(rownames(dt_model_pruned$frame)))),
    nrow(dt_model_pruned$frame) - sum(dt_model_pruned$frame$var == "<leaf>")
  )
)

# ===== COMPLEXITY PARAMETER TABLE =====
cp_table <- as.data.frame(dt_model_pruned$cptable)
cp_table <- cbind(Split = 0:(nrow(cp_table)-1), cp_table)
rownames(cp_table) <- NULL

# ===== PREPARE RESULTS FOR EXCEL =====

# 1. Model Information
model_info <- data.frame(
  Information = c("Method", "Number of Models", "Number of Observations", 
                  "Number of Features", "Overall Accuracy (%)", 
                  "Optimal CP", "Number of Terminal Nodes", "Tree Depth"),
  Value = c("Decision Tree (Pruned)",
            length(levels(analysis_data$Model)),
            nrow(analysis_data),
            length(features),
            round(accuracy * 100, 2),
            round(optimal_cp, 6),
            tree_info$Value[1],
            tree_info$Value[2])
)

# 2. Confusion Matrix
confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted_Model = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

# ===== SAVE TREE PLOT =====
# Save tree visualization as PDF
pdf(here("Results", "OP_feature_analysis", "Decision_Tree_Plot.pdf"), width = 16, height = 10)
rpart.plot(dt_model_pruned, 
           type = 4,                    # Type of plot
           extra = 2,                   # Display number of observations
           under = TRUE,                # Put extra info under the box
           fallen.leaves = TRUE,        # Align leaves at bottom
           main = "Decision Tree for Model Classification",
           cex = 0.6)                   # Text size
dev.off()

cat("\nDecision tree plot saved as PDF\n")

# ===== SAVE ALL RESULTS TO EXCEL =====
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

all_results <- list(
  "Model_Info" = model_info,
  "Accuracy_Summary" = accuracy_summary,
  "Confusion_Matrix" = confusion_df,
  "Variable_Importance" = importance_df,
  "Tree_Structure" = tree_info,
  "CP_Table" = cp_table
)

write_xlsx(all_results, 
           path = here("Results", "OP_feature_analysis", "Decision_Tree_Results.xlsx"))

cat("\n========================================\n")
cat("DECISION TREE RESULTS SAVED\n")
cat("========================================\n")
cat("Excel File:", here("Results", "OP_feature_analysis", "Decision_Tree_Results.xlsx"), "\n")
cat("Tree Plot:", here("Results", "OP_feature_analysis", "Decision_Tree_Plot.pdf"), "\n")
cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")
cat("Number of Terminal Nodes:", tree_info$Value[1], "\n")
cat("Tree Depth:", tree_info$Value[2], "\n")
if (nrow(importance_df) > 0) {
  cat("Top 3 Important Features:\n")
  print(head(importance_df, 3))
}
cat("========================================\n")

# End of Decesion Tree Analysis
#--------------------------------------------------------------------------------

# ===== RANDOM FOREST CLASSIFICATION =====

# Install packages if needed
# install.packages("randomForest")
# install.packages("caret")

library(randomForest)
library(dplyr)
library(readxl)
library(here)
library(writexl)
library(caret)

# Read data
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

# Select features
features <- c("H_Shannon", "H_Renyi", "H_Tsallis", "H_Fisher", 
              "C_Shannon", "H_Disequlibrium", "C_Renyi", "C_Tsallis", "C_Fisher")

analysis_data <- df %>%
  select(Model, all_of(features)) %>%
  na.omit()

# Convert Model to factor
analysis_data$Model <- as.factor(analysis_data$Model)

cat("Number of observations:", nrow(analysis_data), "\n")
cat("Number of models:", length(levels(analysis_data$Model)), "\n")
cat("Number of features:", length(features), "\n\n")

# ===== FIT RANDOM FOREST MODEL =====
set.seed(123)  # For reproducibility

rf_model <- randomForest(Model ~ ., 
                         data = analysis_data,
                         ntree = 500,           # Number of trees
                         mtry = 3,              # Number of variables tried at each split
                         importance = TRUE,     # Calculate variable importance
                         proximity = FALSE,
                         keep.forest = TRUE)

# Print model summary
print(rf_model)

# ===== PREDICTIONS =====
predictions <- predict(rf_model, type = "class")

# ===== CONFUSION MATRIX =====
confusion_matrix <- table(Predicted = predictions, Actual = analysis_data$Model)
print(confusion_matrix)

# ===== ACCURACY =====
# Overall accuracy
accuracy <- mean(predictions == analysis_data$Model)
cat("\nOverall Accuracy:", round(accuracy * 100, 2), "%\n")

# OOB (Out-of-Bag) Error Rate
oob_error <- rf_model$err.rate[500, "OOB"]
oob_accuracy <- 1 - oob_error
cat("OOB Accuracy:", round(oob_accuracy * 100, 2), "%\n\n")

# ===== FEATURE IMPORTANCE =====
# Get importance measures
importance_measures <- importance(rf_model)

# Mean Decrease in Accuracy
importance_accuracy <- data.frame(
  Feature = rownames(importance_measures),
  MeanDecreaseAccuracy = importance_measures[, "MeanDecreaseAccuracy"]
)
importance_accuracy <- importance_accuracy[order(-importance_accuracy$MeanDecreaseAccuracy), ]
rownames(importance_accuracy) <- NULL

# Mean Decrease in Gini
importance_gini <- data.frame(
  Feature = rownames(importance_measures),
  MeanDecreaseGini = importance_measures[, "MeanDecreaseGini"]
)
importance_gini <- importance_gini[order(-importance_gini$MeanDecreaseGini), ]
rownames(importance_gini) <- NULL

cat("=== FEATURE IMPORTANCE (by Mean Decrease Accuracy) ===\n")
print(importance_accuracy)
cat("\n")

# ===== CLASS-SPECIFIC ACCURACY =====
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

# ===== ERROR RATE BY MODEL =====
error_rates <- rf_model$confusion[, "class.error"]
error_rate_df <- data.frame(
  Model = names(error_rates),
  Error_Rate = round(error_rates * 100, 2),
  Accuracy_Percent = round((1 - error_rates) * 100, 2)
)
error_rate_df <- error_rate_df[order(error_rate_df$Error_Rate), ]
rownames(error_rate_df) <- NULL

# ===== PREPARE RESULTS FOR EXCEL =====

# 1. Model Information
model_info <- data.frame(
  Information = c("Method", "Number of Trees", "Variables per Split (mtry)", 
                  "Number of Models", "Number of Observations", 
                  "Number of Features", "Overall Accuracy (%)", 
                  "OOB Accuracy (%)", "OOB Error Rate (%)"),
  Value = c("Random Forest",
            rf_model$ntree,
            rf_model$mtry,
            length(levels(analysis_data$Model)),
            nrow(analysis_data),
            length(features),
            round(accuracy * 100, 2),
            round(oob_accuracy * 100, 2),
            round(oob_error * 100, 2))
)

# 2. Confusion Matrix
confusion_df <- as.data.frame.matrix(confusion_matrix)
confusion_df <- cbind(Predicted_Model = rownames(confusion_df), confusion_df)
rownames(confusion_df) <- NULL

# 3. Variable Importance - Combined
variable_importance <- data.frame(
  Feature = importance_accuracy$Feature,
  MeanDecreaseAccuracy = round(importance_accuracy$MeanDecreaseAccuracy, 4),
  MeanDecreaseGini = round(importance_gini$MeanDecreaseGini, 4),
  Rank_by_Accuracy = 1:nrow(importance_accuracy)
)

# 4. OOB Error Rate Over Trees
error_evolution <- data.frame(
  Tree = 1:rf_model$ntree,
  OOB_Error = rf_model$err.rate[, "OOB"] * 100
)

# ===== SAVE ALL RESULTS TO EXCEL =====
dir.create(here("Results", "OP_feature_analysis"), recursive = TRUE, showWarnings = FALSE)

all_results <- list(
  "Model_Info" = model_info,
  "Accuracy_Summary" = accuracy_summary,
  "Error_Rate_by_Model" = error_rate_df,
  "Confusion_Matrix" = confusion_df,
  "Variable_Importance" = variable_importance,
  "Error_Evolution" = error_evolution
)

write_xlsx(all_results, 
           path = here("Results", "OP_feature_analysis", "Random_Forest_Results.xlsx"))

cat("\n========================================\n")
cat("RANDOM FOREST RESULTS SAVED\n")
cat("========================================\n")
cat("File:", here("Results", "OP_feature_analysis", "Random_Forest_Results.xlsx"), "\n")
cat("Overall Accuracy:", round(accuracy * 100, 2), "%\n")
cat("OOB Accuracy:", round(oob_accuracy * 100, 2), "%\n")
cat("Top 3 Important Features:\n")
print(head(importance_accuracy, 3))
cat("========================================\n")
