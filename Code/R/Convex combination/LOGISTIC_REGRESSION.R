# Install packages if not already
install.packages(c("nnet", "car", "readxl", "caret"))

# Load libraries
library(nnet)
library(car)
library(VGAM)
library(readxl)
library(dplyr)
library(here)
library(caret)

D <- 4
n <- 10000


hc_data_path <- here("Data", "Convex_combination", paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))

output_dir <- here("Results", "Convex_combination", "Multinom_Detailed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Read Excel file
df <- read_excel(hc_data_path)

# 2. Convert Model to a factor 
df$Model <- as.factor(df$Model)

# 3. Fit multinomial logistic regression
model <- multinom(Model ~ 
                    H_Shannon + H_Renyi + H_Tsallis + H_Fisher +
                    C_Shannon + C_Renyi + C_Tsallis + C_Fisher +
                    Disequilibrium +
                    Var_H_Shannon + Var_H_Renyi + Var_H_Tsallis +
                    Var_H_Fisher + Var_C_Shannon,
                  data = df)

# 4. Check multicollinearity (VIF)
lm_test <- glm(
  H_Shannon ~ 
    H_Renyi + H_Tsallis + H_Fisher +
    C_Shannon + C_Renyi + C_Tsallis + C_Fisher +
    Disequilibrium +
    Var_H_Shannon + Var_H_Renyi + Var_H_Tsallis +
    Var_H_Fisher + Var_C_Shannon,
  data = df
)

vif(lm_test)   

# 5. Predictions & confusion matrix
pred <- predict(model, df)
confusionMatrix(as.factor(pred), df$Model)

#------------------------------------------------------------------------------------------
# The above code gave me extreamly large VIF results. It means we can cannot use that model.Therefore, 
# I decided to reduce the number of  predictors as those are severely multicollinear.
# They are nearly perfectly correlated with each other. 
# Multicollinearity causes:
# Unstable model coefficients
# Very large standard errors
# Impossible to interpret model effects
# Very poor generalization
# A multinomial model cannot work properly with this amount of collinearity.


library(readxl)
library(caret)
library(here)

# -------------------------------------
# 1. Load your dataset using your path
# -------------------------------------

D <- 4

hc_data_path <- here("Data", "Convex_combination",
                     paste0("D", D, "_Data"),
                     paste0("HC_Results_D", D, ".xlsx"))

output_dir <- here("Results", "Convex_combination", "Multinom_Detailed")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_excel(hc_data_path)

# Convert Model to factor
df$Model <- as.factor(df$Model)

# -------------------------------------
# 2. Choose a small set of predictors
#    (Option B: low/mild collinearity)
# -------------------------------------
# Recommended stable subset:
#   - H_Shannon     : entropy
#   - C_Fisher      : complexity
#   - Var_H_Tsallis : variability
#   - Disequilibrium (optional, often useful)

predictors_formula <- Model ~ H_Shannon + C_Fisher + Var_H_Tsallis + Disequilibrium

# -------------------------------------
# 3. Fit multinomial logistic regression
# -------------------------------------
model <- multinom(predictors_formula, data = df)

summary(model)

# -------------------------------------
# 4. Evaluate model performance
# -------------------------------------
pred <- predict(model, df)
cm <- confusionMatrix(as.factor(pred), df$Model)
cm

# Save confusion matrix if needed
write.csv(as.data.frame(cm$table),
          file = file.path(output_dir, "Confusion_Matrix.csv"),
          row.names = FALSE)
