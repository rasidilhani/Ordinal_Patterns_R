install.packages("nnet")      # For multinomial logistic regression
install.packages("dplyr")     
install.packages("readxl")    
install.packages("here")      

# Load packages
library(nnet)
library(dplyr)
library(readxl)
library(here)

here()

# read data from directory
df <- read_excel(here("Data", "OP_Features n1000_n5000", "D3 results_tau_1", "D3_n1000_n5000_tau_1.xlsx"), 
                 sheet = "D3n1000")

str(df)
head(df)

# Check unique models
unique(df$Model)

model_list <- split(df, df$Model)
for (model_name in names(model_list)) {
  model_data <- model_list[[model_name]]
  
  
  
}