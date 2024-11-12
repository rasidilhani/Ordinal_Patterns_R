# Required libraries
library(rprojroot)
library(readr)
library(dplyr)
library(tidyverse)

library(StatOrdPattHxC)

data("LinfLsup")

# Define the base path of the project using rprojroot
base_path <- find_root(has_file("README.md"))
data_path <- file.path(base_path, "Data", "csv")

# create Data/csv if not exists. Please copy all the csv files here in this directory
if(!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
}

csv_save_path = file.path(data_path, "HC_Points")

if(!dir.exists(csv_save_path)) {
  dir.create(csv_save_path, recursive = TRUE)
}

dataset_type = "Normal_baseline"

# Build the relative path to the Excel file
file_path <- file.path(base_path, "Data", "csv", sprintf("%s_data.csv", dataset_type))

read_data_file <- function(file_path){
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }
  
  # Read the CSV file
  data_normal <- read_csv(file_path, show_col_types = FALSE)
  
  # Select specific columns for analysis
  selected_columns <- data_normal %>% 
    select(Fault_diameter, Motor_load, RPM, DE_time, FE_time)
  
  # Example: Print the first few rows of the selected columns
  #print(head(selected_columns))
  return (selected_columns)
  
}

df <- read_data_file(file_path)

head(df)

df_keys <- list("DE_time", "FE_time")
#df_keys <- list("DE_time")

# check for NA values for different Fields
sum(is.na(df$DE_time))
sum(is.na(df$FE_time))
#sum(is.na(df$BA_time))

# remove na values
#df = df[rowSums(is.na(df)) == 0, ] 

#head(df)

# ML = 2
# D = 3

for(ML in 0:3){
  for(D in 5:6){
    
    de_time_data <- df[df$Motor_load == ML, ]$DE_time
    fe_time_data <- df[df$Motor_load == ML, ]$FE_time
    
    # Matrix that stores the variances
    Variances <- matrix(nrow=2, ncol=1)
    
    #chunk_limit = length(de_time_data) * 0.1 # 10% from the data length
    
    Variances[1,1] <- sigma2q(de_time_data, emb = D, ent = "S")
    Variances[2,1] <- sigma2q(fe_time_data, emb = D, ent = "S")
    
    Variances[Variances<0] <- 0
    
    x1sub <- de_time_data
    x2sub <- fe_time_data
    
    ShannonEntropies <- c(
      HShannon(OPprob(x1sub, emb=D)),
      HShannon(OPprob(x2sub, emb=D))
    ) 
    
    StatisticalComplexities <- c(
      StatComplexity(OPprob(x1sub, emb=D)),
      StatComplexity(OPprob(x2sub, emb=D))
    )
    
    alpha <- 0.05
    StandardDeviations <- sqrt(Variances[,1])
    SemiLength <- StandardDeviations/sqrt(length(x1sub)-D)*qnorm(1-alpha/2) 
    # The three time series have the same length, but they could be different
    
    features_c <- c("DE time", "FE time")
    
    
    HCPoints <- data.frame(H=ShannonEntropies,
                           C=StatisticalComplexities,
                           STD=StandardDeviations,
                           SemiLength=SemiLength,
                           Series=as.factor(features_c))
    
    write.csv(HCPoints, file.path(csv_save_path, sprintf("HCPoints_%s_ML_%s_D_%s.csv", dataset_type, ML, D)))
    
  }
}
#print("done!!")