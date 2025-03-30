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

number_of_rows <- nrow(df)

data_range <- 10000

number_of_sets <- ceiling(number_of_rows/data_range)
#print(number_of_sets)
#print(ceiling(number_of_sets))

inner_bound <- 1
outer_bound <- data_range

df_keys <- list("DE_time", "FE_time")

de_time_variance = list()
fe_time_variance = list()

idx = 1

for(i in 1:number_of_sets) {
  df[inner_bound:outer_bound, ]
  # process data here
  selected_df <- df[inner_bound:outer_bound, ]
  
  #print(nrow(selected_df))
  #print("-------")
  
  # apply aggregation function
  Normal_S_DE_3 <-aggregate(selected_df,
                            DE_time~ Motor_load,
                            FUN = sigma2q, emb=3)
  de_time_variance[idx] = Normal_S_DE_3
  print(Normal_S_DE_3)

  Normal_S_FE_3 <-aggregate(selected_df,
                            FE_time~ Motor_load,
                            FUN = sigma2q, emb=3)
  fe_time_variance[idx] = Normal_S_FE_3
  print(Normal_S_FE_3)
  
  idx = idx + 1
  
  inner_bound <- outer_bound + 1
  outer_bound <- outer_bound + data_range
}

print(de_time_variance)
print(fe_time_variance)
