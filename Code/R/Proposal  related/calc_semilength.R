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

#dataset_type = "12k_drive_end_bearing_fault"
dataset_type = "48k_drive_end_bearing_fault"

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
    #select(Fault_diameter, Motor_load, RPM, DE_time, FE_time, BA_time)
    select(Fault_diameter, Motor_load, RPM, DE_time, FE_time)
  
  # Example: Print the first few rows of the selected columns
  #print(head(selected_columns))
  return (selected_columns)
  
}

df <- read_data_file(file_path)

#df_keys <- list("DE_time", "FE_time", "BA_time")
#df_keys <- list("DE_time")

# check for NA values for different Fields
sum(is.na(df$DE_time))
sum(is.na(df$FE_time))
#sum(is.na(df$BA_time))

# remove na values
df = df[rowSums(is.na(df)) == 0, ] 

head(df)

# calculate semilength
ML = 2
D = 5

#for(ML in 0:3){
#  for(D in 3:6){

de_time_data <- df[df$Motor_load == ML, ]$DE_time
fe_time_data <- df[df$Motor_load == ML, ]$FE_time
#ba_time_data <- df[df$Motor_load == ML, ]$BA_time

x1sub <- fe_time_data

alpha <- 0.05
StandardDeviations <- 0.371824636


SemiLength <- StandardDeviations/sqrt(length(x1sub)-D)*qnorm(1-alpha/2) 

print(sprintf("ML: %s, D: %s, %s", ML, D, SemiLength))
