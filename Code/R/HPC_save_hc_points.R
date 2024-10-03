# Required libraries
library(rprojroot)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_set(theme_clean()+theme(legend.position = "top"))
library(StatOrdPattHxC)

data("LinfLsup")

# Define the base path of the project using rprojroot
base_path <- find_root(has_file("README.md"))
plot_path <- file.path(base_path, "Plots")
data_path <- file.path(base_path, "Data", "csv")

# create Plots path if not exists
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# create Data/csv if not exists. Please copy all the csv files here in this directory
if(!dir.exists(data_path)) {
  dir.create(data_path, recursive = TRUE)
}

#### uncomment following lines to change the data file #######
dataset_type = "Normal_baseline"
#dataset_type = "12k_drive_end_bearing_fault"
#dataset_type = "12kFan_end_bearing_fault"
#dataset_type = "48k_drive_end_bearing_fault"

# Build the relative path to the Excel file
file_path <- file.path(base_path, "Data", "csv", sprintf("%s_data.csv", dataset_type))

directory_path_fault_analysis = file.path(plot_path, "fault_analysis")

if(!dir.exists(directory_path_fault_analysis)) {
  dir.create(directory_path_fault_analysis)
}


# Check if the file exists
if (!file.exists(file_path)) {
  stop("The file does not exist: ", file_path)
}

# Read the CSV file
data_normal <- read_csv(file_path, show_col_types = FALSE)

# Select specific columns for analysis
df <- data_normal %>% 
  select(Fault_diameter, Motor_load, RPM, DE_time, FE_time, BA_time)

#df <- read_data_file(file_path)

df_keys <- list("DE_time", "FE_time", "BA_time")
plot_colors <- list("green", "orange", "purple")
#df_keys <- list("DE_time")

# check for NA values for different Fields
#sum(is.na(df$DE_time))
#sum(is.na(df$FE_time))
#sum(is.na(df$BA_time))

# remove na values
#df = df[rowSums(is.na(df)) == 0, ] 

#head(df)

print("done!!")