
# Install necessary packages if not already installed
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

# Load the packages
library(readxl)
library(dplyr)
library(openxlsx)

# Specify the path to your Excel file
directory <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/XLS"
file_path <- file.path(directory, "Normal Baseline Data.xlsx")

# List of sheet names
sheet_names <- c("Normal 0", "Normal 1", "Normal 2", "Normal 3")

# Function to read a sheet and select the necessary columns
read_sheet <- function(sheet_name) {
  data <- read_excel(file_path, sheet = sheet_name)
  print(names(data))  # Print names to check
  
  # Use correct column names based on what is printed
  data %>%
    select(Machine_Type = `Machine Type`, DE_Time = `DE_Time`, FE_Time = `FE_Time`)
}

# Read all sheets and concatenate them into one data frame
combined_data <- bind_rows(lapply(sheet_names, read_sheet))

# Write the combined data to a new sheet in the same Excel file
write.xlsx(combined_data, file_path, sheetName = "Combined", append = TRUE)

