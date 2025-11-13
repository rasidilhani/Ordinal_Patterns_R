# Description: This script saves time series data for different ARMA models categorized by cases.
# It reads central point data and corresponding time series data, then compiles them into a new
# Excel workbook with separate sheets for each case, including the emblematic row information.

# --- Required Libraries ---
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)

# --- File Paths ---
central_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/CentralPoints.xlsx"
ts_path      <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/TimeSeries_Data_all_Models.xlsx"
save_path    <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/ARMA Time series results/TimeSeries_by_Case.xlsx"

# --- Cases Definition ---
cases <- list(
  Case1 = c("AR1_M1","AR1_M2","AR2_M1","MA1_M1","MA1_M2","MA2_M1","ARMA11_M1","ARMA11_M2","ARMA22_M1"),
  Case2 = c("AR1_M3","AR1_M4","AR2_M4","MA1_M3","MA1_M4","MA2_M4","ARMA11_M3","ARMA11_M4","ARMA22_M4"),
  Case3 = c("AR2_M2","AR2_M3","MA2_M2","MA2_M3","ARMA22_M2","ARMA22_M3")
)
sample_sizes <- c(500, 1000)

# --- Function to get time series for a given central point ---
get_ts <- function(model_name, n_val, emblematic_row) {
  sheet_name <- paste0(model_name, "_n", n_val)
  # Check if sheet exists
  sheets <- excel_sheets(ts_path)
  if(!(sheet_name %in% sheets)) return(NA)
  
  ts_data <- read_excel(ts_path, sheet = sheet_name)
  
  ts_row <- ts_data %>% filter(Rep == emblematic_row)
  if(nrow(ts_row) == 0) return(NA)
  
  # Return values as comma-separated string
  paste(ts_row$Value, collapse = ",")
}

# --- Create Workbook ---
wb <- createWorkbook()

for(case_name in names(cases)){
  df_list <- list()
  
  for(n_val in sample_sizes){
    # Load central points
    central_df <- read_excel(central_path, sheet = paste0("CentralPoints_n", n_val)) %>%
      filter(Model %in% cases[[case_name]]) %>%
      mutate(Model = as.character(Model))
    
    # Add TimeSeries column
    central_df <- central_df %>%
      rowwise() %>%
      mutate(TimeSeries = get_ts(Model, n_val, Emblematic_Row)) %>%
      ungroup()
    
    # Add column for n value
    central_df <- central_df %>%
      mutate(n = n_val)
    
    df_list[[paste0("n", n_val)]] <- central_df
  }
  
  # Combine both n values for the case
  df_case <- bind_rows(df_list)
  
  # Add sheet to workbook
  addWorksheet(wb, case_name)
  writeData(wb, sheet = case_name, df_case)
}

# --- Save Workbook ---
saveWorkbook(wb, save_path, overwrite = TRUE)

message("✅ Time series data saved by case, including Emblematic_Row and NA where missing.")
