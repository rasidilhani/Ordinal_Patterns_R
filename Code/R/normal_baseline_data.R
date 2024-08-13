# Load required packages
library(statcomp)
library(R.matlab)

# Step 2: Read MATLAB files from the folder
folder_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Data/mat/Normal_Baseline_Data"
mat_files <- list.files(path = folder_path, pattern = "\\.mat$", full.names = TRUE)

# Print files to check path and names
print(mat_files)

# Step 3: Create a list or dataframe to store time series data
time_series_data <- list()

# Load each MATLAB file and extract BE_time and FE_time
for (file in mat_files) {
  mat_data <- readMat(file)
  
  # Print variable names to check content
  print(names(mat_data))
  
  # Check if BE_time and FE_time exist in the file
  if (!("BE_time" %in% names(mat_data)) || !("FE_time" %in% names(mat_data))) {
    warning(paste("Variables BE_time or FE_time not found in file:", file))
    next
  }
  
  be_time <- mat_data$BE_time
  fe_time <- mat_data$FE_time
  
  # Print the first few elements to inspect
  print(head(be_time))
  print(head(fe_time))
  
  # Store the data
  time_series_data[[file]] <- list(BE_time = be_time, FE_time = fe_time)
}

# Initialize result storage
results <- list()

# Define embedding dimension
ndemb <- 3

# Step 4 to Step 9: Process each time series
for (file in names(time_series_data)) {
  be_time <- time_series_data[[file]]$BE_time
  fe_time <- time_series_data[[file]]$FE_time
  
  # Handle missing values by removing them
  be_time <- na.omit(be_time)
  fe_time <- na.omit(fe_time)
  
  # Ensure the data is not empty after removing NA
  if (length(be_time) < 1 || length(fe_time) < 1) {
    warning(paste("Insufficient data after NA removal in file:", file))
    next
  }
  
  # Step 4: Calculate ordinal pattern distribution
  be_opd <- ordinal_pattern_distribution(be_time, ndemb = ndemb)
  fe_opd <- ordinal_pattern_distribution(fe_time, ndemb = ndemb)
  
  # Step 5: Calculate permutation entropy
  tryCatch({
    be_pe <- permutation_entropy(be_opd)
    fe_pe <- permutation_entropy(fe_opd)
  }, error = function(e) {
    warning(paste("Error in permutation_entropy for file:", file, "Error message:", e$message))
    be_pe <- NA
    fe_pe <- NA
  })
  
  # Step 6: Calculate Jensen-Shannon entropy using statcomp
  tryCatch({
    be_js_entropy <- jensen_shannon_entropy(be_opd)
    fe_js_entropy <- jensen_shannon_entropy(fe_opd)
  }, error = function(e) {
    warning(paste("Error in jensen_shannon_entropy for file:", file, "Error message:", e$message))
    be_js_entropy <- NA
    fe_js_entropy <- NA
  })
  
  # Step 7: Calculate complexity plane
  be_histogram <- nbitflips(be_time)
  fe_histogram <- nbitflips(fe_time)
  
  # Step 8: Get complexity plane
  be_cp <- complexity_plane(be_histogram)
  fe_cp <- complexity_plane(fe_histogram)
  
  # Step 9: Compute limit curves
  be_limit_curve <- limit_curves(be_cp)
  fe_limit_curve <- limit_curves(fe_cp)
  
  # Store results
  results[[file]] <- list(
    BE_time = list(
      Ordinal_Pattern_Distribution = be_opd,
      Permutation_Entropy = be_pe,
      Jensen_Shannon_Entropy = be_js_entropy,
      Complexity_Plane = be_cp,
      Limit_Curve = be_limit_curve
    ),
    FE_time = list(
      Ordinal_Pattern_Distribution = fe_opd,
      Permutation_Entropy = fe_pe,
      Jensen_Shannon_Entropy = fe_js_entropy,
      Complexity_Plane = fe_cp,
      Limit_Curve = fe_limit_curve
    )
  )
}

# Example plot for one of the results
if (!is.na(results[[1]]$BE_time$Limit_Curve)) {
  plot(results[[1]]$BE_time$Limit_Curve, main = "Limit Curve for BE_time")
} else {
  warning("No valid limit curve data to plot for BE_time.")
}