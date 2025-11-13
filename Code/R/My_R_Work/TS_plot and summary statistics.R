################################################################################
#AR(2) and AR(1) model time series plots and summary statistics
################################################################################    
library(ggplot2)
library(dplyr)
library(StatOrdPattHxC)
set.seed(1234567890)
N <- c(500, 1000)
ar2_models <- list(
  AR2_M1 = list(ar = c(0.1, 0.8), type = "M1"),
  AR2_M2 = list(ar = c(-0.8, 0.1), type = "M2"),
  AR2_M3 = list(ar = c(0.1, -0.8), type = "M3"),
  AR2_M4 = list(ar = c(-0.8, -0.1), type = "M4")
)
ar1_models <- list(
  AR1_M1 = list(ar = c(0.8), type = "M1"),
  AR1_M2 = list(ar = c(0.1), type = "M2"), 
  AR1_M3 = list(ar = c(-0.8), type = "M3"),
  AR1_M4 = list(ar = c(-0.1), type = "M4")
)
# Function to generate and summarize AR(2) series
generate_and_summarize_ar2 <- function(model_info, n) {
  ts_data <- arima.sim(model = list(ar = model_info$ar), n = n)
  ts_df <- data.frame(Index = 1:length(ts_data), Value = as.numeric(ts_data), Model = model_info$type)
  summary_stats <- data.frame(
    Model = model_info$type,
    Mean = mean(ts_data),
    SD = sd(ts_data),
    Min = min(ts_data),
    Max = max(ts_data)
  )
  return(list(ts_df = ts_df, summary_stats = summary_stats))
}
# Function to generate and summarize AR(1) series
generate_and_summarize_ar1 <- function(model_info, n) {
  ts_data <- arima.sim(model = list(ar = model_info$ar), n = n)
  ts_df <- data.frame(Index = 1:length(ts_data), Value = as.numeric(ts_data), Model = model_info$type)
  summary_stats <- data.frame(
    Model = model_info$type,
    Mean = mean(ts_data),
    SD = sd(ts_data),
    Min = min(ts_data),
    Max = max(ts_data)
  )
  return(list(ts_df = ts_df, summary_stats = summary_stats))
}
# Generate data for each AR(2) model
results_ar2 <- lapply(ar2_models, generate_and_summarize_ar2, n = 1000)
# Generate data for each AR(1) model
results_ar1 <- lapply(ar1_models, generate_and_summarize_ar1, n = 1000)
# Combine time series and summary tables for AR(2)
ts_dfs_ar2 <- do.call(rbind, lapply(results_ar2, function(x) x$ts_df))
stats_all_ar2 <- do.call(rbind, lapply(results_ar2, function(x) x$summary_stats))
# Combine time series and summary tables for AR(1)
ts_dfs_ar1 <- do.call(rbind, lapply(results_ar1, function(x) x$ts_df))
stats_all_ar1 <- do.call(rbind, lapply(results_ar1, function(x) x$summary_stats))
# Time series plot for AR(2)
p_ar2 <- ggplot(ts_dfs_ar2, aes(x = Index, y = Value, color = Model)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Model, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(title = "Simulated AR(2) Time Series (M1–M4)",
       x = "Time Index",
       y = "Series Value")
# Time series plot for AR(1)
p_ar1 <- ggplot(ts_dfs_ar1, aes(x = Index, y = Value, color = Model)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Model, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 12, base_family = "serif") +
  labs(title = "Simulated AR(1) Time Series (M1–M4)",
       x = "Time Index",
       y = "Series Value")
# Summary statistics for AR(2)
print(stats_all_ar2)
# Summary statistics for AR(1)
print(stats_all_ar1)
# Save summary statistics to CSV
write.csv(stats_all_ar2, "ar2_summary_statistics.csv", row.names = FALSE)
write.csv(stats_all_ar1, "ar1_summary_statistics.csv", row.names = FALSE)
# Save time series data to CSV
write.csv(ts_dfs_ar2, "ar2_time_series_data.csv", row.names = FALSE)
write.csv(ts_dfs_ar1, "ar1_time_series_data.csv", row.names = FALSE)
# Save the plots
ggsave("ar2_time_series_plot.pdf", plot = p_ar2, width = 10, height = 6)
ggsave("ar1_time_series_plot.pdf", plot = p_ar1, width = 10, height = 6)
# End of code 

################################################################################
#MA(1) and MA(2) model time series plots and summary statistics   
################################################################################
library(ggplot2)
library(dplyr)
set.seed(1234567890)
N <- c(500, 1000)
ma2_models <- list(
  MA2_M1 = list(ma = c(0.1, 0.8), type = "M1"),
  MA2_M2 = list(ma = c(-0.8, 0.1), type = "M2"),
  MA2_M3 = list(ma = c(0.1, -0.8), type = "M3"),
  MA2_M4 = list(ma = c(-0.8, -0.1), type = "M4")
)
ma1_models <- list(
  MA1_M1 = list(ma = c(0.8), type = "M1"),
  MA1_M2 = list(ma = c(0.1), type = "M2"), 
  MA1_M3 = list(ma = c(-0.8), type = "M3"),
  MA1_M4 = list(ma = c(-0.1), type = "M4")
)
# Function to generate and summarize MA(2) series
generate_and_summarize_ma2 <- function(model_info, n) {
  ts_data <- arima.sim(model = list(ma = model_info$ma), n = n)
  ts_df <- data.frame(Index = 1:length(ts_data), Value = as.numeric(ts_data), Model = model_info$type)
  summary_stats <- data.frame(
    Model = model_info$type,
    Mean = mean(ts_data),
    SD = sd(ts_data),
    Min = min(ts_data),
    Max = max(ts_data)
  )
  return(list(ts_df = ts_df, summary_stats = summary_stats))
}
# Function to generate and summarize MA(1) series
generate_and_summarize_ma1 <- function(model_info, n) {
  ts_data <- arima.sim(model = list(ma = model_info$ma), n = n)
  ts_df <- data.frame(Index = 1:length(ts_data), Value = as.numeric(ts_data), Model = model_info$type)
  summary_stats <- data.frame(
    Model = model_info$type,
    Mean = mean(ts_data),
    SD = sd(ts_data),
    Min = min(ts_data),
    Max = max(ts_data)
  )
  return(list(ts_df = ts_df, summary_stats = summary_stats))
}
# Generate data for each MA(2) model
results_ma2 <- lapply(ma2_models, generate_and_summarize_ma2, n = 1000)
# Generate data for each MA(1) model
results_ma1 <- lapply(ma1_models, generate_and_summarize_ma1, n = 1000)
# Combine time series and summary tables for MA(2)
ts_dfs_ma2 <- do.call(rbind, lapply(results_ma2, function(x) x$ts_df))
stats_all_ma2 <- do.call(rbind, lapply(results_ma2, function(x) x$summary_stats))
# Combine time series and summary tables for MA(1)
ts_dfs_ma1 <- do.call(rbind, lapply(results_ma1, function(x) x$ts_df))
stats_all_ma1 <- do.call(rbind, lapply(results_ma1, function(x) x$summary_stats))
# Time series plot for MA(2)
p_ma2 <- ggplot(ts_dfs_ma2, aes(x = Index, y = Value, color = Model)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Model, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(title = "Simulated MA(2) Time Series (M1–M4)",
       x = "Time Index",
       y = "Series Value")
# Time series plot for MA(1)
p_ma1 <- ggplot(ts_dfs_ma1, aes(x = Index, y = Value, color = Model)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Model, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 12, base_family = "serif") +
  labs(title = "Simulated MA(1) Time Series (M1–M4)",
       x = "Time Index",
       y = "Series Value")
# Summary statistics for MA(2)
print(stats_all_ma2)
# Summary statistics for MA(1)
print(stats_all_ma1)
# Save summary statistics to CSV
write.csv(stats_all_ma2, "ma2_summary_statistics.csv", row.names = FALSE)
write.csv(stats_all_ma1, "ma1_summary_statistics.csv", row.names = FALSE)
# Save time series data to CSV
write.csv(ts_dfs_ma2, "ma2_time_series_data.csv", row.names = FALSE)
write.csv(ts_dfs_ma1, "ma1_time_series_data.csv", row.names = FALSE)
# Save the plots
ggsave("ma2_time_series_plot.pdf", plot = p_ma2, width = 10, height = 6)
ggsave("ma1_time_series_plot.pdf", plot = p_ma1, width = 10, height = 6)
# End of code

################################################################################
#ARMA(1,1) and ARMA(2,2) model time series plots and summary statistics   
################################################################################
library(ggplot2)
library(dplyr)
set.seed(1234567890)
N <- c(500, 1000)
arma11_models <- list(
  ARMA11_M1 = list(ar = c(0.8), ma = c(0.8), type = "M1"),
  ARMA11_M2 = list(ar = c(0.1), ma = c(0.1), type = "M2"),
  ARMA11_M3 = list(ar = c(-0.8), ma = c(-0.8), type = "M3"),
  ARMA11_M4 = list(ar = c(-0.1), ma = c(-0.1), type = "M4")
)
arma22_models <- list(
  ARMA22_M1 = list(ar = c(0.1, 0.8), ma = c(0.1, 0.8), type = "M1"),
  ARMA22_M2 = list(ar = c(-0.8, 0.1), ma = c(-0.8, 0.1), type = "M2"),
  ARMA22_M3 = list(ar = c(0.1, -0.8), ma = c(0.1, -0.8), type = "M3"),
  ARMA22_M4 = list(ar = c(-0.8, -0.1), ma = c(-0.8, -0.1), type = "M4")
)
# Function to generate and summarize ARMA(1,1) series
generate_and_summarize_arma11 <- function(model_info, n) {
  ts_data <- arima.sim(model = list(ar = model_info$ar, ma = model_info$ma), n = n)
  ts_df <- data.frame(Index = 1:length(ts_data), Value = as.numeric(ts_data), Model = model_info$type)
  summary_stats <- data.frame(
    Model = model_info$type,
    Mean = mean(ts_data),
    SD = sd(ts_data),
    Min = min(ts_data),
    Max = max(ts_data)
  )
  return(list(ts_df = ts_df, summary_stats = summary_stats))
}
# Function to generate and summarize ARMA(2,2) series
generate_and_summarize_arma22 <- function(model_info, n) {
  ts_data <- arima.sim(model = list(ar = model_info$ar, ma = model_info$ma), n = n)
  ts_df <- data.frame(Index = 1:length(ts_data), Value = as.numeric(ts_data), Model = model_info$type)
  summary_stats <- data.frame(
    Model = model_info$type,
    Mean = mean(ts_data),
    SD = sd(ts_data),
    Min = min(ts_data),
    Max = max(ts_data)
  )
  return(list(ts_df = ts_df, summary_stats = summary_stats))
}
# Generate data for each ARMA(1,1) model
results_arma11 <- lapply(arma11_models, generate_and_summarize_arma11, n = 1000)
# Generate data for each ARMA(2,2) model
results_arma22 <- lapply(arma22_models, generate_and_summarize_arma22, n = 1000)
# Combine time series and summary tables for ARMA(1,1)
ts_dfs_arma11 <- do.call(rbind, lapply(results_arma11, function(x) x$ts_df))
stats_all_arma11 <- do.call(rbind, lapply(results_arma11, function(x) x$summary_stats))
# Combine time series and summary tables for ARMA(2,2)
ts_dfs_arma22 <- do.call(rbind, lapply(results_arma22, function(x) x$ts_df))
stats_all_arma22 <- do.call(rbind, lapply(results_arma22, function(x) x$summary_stats))
# Time series plot for ARMA(1,1)
p_arma11 <- ggplot(ts_dfs_arma11, aes(x = Index, y = Value, color = Model)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Model, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(title = "Simulated ARMA(1,1) Time Series (M1–M4)",
       x = "Time Index",
       y = "Series Value")
# Time series plot for ARMA(2,2)
p_arma22 <- ggplot(ts_dfs_arma22, aes(x = Index, y = Value, color = Model)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ Model, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 12, base_family = "serif") +  
  labs(title = "Simulated ARMA(2,2) Time Series (M1–M4)",
       x = "Time Index",
       y = "Series Value")
# Summary statistics for ARMA(1,1)
print(stats_all_arma11)
# Summary statistics for ARMA(2,2)
print(stats_all_arma22)
# Save summary statistics to CSV
write.csv(stats_all_arma11, "arma11_summary_statistics.csv", row.names = FALSE)
write.csv(stats_all_arma22, "arma22_summary_statistics.csv", row.names = FALSE)
# Save time series data to CSV
write.csv(ts_dfs_arma11, "arma11_time_series_data.csv", row.names = FALSE)
write.csv(ts_dfs_arma22, "arma22_time_series_data.csv", row.names = FALSE)
# Save the plots
ggsave("arma11_time_series_plot.pdf", plot = p_arma11, width = 10, height = 6)
ggsave("arma22_time_series_plot.pdf", plot = p_arma22, width = 10, height = 6)
# End of code 



