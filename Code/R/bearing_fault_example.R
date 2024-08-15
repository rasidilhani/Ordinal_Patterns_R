# install required packages if not installed
if(!requireNamespace("rprojroot", quietly = TRUE)){
  install.packages("rprojroot")
}

# if(!requireNamespace("readxl", quietly = TRUE)){
#   install.packages("readxl")
# }

if(!requireNamespace("reader", quietly = TRUE)){
  install.packages("reader")
}

if(!requireNamespace("dplyr", quietly = TRUE)){
  install.packages("dplyr")
}

# Required libraries
library(rprojroot)
library(statcomp)
library(entropy)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)

# Example: 
#x <- c(1.5, 2.3, 3.1, 2.9, 4.0, 3.5, 2.7, 4.2, 3.8, 5.0)

# Define the base path of the project using rprojroot
base_path <- find_root(has_file("README.md"))

# Build the relative path to the Excel file
file_path <- file.path(base_path, "Data", "csv", "Normal_baseline_data.csv")
plot_path <- file.path(base_path, "Plots")

if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

#Initialize lists to store time series data and results
DE_Time <- list()
FE_Time <- list()
BA_Time <- list()
results <- list()

read_data_file <- function(file_path){
  # Check if the file exists
  if (!file.exists(file_path)) {
    stop("The file does not exist: ", file_path)
  }
  
  # Read the CSV file
  data_normal <- read_csv(file_path, show_col_types = FALSE)
  
  # Select specific columns for analysis
  selected_columns <- data_normal %>% 
    select(Fault_diameter, Motor_load, RPM, DE_time, FE_time, BA_time)
  
  # Example: Print the first few rows of the selected columns
  #print(head(selected_columns))
  return (selected_columns)
  
}

# function for time series plot
# function for time series plot
plot_time_series <- function(dataset, title, xlabel="Index", ylabel, color="#000000") {
  pdf(file=file.path(plot_path, paste(title, ".pdf")), width = 20, height = 5)
  plot(dataset, type='l', main=title, xlab = xlabel, ylab = ylabel, col=color)
  dev.off()
}

df <- read_data_file(file_path)

df_keys <- list("DE_time", "FE_time", "BA_time")
plot_colors <- list("#1B9E77", "#D95F02", "#7570B3")
#df_keys <- list("DE_time")

for(motor_load in 0:3){
  for(i in seq_along(df_keys)) {
    df_key = df_keys[[i]]
    plot_color = plot_colors[[i]]
    
    motor_load_ind <- which(df$Motor_load == motor_load)
    motor_load_df <- df[motor_load_ind, ]
    
    if(sum(!is.na(motor_load_df[[df_key]])) > 0) {
      #print(head(motor_load_df[, df_key]))
      plot_time_series(motor_load_df[[df_key]], sprintf("%s - Motor Load %s", df_key, motor_load), ylabel = df_key, color=plot_color) 
    }
  }
}

# summary stats by Motor Load
by(df, df[, "Motor_load"], summary)


####################### 
# ordinal patterns
#######################
gen_ord_pattern_and_hist <- function(dataset, ndim=3, title, xlabel="Index", ylabel, color="#000000") {
  ord_patt <- ordinal_pattern_time_series(dataset, ndemb=ndim)
  
  png_file_path = file.path(plot_path, "histograms", paste("Embed_dim_", ndim))
  if(!dir.exists(png_file_path)) {
    dir.create(png_file_path)
  }
  
  pdf(file=file.path(png_file_path, paste("Ordinal pattern ", title, ".pdf")), width=20, height=5)
  hist(ord_patt, main=title, xlab = xlabel, ylab = ylabel, col=color)
  dev.off()
}

for(embed_dim in 3:6){
  for(motor_load in 0:3){
    for(i in seq_along(df_keys)) {
      df_key = df_keys[[i]]
      plot_color = plot_colors[[i]]
      
      motor_load_ind <- which(df$Motor_load == motor_load)
      motor_load_df <- df[motor_load_ind, ]
      
      if(sum(!is.na(motor_load_df[[df_key]])) > 0) {
        #print(head(motor_load_df[, df_key]))
        gen_ord_pattern_and_hist(motor_load_df[[df_key]], ndim = embed_dim, sprintf("%s - Motor Load %s", df_key, motor_load), ylabel = df_key, color=plot_color) 
      }
    }
  }
}

##################################### 
# Compute OPD, entropy and complexity
#####################################

opd_list = list()
S_entropy_list = c()
complexity_pe_list = c()
complexity_mpr_list = c()
motor_load_list = c()
embed_dim_list = c()
df_key_list = c()

for(motor_load in 0:3){
  for(embed_dim in 3:6){
    for(i in seq_along(df_keys)) {
      df_key = df_keys[[i]]
      
      motor_load_ind <- which(df$Motor_load == motor_load)
      motor_load_df <- df[motor_load_ind, ]
      
      if(sum(!is.na(motor_load_df[[df_key]])) > 0) {
        #print(head(motor_load_df[, df_key]))
        motor_load_list = c(motor_load_list, motor_load)
        embed_dim_list = c(embed_dim_list, embed_dim)
        df_key_list = c(df_key_list, df_key)
        
        opd <- ordinal_pattern_distribution(x=motor_load_df[[df_key]], ndemb = embed_dim)
        S_entropy <- permutation_entropy(opd)
        complexity <- global_complexity(opd = opd, ndemb = embed_dim)
        complexity_pe_list <- c(complexity_pe_list, unname(complexity)[1])
        complexity_mpr_list <- c(complexity_mpr_list, unname(complexity)[2])
        
        #opd_list[[length(opd_list)+1]] = opd
        S_entropy_list = c(S_entropy_list, S_entropy)
      }
    }
  }
}

Motor_load <- motor_load_list
Embed_dim <- embed_dim_list
DF_key <- df_key_list
#Opd <- opd_list
S_entropy <- S_entropy_list
P_entropy <- complexity_pe_list
Stat_complexity <- complexity_mpr_list

results_df <- as.data.frame (cbind(Motor_load, Embed_dim, DF_key, S_entropy, P_entropy, Stat_complexity))

results_df

write.csv(results_df, file.path(base_path, "Data", "csv", "Complexity_results1.csv"), row.names = FALSE)

#############################
# Complexity plot
#############################

min_func = list(NULL, NULL,mind3, mind4, mind5, mind6)
max_func = list(NULL, NULL,maxd3, maxd4, maxd5, maxd6)


for(em_dim in 3:6) {
  embed_dim_n_ind = which(results_df$Embed_dim == em_dim)
  embed_dim_n_df <- results_df[embed_dim_n_ind, ]
  embed_dim_n_df$DF_key <- as.factor(embed_dim_n_df$DF_key)
  embed_dim_n_df$Motor_load <- as.factor(embed_dim_n_df$Motor_load)
  embed_dim_n_df$S_entropy <- as.numeric(embed_dim_n_df$S_entropy)
  embed_dim_n_df$Stat_complexity <- as.numeric(embed_dim_n_df$Stat_complexity)
  
  min_limit_data = min_func[em_dim][[1]]
  max_limit_data = max_func[em_dim][[1]]
  
  graph_title = sprintf("Complexity Plane - Embed dimenssion %s", em_dim)

  ggplot() +
    geom_point(data=embed_dim_n_df, aes(x=S_entropy, y=Stat_complexity, colour = DF_key, shape = Motor_load)) +
    geom_line(data = max_limit_data, aes(x=x, y=y), color='red') +
    geom_line(data = min_limit_data, aes(x=x, y=y), color='blue') +
    labs(
      title = graph_title,
      x="Entropy",
      y="Complexity"
    )

  pdf_file_name = file.path(plot_path, "complexity_plane", paste(sprintf("complexity_plane_graph_em_dim_%s", em_dim), ".pdf"))
  ggsave(pdf_file_name, width=20, height = 5, units = "in", dpi = 300)
}


#a = limit_curves(ndemb = 3, fun = 'max')
#print(a)
#plot()

print(dim(c))
#plot(x=cbind(b$x[1:494], c$x[1:494]), y=cbind(b$y[1:494], c$y[1:494]), type = 'l')

comp_plot = ggplot(embed_dim_3_df, aes(x=S_entropy, y=Stat_complexity, colour = DF_key, shape = Motor_load)) + geom_point() + coord_equal(ratio = 1)
comp_plot + geom_line(data = b, aes(x=x, y=y))
comp_plot
embed_dim_3_df



#  + ggplot(data = embed_dim_3_df, aes(x=S_entropy, y=Stat_complexity, colour = DF_key, shape = Motor_load)) + geom_point()

c$x[min(embed_dim_3_df$S_entropy):max(embed_dim_3_df$S_entropy)]

#######################################
#Boundaries for the complexity plane
#############################################



# opd <- ordinal_pattern_distribution(data_normal$DE_time, ndemb = 3)
# opd
# 
# # 
# S_entropy <- permutation_entropy(opd)
# S_entropy
# 
# # # Compute Jensen-Shannon Divergence
# # js_divergence <- jensen_shannon_divergence(data_normal$DE_time, q = "unif")
# # js_divergence
# 
# 
# # Function to calculate Jensen-Shannon divergence
# js_divergence <- function(P, Q) {
#   # Ensure Q is a uniform distribution
#   uniform_Q <- rep(1 / length(Q), length(Q))
# 
#   M <- 0.5 * (P + Q)
# 
#   kl_divergence <- function(P, Q) {
#     sum(P * log2(P / Q), na.rm = TRUE)
#   }
# 
#   jsd <- 0.5 * kl_divergence(P, M) + 0.5 * kl_divergence(Q, M)
#   return(jsd)
# }
# 
# # Example data
# data_normal <- list(
#   DE_time = c(0.1, 0.2, 0.3, 0.4, 0.5)
# )
# 
# # Calculate probability distributions from your data
# probabilities <- function(data) {
#   # Compute frequencies
#   freq <- table(data)
#   # Convert frequencies to probabilities
#   prob <- freq / sum(freq)
#   return(prob)
# }
# 
# # Calculate the Jensen-Shannon divergence
# probabilities_DE_time <- probabilities(data_normal$DE_time)
# uniform_prob <- rep(1 / length(probabilities_DE_time), length(probabilities_DE_time))
# 
# js_divergence_result <- js_divergence(probabilities_DE_time, uniform_prob)
# print(js_divergence_result)
# 
# 
# 
# # Print results
# print(paste("Entropy:", S_entropy))
# print(paste("Jensen-Shannon Divergence:", js_divergence))
# 
# # complexity measure
# complexity <- global_complexity(opd, ndemb = 3)
# print(complexity)



