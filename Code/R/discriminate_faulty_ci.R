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

# Build the relative path to the Excel file
file_path_normal <- file.path(base_path, "Data", "csv", "Normal_baseline_data.csv")
file_path_12k_DE <- file.path(base_path, "Data", "csv", "12k_drive_end_bearing_fault_data.csv")
file_path_12k_FE <- file.path(base_path, "Data", "csv", "12kFan_end_bearing_fault_data.csv")
file_path_48k_DE <- file.path(base_path, "Data", "csv", "48k_drive_end_bearing_fault_data.csv")


plot_path <- file.path(base_path, "Plots")

if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

directory_path_fault_analysis = file.path(plot_path, "fault_analysis")

if(!dir.exists(directory_path_fault_analysis)) {
  dir.create(directory_path_fault_analysis)
}

file_path <- file_path_normal
#file_path <- file_path_12k_DE
#file_path <- file_path_12k_FE
#file_path <- file_path_48k_DE

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

df <- read_data_file(file_path)
#df_12k_de <- read_data_file(file_path_12k_DE)
#df_12k_fe <- read_file(file_path_12k_FE)
#df_48k_de <- read_file(file_path_48k_DE)

df_keys <- list("DE_time", "FE_time", "BA_time")
plot_colors <- list("green", "orange", "purple")
#df_keys <- list("DE_time")

head(df)

# define dataframe to be compared with
df2 = NULL


ML = 0
D = 3

######## Save H, C pints to csv file #########
dataset_type = "12k_DE"

for(ML in 0:3){
  for(D in 3:6){
    plot_title = sprintf("Confidence intervel - Motor load: %s, Embed dim: %s", ML, D)
    
    de_time_data <- df[df$Motor_load == ML, ]$DE_time
    fe_time_data <- df[df$Motor_load == ML, ]$FE_time
    ba_time_data <- df[df$Motor_load == ML, ]$BA_time
    
    # Matrix that stores the variances
    Variances <- matrix(nrow=2, ncol=1)
    
    if(sum(!is.na(ba_time_data)) > 0) {
      Variances <- matrix(nrow=3, ncol=1)
    }
    
    #chunk_limit = length(de_time_data) * 0.1 # 10% from the data length
    
    Variances[1,1] <- sigma2q(de_time_data, emb = D, ent = "S")
    Variances[2,1] <- sigma2q(fe_time_data, emb = D, ent = "S")
    
    if(sum(!is.na(ba_time_data)) > 0) {
      Variances[3,1] <- sigma2q(ba_time_data, emb = D, ent = "S")
    }
    
    Variances[Variances<0] <- 0
    
    x1sub <- de_time_data
    x2sub <- fe_time_data
    x3sub <- ba_time_data
    
    if(sum(!is.na(x3sub)) > 0) {
      ShannonEntropies <- c(
        HShannon(OPprob(x1sub, emb=D)),
        HShannon(OPprob(x2sub, emb=D)),
        HShannon(OPprob(x3sub, emb=D))
      )
      
      StatisticalComplexities <- c(
        StatComplexity(OPprob(x1sub, emb=D)),
        StatComplexity(OPprob(x2sub, emb=D)),
        StatComplexity(OPprob(x3sub, emb=D))
      )
    }else{
      ShannonEntropies <- c(
        HShannon(OPprob(x1sub, emb=D)),
        HShannon(OPprob(x2sub, emb=D))
      ) 
      
      StatisticalComplexities <- c(
        StatComplexity(OPprob(x1sub, emb=D)),
        StatComplexity(OPprob(x2sub, emb=D))
      )
    }
    
    alpha <- 0.05
    StandardDeviations <- sqrt(Variances[,1])
    SemiLength <- StandardDeviations/sqrt(length(x1sub)-D)*qnorm(1-alpha/2) 
    # The three time series have the same length, but they could be different
    
    features_c <- c("DE time", "FE time")
    
    if(sum(!is.na(ba_time_data)) > 0) {
      features_c <- c("DE time", "FE time", "BA_time")
    }
    
    HCPoints <- data.frame(H=ShannonEntropies,
                           C=StatisticalComplexities,
                           STD=StandardDeviations,
                           SemiLength=SemiLength,
                           Series=as.factor(features_c))
    
    write.csv(HCPoints, file.path(base_path, "Data", "csv", sprintf("HCPoints_%s_ML_%s_D_%s.csv", dataset_type, ML, D)))
    
    xlim_threshold = 0.05
    xlim_min_point = min(HCPoints$H)
    xlim_max_point = max(HCPoints$H)
    
    pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("confidence_interval_%s_Embed_dim_%s_Motor_load_%s", dataset_type, D, ML), ".pdf"))
    
    xlim_min_point = xlim_min_point - xlim_threshold
    xlim_max_point = xlim_max_point + xlim_threshold
    
    
    ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
           aes(x=H, y=C)) +
      geom_line() +
      geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
                aes(x=H, y=C)) +
      xlab(expression(italic(H))) +
      ylab(expression(italic(C))) +
      geom_point(data=HCPoints, aes(x=H, y=C, col=Series)) +
      geom_errorbarh(data=HCPoints, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Series, col=Series)) +
      if(is.null(df2)) {
        geom_point(data=HCPoints2, aes(x=H, y=C, col=Series)) +
          geom_errorbarh(data=HCPoints2, aes(xmin=H-SemiLength2, xmax=H+SemiLength2, group=Series, col=Series))
      } +
      ggtitle(plot_title) +
      coord_cartesian(xlim=c(xlim_min_point, xlim_max_point), ylim=c(0, 0.4))
    
    ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)
  }
}

######## End saving H, C data into csv files ###########


######## create graphs using saved H,C points ############
ML = 0
D = 3

dataset_type = "Normal"
#dataset_type = "12k_DE"
#dataset_type = "12k_FE"
#dataset_type = "48k_DE"

for(ML in 0:3){
  for(D in 3:6){
    plot_title = sprintf("Confidence intervel - Motor load: %s, Embed dim: %s", ML, D)
    
    HCPoints_csv = file.path(base_path, "Data", "csv", sprintf("HCPoints_%s_ML_%s_D_%s.csv", dataset_type, ML, D))
    HCPoints = read_csv(HCPoints_csv, show_col_types = FALSE)
    
    HCPoints
    
    SemiLength = max(HCPoints$SemiLength)
    SemiLength
    xlim_threshold = 0
    xlim_min_point = min(HCPoints$H)
    xlim_max_point = max(HCPoints$H)
    
    pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("confidence_interval_%s_Embed_dim_%s_Motor_load_%s", dataset_type, D, ML), ".pdf"))
    
    xlim_min_point = xlim_min_point - SemiLength
    xlim_max_point = xlim_max_point + SemiLength
    
    
    ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
           aes(x=H, y=C)) +
      geom_line() +
      geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
                aes(x=H, y=C)) +
      xlab(expression(italic(H))) +
      ylab(expression(italic(C))) +
      geom_point(data=HCPoints, aes(x=H, y=C, col=Series)) +
      geom_errorbarh(data=HCPoints, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Series, col=Series)) +
      ggtitle(plot_title) +
      coord_cartesian(xlim=c(xlim_min_point, xlim_max_point), ylim=c(0, 0.4))
    
    ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)
  }
}


####### End create graphs using saved H, C points ########## 

xlim_threshold = 0
xlim_min_point = min(HCPoints$H)
xlim_max_point = max(HCPoints$H)

pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("confidence_interval_Embed_dim_%s_Motor_load_%s", D, ML), ".pdf"))

# if compare with other dataframe
if(!is.null(df2)){
  directory_path_fault_analysis <- file.path(plot_path, "fault_analysis", "comparison")
  pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("Normal_with_12k_DE_Embed_dim_%s_Motor_load_%s", D, ML), ".pdf"))
  
  if(!dir.exists(directory_path_fault_analysis)) {
    dir.create(directory_path_fault_analysis)
  }
  
  de_time_data_2 <- df2[df2$Motor_load == ML, ]$DE_time
  fe_time_data_2 <- df2[df2$Motor_load == ML, ]$FE_time
  ba_time_data_2 <- df2[df2$Motor_load == ML, ]$BA_time
  
  Variances2 <- matrix(nrow=2, ncol=1)
  
  if(sum(!is.na(ba_time_data_2)) > 0) {
    Variances2 <- matrix(nrow=3, ncol=1)
  }
  
  Variances2[Variances<0] <- 0
  
  Variances2[1,1] <- sigma2q(de_time_data_2, emb = D, ent = "S")
  Variances2[2,1] <- sigma2q(fe_time_data_2, emb = D, ent = "S")
  
  if(sum(!is.na(ba_time_data_2)) > 0) {
    Variances2[3,1] <- sigma2q(ba_time_data_2, emb = D, ent = "S")
  }
  
  x1sub2 <- de_time_data_2
  x2sub2 <- fe_time_data_2
  x3sub2 <- ba_time_data_2
  
  if(sum(!is.na(x3sub2)) > 0) {
    ShannonEntropies2 <- c(
      HShannon(OPprob(x1sub2, emb=D)),
      HShannon(OPprob(x2sub2, emb=D)),
      HShannon(OPprob(x3sub2, emb=D))
    )
    
    StatisticalComplexities2 <- c(
      StatComplexity(OPprob(x1sub2, emb=D)),
      StatComplexity(OPprob(x2sub2, emb=D)),
      StatComplexity(OPprob(x3sub2, emb=D))
    )
  }else{
    ShannonEntropies2 <- c(
      HShannon(OPprob(x1sub2, emb=D)),
      HShannon(OPprob(x2sub2, emb=D))
    ) 
    
    StatisticalComplexities2 <- c(
      StatComplexity(OPprob(x1sub2, emb=D)),
      StatComplexity(OPprob(x2sub2, emb=D))
    )
  }
  
  StandardDeviations2 <- sqrt(Variances2[,1])
  SemiLength2 <- StandardDeviations2/sqrt(length(x1sub2)-D)*qnorm(1-alpha/2) 
  # The three time series have the same length, but they could be different
  
  features_c2 <- c("DE time", "FE time")
  
  if(sum(!is.na(ba_time_data_2)) > 0) {
    features_c2 <- c("DE time", "FE time", "BA_time")
  }
  
  HCPoints2 <- data.frame(H=ShannonEntropies2,
                         C=StatisticalComplexities2,
                         STD=StandardDeviations2,
                         SemiLength=SemiLength2,
                         Series=as.factor(features_c2))
  
  xlim_min_point = min(HCPoints$H, HCPoints2$H)
  xlim_max_point = max(HCPoints$H, HCPoints2$H)
}

xlim_min_point = xlim_min_point - xlim_threshold
xlim_max_point = xlim_max_point + xlim_threshold


ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
       aes(x=H, y=C)) +
  geom_line() +
  geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
            aes(x=H, y=C)) +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C))) +
  geom_point(data=HCPoints, aes(x=H, y=C, col=Series)) +
  geom_errorbarh(data=HCPoints, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Series, col=Series)) +
  if(is.null(df2)) {
    geom_point(data=HCPoints2, aes(x=H, y=C, col=Series)) +
    geom_errorbarh(data=HCPoints2, aes(xmin=H-SemiLength2, xmax=H+SemiLength2, group=Series, col=Series))
  } +
  ggtitle(plot_title) +
  coord_cartesian(xlim=c(xlim_min_point, xlim_max_point), ylim=c(0, 0.4))

ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)

for(ML in 0:3){
  for(D in 3:6){
    de_time_data <- df[df$Motor_load == ML, ]$DE_time
    fe_time_data <- df[df$Motor_load == ML, ]$FE_time
    ba_time_data <- df[df$Motor_load == ML, ]$BA_time
    
    # Matrix that stores the variances
    Variances <- matrix(nrow=2, ncol=1)
    
    if(sum(!is.na(ba_time_data)) > 0) {
      Variances <- matrix(nrow=3, ncol=1)
    }
    
    #chunk_limit = length(de_time_data) * 0.1 # 10% from the data length
    
    Variances[1,1] <- sigma2q(de_time_data, emb = D, ent = "S")
    Variances[2,1] <- sigma2q(fe_time_data, emb = D, ent = "S")
    
    if(sum(!is.na(ba_time_data)) > 0) {
      Variances[3,1] <- sigma2q(ba_time_data, emb = D, ent = "S")
    }
    
    Variances[Variances<0] <- 0
    
    x1sub <- de_time_data
    x2sub <- fe_time_data
    x3sub <- ba_time_data
    
    if(sum(!is.na(x3sub)) > 0) {
      ShannonEntropies <- c(
        HShannon(OPprob(x1sub, emb=D)),
        HShannon(OPprob(x2sub, emb=D)),
        HShannon(OPprob(x3sub, emb=D))
      )
      
      StatisticalComplexities <- c(
        StatComplexity(OPprob(x1sub, emb=D)),
        StatComplexity(OPprob(x2sub, emb=D)),
        StatComplexity(OPprob(x3sub, emb=D))
      )
    }else{
      ShannonEntropies <- c(
        HShannon(OPprob(x1sub, emb=D)),
        HShannon(OPprob(x2sub, emb=D))
      ) 
      
      StatisticalComplexities <- c(
        StatComplexity(OPprob(x1sub, emb=D)),
        StatComplexity(OPprob(x2sub, emb=D))
      )
    }
    
    alpha <- 0.05
    StandardDeviations <- sqrt(Variances[,1])
    SemiLength <- StandardDeviations/sqrt(length(x1sub)-D)*qnorm(1-alpha/2) 
    # The three time series have the same length, but they could be different
    
    features_c <- c("DE time", "FE time")
    
    if(sum(!is.na(ba_time_data)) > 0) {
      features_c <- c("DE time", "FE time", "BA_time")
    }
    
    HCPoints <- data.frame(H=ShannonEntropies,
                           C=StatisticalComplexities,
                           STD=StandardDeviations,
                           SemiLength=SemiLength,
                           Series=as.factor(features_c))
    
    xlim_min_point = min(HCPoints$H - 0.05)
    xlim_max_point = max(HCPoints$H + 0.05)
    
    pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("confidence_interval_Embed_dim_%s_Motor_load_%s", D, ML), ".pdf"))
    
    ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
           aes(x=H, y=C)) +
      geom_line() +
      geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
                aes(x=H, y=C)) +
      xlab(expression(italic(H))) +
      ylab(expression(italic(C))) +
      geom_point(data=HCPoints, aes(x=H, y=C, col=Series)) +
      geom_errorbarh(data=HCPoints, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Series, col=Series)) +
      coord_cartesian(xlim=c(xlim_min_point, xlim_max_point), ylim=c(0, 0.4))
    
    ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)
  }
}



