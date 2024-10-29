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

directory_path_fault_analysis = file.path(plot_path, "fault_analysis", "final")

if(!dir.exists(directory_path_fault_analysis)) {
  dir.create(directory_path_fault_analysis, recursive = TRUE)
}

dataset_types <- list(
  "Normal_baseline" = "Normal Baseline",
  "12k_drive_end_bearing_fault" = "12K Fan-end",
  "12kFan_end_bearing_fault" = "12K Drive-end",
  "48k_drive_end_bearing_fault" = "48k Drive-end")

hc_df = data.frame()

for(dataset_type in names(dataset_types)) {
  
  ############ Generate Graphs ########################
  
  for(ML in 0:3){
    for(D in 3:4){
      
      HCPoints_csv = file.path(base_path, "Data", "csv", sprintf("HCPoints_%s_ML_%s_D_%s.csv", dataset_type, ML, D))
      HCPoints = read_csv(HCPoints_csv, show_col_types = FALSE)
      
      # add Type column to the dataframe
      HCPoints$Type <- c(dataset_types[[dataset_type]])
      
      # add Dimension column to the dataframe
      HCPoints$D <- c(D)
      
      #print(dataset_types[[dataset_type]])
      #print(HCPoints)
      
      hc_df <- rbind(hc_df, HCPoints[, -1])
      
      
    }
  }
}

print(hc_df, n=20)


############ Generate Graphs ########################

time_serires = list(
  "DE time",
  "FE time",
  "BA_time"
)

plots_list = list()
plt_ind = 1
for (t in time_serires) {
  for (D in 3:4){
    HCPoints_for_graph = hc_df[hc_df$Series == t & hc_df$D == D, ]
    
    SemiLength = max(HCPoints_for_graph$SemiLength)
    SemiLength
    xlim_threshold = 0
    xlim_min_point = min(HCPoints_for_graph$H)
    xlim_max_point = max(HCPoints_for_graph$H)
    
    pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("confidence_interval_%s_Embed_dim_%s", t, D), ".pdf"))
    
    xlim_min_point = xlim_min_point - SemiLength
    xlim_max_point = xlim_max_point + SemiLength
    
    plot_title = sprintf("Time series: %s, Embed dim: %s", t, D)
    
   plt <-  ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
           aes(x=H, y=C)) +
      geom_line() +
      geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
                aes(x=H, y=C)) +
      xlab(expression(italic(H))) +
      ylab(expression(italic(C))) +
      geom_point(data=HCPoints_for_graph, aes(x=H, y=C, col=Type)) +
      geom_errorbarh(data=HCPoints_for_graph, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Type, col=Type)) +
      ggtitle(plot_title) +
      coord_cartesian(xlim=c(xlim_min_point, xlim_max_point), ylim=c(0, 0.3)) +
      theme(panel.border = element_blank(), plot.margin = unit(c(5, 5, 5, 5), "mm"), plot.background = element_rect(color="white", size=2))
    
   plots_list[[plt_ind]] <- plt
   plt_ind = plt_ind + 1
   #ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)
  } 
}

# generate plot grid
ci_plot_path = file.path(plot_path, "fault_analysis", "final")
if(!dir.exists(ci_plot_path)) {
  dir.create(ci_plot_path, recursive = TRUE)
}

grob <- grid.arrange(grobs=plots_list, ncol = 3, nrow = 2, top="Confidence Interval", clip=TRUE)
grob_title <- sprintf("Motor Load %s", motor_load)
ggsave(file.path(ci_plot_path, "confidence_interval.pdf"), grob, width=20, height=10, units = 'in', dpi = 300, bg="white")


for(ML in 0:3){
  for(D in 3:4){
    plot_title = sprintf("Confidence intervel - Embed dim: %s", )
    
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