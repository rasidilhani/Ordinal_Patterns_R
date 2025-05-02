# Required libraries
library(rprojroot)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_set(theme_clean()+theme(legend.position = "top", legend.title = element_text(size = 12, family = 'serif'), legend.text = element_text(size = 12, family = 'serif'), text = element_text(family = 'serif')))
library(StatOrdPattHxC)
library(gridExtra)

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

directory_path_fault_analysis = file.path(plot_path, "fault_analysis", "final", "batches", "complexity_ci")

if(!dir.exists(directory_path_fault_analysis)) {
  dir.create(directory_path_fault_analysis, recursive = TRUE)
}

dataset_types <- list(
  "Normal_baseline" = "Normal Baseline",
  "48k_drive_end_bearing_fault" = "48k Drive-end")

hc_df = data.frame()

batches = list(
  "0" = 255
  #"1" = 1244,
  #"2" = 1268,
  #"3" = 1266
)

for(dataset_type in names(dataset_types)) {
  
  for(ML in names(batches)) {
    number_of_batches = batches[[ML]]
    
    ############ Generate Graphs ########################
    
    for(batch in 1:number_of_batches) {
      #for(ML in 0:3){
        for(D in 3:4){
          #D = 4
          HCPoints_csv = file.path(base_path, "Data", "csv", "HC_Points", "Batches", sprintf("HCPoints_%s_ML_%s_D_%s_Batch_%s.csv", dataset_type, ML, D, batch))
          HCPoints_csv_varc = file.path(base_path, "Data", "csv", "HC_Points", "Batches", "Complexity_CI", sprintf("HCPoints_%s_ML_%s_D_%s_Batch_%s.csv", dataset_type, ML, D, batch))
          
          if(file.exists(HCPoints_csv) & file.exists(HCPoints_csv_varc)) {
            HCPoints = read_csv(HCPoints_csv, show_col_types = FALSE)
            HCPoints_varc = read_csv(HCPoints_csv_varc, show_col_types = FALSE)
            
            # merge with VarC data
            merged_df = merge(HCPoints[, -1], HCPoints_varc[, -1], by="Series")
            
            column_names_lookup <- c(H="H.x", C="C.x", STD="STD.x", SemiLength="SemiLength.x", H_vc="H.y", C_vc="C.y", STD_vc="STD.y", SemiLength_vc="SemiLength.y")
            merged_df = rename(merged_df, all_of(column_names_lookup))
            merged_df = subset(merged_df, select = -c(H_vc, C_vc))
            
            # add Batch column to the dataframe
            merged_df$Batch <- c(batch)
            
            # add Type column to the dataframe
            merged_df$Type <- c(dataset_types[[dataset_type]])
            
            # add Dimension column to the dataframe
            merged_df$D <- c(D)
            
            #print(dataset_types[[dataset_type]])
            #print(HCPoints)
            
            hc_df <- rbind(hc_df, merged_df)
          }
        }
      #}
    }
    
  }
}
hc_df

#print(hc_df, n=50)

############ Generate Graphs ########################

time_serires = list(
  "DE time",
  "FE time"
)

# TODO: add batch
plot_batches = list()
plt_batch_ind = 1
plots_list = list()
plt_ind = 1

for(ML in names(batches)) {
  number_of_batches = batches[[ML]]
  
  for(batch in 1:number_of_batches) {
    plt_ind = 1
    for (t in time_serires) {
      for (D in 3:4){
        #D = 4
        HCPoints_for_graph = hc_df[hc_df$Series == t & hc_df$D == D & hc_df$Batch == batch, ]
        
        SemiLength = max(HCPoints_for_graph$SemiLength)
        #SemiLength
        SemiLength_vc = HCPoints_for_graph$SemiLength_vc
        
        xlim_threshold = 0
        xlim_min_point = min(HCPoints_for_graph$H)
        xlim_max_point = max(HCPoints_for_graph$H)
        
        pdf_file_name <- file.path(directory_path_fault_analysis, paste(sprintf("confidence_interval_%s_Embed_dim_%s_Batch_%s", t, D, batch), ".pdf"))
        
        xlim_min_point = xlim_min_point - SemiLength
        xlim_max_point = xlim_max_point + SemiLength
        
        plot_title = sprintf("%s, Motor Load: %s, Embed dim: %s", t, ML, D)
        
        plt <-  ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
                       aes(x=H, y=C)) +
          geom_line() +
          geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
                    aes(x=H, y=C)) +
          xlab(expression(italic(H))) +
          ylab(expression(italic(C))) +
          geom_point(data=HCPoints_for_graph, aes(x=H, y=C, col=Type)) +
          #geom_errorbar(data=HCPoints_for_graph, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Type, col=Type), width=.02) +
          geom_errorbar(data=HCPoints_for_graph, aes(ymin=C-SemiLength_vc, ymax=C+SemiLength_vc, group=Type, col=Type), width=.01) +
          ggtitle(plot_title) +
          coord_cartesian(xlim=c(xlim_min_point, xlim_max_point)) +
          theme(panel.border = element_blank(), plot.margin = unit(c(5, 5, 5, 5), "mm"), plot.background = element_rect(color="white", size=2), text = element_text(family = 'serif'))
        
        plots_list[[plt_ind]] <- plt
        plt_ind = plt_ind + 1
        #ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)
      } 
    }
    plot_batches[[plt_batch_ind]] <- plots_list
    plt_batch_ind = plt_batch_ind + 1
  }
}

#, ylim=c(0.05, 0.4)

# generate plot grid
#ci_plot_path = file.path(plot_path, "fault_analysis", "final", "Batches")
#if(!dir.exists(ci_plot_path)) {
#  dir.create(ci_plot_path, recursive = TRUE)
#}

batch_ind = 1
for(batch_plots_list in plot_batches) {
  grob <- arrangeGrob(grobs=batch_plots_list, ncol = 4, nrow = 1, clip=TRUE)
  #grob_title <- sprintf("Motor Load %s", motor_load)
  ggsave(file.path(directory_path_fault_analysis, sprintf("confidence_interval_Batch_%s.pdf", batch_ind)), grob, width=18, height=3, units = 'in', dpi = 300, bg="white")
  batch_ind = batch_ind + 1  
}


for(ML in 0:3){
  for(D in 3:4){
    plot_title = sprintf("Confidence intervel - Embed dim: %s", D)
    
    HCPoints_csv = file.path(base_path, "Data", "csv", "HC_Points", sprintf("HCPoints_%s_ML_%s_D_%s.csv", dataset_type, ML, D))
    HCPoints = read_csv(HCPoints_csv, show_col_types = FALSE)
    
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
      geom_errorbarh(data=HCPoints, aes(xmin=H-SemiLength, xmax=H+SemiLength, 
                                        group=Series, col=Series)) +
      ggtitle(plot_title) +
      coord_cartesian(xlim=c(xlim_min_point, xlim_max_point), ylim=c(0, 0.3))
    
    ggsave(pdf_file_name, width=16, height = 10, units = "in", dpi = 300)
  }
}


####### End create graphs using saved H, C points ##########