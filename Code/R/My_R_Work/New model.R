library(ggplot2)
library(readxl)
library(dplyr)
library(StatOrdPattHxC)

data("LinfLsup")

# Read Excel file
setwd("C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Code/R/My_R_Work")
df <- read_excel("Results.xlsx")

# --- Subset for n = 1000 ---
df_1000 <- df %>% filter(n == 1000)

# --- Loop for n = 1000 ---
types <- unique(df_1000$Type)

for(t in types){
  df_sub <- df_1000 %>% filter(Type == t)
  
  p <- ggplot(df_sub, aes(x = Entropy, y = Complexity, color = Model)) +
    geom_point(size = 2) +
    #geom_line(aes(group = Model), linewidth = 0.8) +
    scale_color_manual(values = c("M1" = "red",
                                  "M2" = "blue",
                                  "M3" = "green",
                                  "M4" = "purple")) +
    coord_cartesian(xlim = c(0.85, 1), ylim = c(0, 0.3)) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Entropy–Complexity (Type =", t, ", n=1000)"),
         x = expression(italic(H)),
         y = expression(italic(C)),
         color = "Model")
  
  # Save plot for each type
  ggsave(paste0("plot_", t, "_n1000.pdf"), plot = p, width = 6, height = 4, dpi = 300)
}



###########################################################################
# n=500 seperate graphs
###########################################################################
library(ggplot2)
library(readxl)
library(dplyr)
library(StatOrdPattHxC)

data("LinfLsup")


# Read Excel file
setwd("C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Code/R/My_R_Work")
df <- read_excel("Results.xlsx")

# --- Subset for n = 500 ---
df_500 <- df %>% filter(n == 500)


# --- Loop for n = 500 ---
types <- unique(df_500$Type)

for(t in types){
  df_sub <- df_500 %>% filter(Type == t)
  
  p <- ggplot(df_sub, aes(x = Entropy, y = Complexity, color = Model)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("M1" = "red",
                                  "M2" = "blue",
                                  "M3" = "green",
                                  "M4" = "purple")) +
    coord_cartesian(xlim = c(0.85, 1), ylim = c(0, 0.3)) +
    theme_minimal(base_size = 12) +
    labs(title = paste("Entropy–Complexity (Type =", t, ", n=500)"),
         x = expression(italic(H)),
         y = expression(italic(C)),
         color = "Model")
  print(p)  
  # Save plot for each type
  ggsave(paste0("plot_", t, "_n500.pdf"), plot = p, width = 6, height = 4, dpi = 300)
}