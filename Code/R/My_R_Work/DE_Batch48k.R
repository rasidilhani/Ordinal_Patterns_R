library(dplyr)
library(tidyverse)
library(ggplot2)
library(readxl)

DE_Batch <- read_excel("Data/csv/Batch_48k_Results_DE_Time.xlsx")
View(DE_Batch) 

ggplot(data=DE_Batch, mapping = aes(DE_Batch$H,DE_Batch$C)) + 
  geom_point(color="red", size=1) +
  geom_errorbarh(data=DE_Batch, mapping = aes(xmin=DE_Batch$H-DE_Batch$SemiLength_H,
                                        xmax=DE_Batch$H+DE_Batch$SemiLength_H, height=0.001)
  ) +
  geom_errorbar(data=DE_Batch, mapping=aes(ymin=DE_Batch$C-DE_Batch$SemiLength_C,
                                     ymax=DE_Batch$C+DE_Batch$SemiLength_C, width=0.1)
  ) +
  coord_cartesian(xlim = c(0.295, 0.945), ylim = c(0.025, 0.4275)) +
  labs(title = "Confidence interval for entropy and complexity") + 
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))



library(readxl)
library(ggplot2)
library(dplyr)

# Read the Excel file
DE_Batch <- read_excel("Data/csv/Batch_48k_Results_DE_Time.xlsx")

# ML is a factor
DE_Batch$ML <- as.factor(DE_Batch$ML)

# Confidence interval plot
ggplot(DE_Batch, aes(x = DE_Batch$H, y = DE_Batch$C, color = ML)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = DE_Batch$H - DE_Batch$SemiLength_H,
                     xmax = DE_Batch$H + DE_Batch$SemiLength_H),
                 height = 0.005, alpha = 0.5) +
  geom_errorbar(aes(ymin = DE_Batch$C - DE_Batch$SemiLength_C,
                    ymax = DE_Batch$C + DE_Batch$SemiLength_C),
                width = 0.005, alpha = 0.5) +
  coord_cartesian(xlim = c(0.295, 0.945), ylim = c(0.025, 0.4275)) +
  labs(title = "Confidence Interval for Entropy and Complexity",
       subtitle = "Grouped by Motor Load",
       x = expression(italic(H)),
       y = expression(italic(C)),
       color = "Motor Load") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")
             
