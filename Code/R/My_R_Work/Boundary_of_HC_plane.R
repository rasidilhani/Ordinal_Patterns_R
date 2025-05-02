library(StatOrdPattHxC)
library(ggplot2)
library(dplyr)  
library(tidyverse)

data("LinfLsup")

ggplot(
  data = LinfLsup %>% 
    filter(Dimension %in% c("3","4","5","6")),
  aes(x = H, y = C, 
      color = Dimension,  # Color by dimension
      group = interaction(Side, Dimension))
) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(
    name = "Dimension (D)",
    values = c(
      "3" = "green",  
      "4" = "red",  
      "5" = "black",  
      "6" = "blue"   
    ),
    breaks = c("3", "4", "5", "6")  
  ) +
  labs(x = expression(italic(H)), y = expression(italic(C))) +
  theme_minimal() +
  theme(
    text = element_text(family = "serif", size = 16),
    legend.position = "right",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )



