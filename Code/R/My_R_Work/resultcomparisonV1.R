library(dplyr)
library(ggplot2)
library(readxl)
library(StatOrdPattHxC)

# Read data
FullResults <- read_excel("Data/Resultscomparison.xlsx", sheet = "48k&NormalML0") %>% 
  mutate(Dataset = "FullResults")
Batch48k <- read_excel("Data/Resultscomparison.xlsx", sheet = "Batch results_48k") %>% 
  mutate(Dataset = "Batch48k")
BatchNormal <- read_excel("Data/Resultscomparison.xlsx", sheet = "Batch Results_Normal") %>% 
  mutate(Dataset = "BatchNormal")

# Combine datasets
combined <- bind_rows(FullResults, Batch48k, BatchNormal) %>%
  mutate(across(c(H, C, SemiLength_H), as.numeric)) %>%
  na.omit()

D = 3
data("LinfLsup")

# Create plot
ggplot(data = combined, aes(x = H, y = C, color = Dataset)) + 
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)), 
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)), 
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(aes(size = Dataset)) +
  geom_errorbarh(data = subset(combined, Dataset == "FullResults"), 
                 aes(xmin = H - SemiLength_H, xmax = H + SemiLength_H), height = 0.001) +
  scale_color_manual(values = c("FullResults" = "blue", "Batch48k" = "red", "BatchNormal" = "green")) +
  scale_size_manual(values = c("FullResults" = 2, "Batch48k" = 1, "BatchNormal" = 1)) +
  coord_cartesian(xlim = c(0.55, 0.845), ylim = c(0.125, 0.3275)) +
  labs(title = "Confidence interval for entropy",
       x = expression(italic(H)),
       y = expression(italic(C))) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"))


