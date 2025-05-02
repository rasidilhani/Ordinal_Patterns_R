library(dplyr)
library(ggplot2)
library(readxl)

Normal_48k <- read_excel("Data/Resultscomparison.xlsx", sheet = "48k&NormalML0") %>% 
  mutate(Dataset = "Normal_48k")
Batch48k <- read_excel("Data/Resultscomparison.xlsx", sheet = "Batch results_48k") %>% 
  mutate(Dataset = "Batch48k")
BatchNormal <- read_excel("Data/Resultscomparison.xlsx", sheet = "Batch Results_Normal") %>% 
  mutate(Dataset = "BatchNormal")

combined <- bind_rows(Normal_48k, Batch48k, BatchNormal) %>%
  mutate(across(c(H, C, SemiLength_H), as.numeric)) %>%
  na.omit()

ggplot(data = combined, aes(x = H, y = C, color = Dataset)) + 
  geom_point(aes(size = Dataset)) +
  geom_errorbarh(data = subset(combined, Dataset == "Normal_48k"), 
                 aes(xmin = H - SemiLength_H, xmax = H + SemiLength_H), height = 0.001) +
  scale_color_manual(values = c("Normal_48k" = "blue", "Batch48k" = "red", "BatchNormal" = "green")) +
  scale_size_manual(values = c("Normal_48k" = 2, "Batch48k" = 1, "BatchNormal" = 1)) +
  coord_cartesian(xlim = c(0.55, 0.845), ylim = c(0.125, 0.3275)) +
  labs(title = "Confidence interval for entropy",
       x = expression(italic(H)),
       y = expression(italic(C))) +
  theme_minimal()




