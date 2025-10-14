library(ggplot2)
library(readxl)
library(dplyr)
library(StatOrdPattHxC)

data("LinfLsup")


# Read Excel file
setwd("C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Code/R/My_R_Work")
df <- read_excel("Results.xlsx")

# Filter for n = 500
df_500 <- df %>% filter(n == 500)

p500 <- ggplot(df_500, aes(x = Entropy, y = Complexity, color = Model)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
            aes(x = H, y = C), col = "lightgray", size = 1, inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
            aes(x = H, y = C), col = "lightgray", size = 1, inherit.aes = FALSE) +
  geom_point(size = 2) +
  facet_wrap(~Type, ncol = 2) +   
  scale_color_manual(values = c("M1" = "red",
                                "M2" = "blue",
                                "M3" = "green",
                                "M4" = "purple")) +
  coord_cartesian(xlim = c(0.85, 1), ylim = c(0, 0.3)) +   # <- axis limits
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "bottom",
       #panel.grid = element_blank(),
        strip.text = element_text(face = "bold")) +
  labs(title = "Entropy Complexity plane (n = 500)",
       x = expression(italic(H)),
       y = expression(italic(C)),
       color = "Model")

# Print plots
p500
ggsave("New_model_group_plot_n500.pdf", width = 8, height = 6, dpi = 300)  


###############################################################################################
# Sample size 1000

# Filter for n = 1000
df_1000 <- df %>% filter(n == 1000)

p1000 <- ggplot(df_1000, aes(x = Entropy, y = Complexity, color = Model)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
            aes(x = H, y = C), col = "lightgray", size = 1, inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
            aes(x = H, y = C), col = "lightgray", size = 1, inherit.aes = FALSE) +
  geom_point(size = 2) +
  facet_wrap(~Type, ncol = 2) +   
  scale_color_manual(values = c("M1" = "red",
                                "M2" = "blue",
                                "M3" = "green",
                                "M4" = "purple")) +
  coord_cartesian(xlim = c(0.85, 1), ylim = c(0, 0.3)) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  labs(title = "Entropy–Complexity (n = 1000)",
       x = expression(italic(H)),
       y = expression(italic(C)),
       color = "Model")

# Print plots
p1000

ggsave("New_model_group_plot_n1000.pdf", width = 8, height = 6, dpi = 300)


###############################################################################################
#Plot##############
