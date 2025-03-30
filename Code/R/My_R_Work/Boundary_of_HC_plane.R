library(StatOrdPattHxC)

data("LinfLsup")

# Create plot
ggplot(data = combined, aes(x = H, y = C, color = Boundaries)) + 
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)), 
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)), 
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  labs(title = "Entropy Complexity Plane for Embedding Dimension 3",
       x = expression(italic(H)),
       y = expression(italic(C))) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"))
