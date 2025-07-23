install.packages("magick")
install.packages("gridExtra")
install.packages("pdftools")
library(dplyr)
library(ggplot2)
library(readxl)
library(StatOrdPattHxC)
library(magick)
library(grid)
library(gridExtra)
library(pdftools)

# Read data
TS <- read_excel("Data/10 systems and their point in HC plane.xlsx", 
                                                      sheet = "TS data")
View(TS)
na.omit()
df <- read_excel("Data/10 systems and their point in HC plane.xlsx", 
                                                      sheet = "Entropy&Complexity_D6")
View(df)
names(df) <- trimws(names(df))

TS$Index <- 1:nrow(TS)

# Plot Normal_DE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = Normal_DE)) +
  geom_line(color = "blue") +
  labs(title = "Normal_DE Time Series", x = "Index", y = "Normal_DE") +
  coord_cartesian(xlim = c(0, 250000), ylim = c(-0.2, 0.3))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot Normal_FE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = Normal_FE)) +
  geom_line(color = "red") +
  labs(title = "Normal_FE Time Series", x = "Index", y = "Normal_FE") +
  coord_cartesian(xlim = c(0, 250000), ylim = c(-0.2, 0.3))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot 48K_DE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`48K_DE`)) +
  geom_line(color = "green") +
  labs(title = "48K_DE Time Series", x = "Index", y = "48K_DE") +
  #coord_cartesian(xlim = c(0, 250000), ylim = c(-0.2, 0.3))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot 48K_FE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`48K_FE`)) +
  geom_line(color = "purple") +
  labs(title = "48K_FE Time Series", x = "Index", y = "48K_FE") +
  #coord_cartesian(xlim = c(0, 250000), ylim = c(-0.2, 0.3))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot 12KDrive_DE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`12kDrive_DE`)) +
  geom_line(color = "black") +
  labs(title = "12kDrive_DE Time Series", x = "Index", y = "12kDrive_DE") +
  coord_cartesian(xlim = c(0, 700000))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")


# Plot 12KDrive_FE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`12kDrive_FE`)) +
  geom_line(color = "orange") +
  labs(title = "12kDrive_FE Time Series", x = "Index", y = "12kDrive_FE") +
  coord_cartesian(xlim = c(0, 700000))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot 12KDrive_BA with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`12kDrive_BA`)) +
  geom_line(color = "gray") +
  labs(title = "12kDrive_BA Time Series", x = "Index", y = "12kDrive_BA") +
  coord_cartesian(xlim = c(0, 700000))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot 12KFan_DE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`12kFan_DE`)) +
  geom_line(color = "pink") +
  labs(title = "12kFan_DE Time Series", x = "Index", y = "12kFan_DE") +
  coord_cartesian(xlim = c(0, 700000))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")

# Plot 12KFan_FE with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`12kFan_DE`)) +
  geom_line(color = "brown") +
  labs(title = "12kFan_FE Time Series", x = "Index", y = "12kFan_FE") +
  coord_cartesian(xlim = c(0, 700000))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")



# Plot 12KFan_BA with Index as x-axis
ggplot(data = TS, aes(x = Index, y = TS$`12kFan_DE`)) +
  geom_line(color = "skyblue") +
  labs(title = "12kFan_BA Time Series", x = "Index", y = "12kFan_BA") +
  coord_cartesian(xlim = c(0, 700000))+
  theme_minimal()+
  scale_color_brewer(palette = "Set1")


D = 6
data("LinfLsup")
# Confidence interval plot
df <- df %>%
  mutate(xmin = Entropy - df$SemiLength_H,
         xmax = Entropy + df$SemiLength_H)

#Define your custom colors
custom_colors <- c(
  "Normal_DE"     = "blue",
  "Normal_FE"     = "red",
  "48K_DE"        = "green",
  "48K_FE"        = "purple",
  "12kDrive_DE"   = "black",
  "12kDrive_FE"   = "orange",
  "12kDrive_BA"   = "gray",
  "12kFan_DE"     = "pink",
  "12kFan_FE"     = "brown",
  "12kFan_BA"     = "skyblue"
)

# Plot
ggplot(df, aes(x = Entropy, y = Complexity, color = TS)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)), 
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)), 
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point() +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax),
                 height = 0.005, alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(limits = c(0.25, 1)) +
  labs(title = "Confidence Interval for Entropy, D=6",
       x = expression(italic(H)),
       y = expression(italic(C))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(family = "serif", size = 10),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))



# Set the working path
plot_path <- "C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/Time Series and HC plane"

# Read PDF files and convert to magick images
img_main  <- image_read_pdf(file.path(plot_path, "Confidence Interval.pdf"), density = 300)
img_sub1  <- image_read_pdf(file.path(plot_path, "Normal_DE.pdf"), density = 300)
img_sub2  <- image_read_pdf(file.path(plot_path, "Normal_FE.pdf"), density = 300)
img_sub3  <- image_read_pdf(file.path(plot_path, "48K_DE.pdf"), density = 300)
img_sub4  <- image_read_pdf(file.path(plot_path, "48K_FE.pdf"), density = 300)
img_sub5  <- image_read_pdf(file.path(plot_path, "12kDrive_DE.pdf"), density = 300)
img_sub6  <- image_read_pdf(file.path(plot_path, "12kDrive_FE.pdf"), density = 300)
img_sub7  <- image_read_pdf(file.path(plot_path, "12kDrive_BA.pdf"), density = 300)
img_sub8  <- image_read_pdf(file.path(plot_path, "12kFan_DE.pdf"), density = 300)
img_sub9  <- image_read_pdf(file.path(plot_path, "12kFan_FE.pdf"), density = 300)
img_sub10 <- image_read_pdf(file.path(plot_path, "12kFan_BA.pdf"), density = 300)

# Convert to grobs
g_main <- rasterGrob(img_main)
g1 <- rasterGrob(img_sub1)
g2 <- rasterGrob(img_sub2)
g3 <- rasterGrob(img_sub3)
g4 <- rasterGrob(img_sub4)
g5 <- rasterGrob(img_sub5)
g6 <- rasterGrob(img_sub6)
g7 <- rasterGrob(img_sub7)
g8 <- rasterGrob(img_sub8)
g9 <- rasterGrob(img_sub9)
g10 <- rasterGrob(img_sub10)

# Arrange layout parts
top_row <- arrangeGrob(g1, g2, g3, ncol = 3)
bottom_row <- arrangeGrob(g4, g5, g6, ncol = 3)
right_column <- arrangeGrob(g7, g8, g9, g10, ncol = 1)

center_column <- arrangeGrob(
  top_row,
  g_main,
  bottom_row,
  ncol = 1,
  heights = c(1, 2, 1)
)

# Combine center and right into final plot
final_plot <- grid.arrange(
  center_column,
  right_column,
  ncol = 2,
  widths = c(3, 1)
)

# Save to PDF
ggsave("C:/Users/UserA1/Documents/GitHub/Ordinal_Patterns_R/Plots/combined_asymmetric_layout.pdf",
       plot = final_plot, width = 14, height = 12)











