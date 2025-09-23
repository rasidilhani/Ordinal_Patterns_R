library(ggplot2)
library(readxl)

# Load saved results from Excel
results <- readxl::read_xlsx("arma_results_Shannon.xlsx")

# Thresholds – adjust as needed
high_entropy_thresh <- 0.96   # Example: entropy < 0.98
low_complexity_thresh <- 0.1 # Example: complexity > 0.025

# Subset for error bars: high entropy, low complexity
results_with_ci <- subset(results, Entropy < high_entropy_thresh & Complexity > low_complexity_thresh)
results_no_ci   <- subset(results, !(Entropy < high_entropy_thresh & Complexity > low_complexity_thresh))

# Plot
ggplot() +
  # Points without error bars
  geom_point(data = results_no_ci, aes(x = Entropy, y = Complexity, color = ParameterSet), size = 2) +
  # Points with error bars
  geom_point(data = results_with_ci, aes(x = Entropy, y = Complexity, color = ParameterSet), size = 2) +
  geom_errorbarh(
    data = results_with_ci,
    aes(xmin = Entropy - semiLengthPD, xmax = Entropy + semiLengthPD, y = Complexity),
    height = 0.02
  ) +
  geom_errorbar(
    data = results_with_ci,
    aes(x = Entropy, ymin = Complexity - semiLengthPDC, ymax = Complexity + semiLengthPDC),
    width = 0.02
  ) +
  geom_line(
    data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
    aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE
  ) +
  geom_line(
    data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
    aes(x = H, y = C), color = "black", linetype = "dashed", inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(0.7, 1.0), ylim = c(0, 0.4)) +
  labs(
    title = paste("Entropy and Complexity for ARMA Time Series (D =", D, ", n =", n, ")"),
    x = expression(italic(H)),
    y = expression(italic(C)),
    color = "Parameter Set"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None")
