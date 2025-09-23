# --- Required Packages ---
library(ggplot2)
library(gridExtra)
library(grid)
library(StatOrdPattHxC)
library(ggthemes)
library(ggExtra)
library(writexl)
library(dplyr)

# --- Data and Parameters ---
data("LinfLsup")
set.seed(123456789)
D <- 3
N <- c(500, 1000)
R <- 100

arma11_models <- list(
  ARMA11_M1 = list(ar = c(0.8), ma = c(0.8), type = "ARMA11_M1"),
  ARMA11_M2 = list(ar = c(0.1), ma = c(0.1), type = "ARMA11_M2"),
  ARMA11_M3 = list(ar = c(-0.8), ma = c(-0.8), type = "ARMA11_M3"),
  ARMA11_M4 = list(ar = c(-0.1), ma = c(-0.1), type = "ARMA11_M4")
)

# --- Statistic Calculation Function ---
generate_arma11_data <- function(ar_coef, ma_coef, model_name, n, r_iterations = R) {
  Output <- NULL
  for(r in 1:r_iterations){
    ts_data <- arima.sim(model = list(ar = ar_coef, ma = ma_coef), n = n)
    ProbTS <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbTS)
    Complexity <- StatComplexity(ProbTS)
    Output <- rbind(Output, c(Entropy, Complexity, n, model_name))
  }
  Output <- data.frame(Output, stringsAsFactors = FALSE)
  names(Output) <- c("Entropy", "Complexity", "n", "Model")
  Output$Entropy <- as.numeric(Output$Entropy)
  Output$Complexity <- as.numeric(Output$Complexity)
  Output$n <- as.factor(Output$n)
  return(Output)
}

# --- Dynamic Axis Helper (y starts at 0, x ends at 1) ---
get_axis_range <- function(x, pad_lo = 0.03, pad_hi = 0.03, x_axis = FALSE) {
  rng <- range(x, na.rm = TRUE)
  padding_lo <- pad_lo * diff(rng)
  padding_hi <- pad_hi * diff(rng)
  if (padding_lo == 0) padding_lo <- pad_lo * abs(rng[1])
  if (padding_hi == 0) padding_hi <- pad_hi * abs(rng[2])
  if (x_axis) {
    return(c(rng[1] - padding_lo, 1))
  } else {
    return(c(0, rng[2] + padding_hi))
  }
}

# --- Plot Creation Function (with nice padded axes, no legend per subplot) ---
create_arma11_plot <- function(output_data, title) {
  x_rng <- get_axis_range(output_data$Entropy, pad_lo = 0.02, pad_hi = 0.02, x_axis = TRUE)
  y_rng <- get_axis_range(output_data$Complexity, pad_lo = 0.05, pad_hi = 0.05, x_axis = FALSE)
  p <- ggplot() +
    geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
              aes(x = H, y = C), col = "lightgray", size = 1) +
    geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
              aes(x = H, y = C), col = "lightgray", size = 1) +
    geom_point(data = output_data,
               aes(x = Entropy, y = Complexity, col = n), alpha = 0.7, size = 1.5) +
    xlab(expression(italic(H))) +
    ylab(expression(italic(C))) +
    ggtitle(title) +
    theme_pander() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    guides(color = guide_legend(title = "Sample Size")) +
    xlim(x_rng) +
    ylim(y_rng)
  p_marginal <- ggMarginal(p, groupColour = TRUE, groupFill = FALSE)
  return(p_marginal)
}

# --- Data Generation and Plotting ---
all_arma11_data <- data.frame()
individual_plots <- list()
for(i in 1:length(arma11_models)) {
  model_name <- names(arma11_models)[i]
  ar_coef <- arma11_models[[i]]$ar
  ma_coef <- arma11_models[[i]]$ma
  model_data <- data.frame()
  for(n in N) {
    temp_data <- generate_arma11_data(ar_coef, ma_coef, model_name, n)
    model_data <- rbind(model_data, temp_data)
  }
  individual_plots[[i]] <- create_arma11_plot(
    model_data,
    paste(model_name, "(AR =", ar_coef, ", MA =", ma_coef, ")")
  )
  all_arma11_data <- rbind(all_arma11_data, model_data)
}

# --- Shared Legend Extraction Helper ---
get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

# --- Create Dummy Plot With Legend (grouped data, padded limits) ---
get_plot_bounds <- function(df) {
  x_rng <- get_axis_range(df$Entropy, pad_lo = 0.02, pad_hi = 0.02, x_axis = TRUE)
  y_rng <- get_axis_range(df$Complexity, pad_lo = 0.05, pad_hi = 0.05, x_axis = FALSE)
  list(x = x_rng, y = y_rng)
}
dummy_bounds <- get_plot_bounds(all_arma11_data)

dummy_plot <- ggplot() +
  geom_point(data = all_arma11_data,
             aes(x = Entropy, y = Complexity, col = n), alpha = 0.7, size = 1.5) +
  theme_pander() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Sample Size")) +
  xlim(dummy_bounds$x) +
  ylim(dummy_bounds$y)
shared_legend <- get_only_legend(dummy_plot)

# --- Combine Subplots and Shared Legend ---
combined_plot <- arrangeGrob(
  grobs = individual_plots, ncol = 2, nrow = 2,
  top = "ARMA(1,1) Models - Ordinal Pattern Analysis"
)
grid.arrange(combined_plot, shared_legend, nrow = 2, heights = c(10, 1))

# --- Save Combined Plot to File ---
ggsave(
  "ARMA11_combined_analysis.pdf",
  arrangeGrob(combined_plot, shared_legend, nrow = 2, heights = c(10, 1)),
  width = 12, height = 10, dpi = 300
)