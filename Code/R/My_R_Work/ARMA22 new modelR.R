library(ggplot2)
library(gridExtra)
library(grid)
library(StatOrdPattHxC)
library(ggthemes)
library(ggExtra)
library(writexl)
library(dplyr)

# --- Setup ---
data("LinfLsup")
set.seed(123456789, kind="Mersenne-Twister")
D <- 3
N <- 1000
R <- 100

# --- Define Models ---
arma22_models <- list(
  ARMA22_M1 = list(ar = c(0.3, 0.4), ma = c(0.3, 0.4), type = "M1"),
  ARMA22_M2 = list(ar = c(-0.5, 0.3), ma = c(-0.5, 0.3), type = "M2"),
  ARMA22_M3 = list(ar = c(0.4, -0.25), ma = c(0.4, -0.25), type = "M3"),
  ARMA22_M4 = list(ar = c(0.7, -0.25), ma = c(0.7, -0.25), type = "M4"),
  ARMA22_M5 = list(ar = c(-0.7, 0.25), ma = c(-0.7, 0.25), type = "M5")
)

# --- Generate Data Function ---
generate_arma22_data <- function(ar_coef, ma_coef, model_name, n, r_iterations = R) {
  Output <- NULL
  cat(paste("Generating", model_name, "...\n"))
  
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

# --- Run Analysis ---
results_list <- list()
for(model_name in names(arma22_models)){
  model <- arma22_models[[model_name]]
  res <- generate_arma22_data(model$ar, model$ma, model$type, N, R)
  results_list[[model_name]] <- res
}

# --- Combine and Save ---
all_results <- bind_rows(results_list)
write_xlsx(all_results, "ARMA22_entropy_complexity.xlsx")

cat("\nAnalysis complete! Results saved to ARMA22_entropy_complexity.xlsx\n")
cat(paste("Total observations:", nrow(all_results), "\n"))

# --- Plotting ---
p <- ggplot(all_results, aes(x = Entropy, y = Complexity, color = Model)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"),
            aes(x = H, y = C), col = "lightgray", size = 1, inherit.aes = FALSE) +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"),
            aes(x = H, y = C), col = "lightgray", size = 1, inherit.aes = FALSE) +
  geom_point(size = 2) +
  facet_wrap(~Model, ncol = 2) +   
  scale_color_manual(values = c("M1" = "red",
                                "M2" = "blue",
                                "M3" = "green", 
                                "M4" = "purple",
                                "M5" = "orange")) +
  coord_cartesian(xlim = c(0.85, 1), ylim = c(0, 0.3)) +
  theme_minimal(base_size = 12, base_family = "serif") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold")) +
  labs(title = "Entropy–Complexity Plane for ARMA(2,2) Models",
       x = expression(italic(H)),
       y = expression(italic(C)),
       color = "Model") 
print(p)
ggsave("ARMA22_entropy_complexity_plot.pdf", width = 8, height = 6, dpi = 300)
