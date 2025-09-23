library(StatOrdPattHxC)
library(ggplot2)
library(ggExtra)
library(ggthemes)
library(GGally)

data("LinfLsup")

###############################################################################

set.seed(123456789, kind="Mersenne-Twister")

D <- 3
N <- c(500, 1000)
R <- 100

### Sketch of the experiment (one example)

### ARMA(2,2)

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    #ts_data <- arima.sim(model = list(ar1 = 0.1, ar2 = 0.8, ma1 = 0.1, ma2 = 0.8), n = n)
    ts_data <- arima.sim(model = list(ar = c(0.1, 0.8), ma = c(0.1, 0.8)), n = n)
    ProbARMA <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbARMA)
    Complexity <- StatComplexity(ProbARMA)
    
    Output <- rbind(Output, c(Entropy, Complexity, n))
  }
}

Output <- data.frame(Output)
Output
names(Output) <- c("Entropy", "Complexity", "n")

Output$n <- as.factor(Output$n)
Output.molten <- reshape2::melt(Output)

p <- ggplot() +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_point(data = subset(Output, n=="500" | n=="1000"), 
             aes(x=Entropy, y=Complexity, col=n)
  ) +
  xlim(.90, 1) +
  xlab(expression(italic(H))) +
  ylim(0, 0.15) +
  ylab(expression(italic(C))) +
  theme_pander() +
  theme(legend.position = "bottom")

ggMarginal(p, groupColour=TRUE)

save(Output, file = "arma22_M1.RData")
write_xlsx(Output, path = "arma22_M1.xlsx")
#ggsave("ARMA22_M1.pdf", p, 
#       width = 12, height = 10, dpi = 300)

