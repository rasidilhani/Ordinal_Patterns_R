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

### ARMA(2,2)_M1

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    ts_data <- arima.sim(model = list(ar1 = 0.1, ar2 = 0.8, ma1 = 0.1, ma2 = 0.8), n = n)
    ProbARMA <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbARMA)
    Complexity <- StatComplexity(ProbARMA)
    
    Output <- rbind(Output, c(Entropy, Complexity, n))
  }
}

Output <- data.frame(Output)
names(Output) <- c("Entropy", "Complexity", "n")

Output$n <- as.factor(Output$n)
Output.molten <- reshape2::melt(Output)

p1 <- ggplot() +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_point(data = subset(Output, n=="500" | n=="1000"), 
             aes(x=Entropy, y=Complexity, col=n)
  ) +
  xlim(.99, 1) +
  xlab(expression(italic(H))) +
  ylim(0, 0.01) +
  ylab(expression(italic(C))) +
  ggtitle("ARMA22_M1") +
  theme_pander() +
  theme(legend.position = "bottom")

ggMarginal(p1, groupColour=TRUE)

save(Output, file = "arma22_M1.RData")
write_xlsx(Output, path = "arma22_M1.xlsx")


#####################################################################################################################################
set.seed(123456789, kind="Mersenne-Twister")

D <- 3
N <- c(500, 1000)
R <- 100

### Sketch of the experiment (one example)

### ARMA(2,2)_M2

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    ts_data <- arima.sim(model = list(ar1 = -0.8, ar2 = 0.1, ma1 = -0.8, ma2 = 0.1), n = n)
    ProbARMA <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbARMA)
    Complexity <- StatComplexity(ProbARMA)
    
    Output <- rbind(Output, c(Entropy, Complexity, n))
  }
}

Output <- data.frame(Output)
names(Output) <- c("Entropy", "Complexity", "n")

Output$n <- as.factor(Output$n)
Output.molten <- reshape2::melt(Output)

p2 <- ggplot() +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_point(data = subset(Output, n=="500" | n=="1000"), 
             aes(x=Entropy, y=Complexity, col=n)
  ) +
  xlim(.99, 1) +
  xlab(expression(italic(H))) +
  ylim(0, 0.01) +
  ylab(expression(italic(C))) +
  ggtitle("ARMA22_M2") +
  theme_pander() +
  theme(legend.position = "bottom")

ggMarginal(p2, groupColour=TRUE)

save(Output, file = "arma22_M2.RData")
write_xlsx(Output, path = "arma22_M2.xlsx")



#####################################################################################################################################
set.seed(123456789, kind="Mersenne-Twister")

D <- 3
N <- c(500, 1000)
R <- 100

### Sketch of the experiment (one example)

### ARMA(2,2)_M3

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    ts_data <- arima.sim(model = list(ar1 = 0.1, ar2 = -0.8, ma1 = 0.1, ma2 = -0.8), n = n)
    ProbARMA <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbARMA)
    Complexity <- StatComplexity(ProbARMA)
    
    Output <- rbind(Output, c(Entropy, Complexity, n))
  }
}

Output <- data.frame(Output)
names(Output) <- c("Entropy", "Complexity", "n")

Output$n <- as.factor(Output$n)
Output.molten <- reshape2::melt(Output)

p3 <- ggplot() +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_point(data = subset(Output, n=="500" | n=="1000"), 
             aes(x=Entropy, y=Complexity, col=n)
  ) +
  xlim(.99, 1) +
  xlab(expression(italic(H))) +
  ylim(0, 0.01) +
  ylab(expression(italic(C))) +
  ggtitle("ARMA22_M3") +
  theme_pander() +
  theme(legend.position = "bottom")

ggMarginal(p3, groupColour=TRUE)

save(Output, file = "arma22_M3.RData")
write_xlsx(Output, path = "arma22_M3.xlsx")


###################################################################################################################################

set.seed(123456789, kind="Mersenne-Twister")

D <- 3
N <- c(500, 1000)
R <- 100

### Sketch of the experiment (one example)

### ARMA(2,2)_M4

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    ts_data <- arima.sim(model = list(ar1 = -0.8, ar2 = -0.1, ma1 = -0.8, ma2 = -0.1), n = n)
    ProbARMA <- OPprob(ts_data, emb = D)
    Entropy <- HShannon(ProbARMA)
    Complexity <- StatComplexity(ProbARMA)
    
    Output <- rbind(Output, c(Entropy, Complexity, n))
  }
}

Output <- data.frame(Output)
names(Output) <- c("Entropy", "Complexity", "n")

Output$n <- as.factor(Output$n)
Output.molten <- reshape2::melt(Output)

p4 <- ggplot() +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_point(data = subset(Output, n=="500" | n=="1000"), 
             aes(x=Entropy, y=Complexity, col=n)
  ) +
  xlim(.99, 1) +
  xlab(expression(italic(H))) +
  ylim(0, 0.01) +
  ylab(expression(italic(C))) +
  ggtitle("ARMA22_M4") +
  theme_pander() +
  theme(legend.position = "bottom")

ggMarginal(p4, groupColour=TRUE)

save(Output, file = "arma22_M4.RData")
write_xlsx(Output, path = "arma22_M4.xlsx")

###################################################################################################################################

library(gridExtra)

# Run the simulation and create each plot as shown in your code
# For each model:
# p1 <- ... # for ARMA22_M1, with ggtitle("ARMA22_M1")
# p2 <- ... # for ARMA22_M2, with ggtitle("ARMA22_M2")
# p3 <- ... # for ARMA22_M3, with ggtitle("ARMA22_M3")
# p4 <- ... # for ARMA22_M4, with ggtitle("ARMA22_M4")

# Combine the 4 plots in a 2x2 arrangement
grid.arrange(
  p1, p2,
  p3, p4,
  nrow = 2,
  ncol = 2,
  top = "Entropy–Complexity for ARMA(2,2) Models M1–M4"
)
ggsave("ARMA22_gridplot.png", arrangeGrob(p1, p2, p3, p4, nrow=2, ncol=2), width=11, height=8)

