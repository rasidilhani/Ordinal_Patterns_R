library(StatOrdPattHxC)
library(ggplot2)

D <- 3
n <- 1000

ts_data <- arima.sim(model = list(ar1 = 0.8, ar2 = 0.5, ma1 = 0.8, ma2 = 0.5), n = n)
ProbARMA <- OPprob(ts_data, emb = D)
Entropy <- HShannon(ProbARMA)
Complexity <- StatComplexity(ProbARMA)


#####################################################################################################
library(StatOrdPattHxC)
D <- 3
n <- 100

# Function for a single run
arma_entropy_complexity <- function() {
  ts_data <- arima.sim(model = list(ar1 = 0.8, ar2 = 0.5, ma1 = 0.8, ma2 = 0.5), n = n)
  ProbARMA <- OPprob(ts_data, emb = D)
  H <- HShannon(ProbARMA)
  C <- StatComplexity(ProbARMA)
  c(Entropy = H, Complexity = C)
}

set.seed(123)  # for reproducibility!
results_matrix <- replicate(300, arma_entropy_complexity())

# Convert to a data frame for easy use (transpose if needed)
results_df <- data.frame(t(results_matrix))
head(results_df)

=======
library(ggthemes)
library(GGally)

set.seed(123456789, kind="Mersenne-Twister")

D <- 3
N <- c(100, 300, 500, 1000)
R <- 100

data("LinfLsup")

### ARMA(2,2)

Output <- NULL

for(n in N){
  for(r in 1:R){
  
    ts_data <- arima.sim(model = list(ar1 = 0.8, ar2 = 0.5, ma1 = 0.8, ma2 = 0.5), n = n)
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

ggplot(Output, aes(x=Entropy, y=Complexity, col=n)) +
  geom_point() 

### ARMA(1,1)

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    ts_data <- arima.sim(model = list(ar1 = 0.8, ma1 = 0.8), n = n)
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

ggplot(subset(Output, n=="500" | n=="1000"), aes(x=Entropy, y=Complexity, col=n)) +
  geom_point() +
  theme_pander()  

### AR(2)

Output <- NULL

for(n in N){
  for(r in 1:R){
    
    ts_data <- arima.sim(model = list(ar1 = -0.8, ar2=-.15), n = n)
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

ggplot() +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == "3"), 
            aes(x = H, y = C), col="lightgray") +
  geom_jitter(data = subset(Output, n=="500" | n=="1000"), 
              aes(x=Entropy, y=Complexity, col=n),
              width=0, height=.0005
                ) +
  xlim(.99, 1) +
  ylim(0, 0.01) +
  theme_pander()  

ggMarginal(data = subset(Output, n=="500" | n=="1000"),
  ggplot() +
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
    theme_pander() +
    theme(legend.position = "bottom"),
  groupColour=TRUE
)

