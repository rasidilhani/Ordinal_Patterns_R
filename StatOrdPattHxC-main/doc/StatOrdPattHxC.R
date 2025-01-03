## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi = 700, dev = "svg",
  fig.width=7, fig.height=4,
  warning = FALSE, message = FALSE
)

## ----setup, echo=FALSE--------------------------------------------------------
library(ggplot2)
library(ggthemes)
library(reshape2)
library(tibble)
library(dplyr)
library(prodlim)
library(StatOrdPattHxC)

## ----ThreeTimeSeries, out.width="90%", echo=TRUE, message=FALSE, warning=FALSE----
set.seed(1234567890, kind="Mersenne-Twister")

x1 <- mov.av(rnorm(5000), order = 31) # smoothed white noise with moving average
x2 <- mov.av(rnorm(5000), order = 31) # smoothed white noise with moving average
wn <- rnorm(4970) # white noise, same size as the vectors above

## ----Internal-ProcessingThreeTimeSeries, echo=FALSE, message=FALSE, warning=FALSE----
df <- data.frame(t=1:4970, wn, x1, x2)
df.molten <- reshape2::melt(df, measure.vars = 2:4)
df.molten$Type <- as.factor(rep(c("White Noise", "Moving Average 1", "Moving Average 2"), each=4970))

ggplot(subset(df.molten, t>1000 & t<1300), aes(x=t, y=value, col=Type, group=Type)) +
  geom_line() +
  ggtitle("300 observations of three time series") +
  xlab("Time") +
  ylab("Observations") +
  theme_tufte() +
  theme(legend.position = "top")

## ----OPsHists, out.width="90%"------------------------------------------------

OPprobs.wn <- OPprob(wn, emb=4)
OPprobs.x1 <- OPprob(x1, emb=4)
OPprobs.x2 <- OPprob(x2, emb=4)

## ----Internal-OPsHists, echo=FALSE, messages=FALSE, warning=FALSE-------------
OPprobabilities <- reshape2::melt(data.frame(index=1:factorial(4), 
                                             OPprobs.wn, OPprobs.x1, OPprobs.x2),
                                  measure.vars = 2:4)
## Reordering OPprobabilities$variable
OPprobabilities$variable <- factor(OPprobabilities$variable,
  levels = c("OPprobs.x1", "OPprobs.x2", "OPprobs.wn"))

ggplot(OPprobabilities, aes(x=index, y=value)) +
  geom_col(aes(fill=variable), col="black") +
  xlab("Pattern") +
  ylab("Proportion of patterns") +
  facet_wrap(~variable) +
  theme_tufte() +
  theme(legend.position = "none")

## ----Entropies----------------------------------------------------------------
# White noise
HShannon(OPprobs.wn)
HTsallis(OPprobs.wn, beta=2)
HRenyi(OPprobs.wn, beta=1.8)
HFisher(OPprobs.wn)

# Moving averages
HShannon(OPprobs.x1)
HTsallis(OPprobs.x1, beta=2)
HRenyi(OPprobs.x1, beta=1.8)
HFisher(OPprobs.x1)

## ----AsymptoticVariances, echo=TRUE-------------------------------------------
# Matrix that stores the variances
Variances <- matrix(nrow=3, ncol=4)

Variances[1,1] <- sigma2q(wn[1:1000], emb = 4, ent = "S")
Variances[2,1] <- sigma2q(x1[1:1000], emb = 4, ent = "S")
Variances[3,1] <- sigma2q(x2[1:1000], emb = 4, ent = "S")

Variances[1,2] <- sigma2q(wn[1:1000], emb = 4, ent = "T", beta = 1.5)
Variances[2,2] <- sigma2q(x1[1:1000], emb = 4, ent = "T", beta = 1.5)
Variances[3,2] <- sigma2q(x2[1:1000], emb = 4, ent = "T", beta = 1.5)

Variances[1,3] <- sigma2q(wn[1:1000], emb = 4, ent = "R", beta = 1.5)
Variances[2,3] <- sigma2q(x1[1:1000], emb = 4, ent = "R", beta = 1.5)
Variances[3,3] <- sigma2q(x2[1:1000], emb = 4, ent = "R", beta = 1.5)

Variances[1,4] <- sigma2q(wn[1:1000], emb = 4, ent = "F")
Variances[2,4] <- sigma2q(x1[1:1000], emb = 4, ent = "F")
Variances[3,4] <- sigma2q(x2[1:1000], emb = 4, ent = "F")

Variances[Variances<0] <- 0
rownames(Variances) <- c("White Noise", "Moving Average", "Moving Average")

knitr::kable(Variances, digits = 4, format.args = list(scientific = TRUE),
             col.names = c("Shannon", "Tsallis $\\beta=1.5$", "Rényi $\\beta=1.5$", "Fisher"),
             escape = TRUE,
             caption = "Asymptotic variances")

## ----ApplyTestSamePermutationEntropy------------------------------------------

entropicTest(wn, x1, emb=4, ent="S") # Time series with different structure
entropicTest(x1, x2, emb=4, ent="S") # Time series with same structure

## ----Boundaries, out.width="100%"---------------------------------------------
data("LinfLsup")

ggplot(subset(LinfLsup, Side=="Lower"), 
       aes(x=H, y=C, col=Dimension, group=Dimension)) +
  geom_line() +
  geom_line(data=subset(LinfLsup, Side=="Upper"), 
            aes(x=H, y=C, col=Dimension, group=Dimension)) +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C))) +
  theme_tufte()

## ----PointsWithConfidenceIntervals, out.width="100%", echo=TRUE, message=FALSE----
data("LinfLsup")

D <- 4 # Embedding dimension

wnsub <- wn
x1sub <- x1
x2sub <- x2

ShannonEntropies <- c(
  HShannon(OPprob(wnsub, emb=D)),
  HShannon(OPprob(x1sub, emb=D)),
  HShannon(OPprob(x2sub, emb=D))
)

StatisticalComplexities <- c(
  StatComplexity(OPprob(wnsub, emb=D)),
  StatComplexity(OPprob(x1sub, emb=D)),
  StatComplexity(OPprob(x2sub, emb=D))
  )

alpha <- 0.05
StandardDeviations <- sqrt(Variances[,1])
SemiLength <- StandardDeviations/sqrt(length(wnsub)-D)*qnorm(1-alpha/2) 
  # The three time series have the same length, but they could be different

HCPoints <- data.frame(H=ShannonEntropies,
                       C=StatisticalComplexities,
                       STD=StandardDeviations,
                       SemiLength=SemiLength,
                       Series=as.factor(c("White Noise", "Moving Average 1", "Moving Average 2")))

ggplot(subset(LinfLsup, Side=="Lower" & Dimension==as.character(D)), 
       aes(x=H, y=C)) +
  geom_line() +
  geom_line(data=subset(LinfLsup, Side=="Upper" & Dimension==as.character(D)), 
            aes(x=H, y=C)) +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C))) +
  geom_point(data=HCPoints, aes(x=H, y=C, col=Series)) +
  geom_errorbarh(data=HCPoints, aes(xmin=H-SemiLength, xmax=H+SemiLength, group=Series, col=Series)) +
  coord_cartesian(xlim=c(0.925, 1), ylim=c(0, 0.13)) +
  theme_tufte()

## ----StatComplexityMeanVariance, out.width="100%", echo=TRUE, message=FALSE----

n <- 4970-4

mu.wn <- meanC(OPprobs.wn, n)
sd.wn <- sqrt(varC(OPprobs.wn, n))

mu.x1 <- meanC(OPprobs.x1, n)
sd.x1 <- sqrt(varC(OPprobs.x1, n))

mu.x2 <- meanC(OPprobs.x2, n)
sd.x2 <- sqrt(varC(OPprobs.x2, n))

(Complexities <- data.frame(mu.wn, mu.x1, mu.x2))
(StandardDeviantionsComplexities <- data.frame(sd.wn, sd.x1, sd.x2))
alpha <- 0.05
SemiLength <- StandardDeviantionsComplexities/sqrt(n)*qnorm(1-alpha/2) 

