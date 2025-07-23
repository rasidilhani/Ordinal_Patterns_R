library(StatOrdPattHxC)
library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)

#################################################################################
# Asymptotic Distribution Incorporating dependence patterns

data10 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                                    sheet = "series10")
View(data10)

#Compute Ordinal Patterns
OP_data10 <- OPseq(data10$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data10<-OPprob(data10$DE_time, emb = 3)

# Compute entropy and complexity
H_data10<-HShannon(oprob_data10)
C_data10<-StatComplexity(oprob_data10)

# Compute asymptotic variance
variance_data10 <- sigma2q(data10$DE_time, emb = 3)

#####################################################################################
# number of data in series 100

data100 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                     sheet = "series100")
View(data100)

#Compute Ordinal Patterns
OP_data100 <- OPseq(data100$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data100<-OPprob(data100$DE_time, emb = 3)

# Compute entropy and complexity
H_data100<-HShannon(oprob_data100)
C_data100<-StatComplexity(oprob_data100)

# Compute asymptotic variance for entropy and complexity
variance_data100 <- sigma2q(data100$DE_time, emb = 3)
std_data100 <- sqrt(variance_data100)

semiLength100 <- std_data100/sqrt(100-3)*qnorm(1-0.05/2)

varC100 <- 0.3334
semilengthC100 <- sqrt(varC100)/sqrt(100-3)*qnorm(1-0.05/2)

df100 <- data.frame(H_data100,C_data100,variance_data100, semiLength100, varC100, semilengthC100)

D = 3
data("LinfLsup")

ggplot(df100, mapping = aes(H_data100,C_data100)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
  aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
   aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(df100, mapping = aes(xmin=H_data100-semiLength100,
                                        xmax=H_data100+semiLength100)
  ) +
  geom_errorbar(df100, mapping=aes(ymin=C_data100-semilengthC100,
                                     ymax=C_data100+semilengthC100, width=0.1)
  ) +
  coord_cartesian(xlim = c(0.45, 0.95), ylim = c(0.1, 0.3)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))



###################################################################################
# number of data in series 1000

data1000 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                      sheet = "series1000")
View(data1000)

#Compute Ordinal Patterns
OP_data1000 <- OPseq(data1000$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data1000<-OPprob(data1000$DE_time, emb = 3)

# Compute entropy and complexity
H_data1000<-HShannon(oprob_data1000)
C_data1000<-StatComplexity(oprob_data1000)

# Compute asymptotic variance
variance_data1000 <- sigma2q(data1000$DE_time, emb = 3)
std_data1000 <- sqrt(variance_data1000)

semiLength1000 <- std_data1000/sqrt(1000-3)*qnorm(1-0.05/2)

varC1000 <- 0.0415
semilengthC1000 <- sqrt(varC1000)/sqrt(1000-3)*qnorm(1-0.05/2)

df1000 <- data.frame(H_data1000,C_data1000,variance_data1000, semiLength1000, varC1000, semilengthC1000)

D = 3
data("LinfLsup")

ggplot(df1000, mapping = aes(H_data1000,C_data1000)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(df1000, mapping = aes(xmin=H_data1000-semiLength1000,
                                      xmax=H_data1000+semiLength1000)
  ) +
  geom_errorbar(df1000, mapping=aes(ymin=C_data1000-semilengthC1000,
                                   ymax=C_data1000+semilengthC1000, width=0.1)
  ) +
  coord_cartesian(xlim = c(0.45, 0.95), ylim = c(0.1, 0.3)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))

###################################################################################
# number of data in series 2000

data2000 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                       sheet = "series2000")
View(data2000)

#Compute Ordinal Patterns
OP_data1000 <- OPseq(data2000$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data2000<-OPprob(data2000$DE_time, emb = 3)

# Compute entropy and complexity
H_data2000<-HShannon(oprob_data2000)
C_data2000<-StatComplexity(oprob_data2000)

# Compute asymptotic variance
variance_data2000 <- sigma2q(data2000$DE_time, emb = 3)
std_data2000 <- sqrt(variance_data2000)

semiLength2000 <- std_data2000/sqrt(2000-3)*qnorm(1-0.05/2)

varC2000 <- 0.034
semilengthC2000 <- sqrt(varC2000)/sqrt(2000-3)*qnorm(1-0.05/2)

df2000 <- data.frame(H_data2000,C_data2000,variance_data2000, semiLength2000, varC2000, semilengthC2000)

ggplot(df2000, mapping = aes(H_data2000,C_data2000)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(df2000, mapping = aes(xmin=H_data2000-semiLength2000,
                                       xmax=H_data2000+semiLength2000)
  ) +
  geom_errorbar(df2000, mapping=aes(ymin=C_data2000-semilengthC2000,
                                     ymax=C_data2000+semilengthC2000, width=0.1)
  ) +
  coord_cartesian(xlim = c(0.45, 0.95), ylim = c(0.1, 0.3)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))




#####################################################################################################
#All in one

# Create a unified data frame
df100 <- data.frame(
  H = H_data100,
  C = C_data100,
  semiLengthH = semiLength100,
  semiLengthC = semilengthC100,
  Group = "n=100"
)
df1000 <- data.frame(
  H = H_data1000,
  C = C_data1000,
  semiLengthH = semiLength1000,
  semiLengthC = semilengthC1000,
  Group = "n=1000"
)
df2000 <- data.frame(
  H = H_data2000,
  C = C_data2000,
  semiLengthH = semiLength2000,
  semiLengthC = semilengthC2000,
  Group = "n=2000"
)
results_all <- rbind(df100, df1000, df2000)


ggplot(results_all, aes(x = H, y = C, color = Group)) +
  # Lower and Upper boundaries
  geom_line(
    data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
    aes(x = H, y = C, color = "Lower Boundary"),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_line(
    data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
    aes(x = H, y = C, color = "Upper Boundary"),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  # Sample size points
  geom_point(size = 3) +
  # Horizontal error bars for H
  geom_errorbarh(aes(xmin = H - semiLengthH, xmax = H + semiLengthH), height = 0.004) +
  # Vertical error bars for C
  geom_errorbar(aes(ymin = C - semiLengthC, ymax = C + semiLengthC), width = 0.008) +
  coord_cartesian(xlim = c(0.25, 0.95), ylim = c(0.02, 0.35)) +
  labs(
    title = "Confidence Interval for Entropy and Complexity (Pattern Dependence)",
    x = expression(italic(H)),
    y = expression(italic(C)),
    color = "Sample Size/Boundary"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("n=100" = "red", "n=1000" = "blue", "n=2000" = "green",
               "Lower Boundary" = "grey30", "Upper Boundary" = "grey60"),
    breaks = c("n=100", "n=1000", "n=2000", "Lower Boundary", "Upper Boundary")
  )



#############################################################################################
#############################################################################################
#############################################################################################
# Asymptotic ditribution under multinomial law

data10 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                     sheet = "series10")
View(data10)
n=8

#Compute Ordinal Patterns
OP_data10 <- OPseq(data10$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data10<-OPprob(data10$DE_time, emb = 3)

# Compute entropy and complexity
H_data10<-HShannon(oprob_data10)
C_data10<-StatComplexity(oprob_data10)

# Compute asymptotic variance for H
varianceH_data10 <- asymptoticVarHShannonMultinomial(oprob_data10, n)
stdH10 <- sqrt(varianceH_data10)

# Compute asymptotic variance for C
meanC_data10 <- meanC(oprob_data10, n)
varC_data10 <- varC(oprob_data10, n)

stdC10 <- sqrt(varC_data10)

HSemiLength10 <- stdH10/sqrt(10-3)*qnorm(1-0.05/2)
CSemiLength10 <- stdC10/sqrt(10-3)*qnorm(1-0.05/2)

dfM10 <- data.frame(H_data10,meanC_data10,HSemiLength10, CSemiLength10)

ggplot(dfM10, mapping = aes(H_data10,C_data10)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(dfM10, mapping = aes(xmin=H_data10-HSemiLength10,
                                       xmax=H_data10+HSemiLength10)
  ) +
  geom_errorbar(data=dfM10, mapping=aes(ymin=meanC_data10-CSemiLength10,
                                     ymax=meanC_data10+CSemiLength10, width=0.1)
   ) +
   #coord_cartesian(xlim = c(0.7295, 0.845), ylim = c(0.320, 0.3275)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))


#####################################################################################
# number of data in series 100

data100 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                      sheet = "series100")
View(data100)
n=98

#Compute Ordinal Patterns
OP_data100 <- OPseq(data100$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data100<-OPprob(data100$DE_time, emb = 3)

# Compute entropy and complexity
H_data100<-HShannon(oprob_data100)
C_data100<-StatComplexity(oprob_data100)

# Compute asymptotic variance for H
varianceH_data100 <- asymptoticVarHShannonMultinomial(oprob_data100, n)
stdH100 <- sqrt(varianceH_data100)

# Compute asymptotic variance for C
meanC_data100 <- meanC(oprob_data100, n)
varC_data100 <- varC(oprob_data100, n)

stdC100 <- sqrt(varC_data100)

HSemiLength100 <- stdH100/sqrt(100-3)*qnorm(1-0.05/2)
CSemiLength100 <- stdC100/sqrt(100-3)*qnorm(1-0.05/2)

dfM100 <- data.frame(H_data100,meanC_data100,HSemiLength100, CSemiLength100)

ggplot(dfM100, mapping = aes(H_data100,C_data100)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(dfM100, mapping = aes(xmin=H_data100-HSemiLength100,
                                      xmax=H_data100+HSemiLength100)
  ) +
  geom_errorbar(data=dfM100, mapping=aes(ymin=meanC_data100-CSemiLength100,
                                        ymax=meanC_data100+CSemiLength100, width=0.1)
  ) +
  #coord_cartesian(xlim = c(0.7295, 0.845), ylim = c(0.320, 0.3275)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))



###################################################################################
# number of data in series 1000

data1000 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                       sheet = "series1000")
View(data1000)
n=998

#Compute Ordinal Patterns
OP_data1000 <- OPseq(data1000$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data1000<-OPprob(data1000$DE_time, emb = 3)

# Compute entropy and complexity
H_data1000<-HShannon(oprob_data1000)
C_data1000<-StatComplexity(oprob_data1000)

# Compute asymptotic variance for H
varianceH_data1000 <- asymptoticVarHShannonMultinomial(oprob_data1000, n)
stdH1000 <- sqrt(varianceH_data1000)

# Compute asymptotic variance for C
meanC_data1000 <- meanC(oprob_data1000, n)
varC_data1000 <- varC(oprob_data1000, n)

stdC1000 <- sqrt(varC_data1000)

HSemiLength1000 <- stdH1000/sqrt(1000-3)*qnorm(1-0.05/2)
CSemiLength1000 <- stdC1000/sqrt(1000-3)*qnorm(1-0.05/2)

dfM1000 <- data.frame(H_data1000,meanC_data1000,HSemiLength1000, CSemiLength1000)

ggplot(dfM1000, mapping = aes(H_data1000,C_data1000)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(dfM1000, mapping = aes(xmin=H_data1000-HSemiLength1000,
                                       xmax=H_data1000+HSemiLength1000)
  ) +
  geom_errorbar(data=dfM1000, mapping=aes(ymin=meanC_data1000-CSemiLength1000,
                                         ymax=meanC_data1000+CSemiLength1000, width=0.1)
  ) +
  #coord_cartesian(xlim = c(0.7295, 0.845), ylim = c(0.320, 0.3275)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))

###################################################################################
# number of data in series 2000

data2000 <- read_excel("~/GitHub/Ordinal_Patterns_R/Data/hand_calculation_data.xlsx",
                       sheet = "series2000")
View(data2000)
n=1998

#Compute Ordinal Patterns
OP_data1000 <- OPseq(data2000$DE_time, emb = 3,lag = 1)

#Compute probabilities
oprob_data2000<-OPprob(data2000$DE_time, emb = 3)

# Compute entropy and complexity
H_data2000<-HShannon(oprob_data2000)
C_data2000<-StatComplexity(oprob_data2000)

# Compute asymptotic variance for H
varianceH_data2000 <- asymptoticVarHShannonMultinomial(oprob_data2000, n)
stdH2000 <- sqrt(varianceH_data2000)

# Compute asymptotic variance for C
meanC_data2000 <- meanC(oprob_data2000, n)
varC_data2000 <- varC(oprob_data2000, n)

stdC2000 <- sqrt(varC_data2000)

HSemiLength2000 <- stdH2000/sqrt(2000-3)*qnorm(1-0.05/2)
CSemiLength2000 <- stdC2000/sqrt(2000-3)*qnorm(1-0.05/2)

dfM2000 <- data.frame(H_data2000,meanC_data2000,HSemiLength2000, CSemiLength2000)

ggplot(dfM2000, mapping = aes(H_data2000,C_data2000)) +
  geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
            aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red") +
  geom_errorbarh(dfM2000, mapping = aes(xmin=H_data2000-HSemiLength2000,
                                        xmax=H_data2000+HSemiLength2000)
  ) +
  geom_errorbar(data=dfM2000, mapping=aes(ymin=meanC_data2000-CSemiLength2000,
                                          ymax=meanC_data2000+CSemiLength2000, width=0.1)
  ) +
  #coord_cartesian(xlim = c(0.7295, 0.845), ylim = c(0.320, 0.3275)) +
  labs(title = "Confidence interval for entropy and complexity") +
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))


########################################################################################
#All in one

# Construct data frames
#dfM10    <- data.frame(H = H_data10,    C = meanC_data10,    HSemi = HSemiLength10,  CSemi = CSemiLength10,  Group = "n=10")
dfM100   <- data.frame(H = H_data100,   C = meanC_data100,   HSemi = HSemiLength100, CSemi = CSemiLength100, Group = "n=100")
dfM1000  <- data.frame(H = H_data1000,  C = meanC_data1000,  HSemi = HSemiLength1000,CSemi = CSemiLength1000,Group = "n=1000")
dfM2000  <- data.frame(H = H_data2000,  C = meanC_data2000,  HSemi = HSemiLength2000,CSemi = CSemiLength2000,Group = "n=2000")

# Combine all data into one data frame
#df_all <- rbind(dfM10, dfM100, dfM1000, dfM2000)
df_all <- rbind(dfM100, dfM1000, dfM2000)


D <- 3  # chosen embedding dimension

ggplot(df_all, aes(x = H, y = C, color = Group)) +
  # Add Lower and Upper Boundaries
  geom_line(
    data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)),
    aes(x = H, y = C, color = "Lower Boundary"),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_line(
    data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)),
    aes(x = H, y = C, color = "Upper Boundary"),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  # All sample size points
  geom_point(size = 3) +
  # Horizontal error bars (entropy)
  geom_errorbarh(
    aes(xmin = H - HSemi, xmax = H + HSemi),
    height = 0.002
  ) +
  # Vertical error bars (complexity)
  geom_errorbar(
    aes(ymin = C - CSemi, ymax = C + CSemi),
    width = 0.008
  ) +
  coord_cartesian(xlim = c(0.45, 0.95), ylim = c(0.1, 0.3)) +
  labs(
    title = "Confidence Interval for Entropy and Complexity (Multinomial Model )",
    x = expression(italic(H)),
    y = expression(italic(C)),
    color = "Sample Size/Boundary"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("n=100" = "red", "n=1000" = "blue", "n=2000" = "green",
               "Lower Boundary" = "grey30", "Upper Boundary" = "grey60"),
    breaks = c("n=10", "n=100", "n=1000", "n=2000", "Lower Boundary", "Upper Boundary")
  )







