library(dplyr)
library(tidyverse)
library(ggplot2)

H<-0.8008
Q<-0.4052
C<-0.3245
sigmaH<-0.07694
sigmaC<-1.32*10^-6

HSemiLength <- sigmaH/sqrt(10000-3)*qnorm(1-0.05/2) 
VSemiLength <- sigmaC/sqrt(10000-3)*qnorm(1-0.05/2)

df <- data.frame(H,C,sigmaH,sigmaC, HSemiLength, VSemiLength)

#D = 3
#data("LinfLsup")

ggplot(data=df, mapping = aes(H,C)) + 
  #geom_line(data = subset(LinfLsup, Side == "Lower" & Dimension == as.character(D)), 
            #aes(x = H, y = C, color = "Lower Boundary"), linetype = "dashed") +
  #geom_line(data = subset(LinfLsup, Side == "Upper" & Dimension == as.character(D)), 
           # aes(x = H, y = C, color = "Upper Boundary"), linetype = "dashed") +
  geom_point(color="red", size=1) +
  geom_errorbarh(data=df, mapping = aes(xmin=H-HSemiLength,
                                        xmax=H+HSemiLength, height=0.001)
                 ) +
  geom_errorbar(data=df, mapping=aes(ymin=C-VSemiLength,
                                     ymax=C+VSemiLength, width=0.1)
                ) +
  coord_cartesian(xlim = c(0.7295, 0.845), ylim = c(0.320, 0.3275)) +
  labs(title = "Confidence interval for entropy and complexity") + 
  xlab(expression(italic(H))) +
  ylab(expression(italic(C)))
