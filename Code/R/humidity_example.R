library(StatOrdPattHxC)
library(dplyr)
library(tidyverse)
library(ggplot2)

#Data
Mean_of_Relative_Humidity<-c(77.3,81,82.4,81.7,83.6,85.6,84.4,83.1,78.8,79.6,78.2,78.8)

Month<-c("January","February","March","April","May","June","July","August","September","October","November","December")

humidity_frame<-data.frame(Month,Mean_of_Relative_Humidity, stringsAsFactors = T)

humidity_frame$Month <- factor(humidity_frame$Month,
                               levels = c("January", "February", "March", "April", "May", "June",
                                          "July", "August", "September", "October", "November", "December"))

#Compute Ordinal Patterns
ordinal_patterns<-OPseq(humidity_frame$Mean_of_Relative_Humidity,emb=3,lag = 1)
ordinal_patterns_table <- as.data.frame(table(ordinal_patterns))  # Convert to dataframe for plotting

#Compute probabilities
op_probability<-OPprob(humidity_frame$Mean_of_Relative_Humidity,emb = 3)

# Compute entropy and complexity
H_Shannon<-HShannon(op_probability)
C_Complexity<-StatComplexity(op_probability)

# Print results
print(paste("Shannon Entropy:", H_Shannon))
print(paste("Statistical Complexity:", C_Complexity))

# Line Graph for Original Humidity Data 
ggplot(data=humidity_frame, aes(x=Month, y=Mean_of_Relative_Humidity, group=1)) +
  geom_line(color="blue", size=1, alpha=0.8) +
  geom_point(color="red", size=3, alpha=0.8) +
  labs(x="Month", y="Mean Relative Humidity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(family = "serif", size=16),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))


ggsave(file="../../Text/Proposal/humidity graph.pdf", width = 16, height=10, units = "cm") 

# Plot ordinal pattern distribution
ggplot(data=ordinal_patterns_table, aes(x=ordinal_patterns, y=Freq)) +
  geom_bar(stat="identity", fill="orange", alpha=0.8) +
  labs(x="Ordinal Pattern Type",
       y="Frequency") +
  scale_x_discrete(
    breaks=1:6,
    labels=c(expression(pi^1), 
             expression(pi^2), 
             expression(pi^3), 
             expression(pi^4), 
             expression(pi^5), 
             expression(pi^6))
  ) +
  theme_minimal() +
  theme(text = element_text(family = "serif", size=16), panel.grid = element_blank(), axis.line = element_line(color = "black"))

ggsave(file="../../Text/Proposal/frequency histogram.pdf", width = 16, height=10, units = "cm") 
