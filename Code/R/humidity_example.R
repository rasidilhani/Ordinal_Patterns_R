library(StatOrdPattHxC)
library(dplyr)
library(tidyverse)
library(ggplot2)

#Data
Mean_of_Relative_Humidity<-c(77.3,81,82.4,81.7,83.6,85.6,84.4,83.1,78.8,79.6,78.2,78.8)

Month<-c("January","February","March","April","May","June","July","August","September","October","Nvermber","December")

humidity_frame<-data.frame(Month,Mean_of_Relative_Humidity, stringsAsFactors = F)

humidity_frame

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
  geom_line(color="blue", size=1) +
  geom_point(color="red", size=3) +
  labs(title="Mean Monthly Humidity in Wellington",
       x="Month",
       y="Mean Relative Humidity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), panel.grid = element_blank(), axis.line = element_line(color = "black"))  # Rotate x-axis labels

# Plot ordinal pattern distribution
ggplot(data=ordinal_patterns_table, aes(x=ordinal_patterns, y=Freq)) +
  geom_bar(stat="identity", fill="orange", alpha=0.8) +
  labs(title="Histogram of Proportions",
       x=expression(paste("Ordinal Patterns of type", pi,)),
       y="Frequency") +
  theme_minimal()
 
