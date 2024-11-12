# Required libraries
library(rprojroot)
library(readr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_set(theme_clean()+theme(legend.position = "top"))
library(StatOrdPattHxC)

summary(Normal_baseline_data)
attach(Normal_baseline_data)
t <- unlist(Normal_baseline_data)
data <- Normal_baseline_data
aggregate(data, FE_time ~ Motor_load, FUN = sigma2q, emb=6)



                           
                  