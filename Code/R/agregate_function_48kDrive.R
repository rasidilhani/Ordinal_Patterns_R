library(readr)
library(StatOrdPattHxC)

data <- read_csv("../../Data/csv/48k_drive_end_bearing_fault_data.csv")


Drive48_S_DE_6 <-aggregate(data,
                         DE_time~ Motor_load,
                         FUN = sigma2q, emb=6) 
save(file="../../Data/Results/Drive48_S_DE_6.RData", Drive_S_DE_6)

Drive48_S_FE_6 <-aggregate(data,
                         FE_time~ Motor_load,
                         FUN = sigma2q, emb=6) 
save(file="../../Data/Results/Drive48_S_FE_6.RData", Drive_S_DE_6)


                                                                                         