library("StatOrdPattHxC")

data <- read_csv("~/Library/CloudStorage/OneDrive-VictoriaUniversityofWellington-STAFF/Documents/Alunos/Rasika Dilhani/Ordinal_Patterns_R/Data/csv/12k_drive_end_bearing_fault_data.csv")

Drive_S_DE_5 <-aggregate(data,
                   DE_time~ Motor_load,
                   FUN = sigma2q, emb=5) 
save(file="../Results/Drive_S_DE_5.RData", Drive_S_DE_5)

Drive_S_DE_6 <-aggregate(data,
                   DE_time~ Motor_load,
                   FUN = sigma2q, emb=6) 
save(file="../Results/Drive_S_DE_6.RData", Drive_S_DE_6)

Drive_S_BA_5 <-aggregate(data,
                   BA_time~ Motor_load,
                   FUN = sigma2q, emb=5) 

Drive_S_BA_6 <-aggregate(data,
                   BA_time~ Motor_load,
                   FUN = sigma2q, emb=6)

save(file="../Results/Alejandro3.RData", Drive_S_BA_6)
