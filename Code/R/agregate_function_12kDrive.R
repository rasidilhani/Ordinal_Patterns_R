data <- X12k_drive_end_bearing_fault_data

Drive_S_DE_5 <-aggregate(data,
                   DE_time~ Motor_load,
                   FUN = sigma2q, emb=5) 

Drive_S_DE_6 <-aggregate(data,
                   DE_time~ Motor_load,
                   FUN = sigma2q, emb=6) 

Drive_S_BA_5 <-aggregate(data,
                   BA_time~ Motor_load,
                   FUN = sigma2q, emb=5) 

Drive_S_BA_6 <-aggregate(data,
                   BA_time~ Motor_load,
                   FUN = sigma2q, emb=6)
