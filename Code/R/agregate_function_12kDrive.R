data <- X12k_drive_end_bearing_fault_data

Drive_S_FE_5 <-aggregate(data,
                    FE_time~ Motor_load,
                    FUN = sigma2q, emb=5) 

Drive_S_FE_6 <-aggregate(data,
                   FE_time~ Motor_load,
                   FUN = sigma2q, emb=6) 

Drive_S_BE_5 <-aggregate(data,
                   BE_time~ Motor_load,
                   FUN = sigma2q, emb=5) 

Drive_S_BE_6 <-aggregate(data,
                   BE_time~ Motor_load,
                   FUN = sigma2q, emb=6) 

Drive_S_BA_5 <-aggregate(data,
                   BA_time~ Motor_load,
                   FUN = sigma2q, emb=5) 

Drive_S_BA_6 <-aggregate(data,
                   BA_time~ Motor_load,
                   FUN = sigma2q, emb=6)
