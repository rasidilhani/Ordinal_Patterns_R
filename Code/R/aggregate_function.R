summary(Normal_baseline_data)
attach(Normal_baseline_data)
t <- unlist(Normal_baseline_data)

# df$Motor_load<-as.factor(df$Motor_load)
# df$RPM <- as.factor(df$RPM)

aggregate(data=Normal_baseline_data, 
          DE_time ~ Motor_load,
          FUN = mean)

unique(Normal_baseline_data$Motor_load)


summary(X12kFan_end_bearing_fault_data)
attach(X12kFan_end_bearing_fault_data)
a <-unlist(X12kFan_end_bearing_fault_data)
aggregate(data=X12kFan_end_bearing_fault_data,
          BA_time~ Motor_load,
          FUN = sigma2q, emb=6) 