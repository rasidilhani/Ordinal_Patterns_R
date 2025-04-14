library(StatOrdPattHxC)
library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)

normal <- read_excel("Data/case_study.xlsx",  sheet = "Normal")   
                     
fourty_eightk <- read_excel("Data/case_study.xlsx",  sheet = "48k")   

  
#Compute Ordinal Patterns
OP_Normal_DE <- OPseq(normal$DE_time, emb = 3,lag = 1)
OP_Normal_FE <- OPseq(normal$FE_time, emb = 3, lag = 1)

OP_48k_DE <- OPseq(fourty_eightk$DE_time, emb = 3, lag = 1)
OP_48k_FE <- OPseq(fourty_eightk$FE_time, emb = 3, lag = 1)

# Convert to dataframe for plotting
ordinal_patterns_table <- as.data.frame(table(OP_Normal_DE,OP_Normal_FE,OP_48k_DE,OP_48k_FE))  

#Compute probabilities
oprob_Normal_DE<-OPprob(normal$DE_time, emb = 3)
oprob_Normal_FE<-OPprob(normal$FE_time, emb = 3)
oprob_48k_DE<-OPprob(fourty_eightk$DE_time, emb = 3)
oprob_48k_FE<-OPprob(fourty_eightk$FE_time, emb = 3)

# Convert probability to dataframe for plottings 
probability_table <- as.data.frame(table(oprob_Normal_DE,oprob_Normal_FE,oprob_48k_DE,oprob_48k_FE )) 

# Compute entropy and complexity
H_Shannon_Normal_DE<-HShannon(oprob_Normal_DE)
C_Complexity_Normal_DE<-StatComplexity(oprob_Normal_DE)

H_Shannon_Normal_FE<-HShannon(oprob_Normal_FE)
C_Complexity_Normal_FE<-StatComplexity(oprob_Normal_FE)

H_Shannon_48k_DE<-HShannon(oprob_48k_DE)
C_Complexity_48k_DE<-StatComplexity(oprob_48k_DE)

H_Shannon_48k_FE<-HShannon(oprob_48k_FE)
C_Complexity_48k_FE<-StatComplexity(oprob_48k_FE)


# Organize results into a data frame
results_table <- data.frame(
  Dataset = c("Normal_DE", "Normal_FE", "48k_DE", "48k_FE"),
  Shannon_Entropy_H = c(H_Shannon_Normal_DE, H_Shannon_Normal_FE, H_Shannon_48k_DE, H_Shannon_48k_FE),
  Statistical_Complexity_C = c(C_Complexity_Normal_DE, C_Complexity_Normal_FE, C_Complexity_48k_DE, C_Complexity_48k_FE)
)

# Print the organized table
print(results_table) 

###############################################################################################
###############################################################################################

#Batch Results

# Define the process_batch function
process_batch <- function(data, batch_size = 100, emb = 3, lag = 1) {
  n <- length(data)  # Corrected line
  batches <- seq(1, n, by = batch_size)
  
  op_list <- list()
  oprob_list <- list()
  
  for (i in batches) {
    end <- min(i + batch_size - 1, n)
    batch <- data[i:end]  # Corrected line
    
    op <- OPseq(batch, emb = emb, lag = lag)
    oprob <- OPprob(batch, emb = emb)
    
    # Diagnostic code: check and correct probabilities
    if (any(oprob < 0)) {
      warning("Negative probabilities detected; setting to zero")
      oprob[oprob < 0] <- 0  # Ensure no negative probabilities
    }
    
    if (sum(oprob) == 0){
      warning("All probabilities are zero; setting equal probabilities")
      oprob <- rep(1/length(oprob), length(oprob))
    } else {
      oprob <- oprob / sum(oprob) # Normalize probabilities to sum to 1
    }
    
    op_list[[length(op_list) + 1]] <- op
    oprob_list[[length(oprob_list) + 1]] <- oprob
  }
  
  return(list(op = unlist(op_list), oprob = unlist(oprob_list))) #Corrected line
}

# Process data in batches
Normal_DE_results <- process_batch(normal$DE_time)
Normal_FE_results <- process_batch(normal$FE_time)
k48_DE_results <- process_batch(fourty_eightk$DE_time)
k48_FE_results <- process_batch(fourty_eightk$FE_time)

# Extract results
OP_Normal_DE <- Normal_DE_results$op
OP_Normal_FE <- Normal_FE_results$op
OP_48k_DE <- k48_DE_results$op
OP_48k_FE <- k48_FE_results$op

oprob_Normal_DE <- Normal_DE_results$oprob
oprob_Normal_FE <- Normal_FE_results$oprob
oprob_48k_DE <- k48_DE_results$oprob
oprob_48k_FE <- k48_FE_results$oprob

# Convert to dataframe for plotting
ordinal_patterns_table <- as.data.frame(table(OP_Normal_DE, OP_Normal_FE, OP_48k_DE, OP_48k_FE))

# Convert probability to dataframe for plotting
probability_table <- as.data.frame(table(oprob_Normal_DE, oprob_Normal_FE, oprob_48k_DE, oprob_48k_FE))

# Compute entropy and complexity
H_Shannon_Normal_DE <- tryCatch({
  HShannon(oprob_Normal_DE)
}, error = function(e) {
  warning("HShannon failed: ", e$message)
  return(NA)  # Or some other appropriate default
})

C_Complexity_Normal_DE <- tryCatch({
  StatComplexity(oprob_Normal_DE)
}, error = function(e) {
  warning("StatComplexity failed: ", e$message)
  return(NA)
})

H_Shannon_Normal_FE <- tryCatch({
  HShannon(oprob_Normal_FE)
}, error = function(e) {
  warning("HShannon failed: ", e$message)
  return(NA)  # Or some other appropriate default
})

C_Complexity_Normal_FE <- tryCatch({
  StatComplexity(oprob_Normal_FE)
}, error = function(e) {
  warning("StatComplexity failed: ", e$message)
  return(NA)
})

H_Shannon_48k_DE <- tryCatch({
  HShannon(oprob_48k_DE)
}, error = function(e) {
  warning("HShannon failed: ", e$message)
  return(NA)  # Or some other appropriate default
})

C_Complexity_48k_DE <- tryCatch({
  StatComplexity(oprob_48k_DE)
}, error = function(e) {
  warning("StatComplexity failed: ", e$message)
  return(NA)
})

H_Shannon_48k_FE <- tryCatch({
  HShannon(oprob_48k_FE)
}, error = function(e) {
  warning("HShannon failed: ", e$message)
  return(NA)  # Or some other appropriate default
})

C_Complexity_48k_FE <- tryCatch({
  StatComplexity(oprob_48k_FE)
}, error = function(e) {
  warning("StatComplexity failed: ", e$message)
  return(NA)
})

# Organize results into a data frame
results_table_batch <- data.frame(
  Dataset = c("Normal_DE", "Normal_FE", "48k_DE", "48k_FE"),
  Shannon_Entropy_H = c(H_Shannon_Normal_DE, H_Shannon_Normal_FE, H_Shannon_48k_DE, H_Shannon_48k_FE),
  Statistical_Complexity_C = c(C_Complexity_Normal_DE, C_Complexity_Normal_FE, C_Complexity_48k_DE, C_Complexity_48k_FE)
)

# Print results
print(results_table_batch)
