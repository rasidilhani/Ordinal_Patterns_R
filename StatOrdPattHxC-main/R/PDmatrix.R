#' Computes the matrix of partial derivatives
#' a given vector of probabilities without zero entries, and
#' a given entropy
#'
#' @param q vector of probabilities
#' @param ent type of entropy: : S (Shannon), T (Tsallis, \eqn{beta} is required), R (Rényi, \eqn{beta} is required), or F (Fisher)
#' @param beta (only required for Tsallis and Rényi; default value 0.9)

PDmatrix <- function(q, ent, beta = 0.9){
  
  qpos <- q[q>0]
  
  switch(ent,
         
         # Shannon
         "S" = {
           J <- diag(log(qpos) + 1)
         },
         
         # Tsallis
         "T" = {
           J <- diag(1 - beta * qpos^(beta - 1))
         },
         
         # Rényi
         "R" = {
           J <- diag(beta * qpos^(beta - 1))
         },
         
         # Fisher
         "F" = {
           # Number of probabilities
           a <- length(qpos)
           
           # Auxiliar vectors delete the last element (q1) and the first element (q2)
           q1 <- qpos[-a]
           q2 <- qpos[-1]
           
           J <- matrix(0, nrow = a-1, ncol = a)
           
           # Find the diagonal
           ii <- (sqrt(q1) - sqrt(q2)) / sqrt(q1)
           # Find the other elements
           ij <- (sqrt(q2) - sqrt(q1)) / sqrt(q2)
           
           for (i in 1:(a-1)){
             J[i, i] <- ii[i]
             J[i, i+1] <- ij[i]
           }
         }
  )
  
  return(J)
}
