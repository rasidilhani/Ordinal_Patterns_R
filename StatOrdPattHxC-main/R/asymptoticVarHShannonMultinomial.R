#' Function that computes the asymptotic variance of the normalized Shannon 
#' entropy under the Multinomial model 
#' 
#' @param p a probability function
#' @param n the number of independent trials
#' @returns A scalar
#' 
#' @export

asymptoticVarHShannonMultinomial <- function(p, n){
  
  k <- length(p)
  ppos <- p[p > 0]
  
  if (length(ppos) > 1){
    PartialDer <- diag(log(ppos) + 1)
  } else {
    PartialDer <- as.matrix(log(ppos) + 1)
  }
  
  S <- SigmaMultinomial(ppos)
  
  SDelta <- PartialDer %*% S %*% t(PartialDer)
  
  sig <- sum(SDelta) / (log(k))^2
  return(sig / n)
}