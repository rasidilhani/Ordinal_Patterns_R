#' Function that computes the asymptotic variance of the Jensen-Shannon divergence
#' under the Multinomial model 
#' 
#' @param p a probability function
#' @param n the number of independent trials
#' @returns A scalar
#' 
#' @export

asymptoticVarJSDivergenceMultinomial <- function(p, n){
  
  k <- length(p)
  ppos <- p[p > 0]
  
  if (length(ppos) > 1){
    PartialDer <- diag(log(ppos) - log((ppos + 1/k)/2))
  } else {
    PartialDer <- as.matrix(log(ppos) - log((ppos + 1/k)/2))
  }
  
  S <- SigmaMultinomial(ppos)
  
  SDelta <- PartialDer %*% S %*% t(PartialDer)
  
  sig <- sum(SDelta) / 4
  return(sig / n)
}
