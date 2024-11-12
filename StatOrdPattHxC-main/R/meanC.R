#' Function that computes the mean of the Normal distribution that approximates
#' the asymptotic probability density function of the
#' statistical complexity under the Multinomial model 
#' 
#' @param p a probability function
#' @param n the number of independent trials
#' @returns A scalar
#' 
#' @export

meanC <- function(p, n){
  
  k <- length(p)
  
  ro <- correlationH_JSDivMultinomial(p, n)
  Q0 <- -2 / ((k + 1) / k * log(k + 1) - 2 * log(2*k) + log(k))
  
  mH <- HShannon(p)
  sH <- sqrt(asymptoticVarHShannonMultinomial(p, n))
  mDiv <- asymptoticMeanJSDivergenceMultinomial(p)
  sDiv <- sqrt(asymptoticVarJSDivergenceMultinomial(p, n))
  
  deltaH <- mH / sH
  deltaQ <- mDiv / sDiv
  
  return(Q0 * sH * sDiv * (deltaH * deltaQ + ro))
  
}