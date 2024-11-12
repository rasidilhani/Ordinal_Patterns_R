#' Function that computes the variance of the Normal distribution that approximates
#' the asymptotic probability density function of the
#' statistical complexity under the Multinomial model 
#' 
#' @param p a probability function
#' @param n the number of independent trials
#' @returns A scalar
#' 
#' @export

varC <- function(p, n){
  
  k <- length(p)
  
  ro <- correlationH_JSDivMultinomial(p, n)
  Q0 <- -2 / ((k + 1) / k * log(k + 1) - 2 * log(2*k) + log(k))
  
  mH <- HShannon(p)
  vH <- asymptoticVarHShannonMultinomial(p, n)
  mDiv <- asymptoticMeanJSDivergenceMultinomial(p)
  vDiv <- asymptoticVarJSDivergenceMultinomial(p, n)
  
  deltaH <- mH / sqrt(vH)
  deltaQ <- mDiv / sqrt(vDiv)
  
  return(Q0^2 * vH * vDiv * (deltaH^2 + deltaQ^2 + 2 * ro * deltaH * deltaQ + 
                               1 + ro^2))
  
}