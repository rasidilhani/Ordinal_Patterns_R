#' Function that computes the correlation coefficient
#' between the normalized Shannon entropy and the Jensen-Shannon divergence
#' under the Multinomial model 
#' 
#' @param p a probability function
#' @param n the number of independent trials
#' @returns A scalar
#' 
#' @export

correlationH_JSDivMultinomial <- function(p, n){
  
  mH <- HShannon(p)
  sH <- sqrt(asymptoticVarHShannonMultinomial(p, n))
  mDiv <- asymptoticMeanJSDivergenceMultinomial(p)
  sDiv <- sqrt(asymptoticVarJSDivergenceMultinomial(p, n))
  mpos <- asymptoticMeanH_JSDivMultinomial(p, "s")
  mneg <- asymptoticMeanH_JSDivMultinomial(p, "d")
  s2pos <- asymptoticVarH_JSDivMultinomial(p, n, "s")
  s2neg <- asymptoticVarH_JSDivMultinomial(p, n, "d")
  
  return((s2pos + mpos^2 - s2neg - mneg^2 - 4 * mDiv * mH) / (4 * sDiv * sH))
  
}