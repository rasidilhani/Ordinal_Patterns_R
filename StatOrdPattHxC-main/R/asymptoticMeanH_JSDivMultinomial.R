#' Function that computes the asymptotic mean of the sum and difference 
#' between the normalized Shannon entropy and the Jensen-Shannon divergence
#' under the Multinomial model 
#' 
#' @param p a probability function
#' @param operation "s" for the sum and "d" for de difference
#' @returns A scalar
#' 
#' @export

asymptoticMeanH_JSDivMultinomial <- function(p, operation){ 
  
  k <- length(p)
  ppos <- p[p > 0]
  
  switch(operation,
         
         # sum
         "s" = {
           const <- (log(k) - 2) / log(k)
         },
         
         # difference
         "d" = {
           const <- (log(k) + 2) / log(k)
         })
  
  return((sum(const * ppos * log(ppos)) - 
            sum((p + 1/k) * log((p + 1/k) / 2)) - log(k))/2)
  
}
