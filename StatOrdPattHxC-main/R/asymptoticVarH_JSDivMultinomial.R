#' Function that computes the asymptotic variance of the sum and difference
#' between the normalized Shannon entropy and the Jensen-Shannon divergence
#' under the Multinomial model
#'
#' @param p a probability function
#' @param n the number of independent trials
#' @param operation "s" for the sum and "d" for de difference
#' @returns A scalar
#'
#' @keywords internal

asymptoticVarH_JSDivMultinomial <- function(p, n, operation){

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

  if (length(ppos) > 1){
    PartialDer <- diag(const * (1 + log(ppos)) - 1 - log((ppos + 1/k) / 2))
  } else {
    PartialDer <- as.matrix(const * (1 + log(ppos)) - 1 - log((ppos + 1/k) / 2))
  }

  S <- SigmaMultinomial(ppos)

  SDelta <- PartialDer %*% S %*% t(PartialDer)

  sig <- sum(SDelta) / 4
  return(sig / n)
}
