#' Function that computes the asymptotic mean of the Jensen-Shannon divergence
#' under the Multinomial model
#'
#' @param p a probability function
#' @returns The asymptotic mean (a scalar)
#'
#' @keywords internal

asymptoticMeanJSDivergenceMultinomial <- function(p){

  k <- length(p)
  ppos <- p[p > 0]

  return((sum(ppos * log(ppos)) -
            sum((p + 1/k) * log((p + 1/k) / 2)) - log(k))/2)

}
