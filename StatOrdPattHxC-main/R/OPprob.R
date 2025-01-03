#' Function that finds the Ordinal Patterns probabilities in a time series for a given embedding dimension
#'
#' @usage OPprob(TS, emb)
#' @param TS time series of length \eqn{n-m+1}
#' @param emb embedding dimension \eqn{m}
#' @returns a sequence of \eqn{n} patterns
#'
#' @name OPprob
#'
#' @import tibble
#' @import dplyr
#' @import prodlim
#'
#' @export
#'
#' @examples
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000) # white noise
#' y <- mov.av(x, order=11) # smoothed with moving averages
#' OPprob(x, emb=4)
#' OPprob(y, emb=4)

utils::globalVariables("OP")

OPprob <- function(TS, emb){

  op <- tibble::tibble(
    OP = factor(OPseq(TS, emb), levels = 1:factorial(emb))
  )

  fr <- op %>% count(OP, .drop = FALSE)
  probs <- fr$n / sum(fr$n)
  return(probs)
}


