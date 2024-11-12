#' Function that computes the normalized Shannon entropy (using natural logarithms) of a probability function
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)},
#' its normalized Shannon entropy is
#' \deqn{S^{(S)}[\bm p] = -\frac{1}{\ln k}\sum_{i=1}^k p_i \ln p_i.}
#' We adopt the convention \eqn{0\ln 0 = 0}.
#'
#' @usage HShannon(p)
#' @param p a probability function
#' @returns A value in \eqn{[0,1]}
#'
#' @export
#'
#' @examples
#' # Time series
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000) # white noise
#' y <- mov.av(x, order=50) # smoothed with moving averages
#' # Ordinal Patterns
#' op.wn.4 <- OPprob(x, emb=4)
#' op.ma.4 <- OPprob(y, emb=4)
#' # Shannon entropies
#' HShannon(op.wn.4)
#' HShannon(op.ma.4)

HShannon <- function(p){

  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){

    prob <- p[p > 0]
    N <- length(p)
    H <- -sum(prob * log(prob)) / log(N)

    return(H)
  } else {
    message("ERROR: Not a valid probability function")
    return(NULL)
  }
}

