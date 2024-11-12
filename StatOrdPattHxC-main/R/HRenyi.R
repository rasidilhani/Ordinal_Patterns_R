#' Function that computes normalized Rényi entropy of order beta
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)}, its Rényi entropy of order \eqn{\beta \in \mathbb{R}^+\setminus\{1\}} is
#' \deqn{S^{(R_{\beta})}[\bm p] = \frac{1}{1-\beta} \log \sum_{i=1}^k p_i^\beta.}
#' We adopt the convention \eqn{0\ln 0 = 0}.
#'
#' @param p a probability function
#' @param beta the parameter beta for the Rényi entropy; default 1.5
#' @usage HRenyi(p, beta)
#' @returns a value in \eqn{[0,1]}
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
#' # Rényi entropies
#' HRenyi(op.wn.4, beta=1.5)
#' HRenyi(op.ma.4, beta=1.5)

HRenyi <- function(p, beta=1.5){

  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){

    N <- length(p)
    H <- log(sum(p^beta)) / ((1-beta) * log(N))

    return(H)} else {
      message("ERROR: Not a valid probability function")
      return(NULL)
    }
}


