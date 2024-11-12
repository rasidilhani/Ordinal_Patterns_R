#' Function that computes the normalized Fisher information measure of a probability function
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)},
#' its Fisher information measure is
#' \deqn{H^{(F)}[\bm p] = 4  \sum_{i=1}^{k-1} (\sqrt{p_{i+1}} - \sqrt{p_i})^2.}
#'
#' @usage HFisher(p)
#' @param p a probability function
#' @returns a value in \eqn{[0,1]}
#'
#' @export
#'
#' @examples
#' # Time series
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000) # white noise
#' y <- mov.av(x, order=51) # smoothed with moving averages
#' # Ordinal Patterns
#' op.wn.4 <- OPprob(x, emb=4)
#' op.ma.4 <- OPprob(y, emb=4)
#' # Fisher information measures
#' HFisher(op.wn.4)
#' HFisher(op.ma.4)

HFisher <- function(p){
  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){
    p1 <- p[-1]
    p2 <- p[-length(p)]

    H <- 4 * sum((sqrt(p1) - sqrt(p2))^2)

    return(H)} else {
      message("ERROR: Not a valid probability function")
      return(NULL)
    }
}


