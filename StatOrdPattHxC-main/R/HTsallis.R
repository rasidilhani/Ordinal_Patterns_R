#' Function that computes normalized Tsallis entropy of index beta of a probability function
#'
#' Given a probability function \eqn{\bm p=(p_1,p_2,\dots,p_k)}, its
#' Tsallis entropy with index \eqn{\beta \in \mathbb{R}\setminus\{1\}}$ is
#' \deqn{
#' S^{(T_{\beta})}[\bm p] = \sum_{i=1}^k \frac{p_i - p_i^\beta}{\beta-1}.}
#'
#' @param p a probability function
#' @param beta the parameter beta for the Tsallis entropy; default 1.5
#' @usage HTsallis(p, beta)
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
#' # Tsallis entropies
#' HTsallis(op.wn.4, beta=1.5)
#' HTsallis(op.ma.4, beta=1.5)

HTsallis <- function(p, beta=1.5){
  if(length(p) >= 2 & min(p) >= 0 & sum(p) <= (1+.Machine$double.eps)){
  N <- length(p)
  H <- sum(p - p^beta) / (1 - N^(1-beta))

  return(H)} else {
    message("ERROR: Not a valid probability function")
    return(NULL)
  }
}

