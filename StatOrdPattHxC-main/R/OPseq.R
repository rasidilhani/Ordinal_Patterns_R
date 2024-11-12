#' Function that finds the OP sequence of a time series, given the
#' embedding dimension and lag
#'
#' @param TS time series
#' @param emb embedding dimension
#' @param lag time lag (default value: 1)
#'
#' @importFrom stats rnorm
#' @export
#'
#' @examples
#' # Generate a time series and compute its ordinal patterns
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000) # white noise
#' OPseq(x, emb=4, lag=1)
#'


OPseq <- function(TS, emb, lag=1){
  # number of OP for lag = 1
  el <- length(TS) - emb + 1

  # OP sequence for lag = 1
  seqOP <- vector()
  for (i in 1:el){
    seqOP[i] <- pi_i(ind_pos(TS[i:(i + emb - 1)]))
  }

  # OP sequence for the given lag
  return(seqOP[seq(1, el, lag)])
}

