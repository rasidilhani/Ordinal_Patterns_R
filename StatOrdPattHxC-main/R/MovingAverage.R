#' Moving average
#'
#' Computes the moving average of a time series over a specified window size
#'
#' @details
#' Used only for the examples
#'
#' @export
#'
#' @usage mov.av(x, order)
#' @param x time series
#' @param order window size (default: 1, no filter)
#' @returns the filtered time series, smaller that the input (only values are returned)
#'
#' @importFrom stats filter
#'
#' @examples
#' set.seed(1234567890, kind="Mersenne-Twister")
#' x <- rnorm(1000)
#' x.ma <- mov.av(x, order=11)
#'

mov.av <- function(x, order=1){
  full <- unclass(stats::filter(x, rep(1 / order, order), sides = 2))
  return(full[!is.na(full)])
  }
