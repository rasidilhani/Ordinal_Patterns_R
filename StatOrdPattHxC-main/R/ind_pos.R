#' Function that assigns the corresponding OP to a partition (word)
#'
#' @description
#' Internal function. Candidate for C++ implementation.
#'
#'
#' @param part a partition of a time series
#'

ind_pos <- function(part){
  v <- sort(unique(part))
  count <- 0
  pos <- vector()
  for (i in 1:length(v)){
    id <- which(part == v[i])
    for (j in 1:length(id)){
      count <- count + 1
      pos[id[j]] <- count
    }
  }
  return(pos-1)
}


