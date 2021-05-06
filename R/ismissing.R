#' function to indicate missing values 
#'
#' @param x data
#'
#' @return matrix same size as x, NA's are 1, -Inf's are 2, observed data are 0
#' @export
#'

ismissing <- function(x) {
  if(is.na(x)) return(1) # MAR
  else if(x==-Inf) return(2) # lod
  else return(0) 
  } # observed
