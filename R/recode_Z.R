#' Relabel z to 1:length(z) while retaining order of z
#'
#' @param z list of states
#'
#' @export 
#' @return recoded z

recode_Z <- function(z){
  levels <- sort(unique(z))
  for (i in 1:length(levels)){
    z[which(z == levels[i])] <- i
  }
  return(z)
}