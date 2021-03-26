#' Update W for no covariates
#'
#' @param t t
#' @param l l
#' @param i i
#' @param alpha.0k alpha0 
#' @param alpha.jk alpha.jk
#' @param z states
#'
#' @importFrom truncnorm rtruncnorm 
#'
#' @return w 
#' @export
#'
upWnone <- function(t,l,i, alpha.0k, alpha.jk, z){
  if(t==1){
    trunc.mean <- as.numeric(alpha.0k[l])
  }else{
    trunc.mean <- as.numeric(alpha.jk[[z[[i]][t-1]]][l])
  }
  return((z[[i]][t]==l)*rtruncnorm(1, a = 0, b = Inf, mean = trunc.mean, sd = 1) +
           (z[[i]][t]>l)*rtruncnorm(1, a = -Inf, b = 0, mean = trunc.mean, sd = 1))
}
