#' Update W with or without repeated measures 
#'
#' @param t time point
#' @param l state
#' @param i time series
#' @param alpha.0k alpha0k
#' @param X covariates
#' @param beta.k overall beta
#' @param beta.sk subject specific beta
#' @param alpha.jk alphajk
#' @param z hidden states 
#'
#' @return sampled value for W 
#' @export
#'
updateW <- function(t,l,i, alpha.0k, X, beta.k, beta.sk = NULL, alpha.jk, z){
  
  if(!is.null(beta.sk)){
    
    if(t==1){
      trunc.mean <- as.numeric(alpha.0k[l] + crossprod(X[t,],beta.k[[l]]) + crossprod(X[t,],beta.sk[[l]])) 
    }else{
      trunc.mean <- as.numeric(alpha.jk[[z[[i]][t-1]]][l] + crossprod(X[t,],beta.k[[l]]) + crossprod(X[t,],beta.sk[[l]]))
    }
    return((z[[i]][t]==l)*rtruncnorm(1, a = 0, b = Inf, mean = trunc.mean, sd = 1) +
             (z[[i]][t]>l)*rtruncnorm(1, a = -Inf, b = 0, mean = trunc.mean, sd = 1))
    
  }else{
    if(t==1){
      trunc.mean <- as.numeric(alpha.0k[l] + crossprod(X[t,],beta.sk[[l]])) 
    }else{
      trunc.mean <- as.numeric(alpha.jk[[z[[i]][t-1]]][l] + crossprod(X[t,],beta.sk[[l]]))
    }
    return((z[[i]][t]==l)*rtruncnorm(1, a = 0, b = Inf, mean = trunc.mean, sd = 1) +
             (z[[i]][t]>l)*rtruncnorm(1, a = -Inf, b = 0, mean = trunc.mean, sd = 1))
  }
  
}
