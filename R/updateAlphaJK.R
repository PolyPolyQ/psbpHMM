#' Update alpha_jk
#'
#' @param j previous state
#' @param k current state
#' @param n number of time series
#' @param t.max length of time series
#' @param z hidden states
#' @param vinv.alpha prior
#' @param sig2inv.alpha prior
#' @param w.z w
#' @param X covariates
#' @param beta.k overall beta
#' @param beta.sk subject specific beta
#' @param m.alpha prior
#' @param mu.alpha prior
#'
#' @return sampled value for alpha_jk

updateAlphaJK <- function(j, k, n, t.max, z, vinv.alpha, sig2inv.alpha, w.z, X, beta.k, beta.sk = NULL, 
                             m.alpha, mu.alpha){
  
  itimes <- lapply(1:n, FUN = function(i) {
    which(sapply(2:t.max, FUN = function(t) (z[[i]][t]>=k & z[[i]][t-1] == j))) + 1 # because I start at 2
  })
  njk.tilde <- length(unlist(itimes))
  if(njk.tilde > 0){
    if(j==k){
      sig2.jk <- 1/(vinv.alpha + njk.tilde) 
    }else{
      sig2.jk <- 1/(sig2inv.alpha + njk.tilde) 
    }
    
    sum.w <- sum(sapply(1:n, FUN = function(i){
      if(any(itimes[[i]] > 0)){
        # sum over t for this i
        return(sum(sapply(itimes[[i]], 
                          FUN = function(t) wMinusb(i=i, t=t, k=k, w.z=w.z, beta.k = beta.k, beta.sk=beta.sk[[i]], X=X[[i]]))))
      }else{
        return(0)
      }
    }))
    
    if(j==k){
      mu.jk <- sig2.jk*(m.alpha*vinv.alpha + sum.w) # when j=k use a different prior mean 
    }else{
      mu.jk <- sig2.jk*(mu.alpha*sig2inv.alpha + sum.w)
    }
    return(rnorm(1, mu.jk, sqrt(sig2.jk)))
  }else{
    # update from prior 
    if(j==k){
      return(rnorm(1, m.alpha, sqrt(1/vinv.alpha)))
    }else{
      return(rnorm(1, mu.alpha, sqrt(1/sig2inv.alpha)))
    }
  }
}
