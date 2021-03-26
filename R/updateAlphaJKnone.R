#' update alpha.jk for no covariates
#'
#' @param j j
#' @param k k
#' @param n n
#' @param t.max t.max
#' @param z z
#' @param vinv.alpha vinv.alpha 
#' @param sig2inv.alpha sig2inv.alpha
#' @param w.z w.z
#' @param m.alpha m.alpha
#' @param mu.alpha mu.alpha
#'
#' @return alpha.jk
#' @export
#'
updateAlphaJKnone <- function(j, k, n, t.max, z, vinv.alpha, sig2inv.alpha, w.z, m.alpha, mu.alpha){
  
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
        return(sum(sapply(itimes[[i]], FUN = function(t)  w.z[[i]][[t]][k])))
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
