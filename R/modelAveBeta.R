#' calculate model averaged estimates of harmonic trend (X %*% beta)
#'
#' @param fit object of type "ihmm"
#' @param zbest1 results from bestClusteriHMM
#' @param X covariates 
#'
#' @importFrom stats quantile
#' @return model averaged estimates of X %*% beta, state-specific covariate effect
#' @export
#'
modelAveBeta <- function(fit, zbest1, X){
  # model averaged estimates of the harmonic trend (beta only)
  K <- max(unlist(zbest1)) # number of best clusters 
  niter <- length(fit$beta.save);niter
  betaXma <- list()
  betaXupr <- list()
  betaXlwr <- list()
  for(k in 1:K){
    whoBest <- which(unlist(zbest1)==k)
    beta.ave <- matrix(NA, 4, niter)
    beta.x.ave <- matrix(NA, nrow(X), niter)
    for(s in 1:niter){
      these.states <- unlist(fit$z.save[[s]])[whoBest] # average beta over these states
      beta.matrix <- matrix(0, 4, 1) # place holder 
      for(j in 1:length(these.states)){
        beta.matrix <- cbind(beta.matrix, unlist(fit$beta.save[[s]][these.states[j]]))
      }
      beta.ave[,s] <- apply(matrix(beta.matrix[,-1], nrow = 4), 1, mean)
      beta.x.ave[,s] <- as.numeric(beta.ave[,s]%*%t(X))
      
    }

    # distribution of X%*%beta 
    betaXma[[k]] <- apply(beta.x.ave, 1, mean)
    betaXlwr[[k]] <- apply(beta.x.ave, 1, FUN = function(x) quantile(x, .05))
    betaXupr[[k]] <- apply(beta.x.ave, 1, FUN = function(x) quantile(x, .95))
    
  }
  return(list(betaXma = betaXma, betaXupr = betaXupr, betaXlwr = betaXlwr))
}
