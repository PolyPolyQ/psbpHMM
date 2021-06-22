#' Simulate data for psbpHMM package 
#'
#' Simulate multiple time series data from a fixed-state MVN hidden markov model ∂
#'
#' @param n number of time series
#' @param t.max length of each time series
#' @param K number of clusters to be shared across all time series 
#' @param p number of exposures 
#' @param tempTrend logical; if true, all time series have a common temporal trend
#' @param lodRemove logical; if true, remove data emulating observations below the limit of detection
#' @param marRemove logical; if true, remove data emulating observations missing at random 
#' @param lodmis percent of data to remove as below limit of detection 
#' @param marmis percent of data to remove as missing at random 
#' @param sf parameter controlling level of variation in state-specific covariance matrices
#'
#' @return a list with components
#' \itemize{
#'     \item y: scaled simulated data with missing observations set to NA for MAR or -Inf for below LOD 
#'     \item y.complete: scaled complete data 
#'     \item z.true: true hidden state trajectories for each time series
#'     \item K: true number of hiddens states
#'     \item mu.true: scaled true state-specific means
#'     \item Sigma.true: unscaled true state-specific covariance matrices
#'     \item lod: limit of detection
#' }     
#' @export
#'
#' @examples
simdat <- function(n, t.max=24, K = 6, p = 3, tempTrend = TRUE, lodRemove = FALSE, marRemove = FALSE,
                         lodmis = 0, marmis = 0, sf = 0.1){
  
  #K <- 6
  if(tempTrend){
    z.true <- list()
    for(i in 1:n){
      probs <- rdirichlet(1, rep(10, K));probs # different for each person, but from same dist and same order
      nums <- rmultinom(1, t.max, prob = probs); nums; sum(nums)
      num1 <- nums[1]; num1
      numFirsthalf <- floor(nums[1]/2); numFirsthalf
      numSecondhalf <- num1 - numFirsthalf; numSecondhalf
      nums[1] <- numFirsthalf; nums
      nums <- c(nums, numSecondhalf); nums; sum(nums)
      states <- c(1:K, 1)
      z <- unlist(lapply(1:(K+1), FUN = function(j) rep(states[j], nums[j])))
      #idx <- sample(c(0,1,2), 1)
      z.true[[i]] <- z 
    }
  }else{
    # set states: this needs to be a temporal trend, but not shared among people 
    z.true <- list()
    for(i in 1:n){
      probs <- rdirichlet(1, rep(10, K))
      ord <- sample(1:K, size = K, replace = FALSE) # order
      samps <- sample(1:K, t.max, replace = TRUE, prob = probs) # count
      reps <- sapply(1:K, FUN = function(i) length(which(samps==i))) # num 
      states <- c(unlist(lapply(1:K, FUN = function(k) rep(ord[k], reps[k])))); states
      half1 <- ceiling(0.5*length(which(states==states[1])))
      z <- c(states[-(1:half1)], rep(states[1], half1))
      #idx <- sample(c(0,1,2), 1)
      z.true[[i]] <- z 
    }
  }
  
  #p <- 3 # number of exposures
  
  # different covariances for each cluster 
  Sigma.true <- list()
  D <- sf*diag(p)
  for(k in 1:K){
    L <- diag(p)
    lowerTriangle(L) <- rnorm(p, 0, .5) 
    Sigma.true[[k]] <- solve(L)%*%D%*%t(solve(L))
  }
  
  # different means for each cluster 
  corMean <- matrix(c(1, 0.7, 0.4, 0.7, 1, -0.2, 0.4, -0.2, 1), ncol = 3); is.positive.definite(corMean)
  mu.true <- rmvnorm(K, mean = rep(0, p), sigma = corMean); mu.true
  
  y <- list()
  for(i in 1:n){
    y[[i]] <- matrix(0, t.max, p) # p exposures at t time points 
    for(k in 1:K){
      if(any(z.true[[i]]==k)){
        which.zk <- which(z.true[[i]] == k)
        y[[i]][which.zk,] <- rmvnorm(length(which.zk), mu.true[k,], Sigma.true[[k]])
      }
    }
  }
  
  # standardize data 
  yMat <- numeric()
  for(i in 1:n){
    yMat <- rbind(yMat, y[[i]])
  }
  
  # scale mu for MSE 
  meanY <- apply(yMat, 2, mean)
  sdY <- apply(yMat, 2, sd)
  mu.true.sc <- cbind((mu.true[,1] - meanY[1])/sdY[1],
                      (mu.true[,2] - meanY[2])/sdY[2],
                      (mu.true[,3] - meanY[3])/sdY[3])
  
  # i don't use sigma again so not worth scaling it 
  ys <- apply(yMat, 2, scale)
  
  # relist 
  splits <- seq(1,n*t.max, t.max)
  y.list <- lapply(1:length(splits), FUN = function(i) ys[splits[i]:(splits[i]+t.max-1),])
  y <- y.list
  y.complete <- y # list of complete data 
  
  
  ymat <- NULL
  for(i in 1:n){
    ymat <- rbind(ymat, y[[i]])
  }
  
  # chunk removal
  if(marRemove){
    # chunks of size 1 to 10
    misnum <- ceiling(marmis*t.max*p)
    numChunks <- marmis*200; numChunks
    chunks <- rmultinom(1, size = misnum, prob = rep(1/numChunks, numChunks)) # generate size of chunks 
    for(i in 1:n){
      for(k in 1:numChunks){
        j <- sample(1:p, 1, prob = rep(1/p,p)) # randomly choose an exposure
        repeat{
          timePoint <- sample(2:(t.max-chunks[k,]),1, prob = rep(1/(t.max-chunks[k,]), (t.max-chunks[k,]-1)) ) # choose a time point to start at
          if(!anyNA(y[[i]][,j][timePoint:(timePoint+chunks[k,]-1)])){
            break
          }
        }
        y[[i]][,j][timePoint:(timePoint+chunks[k,]-1)] <- NA
      }
    }
  }
  
  
  if(lodRemove){
    lod = list()
    lod.old = apply(ymat, 2, FUN = function(x) quantile(x, lodmis, na.rm = TRUE))
    
    for(i in 1:n){
      lod[[i]] = lod.old
      #lod[[i]] = apply(y[[i]], 2, FUN = function(x) quantile(x, lodmis, na.rm = TRUE))
      for(j in 1:p){
        y[[i]][which(y[[i]][,j] <= lod[[i]][j]),j] <- -Inf
      }
    }
  }else{
    lod <- NULL
  }
  
  list1 <- list(y = y, y.complete = y.complete,
                z.true = z.true, K = K, mu.true = mu.true.sc, Sigma.true = Sigma.true,
                lod = lod)
  return(list1)
}
