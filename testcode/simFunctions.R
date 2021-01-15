
# simulate data 
simdatsimple <- function(n, t.max=24, tempTrend = TRUE, lodRemove = FALSE, marRemove = FALSE,
                         lodmis = 0, marmis = 0, sf = 1){
  
  K <- 6
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
  
  p <- 3 # number of exposures
  
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
    lod <- apply(ymat, 2, FUN = function(x) quantile(x, lodmis))
    for(i in 1:n){
      for(j in 1:p){
        y[[i]][which(y[[i]][,j] <= lod[j]),j] <- -Inf
      }
    }
  }else{
    lod <- NULL
  }
  
  # # check
  # ycomp <- rbind(y[[1]],y[[2]])
  # length(which(ycomp[,1]<lod[1]))
  # length(which(ycomp[,2]<lod[2]))
  # length(which(ycomp[,3]<lod[3]))
  # 
  # #check how much missing
  # for(i in 1:n){
  #   print(length(which(is.na(unlist(y[[i]])))))
  #   print(length(which(unlist(y[[i]])==-Inf)))
  # }
  
  length(which(unlist(y)==-Inf))
  
  ycomp <- rbind(y.complete[[1]], y.complete[[2]])
  ymiss <- rbind(y[[1]], y[[2]])
  
  length(which(ycomp[,1]<lod[1]))
  length(which(ymiss[,1]==-Inf))
  
  length(which(ycomp[,2]<lod[2]))
  length(which(ymiss[,2]==-Inf))
  
  length(which(ycomp[,3]<lod[3]))
  length(which(ymiss[,3]==-Inf))
  
  
  list1 <- list(y = y, y.complete = y.complete,
                z.true = z.true, K = K, mu.true = mu.true.sc, Sigma.true = Sigma.true,
                lod = lod)
  return(list1)
}


###############################################################
### function for unlisting results from independent methods ###
###############################################################

unlistMethod <- function(fit1, marmiss=NULL, lodmiss=NULL, y, ycomplete){
  
  n = length(fit1)
  t.max = 288
  
  ham1 <- unlist(mclapply(1:n, FUN = function(i){
    if(!anyNA(fit1[[i]]$hamming)) return(fit1[[i]]$hamming)
    else return(NULL)
  }))
  
  mu.sse1 <- mclapply(1:n, FUN = function(i){
    if(!anyNA(fit1[[i]]$mu.sse)) return(fit1[[i]]$mu.sse)
    else return(NULL)
  })
  muMat <- numeric()
  for(i in 1:n){
    muMat <- cbind(muMat, mu.sse1[[i]])
  }
  mu.mse1 <- apply(muMat,1, sum)/(n*t.max)
  
  if(!is.null(marmiss)){
    # mar MSE #
    mar.sse1 <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$mar.sse)) return(fit1[[i]]$mar.sse)
      else return(NULL)
    })
    marMat <- numeric()
    for(i in 1:n){
      marMat <- cbind(marMat, mar.sse1[[i]])
    }
    mar.mse1 <- apply(marMat,1,sum)/marmiss # sum the rows to get total sse, then divide by total mar data
    
    # posterior predictive mean imputes 
    # ppmarMean
    ppmar.sse <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$mar.sse)) {
        truemar <- unlist(ycomplete[[i]])[which(is.na(unlist(y[[i]])))]
        ymarmean <- apply(matrix(fit1[[i]]$ymar, nrow = len.imp), 2, mean)
        ppmarSum <- sum((ymarmean - truemar)^2)
      }else return(NULL)
    })
    ppMarMean <- sum(unlist(ppmar.sse))/marmiss
    
    # mar bias #
    mar.sbias1 <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$mar.sum.bias)) return(fit1[[i]]$mar.sum.bias)
      else return(NULL)
    })
    marMatb <- numeric()
    for(i in 1:n){
      marMatb <- cbind(marMatb, mar.sbias1[[i]])
    }
    mar.bias1 <- apply(marMatb,1,sum)/marmiss # sum the rows to get total sse, then divide by total mar data
  }else{
    mar.mse1 = NULL
    mar.bias1 = NULL
  }
  
  if(!is.null(lodmiss)){
    # lod MSE #
    lod.sse1 <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$lod.sse)) return(fit1[[i]]$lod.sse)
      else return(NULL)
    })
    lodMat <- numeric()
    for(i in 1:n){
      lodMat <- cbind(lodMat, lod.sse1[[i]])
    }
    lod.mse1 <- apply(lodMat,1,sum)/lodmiss # sum the rows to get total sse, then divide by total mar data
    
    # pplodMean
    pplod.sse <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$lod.sse)) {
        truelod <- unlist(ycomplete[[i]])[which(unlist(y[[i]])==-Inf)]
        ylodmean <- apply(matrix(fit1[[i]]$ylod, nrow = len.imp), 2, mean)
        pplodSum <- sum((ylodmean - truelod)^2)
      }else return(NULL)
    })
    ppLodMean <- sum(unlist(pplod.sse))/lodmiss
    
    # lod bias #
    lod.sbias1 <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$lod.sum.bias)) return(fit1[[i]]$lod.sum.bias)
      else return(NULL)
    })
    lodMatb <- numeric()
    for(i in 1:n){
      lodMatb <- cbind(lodMatb, lod.sbias1[[i]])
    }
    lod.bias1 <- apply(lodMatb,1,sum)/lodmiss # sum the rows to get total sse, then divide by total mar data
  }else{
    lod.mse1 = NULL
    lod.bias1 = NULL
    ppLodMean = NULL
    ppMarMean = NULL
  }
  list1 = list(ham = ham1, mu.mse = mu.mse1, ppMarMean = ppMarMean, mar.mse = mar.mse1, 
               ppLodMean = ppLodMean, lod.mse = lod.mse1, mar.bias = mar.bias1, lod.bias = lod.bias1)
  
  
}


###############################################################
### function for unlisting a partial share methods in FCCS ###
###############################################################

unlistPartialShare <- function(fit1s, fit2s, fit3s, fit4s, marmiss=NULL, lodmiss=NULL){
  
  n1 = length(fit1s) # indep
  # n2 = length(fit2s) # shared 2
  # n3 = length(fit3s) # shared 2
  # n4 = length(fit4s) # shared 2
  t.max = 288
  
  # mar MSE #
  mar.sse1 <- mclapply(1:n1, FUN = function(i){
    if(!anyNA(fit1s[[i]]$mar.sse)) return(fit1s[[i]]$mar.sse)
    else return(NULL)
  })
  marMat <- numeric()
  for(i in 1:n1){
    marMat <- cbind(marMat, mar.sse1[[i]])
  }
  marMat <- cbind(marMat, fit2s$mar.sse, fit3s$mar.sse, fit4s$mar.sse)
  mar.mse1 <- apply(marMat,1,sum)/marmiss # sum the rows to get total sse, then divide by total mar data
  
  # mar bias #
  mar.sbias1 <- mclapply(1:n1, FUN = function(i){
    if(!anyNA(fit1s[[i]]$mar.sum.bias)) return(fit1s[[i]]$mar.sum.bias)
    else return(NULL)
  })
  marMatb <- numeric()
  for(i in 1:n1){
    marMatb <- cbind(marMatb, mar.sbias1[[i]])
  }
  marMatb <- cbind(marMatb, fit2s$mar.sum.bias, fit3s$mar.sum.bias, fit4s$mar.sum.bias)
  mar.bias1 <- apply(marMatb,1,sum)/marmiss # sum the rows to get total sse, then divide by total mar data
  
  # lod MSE #
  lod.sse1 <- mclapply(1:n1, FUN = function(i){
    if(!anyNA(fit1s[[i]]$lod.sse)) return(fit1s[[i]]$lod.sse)
    else return(NULL)
  })
  lodMat <- numeric()
  for(i in 1:n1){
    lodMat <- cbind(lodMat, lod.sse1[[i]])
  }
  lodMat <- cbind(lodMat, fit2s$lod.sse, fit3s$lod.sse, fit4s$lod.sse)
  lod.mse1 <- apply(lodMat,1,sum)/lodmiss 
  
  # lod bias #
  lod.sbias1 <- mclapply(1:n1, FUN = function(i){
    if(!anyNA(fit1s[[i]]$lod.sum.bias)) return(fit1s[[i]]$lod.sum.bias)
    else return(NULL)
  })
  lodMatb <- numeric()
  for(i in 1:n1){
    lodMatb <- cbind(lodMatb, lod.sbias1[[i]])
  }
  lodMatb <- cbind(lodMatb, fit2s$lod.sum.bias, fit3s$lod.sum.bias, fit4s$lod.sum.bias)
  lod.bias1 <- apply(lodMatb,1,sum)/lodmiss # sum the rows to get total sse, then divide by total mar data
  
  list1 = list(ham = NULL, mu.mse = NULL, mar.mse = mar.mse1, 
               lod.mse = lod.mse1, mar.bias = mar.bias1, lod.bias = lod.bias1)
  
  
}

