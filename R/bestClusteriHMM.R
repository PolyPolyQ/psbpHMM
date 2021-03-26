#' Identify best clustering of states in PSBP-iHMM model 
#'
#' @param ihmm object of type "ihmm"
#'
#' @return best clustering of z
#' @export
#'
#' 


## this code assumes each time series has equal length ## 

bestClusteriHMM <- function(ihmm){
  
  # unlist Z for each iteration
  # then relist Z for each person 
  
  Z_keep <- ihmm$z.save
  K <- max(unlist(ihmm$K.save))
  
  Z_newlist <- list()
  for(s in 1:length(Z_keep)){
    Z_newlist[[s]] <- unlist(Z_keep[[s]])
  }
  
  n <- length(Z_newlist[[1]]) # how many time points to cluster 
  niter <- length(Z_newlist)
  Z_mat <- matrix(unlist(Z_newlist), ncol = n, nrow = niter, byrow = TRUE)
  
  #Z_mat[1,] == Z_newlist[[1]]

  # calculate probability matrix P (symmetric)
  # Rcpp
  Prob <- createMat(n=n, niter = niter, Zmat = Z_mat)
  
  # Prob <- matrix(NA, n, n)
  # for (i in 1:n){
  #   print(i)
  #   for (j in 1:i){
  #     Prob[i, j] <- Prob[j, i] <- mean(Z_mat[,i]==Z_mat[,j])
  #   }
  # }
  
  # create similarity matrix S for each clustering in Z_keep
  # 1 in the i,j location if i and j were clustered together at that iteration 
  # calculate squared distance from S to P
  lstest = sapply(1:niter, FUN = function(s){
    S <- matrix(0, n, n)
    for (k in 1:K){
      S[which(Z_mat[s,] == k), which(Z_mat[s,] == k)] <- 1
    }
    return(sum( (S - Prob)^2 ))
  })
  # select Z_keep that corresponds to the S that minimizes LS distance to P
  best <- which.min(lstest)
  
  Zbest <- Z_mat[best,]
  # relist by person 
  
  numSub <- length(Z_keep[[1]]); numSub
  t.max <- n/numSub; t.max
  
  splits <- seq(1,numSub*t.max, t.max)
  z.list <- lapply(1:length(splits), FUN = function(i) Zbest[splits[i]:(splits[i]+t.max-1)])
  # z.list is a list of best clusters for each person 

  return(z.list)
  
  
}
