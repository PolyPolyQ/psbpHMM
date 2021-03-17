
rm(list=ls())
gc()

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

# getwd()
# sourceCpp("src/tester.cpp")
compileAttributes()
devtools::build()
devtools::install()
library(psbpHMM)

createMat(4,5)
dim(Z_mat)

t1 <- c(1,1,2,2,2,1,1,1)
t2 <- c(1,1,1,1,1,1,1,1)
t3 <- c(2,2,2,2,2,1,1,1)

Ztest = cbind(t1, t2, t3); Ztest
createMat(n=3, niter = 8, Zmat = Ztest)

vectorMean(t1, t2)

# create similarity matrix S for each clustering in Z_keep
# 1 in the i,j location if i and j were clustered together at that iteration 
# calculate squared distance from S to P
LSdist <- rep(NA, niter)
for (s in 1:niter){
  print(s)
  S <- matrix(0, n, n)
  for (k in 1:K){
    S[which(Z_mat[s,] == k), which(Z_mat[s,] == k)] <- 1
  }
  LSdist[s] <- sum( (S - Prob)^2 )
}





a <- c(1,1,2,3,4,4,4)
b <- c(1,1,3,2,4,4,5)

vectorMean(a,b)
