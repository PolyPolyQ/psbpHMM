### practice c++ functions
library(Rcpp)
library(RcppArmadillo)
compileAttributes()
devtools::build()
devtools::install()
library(psbpHMM)
library(microbenchmark)

### testing changes to update pi ###
ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)

### rcpp functions
testpi_rm = updatePi_rm(beta = beta.k, beta_sk = beta.sk, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)
testpi_original = updatePi(beta = beta.k, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)

### r functions
oldpi_rm = mclapply(1:n, FUN = function(i){
  mclapply(1:t.max, FUN = function(t){
    uppi(t, Xi = X[[i]], beta0 = beta.k, beta.ik = beta.sk[[i]])
  })
})
oldpi_original = mclapply(1:t.max, FUN = uppi_o)


f1 = oldpi_rm[[5]][[2]]; head(f1)
f2 = testpi_rm[[5]][[2]]; head(f2)


length(testpi_rm)

g1 = unlist(oldpi_rm)
g2 = unlist(testpi_rm)
length(g2)


range(f1[,31]-f2[,31])

which(abs(unlist(testpi_rm)-unlist(oldpi_rm))>0.99)


length(unlist(oldpi_rm))


anyNA(unlist(testpi_rm))
anyNA(unlist(oldpi_rm))


range(unlist(testpi_rm)-unlist(oldpi_rm))
range(unlist(testpi_original)-unlist(oldpi_original))



xx=runif(100)
microbenchmark(R = sqrt(xx), Rcpp = xx^(0.5))
n=5

# with repeated measures update pi is faster 
microbenchmark(Rcpp = updatePi_rm(beta = beta.k, beta_sk = beta.sk, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max),
               R = mclapply(1:n, FUN = function(i){
                 mclapply(1:t.max, FUN = function(t){
                   uppi(t, X[[i]], beta.k, beta.sk[[i]])
                 })
               }), times = 10)



range(unlist(testold)-unlist(testpi))


fun1 <- function(Xi, beta0, beta.ik){
  first1 = sapply(1:K, FUN = function(k) pnorm(alpha.0k[k] + crossprod(Xi[1,], beta0[[k]]) + crossprod(Xi[1,],beta.ik[[k]])))
  second1 = sapply(1:(K-1), FUN = function(k) 1-pnorm(alpha.0k[k] + crossprod(Xi[1,],beta0[[k]]) + crossprod(Xi[1,],beta.ik[[k]])))
  prod1 = c(1, cumprod(second1))
  return(c(first1*prod1, 1 - sum(first1*prod1)))
}

# takes in t, Xi = X[[i]], beta.ik = beta.k[[i]]
fun2 <- function(t, Xi, beta0, beta.ik){
  t(sapply(1:(K), FUN = function(j){
    first = sapply(1:K, FUN = function(k) pnorm(alpha.jk[[j]][k] + crossprod(Xi[t,],beta0[[k]]) + crossprod(Xi[t,],beta.ik[[k]])))
    second = sapply(1:(K-1), FUN = function(k) 1-pnorm(alpha.jk[[j]][k] + crossprod(Xi[t,],beta0[[k]]) + crossprod(Xi[t,],beta.ik[[k]])))
    prod2 = c(1, cumprod(second))
    return(c(first*prod2, 1 - sum(first*prod2))) 
  }))
}

# takes in t, Xi = X[[i]], and beta.ik = beta.k[[i]]
uppi <- function(t, Xi, beta0, beta.ik){
  if(t == 1) return(fun1(Xi, beta0, beta.ik))
  else return(fun2(t, Xi, beta0, beta.ik))
}



Xind = X[[i]]
# go into K states, then sum to 1 for K+1 
fun1_o <- function(){
  first1 = sapply(1:K, FUN = function(k) pnorm(alpha.0k[k] + crossprod(Xind[1,],beta.k[[k]])))
  second1 = sapply(1:(K-1), FUN = function(k) 1-pnorm(alpha.0k[k] + crossprod(Xind[1,],beta.k[[k]])))
  prod1 = c(1, cumprod(second1))
  return(c(first1*prod1, 1 - sum(first1*prod1)))
}

# K only 
fun2_o <- function(t){
  t(sapply(1:(K), FUN = function(j){
    first = sapply(1:K, FUN = function(k) pnorm(alpha.jk[[j]][k] + crossprod(Xind[t,],beta.k[[k]])))
    second = sapply(1:(K-1), FUN = function(k) 1-pnorm(alpha.jk[[j]][k] + crossprod(Xind[t,],beta.k[[k]])))
    prod2 = c(1, cumprod(second))
    return(c(first*prod2, 1 - sum(first*prod2))) 
  }))
}

uppi_o <- function(t){
  if(t == 1) return(fun1_o())
  else return(fun2_o(t))
}




