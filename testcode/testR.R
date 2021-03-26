
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

X.save = X
X = X[[1]]
X = X.save
t.max = 20

### this is slow: C++
alpha.jk <- lapply(1:(K), FUN = function(j){
  unlist(lapply(1:(K), FUN = function(k){
    updateAlphaJK(j=j, k=k, n=n, t.max=t.max, z=z, vinv.alpha=vinv.alpha,
                  sig2inv.alpha = sig2inv.alpha, w.z = w.z, X = X, beta.k = beta.k,
                  beta.sk = beta.sk, m.alpha = m.alpha, mu.alpha = priors$mu.alpha)
  }))
})

# for each t, w.z gives me a VECTOR based on the previous time point and values up to the current time point
# so each w.z[[t]] should be a VECTOR of length z_t

### this is slow: C++
w.z <- list()
for(i in 1:n){
  w.z[[i]] <- list()
  for(t in 1:t.max){
    w.z[[i]][[t]] <- sapply(1:z[[i]][t], FUN = function(l) updateW(t=t, l=l, i=i, alpha.0k=alpha.0k, X=X[[i]], 
                                                                   beta.k=beta.k, beta.sk = beta.sk[[i]], 
                                                                   alpha.jk=alpha.jk, z=z))
  }
}