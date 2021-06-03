library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

compileAttributes()
devtools::build()
devtools::install()
library(psbpHMM)

#################### test rcpp code #######################


for(k in 1:K){
  sigK = Sigma[[k]]
  muK = mu[[k]]
  yit = y[[i]][t,]
  d = p 
  temp3 = det(sigK)^(-0.5)
  
  (2*pi)^(-d/2) * temp3 * exp(-0.5 * (yit - muK) %*% solve(sigK) %*% t(yit - muK))
}

testupZ(stateList = state.list, y = y, mu = mu, Sigma = Sigma, logStuff = log.stuff, 
         nudf = nu.df, detRstar = detR.star, piz = pi.z, u = u, tmax = t.max, K = K, n = n, d = p)
sample(1:11, 1, 0, tprobs)

upZ(stateList = state.list, y = y, mu = mu, Sigma = Sigma, logStuff = log.stuff, 
         nudf = nu.df, detRstar = detR.star, piz = pi.z, u = u, tmax = t.max, K = K, n = n, d = p)

############################################################


detzAll <- mclapply(1:n, FUN = function(i){
  detZminus1(i = i, state.list.i = state.list[[i]], pi.z = pi.z[[i]], u.i = u[[i]], t.max = t.max)
})

################
### Sample Z ### 
################

cholSigma <- lapply(1:K, FUN = function(k) chol(Sigma[[k]]))
K <- length(unique(unlist(z))); K

ztest = list()
for(i in 1:50){
 ztest[[i]] = updateZold(i=i, y.i=y[[i]], mu=mu, cholSigma=cholSigma, detR.stari=detR.star[[i]],
             nu.df=nu.df, K=K, log.stuff=log.stuff, t.max=t.max, 
             detz = detzAll[[i]], state.list.i = state.list[[i]])  
}

anyNA(unlist(ztest))


y.i=y[[i]]
detR.stari=detR.star[[i]]
detz = detzAll[[i]]
state.list.i = state.list[[i]]


z <- mclapply(1:n, FUN = function(i) {
  updateZold(i=i, y.i=y[[i]], mu=mu, cholSigma=cholSigma, detR.stari=detR.star[[i]],
          nu.df=nu.df, K=K, log.stuff=log.stuff, t.max=t.max, 
          detz = detzAll[[i]], state.list.i = state.list[[i]])
})
