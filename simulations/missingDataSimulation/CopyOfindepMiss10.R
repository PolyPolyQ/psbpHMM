# 10% percent missing simulation INDEP
# rm(list=ls())
# gc()

#################
### Libraries ###
#################

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

library(markovPSBP)
library(psbpHMM)

library(gdata)
library(invgamma)
library(gtools)
library(mvtnorm)
library(matrixcalc)
library(tmvmixnorm)
library(parallel)
library(mvnfast)
library(rlist)
library(truncnorm)

## Summit
simnum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
source("/projects/lvheck@colostate.edu/markovPSBP/simulations/simFunctions.R")

## Mac 
# simnum <- ceiling(runif(1, 0, 1000));simnum
# source("simulations/simFunctions.R")

#############
### Setup ###
#############

# setup
set.seed(22*simnum)
n <- 20
t.max <- 288
lodmis <- 0.05
marmis <- 0.05
sf = 0.1
dat1 <- simdatsimple(n = n, t.max = t.max, tempTrend = TRUE, lodRemove = TRUE, marRemove = TRUE,
                     lodmis = lodmis, marmis = marmis, sf = sf)

# ## look at data ##
# par(mfrow = c(2,2))
# plot(1:t.max, dat1$y.complete[[1]][,1], col = dat1$z.true[[1]], pch = 19)
# plot(1:t.max, dat1$y.complete[[1]][,2], col = dat1$z.true[[1]], pch = 19)
# plot(1:t.max, dat1$y.complete[[3]][,2], col = dat1$z.true[[3]], pch = 19)
# plot(1:t.max, dat1$y.complete[[3]][,3], col = dat1$z.true[[3]], pch = 19)


lodmiss <- sum(unlist(lapply(1:n, FUN = function(i){
  lmiss <- length(which(dat1$y[[i]]==-Inf))
  return(lmiss)
})))

marmiss <- sum(unlist(lapply(1:n, FUN = function(i){
  mmiss <- length(which(is.na(dat1$y[[i]])))
  return(mmiss)
})))

########################
### Set up time data ###
########################

transT <- seq(1:t.max)/t.max*2*pi
X <- cbind(sin(transT), cos(transT), sin(2*transT), cos(2*transT))
# X a list with a matrix for each i 
X1 = list()
for(i in 1:n){
  X1[[i]] = X
}
X = X1
rmlist = NULL
q <- ncol(X[[1]])
p <- 3 

#################################
### Set Priors and Parameters ###
#################################

niter = 100
nburn = 25
len.imp = 10

priors = list(nu = p+2, R = diag(p), m0 = 0, v0 = 1)
K.start = NULL
missing = TRUE 
tau2 = .25 
a.tune = 10
b.tune = 1
resK = TRUE
eta.star = 3
y = dat1$y
ycomplete = dat1$y.complete
z.true = dat1$z.true
lod = dat1$lod
mu.true = dat1$mu.true
rmlist = NULL


######################
### Fit Full Model ###
######################

# see if you need resK = FALSE for these # 
###########

start <- Sys.time()
fit1 <- mclapply(1:n, FUN = function(i) {
  mciHMM(niter=niter, nburn=nburn, y=y[[i]], rmlist=rmlist, ycomplete=ycomplete[[i]], X=list(X[[i]]),
         priors=priors, K.start=K.start, z.true=z.true[[i]], lod=lod,
         mu.true=mu.true, missing = missing, 
         tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
         resK = resK, eta.star = eta.star, len.imp = len.imp)
})

fit1nox <- mclapply(1:n, FUN = function(i) {
  fitMarkovNone(niter=niter, nburn=nburn, y=dat1$y[[i]], ycomplete=dat1$y.complete[[i]],
                priors=priors, K.start=K.start, z.true=dat1$z.true[[i]], lod=dat1$lod,
                mu.true=dat1$mu.true, SigmaPrior = SigmaPrior, 
                algorithm = algorithm, tau2 = tau2, a.tune = a.tune, b.tune = b.tune, 
                resK = TRUE, eta.star = 2, len.imp = len.imp)
})
end <- Sys.time()
end-start

###############
### Results ###
###############

### get coverage ### 

# cyc 
marCovList <- mclapply(1:n, FUN = function(i){
  marcov1 <- cbind(t(apply(fit1[[i]]$ymar, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
        as.vector(unlist(dat1$y.complete[[i]]))[which(is.na(as.numeric(unlist(dat1$y[[i]]))))])
  colnames(marcov1) <- c("mean", "lwr", "upr", "true")
  return((marcov1[,"true"] > marcov1[,"lwr"]) & (marcov1[,"true"] < marcov1[,"upr"]))
})
marcov1 <- mean(unlist(marCovList))

lodCovList <- mclapply(1:n, FUN = function(i){
  lodcov1 <- cbind(t(apply(fit1[[i]]$ylod, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                   as.vector(unlist(dat1$y.complete[[i]]))[which((as.numeric(unlist(dat1$y[[i]])))==-Inf)])
  colnames(lodcov1) <- c("mean", "lwr", "upr", "true")
  return((lodcov1[,"true"] > lodcov1[,"lwr"]) & (lodcov1[,"true"] < lodcov1[,"upr"]))
})
lodcov1 <- mean(unlist(lodCovList))

# base 
marCovList2 <- mclapply(1:n, FUN = function(i){
  marcov1 <- cbind(t(apply(fit1nox[[i]]$ymar, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                   as.vector(unlist(dat1$y.complete[[i]]))[which(is.na(as.numeric(unlist(dat1$y[[i]]))))])
  colnames(marcov1) <- c("mean", "lwr", "upr", "true")
  return((marcov1[,"true"] > marcov1[,"lwr"]) & (marcov1[,"true"] < marcov1[,"upr"]))
})
marcov2 <- mean(unlist(marCovList2))

lodCovList2 <- mclapply(1:n, FUN = function(i){
  lodcov1 <- cbind(t(apply(fit1nox[[i]]$ylod, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                   as.vector(unlist(dat1$y.complete[[i]]))[which((as.numeric(unlist(dat1$y[[i]])))==-Inf)])
  colnames(lodcov1) <- c("mean", "lwr", "upr", "true")
  return((lodcov1[,"true"] > lodcov1[,"lwr"]) & (lodcov1[,"true"] < lodcov1[,"upr"]))
})
lodcov2 <- mean(unlist(lodCovList2))

## other stuff 
fit1cyc <- unlistMethod(fit1, marmiss, lodmiss, y = dat1$y, ycomplete = dat1$y.complete)
fit1base <- unlistMethod(fit1nox, marmiss, lodmiss, y = dat1$y, ycomplete = dat1$y.complete)

## Get the MSE for posterior predictive mean  
ppMarMean1 <- fit1cyc$ppMarMean
ppLodMean1 <- fit1cyc$ppLodMean

ppMarMean2 <- fit1base$ppMarMean
ppLodMean2 <- fit1base$ppLodMean

# hamming distance for best cluster 
bc1 <- lapply(1:n, FUN = function(i) bestClusteriHMM(fit1[[i]]))
ham.bc1 <- sum(unlist(lapply(1:n, FUN = function(i) {
  hamdist(unlist(dat1$z.true[[i]]), unlist(bc1[[i]]))
})))/(n*t.max)  # proportion of misplaced states

bc2 <- lapply(1:n, FUN = function(i) bestClusteriHMM(fit1nox[[i]]))
ham.bc2 <- sum(unlist(lapply(1:n, FUN = function(i) {
  hamdist(unlist(dat1$z.true[[i]]), unlist(bc2[[i]]))
})))/(n*t.max)  # proportion of misplaced states


hams <- cbind(fit1cyc$ham, fit1base$ham)
mse <- cbind(fit1cyc$mu.mse, fit1base$mu.mse)
mar <- cbind(fit1cyc$mar.mse, fit1base$mar.mse)
lod <- cbind(fit1cyc$lod.mse, fit1base$lod.mse)
mar.bias <- cbind(fit1cyc$mar.bias, fit1base$mar.bias)
lod.bias <- cbind(fit1cyc$lod.bias, fit1base$lod.bias)

ham.bc <- c(ham.bc1, ham.bc2)
hamming <- apply(hams, 2, mean)
mu.mse <- apply(mse, 2, mean)
ppMarMean <- c(ppMarMean1, ppMarMean2)
marMean <- apply(mar, 2, mean)
ppLodMean <- c(ppLodMean1, ppLodMean2)
lodMean <- apply(lod, 2, mean)
marbiasMean <- apply(mar.bias, 2, mean)
lodbiasMean <- apply(lod.bias, 2, mean)
marCoverage <- c(marcov1, marcov2)
lodCoverage <- c(lodcov1, lodcov2)

missingResults <- data.frame(ham.bc, hamming, mu.mse, ppMarMean, marMean, ppLodMean, lodMean, marbiasMean, lodbiasMean, marCoverage, lodCoverage)
rownames(missingResults) <- c("full10", "nox10")
colnames(missingResults) <- NULL
mri <- missingResults; mri

#cat results/ind10miss*.csv > combined_results/miss10IndResults.txt

list.save(dat1, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/indmiss10data", simnum, ".rds"))
write.table(missingResults, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/ind10miss", simnum, ".csv"), 
            row.names = TRUE, col.names = FALSE, sep = ",")



