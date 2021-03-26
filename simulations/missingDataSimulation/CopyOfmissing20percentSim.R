# 20% percent missing simulation 
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

set.seed(22*simnum)
n <- 20
t.max <- 288
lodmis <- 0.10
marmis <- 0.10
sf = .1
tempTrend = TRUE
lodRemove = TRUE
marRemove = TRUE
dat1 <- simdatsimple(n = n, t.max = t.max, tempTrend = TRUE, lodRemove = TRUE, marRemove = TRUE,
                     lodmis = lodmis, marmis = marmis, sf = sf)

## look at data ##
# par(mfrow = c(2,2))
# plot(1:t.max, dat1$y.complete[[1]][,1], col = dat1$z.true[[1]], pch = 19)
# plot(1:t.max, dat1$y.complete[[1]][,2], col = dat1$z.true[[1]], pch = 19)
# plot(1:t.max, dat1$y.complete[[3]][,2], col = dat1$z.true[[3]], pch = 19)
# plot(1:t.max, dat1$y.complete[[3]][,3], col = dat1$z.true[[3]], pch = 19)

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

# fit models with and without covariates #
start <- Sys.time()
# DPMM
fitdpmm <- fitDPMM(niter = niter, nburn = nburn, y=dat1$y, ycomplete = dat1$y.complete,
                   priors = priors, K.start = NULL, z.true = dat1$z.true, lod = dat1$lod, 
                   mu.true = dat1$mu.true, SigmaPrior = "wish", algorithm = "MH", 
                   tau2 = .25, a.tune = 10, b.tune = 2, resK = TRUE, eta.star = 6, len.imp = len.imp)
# Cyclical MVN Rcpp
fit1cyc <- mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
                  priors=priors, K.start=K.start, z.true=z.true, lod=lod,
                  mu.true=mu.true, missing = missing, 
                  tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
                  resK = resK, eta.star = eta.star, len.imp = len.imp)
# No covariates MVN
fit1base <- fitMarkovNone(niter = niter, nburn = nburn, y=dat1$y, ycomplete = dat1$y.complete,
                          priors = priors, K.start = NULL, z.true = dat1$z.true, lod = dat1$lod, 
                          mu.true = dat1$mu.true, SigmaPrior = SigmaPrior, algorithm = algorithm, 
                          tau2 = tau2, a.tune = a.tune, b.tune = b.tune, 
                          resK = TRUE, eta.star = 3, len.imp = len.imp)
end <- Sys.time()
end-start

###############
### Results ###
###############

# check the MH ratio 

## coverage for imputations ##
# fit1cyc
marcov1 <- cbind(t(apply(fit1cyc$ymar, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                 as.vector(unlist(dat1$y.complete))[which(is.na(as.numeric(unlist(dat1$y))))])
colnames(marcov1) <- c("mean", "lwr", "upr", "true")
marCoverage1 <- mean((marcov1[,"true"] > marcov1[,"lwr"]) & (marcov1[,"true"] < marcov1[,"upr"]))

lodcov1 <- cbind(t(apply(fit1cyc$ylod, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                 as.vector(unlist(dat1$y.complete))[which((as.numeric(unlist(dat1$y))==-Inf))])
colnames(lodcov1) <- c("mean", "lwr", "upr", "true")
lodCoverage1 <- mean((lodcov1[,"true"] > lodcov1[,"lwr"]) & (lodcov1[,"true"] < lodcov1[,"upr"]))
# fit1base
marcov2 <- cbind(t(apply(fit1base$ymar, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                 as.vector(unlist(dat1$y.complete))[which(is.na(as.numeric(unlist(dat1$y))))])
colnames(marcov2) <- c("mean", "lwr", "upr", "true")
marCoverage2 <- mean((marcov2[,"true"] > marcov2[,"lwr"]) & (marcov2[,"true"] < marcov2[,"upr"]))

lodcov2 <- cbind(t(apply(fit1base$ylod, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                 as.vector(unlist(dat1$y.complete))[which((as.numeric(unlist(dat1$y))==-Inf))])
colnames(lodcov2) <- c("mean", "lwr", "upr", "true")
lodCoverage2 <- mean((lodcov2[,"true"] > lodcov2[,"lwr"]) & (lodcov2[,"true"] < lodcov2[,"upr"]))
# dpmm 
marcov3 <- cbind(t(apply(fitdpmm$ymar, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                 as.vector(unlist(dat1$y.complete))[which(is.na(as.numeric(unlist(dat1$y))))])
colnames(marcov3) <- c("mean", "lwr", "upr", "true")
marCoverage3 <- mean((marcov3[,"true"] > marcov3[,"lwr"]) & (marcov3[,"true"] < marcov3[,"upr"]))

lodcov3 <- cbind(t(apply(fitdpmm$ylod, 2, FUN = function(x) c(mean(x), quantile(x, .025), quantile(x, 0.975)))),
                 as.vector(unlist(dat1$y.complete))[which((as.numeric(unlist(dat1$y))==-Inf))])
colnames(lodcov3) <- c("mean", "lwr", "upr", "true")
lodCoverage3 <- mean((lodcov3[,"true"] > lodcov3[,"lwr"]) & (lodcov3[,"true"] < lodcov3[,"upr"]))

# posterior predictive mean MSE for imputations
truemar <- unlist(dat1$y.complete)[which(is.na(unlist(dat1$y)))]
truelod <- unlist(dat1$y.complete)[which(unlist(dat1$y)==-Inf)]
# fit1cyc
ymarmean1 <- apply(matrix(fit1cyc$ymar, nrow = len.imp), 2, mean)
ylodmean1 <- apply(matrix(fit1cyc$ylod, nrow = len.imp), 2, mean)
ppmarMean1 <- mean((ymarmean1 - truemar)^2)
pplodMean1 <- mean((ylodmean1 - truelod)^2)
# fit1base
ymarmean2 <- apply(matrix(fit1base$ymar, nrow = len.imp), 2, mean)
ylodmean2 <- apply(matrix(fit1base$ylod, nrow = len.imp), 2, mean)
ppmarMean2 <- mean((ymarmean2 - truemar)^2)
pplodMean2 <- mean((ylodmean2 - truelod)^2)
# fitdpmm
ymarmean3 <- apply(matrix(fitdpmm$ymar, nrow = len.imp), 2, mean)
ylodmean3 <- apply(matrix(fitdpmm$ylod, nrow = len.imp), 2, mean)
ppmarMean3 <- mean((ymarmean3 - truemar)^2)
pplodMean3 <- mean((ylodmean3 - truelod)^2)

# all 
hams <- cbind(fit1cyc$hamming, fit1base$hamming, fitdpmm$hamming)
mse <- cbind(fit1cyc$mu.mse, fit1base$mu.mse, fitdpmm$mu.mse)
mar.mse <- cbind(fit1cyc$mar.mse, fit1base$mar.mse, fitdpmm$mar.mse)
lod.mse <- cbind(fit1cyc$lod.mse, fit1base$lod.mse, fitdpmm$lod.mse)
mar.bias <- cbind(fit1cyc$mar.bias, fit1base$mar.bias, fitdpmm$mar.bias)
lod.bias <- cbind(fit1cyc$lod.bias, fit1base$lod.bias, fitdpmm$lod.bias)


# hamming for best cluster 
bc1 <- bestClusteriHMM(fit1cyc)
ham.bc1 <- hamdist(unlist(dat1$z.true), unlist(bc1))/(n*t.max) # proportion of misplaced states
bc2 <- bestClusteriHMM(fit1base)
ham.bc2 <- hamdist(unlist(dat1$z.true), unlist(bc2))/(n*t.max) # proportion of misplaced states
bc3 <- bestClusteriHMM(fitdpmm)
ham.bc3 <- hamdist(unlist(dat1$z.true), unlist(bc3))/(n*t.max) # proportion of misplaced states

# means 
ppmarMean <- c(ppmarMean1, ppmarMean2, ppmarMean3)
pplodMean <- c(pplodMean1, pplodMean2, pplodMean3)
ham.bc <- c(ham.bc1, ham.bc2, ham.bc3)
hamming <- c(apply(hams, 2, mean))
mu.mse <- c(apply(mse, 2, mean))
marMean <- c(apply(mar.mse, 2, mean))
lodMean <- c(apply(lod.mse, 2, mean))
marbiasMean <- c(apply(mar.bias, 2, mean))
lodbiasMean <- c(apply(lod.bias, 2, mean))
marCoverage <- c(marCoverage1, marCoverage2, marCoverage3)
lodCoverage <- c(lodCoverage1, lodCoverage2, lodCoverage3)

missingResults <- data.frame(ham.bc, hamming, mu.mse, ppmarMean, marMean, pplodMean, lodMean, marbiasMean, lodbiasMean, marCoverage, lodCoverage)
rownames(missingResults) <- c("full20", "nox20", "dpmm20")
colnames(missingResults) <- NULL

#cat results/miss20percent*.csv > combined_results/miss20percentResults.txt

list.save(dat1, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/miss20data", simnum, ".rds"))
write.table(missingResults, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/miss20percent", simnum, ".csv"), 
            row.names = TRUE, col.names = FALSE, sep = ",")


