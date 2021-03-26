### Complete Data simulation script to run on SUMMIT
### No shared states in design
### INDEP methods 
### Lauren Hoskovec

# rm(list = ls())
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

#############
### Setup ###
#############

## Summit
simnum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
source("/projects/lvheck@colostate.edu/markovPSBP/simulations/simFunctions.R")

## Mac 
# simnum <- ceiling(runif(1, 0, 1000));simnum
# source("simulations/simFunctions.R")

#####################
### Simulate Data ###
#####################

# complete data case #
set.seed(222*simnum)
n <- 20
t.max <- 288 # every 5 minutes for 24 hours 
lodmis <- 0
marmis <- 0
sf = 0.1

dat1 <- simdatsimple(n = n, t.max = t.max, tempTrend = FALSE, lodRemove = FALSE, marRemove = FALSE,
                     lodmis = lodmis, marmis = marmis, sf = sf)

########################
### Set up time data ###
########################

transT <- seq(1:t.max)/t.max*2*pi
X <- cbind(sin(transT), cos(transT), sin(2*transT), cos(2*transT))
X1 = list()
for(i in 1:n){
  X1[[i]] = X
}
X = X1 # X is list with a matrix for each i 
rmlist = NULL
q <- ncol(X[[1]])
p <- 3 

#################################
### Set Priors and Parameters ###
#################################

niter = 100
nburn = 25

y = dat1$y
rmlist = NULL
ycomplete = dat1$y.complete
priors = list(bj = rep(10, p))
K.start = NULL
z.true = dat1$z.true
lod = NULL
mu.true = dat1$mu.true
missing = FALSE
tau2 = NULL
a.tune = NULL
b.tune = NULL
resK = FALSE
eta.star = NULL
len.imp = NULL


##################
### Fit Models ###
##################

fitindep1 <- mclapply(1:n, FUN = function(i){
  mciHMM(niter=niter, nburn=nburn, y=y[[i]], ycomplete=ycomplete[[i]], X=list(X[[i]]),
         priors=priors, z.true=z.true[[i]],
         mu.true=mu.true, missing = missing)
})


fitindep1nox <- mclapply(1:n, FUN = function(i){
  fitMarkovNone(niter = niter, nburn = nburn, y=dat1$y[[i]], ycomplete = dat1$y.complete[[i]], 
                priors = priors, K.start = NULL, z.true = dat1$z.true[[i]], lod = dat1$lod, 
                mu.true = dat1$mu.true, SigmaPrior = SigmaPrior, algorithm = algorithm)
})

#################################
### Posterior Summary: Shared ###
#################################

# hamming distance for best cluster 
bc1 <- lapply(1:n, FUN = function(i) bestClusteriHMM(fitindep1[[i]]))
ham.bc1 <- sum(unlist(lapply(1:n, FUN = function(i) {
  hamdist(unlist(dat1$z.true[[i]]), unlist(bc1[[i]]))
})))/(n*t.max)  # proportion of misplaced states

bc2 <- lapply(1:n, FUN = function(i) bestClusteriHMM(fitindep1nox[[i]]))
ham.bc2 <- sum(unlist(lapply(1:n, FUN = function(i) {
  hamdist(unlist(dat1$z.true[[i]]), unlist(bc2[[i]]))
})))/(n*t.max)  # proportion of misplaced states

fitind <- unlistMethod(fitindep1)
fitindbase <- unlistMethod(fitindep1nox)

ham.bc <- c(ham.bc1, ham.bc2)
ham <- c(mean(fitind$ham), mean(fitindbase$ham))
mu.mse <- c(mean(fitind$mu.mse), mean(fitindbase$mu.mse))
df2 <- data.frame(ham.bc, ham, mu.mse)
rownames(df2) <- c("ind-cyclical-PSBP-iHMM", "ind-PSBP-iHMM")

#df2
# cat results/ccNoTrendIndep*.csv > combined_results/ccrNoTrendIndepResults.txt
write.table(df2, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/ccNoTrendIndep", simnum, ".csv"), 
            row.names = TRUE, col.names = FALSE, sep = ",")



