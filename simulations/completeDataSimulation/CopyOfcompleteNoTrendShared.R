### Complete Data simulation script to run on SUMMIT
### No shared states in design
### Shared states methods 
### Lauren Hoskovec

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

fit1 <- mciHMM(niter=niter, nburn=nburn, y=y, ycomplete=ycomplete, X=X,
               priors=priors, z.true=z.true,
               mu.true=mu.true, missing = missing)

fit1nox <- fitMarkovNone(niter = niter, nburn = nburn, y=dat1$y, ycomplete = dat1$y.complete,
                         priors = priors, K.start = NULL, z.true = dat1$z.true, lod = dat1$lod, 
                         mu.true = dat1$mu.true, SigmaPrior = SigmaPrior, algorithm = algorithm)

#################################
### Posterior Summary: Shared ###
#################################

# hamming distance of best clustering 
bc1 <- bestClusteriHMM(fit1)
ham.bc1 <- hamdist(unlist(dat1$z.true), unlist(bc1))/(n*t.max) # proportion of misplaced states
bc2 <- bestClusteriHMM(fit1nox)
ham.bc2 <- hamdist(unlist(dat1$z.true), unlist(bc2))/(n*t.max) # proportion of misplaced states

ham.bc <- c(ham.bc1, ham.bc2)
ham <- c(mean(fit1$hamming), mean(fit1nox$hamming))
mu.mse <- c(mean(fit1$mu.mse), mean(fit1nox$mu.mse))
df <- data.frame(ham.bc, ham, mu.mse)
rownames(df) <- c("cyclical-PSBP-iHMM", "PSBP-iHMM")

#df
# cat results/ccNoTrendShared*.csv > combined_results/ccrNoTrendSharedResults.txt
write.table(df, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/ccNoTrendShared", simnum, ".csv"),
            row.names = TRUE, col.names = FALSE, sep = ",")




