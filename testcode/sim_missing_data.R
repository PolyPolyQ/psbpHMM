rm(list = ls())
gc()

#################
### Libraries ###
#################

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
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

simnum <- sample(1:100, 1); simnum
source("testcode/simFunctions.R")
set.seed(22*simnum)

#####################
### Simulate Data ###
#####################

n = 20 # sampling days
t.max <- 288
lodmis <- 0.05
marmis <- 0.05
sf = 0.1 # go back to 0.1 
tempTrend = TRUE
lodRemove = TRUE
marRemove = TRUE
dat1 <- simdatsimple(n = n, t.max = t.max, tempTrend = TRUE, lodRemove = TRUE, marRemove = TRUE,
                     lodmis = lodmis, marmis = marmis, sf = sf)

transT <- seq(1:t.max)/t.max*2*pi
X <- cbind(sin(transT), cos(transT), sin(2*transT), cos(2*transT))
# X a list with a matrix for each i 
X1 = list()
for(i in 1:n){
  X1[[i]] = X
}
X = X1
q <- ncol(X[[1]])
p <- 3 

## repeated measures test 
# rmlist = c(rep(1,5), rep(2,5), rep(3,6), rep(4,5), rep(5,5), rep(6,5),
#             rep(7,9), rep(8, 5), rep(9,5))


#################################
### Set Priors and Parameters ###
#################################

niter = 50
nburn = 25
y = dat1$y
rmlist = NULL
ycomplete = dat1$y.complete
priors = list(nu = p+2, R = diag(p), m0 = 0, v0 = 1)
K.start = NULL
z.true = dat1$z.true
lod = dat1$lod
mu.true = dat1$mu.true
missing = TRUE 
tau2 = .25 
a.tune = 10
b.tune = 1
resK = TRUE
eta.star = 3
len.imp = 4


st1 = Sys.time()
mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
                   priors=priors, K.start=K.start, z.true=z.true, lod=lod,
                   mu.true=mu.true, missing = missing, 
                   tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
                   resK = resK, eta.star = eta.star, len.imp = len.imp)
en1 = Sys.time()
en1-st1

library(markovPSBP) # this brings in all my other functions yay! 
X1 = X[[1]]

st2 = Sys.time()
fitMarkovSame(niter=niter, nburn=nburn, y=y, ycomplete=ycomplete, X=X1,
       priors=priors, K.start=K.start, z.true=z.true, lod=lod,
       mu.true=mu.true, missing = missing, 
       tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
       resK = resK, eta.star = eta.star, len.imp = len.imp)
en2 = Sys.time()
en2-st2


