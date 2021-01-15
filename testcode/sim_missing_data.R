#################
### Libraries ###
#################

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

library(gdata)
library(invgamma)
library(markovPSBP) # this brings in all my other functions yay! 
library(gtools)
library(mvtnorm)
library(matrixcalc)
library(tmvmixnorm)
library(parallel)
library(mvnfast)
library(rlist)

#############
### Setup ###
#############

simnum <- sample(1:100, 1); simnum
source("testcode/simFunctions.R")
set.seed(22*simnum)

#####################
### Simulate Data ###
#####################

n = 5 # sampling days
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

#################################
### Set Priors and Parameters ###
#################################

priors = list(nu = p+2, R = diag(p), m0 = 0, v0 = 1)

niter = 50
nburn = 25
len.imp = 4
K.start = NULL
missing = TRUE 

tau2 = .25 
a.tune = 10
b.tune = 1

y <- dat1$y
ycomplete <- dat1$y.complete
z.true <- dat1$z.true
lod <- dat1$lod
mu.true <- dat1$mu.true
resK = TRUE
eta.star = 3

#################
### Fit Model ###
#################



