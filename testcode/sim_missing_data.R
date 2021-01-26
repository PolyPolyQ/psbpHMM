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
#library(markovPSBP) # this brings in all my other functions yay! 
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


## repeated measures test 
rmlist = c(rep(1,5), rep(2,5), rep(3,6), rep(4,5), rep(5,5), rep(6,5),
            rep(7,9), rep(8, 5), rep(9,5))
rmlist = NULL

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

fitMarkovRM <- function(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
                        priors=priors, K.start=K.start, z.true=z.true, lod=lod,
                        mu.true=mu.true, missing = missing, 
                        tau2 = tau2, a.tune = a.tunee, b.tune = b.tune,
                        resK = resK, eta.star = eta.star, len.imp = len.imp)




## need to bring the data and FCCS repeaated measures script over here and 
## run the 50 sample data analysis 
## then try to make some functions faster so you can do some sensitivity analyses
## also need to redo the validation study simulations because you found one bug
## in the repeated measures function




