## data analysis script for 50 samples 
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

####################
### Read in Data ###
####################

# Summit 
# yMatScaled <- as.matrix(read.table("/projects/lvheck@colostate.edu/markovPSBP/simulations/fccsDataMatrix_50samples.csv"),
#           row.names = NULL, col.names = NULL)

# Mac
yMatScaled <- as.matrix(read.table("~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataMatrix_50samples.csv"),
                        row.names = NULL, col.names = NULL)
colnames(yMatScaled) <- NULL
lod <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/limitOfDetection.rds")
lod <- as.numeric(lod); lod

############################
### Set up exposure data ###
############################

N <- 50
t.max <- 288
splits <- seq(1, nrow(yMatScaled)-287, length.out = N)
y <- lapply(1:length(splits), FUN = function(t) yMatScaled[splits[t]:(splits[t]+t.max-1), ])

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

niter = 50
nburn = 25
ycomplete = NULL
priors = list(nu = p+2, R = diag(p), m0 = 0, v0 = 1)
K.start = NULL
z.true = NULL
mu.true = NULL
missing = TRUE 
tau2 = .25 
a.tune = 10
b.tune = 1
resK = TRUE
eta.star = 3
len.imp = 4

#################################################
### FCCS data analysis without random effects ###
#################################################

rmlist = NULL

st1 = Sys.time()
fit_rcpp <- mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
       priors=priors, K.start=K.start, z.true=z.true, lod=lod,
       mu.true=mu.true, missing = missing, 
       tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
       resK = resK, eta.star = eta.star, len.imp = len.imp)
en1 = Sys.time()
en1-st1

#################################################
### FCCS data analysis WITH random effects ###
#################################################

rmlist = c(rep(1,5), rep(2,5), rep(3,6), rep(4,5), rep(5,5), rep(6,5),
           rep(7,9), rep(8, 5), rep(9,5))

st2 = Sys.time()
fit_rcpp <- mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
                   priors=priors, K.start=K.start, z.true=z.true, lod=lod,
                   mu.true=mu.true, missing = missing, 
                   tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
                   resK = resK, eta.star = eta.star, len.imp = len.imp)
en2 = Sys.time()
en2-st2


## R code: repeated measures, 50 people, 50 iterations took 1 hour 




#################### end script #########################


# save model fit 
fit.save = list(fit = fitMI1)
list.save(fit.save, file = "~/Documents/Lauren/Rpackages/markovPSBP/simulations/FCCSstudyscripts/fitSave_n50_3000.rds")
fit.save <- list.load(file = "~/Documents/Lauren/Rpackages/markovPSBP/simulations/FCCSstudyscripts/fitSave_n50_3000.rds")
fit2 <- fit.save$fit

############
### MISC ###
############

# bad mixing!
fitMI1$MH.arate; fitMI1$MH.lamrate
fit2_rm$MH.arate; fit2_rm$MH.lamrate

####################
### Save Results ###
####################

## best clustering for inference on states ##
zbest1 <- bestClusteriHMM(fitMI1)
## model averaged estimates of harmonic trend ##
beta.save1 <- modelAveBeta(fit = fitMI1, zbest1 = zbest1, X = X)
beta.mean1 <- beta.save1$betaXma
beta.upr1 <- beta.save1$betaXupr
beta.lwr1 <- beta.save1$betaXlwr
## model averaged estimates of mu 
mu.save <- modelAveMu(fit = fitMI1, zbest1 = zbest1)
mu.mean <- mu.save$mu_ma
mu.lwr <- mu.save$mu_lwr
mu.upr <- mu.save$mu_upr
## list of results

which(unlist(fitMI1$mismat) == 2)
fitMI1$ylod


# what else do I need to save? 
MIresults1 <- list(mse.bias = cbind(fitMI1$mar.mse, fitMI1$lod.mse, fitMI1$mar.bias, fitMI1$lod.bias),
                   zbest = zbest1, 
                   mu.mean = mu.mean, mu.lwr = mu.lwr, mu.upr = mu.upr,
                   beta.mean = beta.mean1, beta.upr = beta.upr1, beta.lwr = beta.lwr1,
                   lod.imputes = fitMI1$ylod, mar.imputes = fitMI1$ymar, ycomplete = fitMI1$ycomplete,
                   mismat = fitMI1$mismat, MH1 = fitMI1$MH.arate, MH2 = fitMI1$MH.lamrate)
# Mac
list.save(MIresults1, file = "~/Documents/Lauren/Rpackages/markovPSBP/simulations/FCCSstudyscripts/cyclicalFCCS_50samps.rds")
# Summit 
#list.save(MIresults1, file = "/projects/lvheck@colostate.edu/markovPSBP/simulations/results/cyclicalFCCS_50samps.rds")
## analyze with simulations/combined_results/ScriptsToAnalyzeResults/fccs_dataAnalysisResults_50samps.R

######################################
### Save Results Repeated Measures ###
######################################

## best clustering for inference on states ##
zbest2 <- bestClusteriHMM(fit2_rm)
## model averaged estimates of harmonic trend ##
beta.save2 <- modelAveBeta(fit = fit2_rm, zbest1 = zbest2, X = X)
beta.mean2 <- beta.save2$betaXma
beta.upr2 <- beta.save2$betaXupr
beta.lwr2 <- beta.save2$betaXlwr
## model averaged estimates of mu 
mu.save <- modelAveMu(fit = fit2_rm, zbest1 = zbest2)
mu.mean <- mu.save$mu_ma
mu.lwr <- mu.save$mu_lwr
mu.upr <- mu.save$mu_upr
## list of results
rm_50_results2 <- list(mse.bias = cbind(fit2_rm$mar.mse, fit2_rm$lod.mse, fit2_rm$mar.bias, fit2_rm$lod.bias),
                       zbest = zbest2, 
                       mu.mean = mu.mean, mu.lwr = mu.lwr, mu.upr = mu.upr,
                       lod.imputes = fit2_rm$ylod, mar.imputes = fit2_rm$ymar, ycomplete = fit2_rm$ycomplete,
                       mismat = fit2_rm$mismat, MH1 = fit2_rm$MH.arate, MH2 = fit2_rm$MH.lamrate)
# Mac
list.save(rm_50_results2, file = "~/Documents/Lauren/Rpackages/markovPSBP/simulations/FCCSstudyscripts/rm_FCCS_50samps.rds")

# Summit 
#list.save(MIresults2, file = "/projects/lvheck@colostate.edu/markovPSBP/simulations/results/rm_FCCS_50samps.rds")
## analyze with simulations/combined_results/ScriptsToAnalyzeResults/fccs_dataAnalysisResults_50samps.R

