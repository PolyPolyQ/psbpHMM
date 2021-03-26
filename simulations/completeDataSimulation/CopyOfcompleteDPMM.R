### Complete Data simulation script to run on SUMMIT
### Both designs
### DPMM 
### Lauren Hoskovec

# rm(list=ls())
# gc()

#################
### Libraries ###
#################

library(gdata)
library(markovPSBP)
library(gtools)
library(mvtnorm)
library(matrixcalc)
library(tmvmixnorm)
library(parallel)
library(mvnfast)
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

dat1 <- simdatsimple(n = n, t.max = t.max, tempTrend = TRUE, lodRemove = FALSE, marRemove = FALSE,
                     lodmis = lodmis, marmis = marmis, sf = sf)
dat2 <- simdatsimple(n = n, t.max = t.max, tempTrend = FALSE, lodRemove = FALSE, marRemove = FALSE,
                     lodmis = lodmis, marmis = marmis, sf = sf)

p=3
priors = list(bj = rep(10, p))
K.start <- NULL
SigmaPrior = "non-informative" 
algorithm = "Gibbs"
H = "jointNIW"

niter = 5000
nburn = 2500


niter = 500
nburn = 250

##################
### Fit Models ###
##################

st = Sys.time()
# shared trends
fitdpmm1 <- fitDPMM(niter = niter, nburn = nburn, y=dat1$y, ycomplete = dat1$y.complete,
                   priors = priors, K.start = NULL, z.true = dat1$z.true, lod = dat1$lod, 
                   mu.true = dat1$mu.true, SigmaPrior = "wish", algorithm = "MH", 
                   tau2 = 0.3, a.tune = 10, b.tune = 2)

# distinct trends
fitdpmm2 <- fitDPMM(niter = niter, nburn = nburn, y=dat2$y, ycomplete = dat2$y.complete,
                   priors = priors, K.start = NULL, z.true = dat2$z.true, lod = dat2$lod, 
                   mu.true = dat2$mu.true, SigmaPrior = "wish", algorithm = "MH", 
                   tau2 = 0.3, a.tune = 10, b.tune = 2)
en = Sys.time()
en - st

#################################
### Posterior Summary: Shared ###
#################################

# hamming distance of best clustering 
st = Sys.time()
bc1 <- bestClusteriHMM(fitdpmm1)
ham.bc1 <- hamdist(unlist(dat1$z.true), unlist(bc1))/(n*t.max) # proportion of misplaced states
bc2 <- bestClusteriHMM(fitdpmm2)
ham.bc2 <- hamdist(unlist(dat2$z.true), unlist(bc2))/(n*t.max) # proportion of misplaced states
en = Sys.time()


ham.bc <- c(ham.bc1, ham.bc2)
ham <- c(mean(fitdpmm1$hamming), mean(fitdpmm2$hamming))
mu.mse <- c(mean(fitdpmm1$mu.mse), mean(fitdpmm2$mu.mse))
df <- data.frame(ham.bc, ham, mu.mse)
rownames(df) <- c("shared-multi-DPMM", "distinct-multi-DPMM")

#df
#cat results/completeDPMM*.csv > combined_results/ccDPMM.txt

write.table(df, file = paste0("/projects/lvheck@colostate.edu/markovPSBP/simulations/results/completeDPMM", simnum, ".csv"), 
            row.names = TRUE, col.names = FALSE, sep = ",")


#cat results/completeDPMM*.csv > combined_results/ccDPMM.txt
