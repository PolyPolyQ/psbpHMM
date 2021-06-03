### Complete Data simulation script to run on SUMMIT
### Shared states in design 
### Shared states methods 
### Lauren Hoskovec

rm(list=ls())
gc()

#################
### Libraries ###
#################

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)

# compileAttributes()
# devtools::build()
# devtools::install()
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

# # ## Summit
# simnum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# source("/projects/lvheck@colostate.edu/markovPSBP/simulations/simFunctions.R")

## Mac 
simnum <- ceiling(runif(1, 0, 1000));simnum
source("simulations/simFunctions.R")

#####################
### Simulate Data ###
#####################

# complete data case #
set.seed(222*simnum)
n <- 5
t.max <- 288 # every 5 minutes for 24 hours 
lodmis <- 0
marmis <- 0
sf = 0.1
# shared states AND shared trends among time series
dat1 <- simdatsimple(n = n, t.max = t.max, tempTrend = TRUE, lodRemove = FALSE, marRemove = FALSE,
                     lodmis = lodmis, marmis = marmis, sf = sf)

######################
### Visualize Data ###
######################

par(mfrow = c(2,2))
plot(1:t.max, dat1$y.complete[[1]][,1], col= dat1$z.true[[1]], pch = 19)
plot(1:t.max, dat1$y.complete[[2]][,1], col= dat1$z.true[[2]], pch = 19)
plot(1:t.max, dat1$y.complete[[3]][,1], col= dat1$z.true[[3]], pch = 19)
plot(1:t.max, dat1$y.complete[[4]][,1], col= dat1$z.true[[4]], pch = 19)

### color lines person 1 ###
par(mfrow = c(2,1))
t1 <- min(which(dat1$z.true[[1]]==2))-1
t2 <- min(which(dat1$z.true[[1]]==3))-1
t3 <- min(which(dat1$z.true[[1]]==4))-1
t4 <- min(which(dat1$z.true[[1]]==5))-1
t5 <- min(which(dat1$z.true[[1]]==6))-1
t6 <- max(which(dat1$z.true[[1]]==6))
plot(1:t1, dat1$y.complete[[1]][1:t1,1], type = "l", pch = 19, xlim = c(0, 288),
     ylim = c(-3,3), xlab = "time", ylab = "exposure 1", axes = FALSE,
     main = "subject 1", lwd = 2)
axis(1, at=seq(0, 288, 48), labels=seq(0, 288, 48))
axis(2, at=seq(-3,3,1), labels = seq(-3,3,1), las = 1)
lines((t1+1):t2, dat1$y.complete[[1]][(t1+1):t2,1], type = "l", col = 2, lwd = 2)
lines((t2+1):t3, dat1$y.complete[[1]][(t2+1):t3,1], type = "l", col = 3, lwd = 2)
lines((t3+1):t4, dat1$y.complete[[1]][(t3+1):t4,1], type = "l", col = 4, lwd = 2)
lines((t4+1):t5, dat1$y.complete[[1]][(t4+1):t5,1], type = "l", col = 5, lwd = 2)
lines((t5+1):t6, dat1$y.complete[[1]][(t5+1):t6,1], type = "l", col = 6, lwd = 2)
lines((t6+1):288, dat1$y.complete[[1]][(t6+1):288,1], type = "l", col = 1, lwd = 2)

### color lines person 2 ###
ts1 <- min(which(dat1$z.true[[2]]==2))-1
ts2 <- min(which(dat1$z.true[[2]]==3))-1
ts3 <- min(which(dat1$z.true[[2]]==4))-1
ts4 <- min(which(dat1$z.true[[2]]==5))-1
ts5 <- min(which(dat1$z.true[[2]]==6))-1
ts6 <- max(which(dat1$z.true[[2]]==6))

plot(1:ts1, dat1$y.complete[[2]][1:ts1,1], type = "l", pch = 19, xlim = c(0, 288),
     ylim = c(-3,3), xlab = "time", ylab = "exposure 1", axes = FALSE,
     main = "subject 2", lwd = 2)
axis(1, at=seq(0, 288, 48), labels=seq(0, 288, 48))
axis(2, at=seq(-3,3,1), labels = seq(-3,3,1), las = 1)
lines((ts1+1):ts2, dat1$y.complete[[2]][(ts1+1):ts2,1], type = "l", col = 2, lwd = 2)
lines((ts2+1):ts3, dat1$y.complete[[2]][(ts2+1):ts3,1], type = "l", col = 3, lwd = 2)
lines((ts3+1):ts4, dat1$y.complete[[2]][(ts3+1):ts4,1], type = "l", col = 4, lwd = 2)
lines((ts4+1):ts5, dat1$y.complete[[2]][(ts4+1):ts5,1], type = "l", col = 5, lwd = 2)
lines((ts5+1):ts6, dat1$y.complete[[2]][(ts5+1):ts6,1], type = "l", col = 6, lwd = 2)
lines((ts6+1):288, dat1$y.complete[[2]][(ts6+1):288,1], type = "l", col = 1, lwd = 2)

# plot(1:t.max, dat1$y.complete[[2]][,1], type = "l", pch = 19)
# plot(1:t.max, dat1$y.complete[[3]][,1], type = "l", pch = 19)
# plot(1:t.max, dat1$y.complete[[4]][,1], type = "l", pch = 19)

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

niter = 1000
nburn = 500

y = dat1$y
rmlist = NULL
ycomplete = dat1$y.complete
priors = list(bj = rep(10, p))
K.start = 12
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

# fit model with harmonic trend
st = Sys.time()
fit1 <- mciHMM(niter=niter, nburn=nburn, y=y, ycomplete=ycomplete, X=X,
                           priors=priors, z.true=z.true, K.start = K.start,
                           mu.true=mu.true, missing = missing)
en = Sys.time()
en - st


X = matrix(X[[1]], ncol = q); X
head(X)
dim(X)
# fit model without harmonic trend 
st = Sys.time()
fit1nox <- miHMM(niter=niter, nburn=nburn, y=y, ycomplete=ycomplete,
        priors=priors, z.true=z.true,
        mu.true=mu.true, missing = missing)
en = Sys.time()
en - st

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
# cat results/completeShared.csv > combined_results/completeSharedResults.txt
write.table(df, file = paste0("/projects/lvheck@colostate.edu/psbpHMM/simulations/results/completeShared", simnum, ".csv"), 
            row.names = TRUE, col.names = FALSE, sep = ",")


