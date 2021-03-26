## data analysis script for 50 samples 
rm(list = ls())
gc()

#################
### Libraries ###
#################

library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)


library(mvnfast)
#compileAttributes()
#devtools::build()
#devtools::install()
library(psbpHMM)

library(gdata)
library(invgamma)
library(gtools)
library(mvtnorm)
library(matrixcalc)
library(tmvmixnorm)
library(parallel)

library(rlist)
library(truncnorm)
library(ggplot2)

####################
### Read in Data ###
####################

# Summit 
# yMatScaled <- as.matrix(read.table("/projects/lvheck@colostate.edu/markovPSBP/simulations/fccsDataMatrix_50samples.csv"),
#           row.names = NULL, col.names = NULL)

# Mac
yMatScaled <- as.matrix(read.table("~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fccsDataMatrix_50samples.csv"),
                        row.names = NULL, col.names = NULL)
colnames(yMatScaled) <- NULL
lod <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/limitOfDetection.rds")
lod <- as.numeric(lod); lod
menv <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/micEnv_50samples.rds")


############################
### Set up exposure data ###
############################

n <- 50
t.max <- 288
splits <- seq(1, nrow(yMatScaled)-287, length.out = n)
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
q <- ncol(X[[1]])
p <- 3 

rmlist = c(rep(1,5), rep(2,5), rep(3,6), rep(4,5), rep(5,5), rep(6,5),
           rep(7,9), rep(8, 5), rep(9,5))

#################################
### Set Priors and Parameters ###
#################################

niter = 200
nburn = 100
len.imp = 50

ycomplete = NULL
priors = list(nu = p+2, R = diag(p), m0 = 0, v0 = 1)
K.start = NULL
z.true = NULL
mu.true = NULL
missing = TRUE 
tau2 = .5 
a.tune = 1
b.tune = 1
resK = TRUE
eta.star = 5

# tuning parameters 
# tau2 = 0.25; a.tune = 10; b.tune = 1; eta.star = 5



#################################################
### FCCS data analysis without random effects ###
#################################################

# rmlist = NULL
# st1 = Sys.time()
# fit_rcpp <- mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
#        priors=priors, K.start=K.start, z.true=z.true, lod=lod,
#        mu.true=mu.true, missing = missing, 
#        tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
#        resK = resK, eta.star = eta.star, len.imp = len.imp)
# en1 = Sys.time()
# en1-st1
# 
# fit_rcpp_save = list(fit1 = fit_rcpp) 
# fit_rcpp$MH.arate; fit_rcpp$MH.lamrate
# zbest = bestClusteriHMM(fit_rcpp); zbest[[1]]
# beta.save <- modelAveBeta(fit = fit_rcpp, zbest1 = zbest, X = X[[1]]) # because all X the same here 
# mu.save <- modelAveMu(fit = fit_rcpp, zbest1 = zbest)
# 
# fit_rcpp_save = list(fit1 = fit_rcpp, zbest = zbest, mu.modelave = mu.save, beta.modelave = beta.save) 
# list.save(fit_rcpp_save, file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_save.rds")
# 
# 
# ## load 
# fit_rcpp_results <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_save.rds")
# fit_rcpp <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_save.rds")$fit1
# 
# beta.save = fit_rcpp_results$beta.modelave
# mu.save = fit_rcpp_results$mu.modelave
# zbest = fit_rcpp_results$zbest

# time for 100 iterations: 10 minutes

#################################################
### FCCS data analysis WITH random effects ###
#################################################

rmlist = c(rep(1,5), rep(2,5), rep(3,6), rep(4,5), rep(5,5), rep(6,5),
           rep(7,9), rep(8, 5), rep(9,5))

st2 = Sys.time()
fit_rcpp_rm <- mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
                   priors=priors, K.start=K.start, z.true=z.true, lod=lod,
                   mu.true=mu.true, missing = missing, 
                   tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
                   resK = resK, eta.star = eta.star, len.imp = len.imp)
en2 = Sys.time()
en2-st2

fit_rcpp_rm$MH.arate; fit_rcpp_rm$MH.lamrate
plot(1:(niter-nburn), fit_rcpp_rm$K.save, pch =19)



fit_rcpp_rm_save = list(fit1 = fit_rcpp_rm) 
zbest_rm = bestClusteriHMM(fit_rcpp_rm); zbest_rm[[1]]
beta.save_rm <- modelAveBeta(fit = fit_rcpp_rm, zbest1 = zbest_rm, X = X[[1]]) # because all X the same here 
mu.save_rm <- modelAveMu(fit = fit_rcpp_rm, zbest1 = zbest_rm)

fit_rcpp_rm_save = list(fit1 = fit_rcpp_rm, zbest = zbest_rm, mu.modelave = mu.save_rm, beta.modelave = beta.save_rm) 
list.save(fit_rcpp_rm_save, file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_rm_save.rds")

## load 
fit_rcpp_rm_results <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_rm_save.rds")
fit_rcpp_rm <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_rm_save.rds")$fit1

beta.save_rm = fit_rcpp_rm_results$beta.modelave
mu.save_rm = fit_rcpp_rm_results$mu.modelave
zbest_rm = fit_rcpp_rm_results$zbest

# time for 5000 iterations: 15 hours 

###################
### imputations ###
###################

fit_rcpp$ymar
t(fit_rcpp$ylod) 
dim(fit_rcpp$ymar)
dim(fit_rcpp$ylod)
unlist(fit_rcpp$mismat) # 1 is MAR, 2 is LOD 

#########################
### Check convergence ###
#########################

# some meausre of convergence  

par(mfrow = c(1,2))
plot(1:(niter-nburn), fit_rcpp_rm$K.save, main = "mc-iHMM-rm", 
     ylab = "Number of states", xlab = "iteration")
plot(1:(niter-nburn), fit_rcpp$K.save,  main = "mc-iHMM", 
     ylab = "Number of states", xlab = "iteration")




#########################
### Figures and Story ###
#########################

# ggplot of two sampling days for same person 
# plot of model averaged exposure means 
# stacked bar chart of microenvironments 
# ??? 

#####################################
### Model Averaged Exposure Means ###
#####################################

#############
## mc-iHMM ##
#############

mmean <- matrix(unlist(mu.save$mu_ma), ncol = 3, byrow = TRUE)
mlwr <- matrix(unlist(mu.save$mu_lwr), ncol = 3, byrow = TRUE)
mupr <- matrix(unlist(mu.save$mu_upr), ncol = 3, byrow = TRUE)

K = length(unique(unlist(zbest)));K

cols = c("blue", "red", "green", "purple", "orange", "black", "gray", "darkblue",
         "yellow", "pink", "darkgreen", "lightgray", "lightgreen", "darkred")

col = sample(cols, K, replace = TRUE)

mdat1 = data.frame(x = (1:K)-0.2, y = mmean[,1], lwr = mlwr[,1], upr = mupr[,1], col = col)
mdat2 = data.frame(x = (1:K), y = mmean[,2], lwr = mlwr[,2], upr = mupr[,2], col = col)
mdat3 = data.frame(x = (1:K)+0.2, y = mmean[,3], lwr = mlwr[,3], upr = mupr[,3], col = col)

mplot1 <- ggplot() + ggtitle("mc-iHMM") + theme(legend.position = "none") + 
  scale_x_discrete(name = "Hidden States", limits = factor(1:K), breaks = 1:K) + 
  scale_y_continuous(name = "Model Averaged Exposure Means") + 
  geom_point(data = mdat1, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 15) + 
  geom_errorbar(data = mdat1, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  geom_point(data = mdat2, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 16) + 
  geom_errorbar(data = mdat2, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  geom_point(data = mdat3, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 17) + 
  geom_errorbar(data = mdat3, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(size = 15, angle = 0, hjust = 1), 
        axis.text.y = element_text(size = 15, angle = 0, hjust = 1), 
        plot.title = element_text(size = 15)) 

mplot1

################
## mc-iHMM-rm ##
################

mmean <- matrix(unlist(mu.save_rm$mu_ma), ncol = 3, byrow = TRUE)
mlwr <- matrix(unlist(mu.save_rm$mu_lwr), ncol = 3, byrow = TRUE)
mupr <- matrix(unlist(mu.save_rm$mu_upr), ncol = 3, byrow = TRUE)

K = length(unique(unlist(zbest_rm)));K

cols = c("blue", "red", "green", "purple", "orange", "black", "gray", "darkblue",
         "yellow", "pink", "darkgreen", "lightgray", "lightgreen", "darkred")

col = sample(cols, K, replace = TRUE)

mdat1 = data.frame(x = (1:K)-0.2, y = mmean[,1], lwr = mlwr[,1], upr = mupr[,1], col = col)
mdat2 = data.frame(x = (1:K), y = mmean[,2], lwr = mlwr[,2], upr = mupr[,2], col = col)
mdat3 = data.frame(x = (1:K)+0.2, y = mmean[,3], lwr = mlwr[,3], upr = mupr[,3], col = col)

mplot2 <- ggplot() + ggtitle("mc-iHMM-rm") + theme(legend.position = "none") + 
  scale_x_discrete(name = "Hidden States", limits = factor(1:K), breaks = 1:K) + 
  scale_y_continuous(name = "Model Averaged Exposure Means") + 
  geom_point(data = mdat1, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 15) + 
  geom_errorbar(data = mdat1, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  geom_point(data = mdat2, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 16) + 
  geom_errorbar(data = mdat2, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  geom_point(data = mdat3, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 17) + 
  geom_errorbar(data = mdat3, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(size = 15, angle = 0, hjust = 1), 
        axis.text.y = element_text(size = 15, angle = 0, hjust = 1), 
        plot.title = element_text(size = 15)) 

mplot2

#########################
### Stacked bar chart ###
#########################

#############
## mc-iHMM ##
#############

K <- length(unique(unlist(zbest))); K
zall <- unlist(zbest) # all hidden states
mall <- unlist(menv) # all microenvironments
which.me <- list()
for(k in 1:K){
  which.me[[k]] <- mall[which(zall==k)] # ME for k = 1
}
hid <- factor(as.numeric(sapply(1:K, FUN = function(k) rep(k, 5)))); hid
ms <- factor(rep(1:5, K))
nums <- as.numeric(sapply(1:K, FUN = function(k){
  c(length(which(which.me[[k]]=="other")),
    length(which(which.me[[k]]=="eateries")),
    length(which(which.me[[k]]=="transit")),
    length(which(which.me[[k]]=="work")),
    length(which(which.me[[k]]=="home")))
}))
dat <- data.frame(hid, ms, nums)
dat
hidar <- sort(table(unlist(zbest)), decreasing = TRUE)
hidar <- as.data.frame(hidar)
hmode <- hidar[which(hidar$Freq>112),]
hmode$Var1
datr <- dat[which(dat$hid%in%hmode$Var1),]
datr
# all hidden states 
bchart1 <- ggplot(data=dat, aes(fill=ms, y=nums, x=hid)) + 
  ggtitle("mc-iHMM") +
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 10)) + 
  geom_bar(position="stack", stat="identity", alpha = 0.9) +
  scale_fill_manual(name = "Microenvironment", values = c("#fff700", "#9582ff", "#ff0004", "#02e0ad", "#00a6ff"),
                    labels = c("other", "eateries", "transit", "work", "home")) +
  xlab("Hidden State") + ylab("Count of Time Points") 
bchart1

################
## mc-iHMM-rm ##
################

K <- length(unique(unlist(zbest_rm))); K
zall <- unlist(zbest_rm) # all hidden states
mall <- unlist(menv) # all microenvironments
which.me <- list()
for(k in 1:K){
  which.me[[k]] <- mall[which(zall==k)] # ME for k = 1
}
hid <- factor(as.numeric(sapply(1:K, FUN = function(k) rep(k, 5)))); hid
ms <- factor(rep(1:5, K))
nums <- as.numeric(sapply(1:K, FUN = function(k){
  c(length(which(which.me[[k]]=="other")),
    length(which(which.me[[k]]=="eateries")),
    length(which(which.me[[k]]=="transit")),
    length(which(which.me[[k]]=="work")),
    length(which(which.me[[k]]=="home")))
}))
dat <- data.frame(hid, ms, nums)
dat
hidar <- sort(table(unlist(zbest_rm)), decreasing = TRUE)
hidar <- as.data.frame(hidar)
hmode <- hidar[which(hidar$Freq>112),]
hmode$Var1
datr <- dat[which(dat$hid%in%hmode$Var1),]
datr
# all hidden states 
bchart2 <- ggplot(data=dat, aes(fill=ms, y=nums, x=hid)) + 
  ggtitle("mc-iHMM-rm") +
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 10)) + 
  geom_bar(position="stack", stat="identity", alpha = 0.9) +
  scale_fill_manual(name = "Microenvironment", values = c("#fff700", "#9582ff", "#ff0004", "#02e0ad", "#00a6ff"),
                    labels = c("other", "eateries", "transit", "work", "home")) +
  xlab("Hidden State") + ylab("Count of Time Points") 
bchart2

#############
### Table ###
#############

library(xtable)

#############
## mc-iHMM ##
#############

K <- length(unique(unlist(zbest))); K
tab1 <- matrix(NA, nrow = K, ncol = 3)
for(k in 1:K){
  # how many sampling days 
  sdays <- sapply(1:n, FUN = function(i){
    length(which(zbest[[i]]==k))
  })
  sdaysNum <- length(which(sdays>0))
  # how many unique people
  ppl <- data.frame(rmlist, sdays)
  jTotal <- length(unique(rmlist))
  sppl <- sum(sapply(1:jTotal, FUN = function(j){
    as.numeric(sum(ppl$sdays[which(ppl$rmlist==j)])>0)
  }))
  # total number of time points
  timeTotal <- length(which(unlist(zbest)==k))
  
  tab1[k,] <- c(sdaysNum, sppl, timeTotal)
}
colnames(tab1) <- c("sampling days", "unique people", "time points")

xtable(tab1, digits = 0)


#############
## mc-iHMM-rm ##
#############

K <- length(unique(unlist(zbest_rm))); K
tab2 <- matrix(NA, nrow = K, ncol = 3)
for(k in 1:K){
  # how many sampling days 
  sdays <- sapply(1:n, FUN = function(i){
    length(which(zbest_rm[[i]]==k))
  })
  sdaysNum <- length(which(sdays>0))
  # how many unique people
  ppl <- data.frame(rmlist, sdays)
  jTotal <- length(unique(rmlist))
  sppl <- sum(sapply(1:jTotal, FUN = function(j){
    as.numeric(sum(ppl$sdays[which(ppl$rmlist==j)])>0)
  }))
  # total number of time points
  timeTotal <- length(which(unlist(zbest_rm)==k))
  
  tab2[k,] <- c(sdaysNum, sppl, timeTotal)
}
colnames(tab2) <- c("sampling days", "unique people", "time points")

xtable(tab2, digits = 0)



