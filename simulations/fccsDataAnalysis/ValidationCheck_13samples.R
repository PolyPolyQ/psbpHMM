# Case study for mc-iHMM on FCCS data 

rm(list = ls())
gc()

#################
### Libraries ###
#################

library(psbpHMM)
library(rlist)

####################
### Read in Data ###
####################

# Mac
yMatScaled <- as.matrix(read.table("~/Documents/Lauren/Rpackages/markovPSBP/simulations/FCCS_cleanedMatrix.csv"),
                       row.names = NULL, col.names = NULL)
colnames(yMatScaled) <- NULL

############################
### Set up exposure data ###
############################

n <- 13
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
rmlist = NULL
q <- ncol(X[[1]])
p <- 3 

#############################
### take out missing data ###
#############################

# 10% missing 
marmis <- .05
lodmis <- .05
ycomplete <- y # save
marRemove <- TRUE
lodRemove <- TRUE
ymat <- NULL

for(i in 1:n){
  ymat <- rbind(ymat, y[[i]])
}
if(lodRemove){
  lod <- apply(ymat, 2, FUN = function(x) quantile(x, lodmis))
  for(i in 1:n){
    for(j in 1:p){
      y[[i]][which(y[[i]][,j] <= lod[j]),j] <- -Inf
    }
  }
}
if(marRemove){
  # chunks of size 1 to 10
  misnum <- ceiling(marmis*t.max*p)
  numChunks <- marmis*200; numChunks
  chunks <- rmultinom(1, size = misnum, prob = rep(1/numChunks, numChunks)) # generate size of chunks 
  for(i in 1:n){
    for(k in 1:numChunks){
      j <- sample(1:p, 1, prob = rep(1/p,p)) # randomly choose an exposure
      repeat{
        timePoint <- sample(2:(t.max-chunks[k,]),1, prob = rep(1/(t.max-chunks[k,]), (t.max-chunks[k,]-1)) ) # choose a time point to start at
        if(!anyNA(y[[i]][,j][timePoint:(timePoint+chunks[k,]-1)])){
          break
        }
      }
      y[[i]][,j][timePoint:(timePoint+chunks[k,]-1)] <- NA
    }
  }
}

########################
### Model Parameters ###
########################

#################################
### Set Priors and Parameters ###
#################################

niter = 5000
nburn = 2500
len.imp = 400

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

#################################################
### FCCS data analysis without random effects ###
#################################################

rmlist = NULL
st1 = Sys.time()
fitMI1 <- mciHMM(niter=niter, nburn=nburn, y=y, rmlist=rmlist, ycomplete=ycomplete, X=X,
                   priors=priors, K.start=K.start, z.true=z.true, lod=lod,
                   mu.true=mu.true, missing = missing, 
                   tau2 = tau2, a.tune = a.tune, b.tune = b.tune,
                   resK = resK, eta.star = eta.star, len.imp = len.imp)
en1 = Sys.time()
en1-st1

####################
### Save results ###
####################

## best clustering for inference on states ##
zbest1 <- bestClusteriHMM(fitMI1)

## model averaged estimates of harmonic trend ##
beta.save1 <- modelAveBeta(fit = fitMI1, zbest1 = zbest1, X = X[[1]])
beta.mean1 <- beta.save1$betaXma
beta.upr1 <- beta.save1$betaXupr
beta.lwr1 <- beta.save1$betaXlwr

## model averaged estimates of mu 
mu.save <- modelAveMu(fit = fitMI1, zbest1 = zbest1)
mu.mean <- mu.save$mu_ma
mu.lwr <- mu.save$mu_lwr
mu.upr <- mu.save$mu_upr

MIresults1 <- list(mse.bias = cbind(fitMI1$mar.mse, fitMI1$lod.mse, fitMI1$mar.bias, fitMI1$lod.bias),
                   zbest = zbest1, 
                   mu.mean = mu.mean, mu.lwr = mu.lwr, mu.upr = mu.upr,
                   beta.mean = beta.mean1, beta.upr = beta.upr1, beta.lwr = beta.lwr1,
                   lod.imputes = fitMI1$ylod, mar.imputes = fitMI1$ymar, ycomplete = fitMI1$ycomplete,
                   mismat = fitMI1$mismat, MH1 = fitMI1$MH.arate, MH2 = fitMI1$MH.lamrate)

list.save(MIresults1, file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/validation13samples.rds")

####################
### View Results ###
####################


menv <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/menv13samples.rds")

####################
## exposure means ## 
####################

fccs1 <- list.load("~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/validation13samples.rds")

mmean <- matrix(unlist(fccs1$mu.mean), ncol = 3, byrow = TRUE)
mlwr <- matrix(unlist(fccs1$mu.lwr), ncol = 3, byrow = TRUE)
mupr <- matrix(unlist(fccs1$mu.upr), ncol = 3, byrow = TRUE)

K = length(unique(unlist(fccs1$zbest)));K

cols = c("blue", "red", "green", "purple", "orange", "black", "gray", "darkblue",
         "yellow", "pink", "darkgreen", "lightgray", "lightgreen", "darkred")

col = sample(cols, K, replace = TRUE)

mdat1 = data.frame(x = (1:K)-0.2, y = mmean[,1], lwr = mlwr[,1], upr = mupr[,1], col = col)
mdat2 = data.frame(x = (1:K), y = mmean[,2], lwr = mlwr[,2], upr = mupr[,2], col = col)
mdat3 = data.frame(x = (1:K)+0.2, y = mmean[,3], lwr = mlwr[,3], upr = mupr[,3], col = col)

mplot1 <- ggplot() + ggtitle("mc-iHMM model averaged exposure means") + theme(legend.position = "none") + 
  scale_x_discrete(name = "Hidden States", limits = factor(1:K), breaks = 1:K) + 
  scale_y_continuous(name = "Model Averaged Exposure Means") + 
  geom_point(data = mdat1, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 15) + 
  geom_errorbar(data = mdat1, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  geom_point(data = mdat2, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 16) + 
  geom_errorbar(data = mdat2, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  geom_point(data = mdat3, mapping = aes(x = x, y = y, colour = col), size = 4, shape = 17) + 
  geom_errorbar(data = mdat3, mapping = aes(x = x, ymin = lwr, ymax = upr, colour = col), width = .4) +
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 20, angle = 0, hjust = 1), 
        axis.text.y = element_text(size = 15, angle = 0, hjust = 1), 
        plot.title = element_text(size = 15)) 

mplot1

#######################
## stacked bar chart ##
#######################

K <- length(unique(unlist(fccs1$zbest))); K
zall <- unlist(fccs1$zbest) # all hidden states
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
hidar <- sort(table(unlist(fccs1$zbest)), decreasing = TRUE)
hidar <- as.data.frame(hidar)
hmode <- hidar[which(hidar$Freq>112),]
hmode$Var1
datr <- dat[which(dat$hid%in%hmode$Var1),]
datr
# all hidden states 
bchart1 <- ggplot(data=dat, aes(fill=ms, y=nums, x=hid)) + 
  ggtitle("mc-iHMM microenvironments and hidden states") +
  theme(text = element_text(size = 10), axis.text.x = element_text(size = 10)) + 
  geom_bar(position="stack", stat="identity", alpha = 0.9) +
  scale_fill_manual(name = "Microenvironment", values = c("#fff700", "#9582ff", "#ff0004", "#02e0ad", "#00a6ff"),
                    labels = c("other", "eateries", "transit", "work", "home")) +
  xlab("Hidden State") + ylab("Count of Time Points") 
bchart1








