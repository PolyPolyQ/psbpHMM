
rm(list = ls())
gc()

#############
### Facet ###
#############
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

library(ggplot2)
library(gridExtra)
library(ggthemes)
library(rlist)


# data 
yMatScaled <- as.matrix(read.table("~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fccsDataMatrix_50samples.csv"),
                        row.names = NULL, col.names = NULL)
colnames(yMatScaled) <- NULL
lod <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/limitOfDetection.rds")
lod <- as.numeric(lod); lod
n <- 50
t.max <- 288
splits <- seq(1, nrow(yMatScaled)-287, length.out = n)
y <- lapply(1:length(splits), FUN = function(t) yMatScaled[splits[t]:(splits[t]+t.max-1), ])
# microenvironments
menv <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/micEnv_50samples.rds")

# results 
fit_rcpp_results <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_save.rds")
fit_rcpp_rm_results <- list.load(file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fit_rcpp_rm_save.rds")

fit_rcpp = fit_rcpp_results$fit1
zbest = fit_rcpp_results$zbest

fit_rcpp_rm = fit_rcpp_rm_results$fit1
zbest_rm = fit_rcpp_rm_results$zbest

rmlist = c(rep(1,5), rep(2,5), rep(3,6), rep(4,5), rep(5,5), rep(6,5),
           rep(7,9), rep(8, 5), rep(9,5))

###################
### Facet Grids ###
###################

########################
## micro-environments ##
########################

# person-days to go in figure 
i1 = 49
i2 = 50

# y limits on figure 
rlow = min(unlist(y)[which(unlist(y)>-Inf)], na.rm=T); rlow
rhi = max(unlist(y), na.rm=T); rhi


# factor the microenvironents 

# microenvironment person i1 
states <- menv[[i1]]; states
xstarts = numeric()
for(t in 1:(length(states)-1)){
  if(states[t]!=states[t+1]) {
    print(t+1)
    xstarts = c(xstarts, t+1)
  }
}
xe = c(xstarts-1, 288); xe
xs = c(1, xstarts); xs

stateCol <- states[xs]; stateCol
rects <- data.frame(xstart = xs, xend = xe, col = stateCol)
rects
micen <- rects
micen$xstart <- rects$xstart/288*24
micen$xend <- micen$xend/288*24
#micen$col <- as.numeric(micen$col)
micen
dim(micen)
me1 <- cbind(rbind(micen, micen, micen), 
             person = paste("i=", i1), 
             j=c(rep("BC",nrow(micen)), rep("CO",nrow(micen)), rep("PM2.5",nrow(micen)))); me1

# microenvironment person i2
states <- menv[[i2]]; states
xstarts = numeric()
for(t in 1:(length(states)-1)){
  if(states[t]!=states[t+1]) {
    print(t+1)
    xstarts = c(xstarts, t+1)
  }
}
xe = c(xstarts-1, 288); xe
xs = c(1, xstarts); xs
stateCol <- states[xs]; stateCol
rects <- data.frame(xstart = xs, xend = xe, col = stateCol)
rects
micen <- rects
micen$xstart <- rects$xstart/288*24
micen$xend <- micen$xend/288*24
#micen$col <- as.numeric(micen$col)
micen
dim(micen)
me2 <- cbind(rbind(micen, micen, micen), 
             person = paste("i=", i2), 
             j=c(rep("BC",nrow(micen)), rep("CO",nrow(micen)), rep("PM2.5",nrow(micen)))); me2

# combined ME's
meboth <- data.frame(rbind(me1, me2)); meboth
menv_cols = rev(c("#fff700", "#00a6ff", "#ff0004", "#9582ff","#02e0ad")); menv_cols
menv_labs = c("home", "other", "transit", "work", "eateries")

#############
## mc-iHMM ##
#############

# x and y data 
ys = c(as.vector(y[[i1]]), as.vector(y[[i2]])); ys
ys[which(ys==-Inf)] = NA
times = rep(((1:288)/288)*24,6); times
js <- rep(c(rep("BC", 288), rep("CO", 288), rep("PM2.5", 288)),2) # which exposure
cols <- factor(c(rep(zbest[[i1]],3), rep(zbest[[i2]],3))) # hidden states
pers <- c(rep(paste("i=", i1), 288*3), rep(paste("i=",i2), 288*3)) # which person

# combined ME's
meboth <- data.frame(rbind(me1, me2)); meboth

# make data frame
dat.test <- data.frame(x = times, j = js, y = ys,
                       col = cols, person = pers)

which_menv = unique(c(menv[[i1]], menv[[i2]])); which_menv
unique(menv[[i1]])
which_menv_cols = menv_cols[which(menv_labs %in% which_menv)]; which_menv_cols
which_menv_labs = menv_labs[which(menv_labs %in% which_menv)]; which_menv_labs

p1 <- ggplot(data = dat.test, aes(x = x, y = y)) + ggtitle("mc-iHMM") +
  scale_x_continuous(name=NULL, limits=c(0, 24), breaks = seq(0,24,4),
                     labels = c("9:00pm", "1:00am", "5:00am", "9:00am", "1:00pm", "5:00pm", "9:00pm")) +
  scale_y_continuous(name=NULL, limits=c(rlow, rhi), breaks = seq(-5,5,2)) +
  theme_bw() +
  geom_rect(data = meboth, aes(x = NULL, y = NULL, xmin = xstart, xmax = xend, ymin = rlow,
                               ymax = rhi, fill = factor(col)), alpha = .4) +
  scale_fill_manual(values = which_menv_cols,
                    labels = which_menv_labs) +
  geom_point(aes(x = x, y = y, colour = col), show.legend = FALSE) +
  theme(legend.position = "bottom") + labs(fill = "Microenvironment") +
  facet_grid(rows = vars(j), cols = vars(person))
p1


################
## mc-iHMM-rm ##
################

# x and y data 
ys = c(as.vector(y[[i1]]), as.vector(y[[i2]])); ys 
ys[which(ys==-Inf)] = NA
times = rep(((1:288)/288)*24,6); times
js <- rep(c(rep("BC", 288), rep("CO", 288), rep("PM2.5", 288)),2) # which exposure
cols <- factor(c(rep(zbest_rm[[i1]],3), rep(zbest_rm[[i2]],3))) # hidden states
pers <- c(rep(paste("i=", i1), 288*3), rep(paste("i=",i2), 288*3)) # which person

# make data frame 
dat.test <- data.frame(x = times, j = js, y = ys, 
                       col = cols, person = pers)

which_menv = unique(c(menv[[i1]], menv[[i2]])); which_menv

which_menv_cols = menv_cols[which(menv_labs %in% which_menv)]; which_menv_cols
which_menv_labs = menv_labs[which(menv_labs %in% which_menv)]; which_menv_labs

p2 <- ggplot(data = dat.test, aes(x = x, y = y)) + ggtitle("mc-iHMM-rm") +
  scale_x_continuous(name=NULL, limits=c(0, 24), breaks = seq(0,24,4),
                     labels = c("9:00pm", "1:00am", "5:00am", "9:00am", "1:00pm", "5:00pm", "9:00pm")) + 
  scale_y_continuous(name=NULL, limits=c(rlow, rhi), breaks = seq(-5,5,2)) +
  theme_bw() + 
  geom_rect(data = meboth, aes(x = NULL, y = NULL, xmin = xstart, xmax = xend, ymin = rlow,
                               ymax = rhi, fill = factor(col)), alpha = .4) +
  scale_fill_manual(values = which_menv_cols,
                    labels = which_menv_labs) + 
  geom_point(aes(x = x, y = y, colour = col), show.legend = FALSE) +
  theme(legend.position = "bottom") + labs(fill = "Microenvironment") +
  facet_grid(rows = vars(j), cols = vars(person)) 
p2

