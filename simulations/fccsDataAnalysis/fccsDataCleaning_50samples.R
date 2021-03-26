## save cleaned data 
# take out some observations and test imputation method
rm(list = ls())
gc()

# libraries
library(R.matlab)
library(rlist)


# function 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# Fort Collins Commuter Study participant data 
alldata <- readMat(con = "simulations/fccsDataAnalysis/Participant_Final_deidentified.mat")
attributes(alldata)
dim(alldata$Participant)
head(alldata$Participant)
class(alldata$Participant)
alldata$Participant[1,1,46] # 46 participants 

p1 <- alldata$Participant[1,1,1]$Commute[,,5]

# Select only complete day for all three exposures
cc.df.new <- data.frame()
howMuchMiss <- list()
tpm <- list()
for(i in 1:46){
  numDays <- dim(alldata$Participant[1,1,i]$Commute)[3]; numDays
  if(is.null(numDays)) howMuchMiss[[i]] <- NULL
  if(!is.null(numDays)){
    howMuchMiss[[i]] <- numeric(numDays)
    tpm[[i]] <- numeric(numDays)
    for(d in 1:numDays){
      personday <- alldata$Participant[1,1,i]$Commute[,,d]
      
      # MAR
      l1 <- length(which(is.na(personday$BC)))
      l2 <- length(which(is.na(personday$CO)))
      l3 <- length(which(is.na(personday$PM25)))
      # below LOD 
      l4 <- length(which(personday$BC == 0))
      l5 <- length(which(personday$CO == 0))
      l6 <- length(which(personday$PM25 == 0))
      
      howMuchMiss[[i]][d] <- max(c(l1+l4, l2+l5, l3+l6))/8640 # max percent missing over exposures 
      percentMiss <- c(l1+l4, l2+l5, l3+l6)/8640; percentMiss # percent missing for each exposure 
      totalPercentMiss <- (l1+l2+l3+l4+l5+l6)/(8640*3); totalPercentMiss # total percent missing 
      
      tpm[[i]][d] <- totalPercentMiss
      
      marmiss <- (l1+l2+l3)/(8640*3); marmiss
      lodmiss <- (l4+l5+l6)/(8640*3); lodmiss
      
      # 10 % total missing 
      if(totalPercentMiss < 0.10){
        cc.df.new <- rbind(cc.df.new, c(i,d,round(totalPercentMiss, 3), round(marmiss,3), round(lodmiss,3)))
      }
      
    }
  }
}
colnames(cc.df.new) <- c("person", "day", "percentMiss", "marmiss", "lodmiss"); cc.df.new
length(unique(cc.df.new$person))
dim(cc.df.new)

percentMiss <- unlist(tpm)
range(tpm)
range(cc.df.new$percentMiss)

# plot some person days 
idx <- 1
i <- cc.df.new$person[idx]; i
d <- cc.df.new$day[idx]; d
timeLength <- 8640
personday <- alldata$Participant[1,1,i]$Commute[,,d]
plot(1:timeLength, personday$PM25, type = "l", col = "red")
lines(1:timeLength, personday$CO, type = "l", col = "blue")
lines(1:timeLength, personday$BC, type = "l", col = "black")

# less than 10% missing data with at least 5 repeated sampling days 
wy <- unique(cc.df.new$person)[as.numeric(which(table(cc.df.new$person)>=5))]

cc.df.save <- cc.df.new[which(cc.df.new$person%in%wy),]
cc.df.save
dim(cc.df.save) # 50 sampling days 
length(unique(cc.df.save$person)) # 9 unique people with at least 5 sampling days each 

cc.df = cc.df.save

# plot some person days 
idx <- 1
i <- cc.df$person[idx]; i
d <- cc.df$day[idx]; d
timeLength <- 8640
personday <- alldata$Participant[1,1,i]$Commute[,,d]
plot(1:timeLength, personday$PM25, type = "l", col = "red")
lines(1:timeLength, personday$CO, type = "l", col = "blue")
lines(1:timeLength, personday$BC, type = "l", col = "black")

# other covariates
personday <- alldata$Participant[1,1,i]$Commute[,,d]
plot(1:timeLength, personday$BPM, type = "l") # heartrate?
plot(1:timeLength, personday$activ, type = "l")
plot(1:timeLength, personday$Temp, type = "l")
plot(1:timeLength, personday$RH, type = "l")
plot(1:timeLength, personday$Lu, type = "l")
plot(1:timeLength, personday$acc, type = "l")


# DATA ANALYSIS 
y <- list() # exposures averaged over 5 minutes 
N <- nrow(cc.df) # number of persondays 
me <- list() # FCCS microenvironment

for(idx in 1:N){
  i <- cc.df$person[idx]
  d <- cc.df$day[idx]
  personday <- alldata$Participant[1,1,i]$Commute[,,d]
  pde <- cbind(personday$BC, personday$CO, personday$PM25)
  
  # mode of microenvironment during 5min period 
  me1 <- personday$Activity
  me.temp <- apply(me1, 2, FUN = function(x) {
    sapply(1:288, FUN = function(t){
      getmode(x[((30*t + 1)-30): (30*t)])
    })
  })
  me_names <- as.numeric(me.temp)
  # relable microenvironments
  me_names[which(me_names==0)] <- "other"
  me_names[which(me_names==1)] <- "home"
  me_names[which(me_names==2)] <- "work"
  me_names[which(me_names==3)] <- "transit"
  me_names[which(me_names==5)] <- "transit"
  me_names[which(me_names==6)] <- "other"
  me_names[which(me_names==7)] <- "other"
  me_names[which(me_names==8)] <- "other"
  me_names[which(me_names==9)] <- "transit"
  me_names[which(me_names==11)] <- "transit"
  me_names[which(me_names==12)] <- "other"
  me_names[which(me_names==13)] <- "eateries"
  me_names[which(me_names==14)] <- "eateries"
  me[[idx]] <- me_names
  
  # observed, LOD, or MAR over 5 minute period
  # if over 75% observed then calculate mean 
  # if less than 75% observed, whatever missing data is more 
  trim1 <- sapply(1:288, FUN = function(t){
    fiveMin <- pde[,1][((30*t + 1)-30): (30*t)]
  })
  mar1 <- apply(trim1, 2, FUN = function(x){
    length(which(is.na(x)))/30
  })
  lod1 <- apply(trim1, 2, FUN = function(x){
    length(which(x==0))/30
  })
  
  trim2 <- sapply(1:288, FUN = function(t){
    fiveMin <- pde[,2][((30*t + 1)-30): (30*t)]
  })
  mar2 <- apply(trim2, 2, FUN = function(x){
    length(which(is.na(x)))/30
  })
  lod2 <- apply(trim2, 2, FUN = function(x){
    length(which(x==0))/30
  })
  
  trim3 <- sapply(1:288, FUN = function(t){
    fiveMin <- pde[,3][((30*t + 1)-30): (30*t)]
  })
  mar3 <- apply(trim3, 2, FUN = function(x){
    length(which(is.na(x)))/30
  })
  lod3 <- apply(trim3, 2, FUN = function(x){
    length(which(x==0))/30
  })
  
  obsP <- 0.9 # percent of time period that must be observed for observed data to be calculated 
  
  obsT1 <- which(mar1+lod1<(1-obsP))
  marT1 <- which(mar1 > obsP/2)
  lodT1 <- which(lod1 >= obsP/2)
  finTrim1 <- apply(trim1, 2, FUN = function(x) mean(x, na.rm = TRUE))
  finTrim1[marT1] <- NA
  finTrim1[lodT1] <- -Inf
  
  
  obsT2 <- which(mar2+lod2<(1-obsP))
  marT2 <- which(mar2 > obsP/2)
  lodT2 <- which(lod2 >= obsP/2)
  finTrim2 <- apply(trim2, 2, FUN = function(x) mean(x, na.rm = TRUE))
  finTrim2[marT2] <- NA
  finTrim2[lodT2] <- -Inf
  
  
  obsT3 <- which(mar3+lod3<(1-obsP))
  marT3 <- which(mar3 > obsP/2)
  lodT3 <- which(lod3 >= obsP/2)
  finTrim3 <- apply(trim3, 2, FUN = function(x) mean(x, na.rm = TRUE))
  finTrim3[marT3] <- NA
  finTrim3[lodT3] <- -Inf
  
  y[[idx]] <- cbind(finTrim1, finTrim2, finTrim3)
  
}

t.max <- 288 # length of time series 

# plot the trimmed data 
i = 13
plot(1:t.max, y[[i]][,3], type = "l", lwd = 2)
lines(1:t.max, y[[i]][,2], type = "l", lwd = 2, col = "red")
lines(1:t.max, y[[i]][,1], type = "l", lwd = 2, col = "blue")

ymatrix <- NULL
for(idx in 1:N){
  ymatrix <- rbind(ymatrix, y[[idx]])
}

apply(ymatrix, 2, FUN = function(x){
  length(which(is.na(x)))/length(x)
})

apply(ymatrix, 2, FUN = function(x){
  length(which(x==-Inf))/length(x)
})

#
y.save = ymatrix
y.save[which(ymatrix==-Inf)] <- NA


# lod from Kirsten's paper 
# BC=0.01, CO=0.01, PM=1

msd <- apply(log(y.save), 2, FUN= function(y){
  c(mean(y, na.rm = TRUE), sd(y, na.rm = TRUE))
})

colnames(msd) <- c("BC", "CO", "PM")

msd

bcLod <- (log(0.01)-msd[1,1])/msd[2,1]; bcLod
coLod <- (log(0.01)-msd[1,2])/msd[2,2]; coLod
pmLod <- (log(1)-msd[1,3])/msd[2,3]; pmLod


lods <- data.frame(bcLod, coLod, pmLod)
rownames(lods) = NULL
lods

# log transform and scale the exposures that are not missing 
yMatScaled <- apply(y.save, 2, FUN = function(y){
  scale(log(y))
})

hist(yMatScaled[,1])
hist(yMatScaled[,2])
hist(yMatScaled[,3])

apply(yMatScaled, 2, FUN = function(y){
  min(y, na.rm = TRUE)
})

dim(ymatrix)
which(ymatrix==-Inf)

yMatScaled[which(ymatrix==-Inf)]<--Inf

apply(yMatScaled, 2, FUN = function(x){
  length(which(is.na(x)))/length(x)
})
apply(yMatScaled, 2, FUN = function(x){
  length(which(x==-Inf))/length(x)
})

# save the data to .csv
write.table(yMatScaled, "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/fccsDataMatrix_50samples.csv",
            row.names = FALSE, col.names = FALSE) 

# me is a list of microenvironments for each person ##
list.save(me, file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/micEnv_50samples.rds")
list.save(lods, file = "~/Documents/Lauren/Rpackages/psbpHMM/simulations/fccsDataAnalysis/limitOfDetection.rds")



