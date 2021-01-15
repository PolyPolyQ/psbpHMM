

####################################################################################
### Function 1: mc-iHMM ###
### fits the cyclical model for multple time series ###
### allows for different X_i for each i to permit covariates or different times ###
### all time series are the same length ### 
####################################################################################

### X is a list with a matrix for each i ### 
if(missing){
  algorithm = "MH"
  SigmaPrior = "wishart"
}else{
  algorithm = "Gibbs"
  SigmaPrior = "non-informative"
}

#####################
### Initial Setup ###
#####################

# how many time series and exposures 
if(class(y)=="list"){
  p <- ncol(y[[1]]) # number of exposures 
  n <- length(y) # number of time series
  t.max <- nrow(y[[1]]) # number of time points 
}else if(class(y)=="matrix"){
  p <- ncol(y)
  n <- 1
  t.max <- nrow(y)
  y <- list(y) # make y into a list 
  ycomplete <- list(ycomplete)
  z.true <- list(z.true)
}else if(class(y)=="numeric"){
  p <- 1
  n <- 1
  t.max <- length(y)
  y <- list(matrix(y, ncol = 1))
  ycomplete <- list(matrix(ycomplete, ncol= 1))
  z.true <- list(z.true)
}

# X is a list 
q <- ncol(X[[1]])

##############
### Priors ###
##############

if(missing(priors)) priors <- NULL
# mu
if(is.null(priors$mu0)) priors$mu0 <- matrix(0, p, 1) # mu_k|Sigma_k ~ N(mu0, 1/lambda*Sigma_k) prior mean on p exposures
# alpha
if(is.null(priors$mu.alpha)) priors$mu.alpha <- 0 # fixed mean parameter prior on intercepts alpha.jk
if(is.null(priors$m0)) priors$m0 <- 0 # fixed mean on m.alpha, prior mean for alpha.jj # higher so self-transition prob is higher
if(is.null(priors$v0)) priors$v0 <- 1 # fixed variance on m.alpha, prior mean for alpha.jj
# beta
if(is.null(priors$mu.beta)) priors$mu.beta <- rep(0, q) # basis weight prior mean 
if(is.null(priors$Sigma.beta)) priors$Sigma.beta <- diag(q) # basis weight prior variance
if(!is.null(X)) priors$SigInv.beta <- invMat(priors$Sigma.beta) # cppFunction

# hyperiors on alpha
if(is.null(priors$a1)) priors$a1 <- 1 # shape for sig2inv.alpha
if(is.null(priors$b1)) priors$b1 <- 1 # rate for sig2inv.alpha
if(is.null(priors$a2)) priors$a2 <- 1 # shape for vinv.alpha
if(is.null(priors$b2)) priors$b2 <- 1 # rate for vinv.alpha

if(SigmaPrior == "wishart"){
  # missing data case (MH udpates)
  if(is.null(priors$R)) priors$R <- diag(p) # Sigma_k ~ Inv.Wish(nu, R) hyperparameter for Sigma_k
  if(is.null(priors$nu)) priors$nu <- p+2 # Sigma_k ~ Inv.Wish(nu, R) hyperparameter for Sigma_k, nu > p+1
  nu.df <- priors$nu # 
  R.mat <- priors$R 
}else if(SigmaPrior == "non-informative"){
  # complete data case (Gibbs updates)
  if(is.null(priors$bj)) priors$bj <- rep(1, p) # Huang and Wand advise 10e5
  if(is.null(priors$nu)) priors$nu <- 2 # Huang and Wand advise 2, p+4 for so the variance exists
  nu.df <- priors$nu + p - 1 # Huang and Wand, prior df 
  aj.inv <- rgamma(p, shape = 1/2, rate = 1/(priors$bj^2)) # these are 1/aj 
  # starting value for R.mat, the matrix parameter on Sigma.Inv
  R.mat <- 2*priors$nu*diag(aj.inv) 
  # we model aj.inv with gamma(shape = 1/2, rate = 1/(b_j^2))
}

if(is.null(priors$lambda)) priors$lambda <- 1 

#############################
### Indicate Missing Type ###
#############################

# indicate missing data: obs = 0, mar = 1, lod = 2
mismat <- list()
for(i in 1:n){
  mismat[[i]] <- matrix(sapply(y[[i]], ismissing), ncol = p)
}

mism <- numeric()
for(i in 1:n){
  mism <- rbind(mism, mismat[[i]])
}

# for each i, which time points have any missing data? 
missingTimes <- lapply(1:n, FUN = function(i){
  which(apply(mismat[[i]],1,sum)>0)
})

# for each i, which times points are observed? 
observedTimes <- lapply(1:n, FUN = function(i){
  which(apply(mismat[[i]], 1, sum)==0)
})

###############################################
### Impute Starting Values for Missing Data ###
###############################################

for(i in 1:n){
  if(any(mismat[[i]]==2)){ # lod 
    expLod <- exp(lod)
    numlod <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==2))) # how many LOD
    for(t in which(numlod>0)){ # loop thru LOD data
      whichlod <- which(mismat[[i]][t,]==2) # which exposures are below LOD 
      y[[i]][t,whichlod] <- log(expLod[whichlod]/sqrt(2)) # impute the LOD with the log(LOD/sqrt(2))
    }
  }
  if(any(mismat[[i]]==1)){ # mar
    nummis <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==1))) # how many missing at each time point 
    for(t in which(nummis>0)){ # only loop thru time points with MAR
      whichmis <- which(mismat[[i]][t,]==1) # which exposures are missing at random 
      if(t == 1){
        for(ws in whichmis){
          lastT = max(which(!is.na(y[[i]][,ws])))
          y[[i]][t,ws] <- y[[i]][lastT, ws] # fill in with the last observed value from the end of the time series
        }
      }else{
        y[[i]][t,whichmis] <- y[[i]][t-1, whichmis] # fill in the missing by LVCF
      }
    }
  }
} 

#####################
## Starting Values ##  
#####################

z <- list()
for(i in 1:n){
  K <- K.start
  if(is.null(K)) K <- 12
  z[[i]] <- sample(1:K, t.max, replace = TRUE)
}

mu <- list()
Sigma <- list()
D <- list()
L <- list()
lams <- list()
al <- list()

ymatrix <- NULL
for(i in 1:n){
  ymatrix <- rbind(ymatrix, y[[i]])
}

for(k in 1:K){
  if(algorithm == "MH"){
    # we reparameterize Sigma and model L and D instead 
    vj0 <- sapply(1:p, FUN = function(j) priors$nu + j - p); vj0 # fixed for each k 
    deltaj0 <- rep(1,p); deltaj0 # fixed for each k 
    lams[[k]] <- 1/rgamma(3, vj0, rate = deltaj0)
    D[[k]] <- diag(lams[[k]])
    al.list <- list()
    for(j in 2:p){
      al.list[[j-1]] <- rnorm(j-1, 0, lams[[k]][j])
    }
    al[[k]] <- unlist(al.list) # for j = 2 to p
    which.lams <- unlist(sapply(2:p, FUN = function(j) rep(j,j-1))) # which lams to use for each al 
    L[[k]] <- diag(p)
    lowerTriangle(L[[k]]) <- al[[k]]
    
    # make a C++ function for solve 
    #Sigma[[k]] <- solve(L[[k]])%*%D[[k]]%*%t(solve(L[[k]])) 
    Sigma[[k]] <- mhDecomp(L[[k]], D[[k]]) # cppFunction
    mu[[k]] <- mvnfast::rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
  }else{
    Sigma[[k]] <- invMat(matrix(rWishart(1, df = nu.df, Sigma = invMat(R.mat)),p,p))
    mu[[k]] <- rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
  }
}


# K only 
alpha.0k <- rep(0, K) # initial state intercept
alpha.jk <- list() # state intercepts
m.alpha <- priors$m0 # mean on alpha.jj
sig2inv.alpha <- priors$a1/priors$a2 # precision on alpha.jk
vinv.alpha <- priors$b1/priors$b2 # precision on alpha.jj

# K only 
beta.k <- list() # regression coefficients
for(k in 1:(K)){
  alpha.jk[[k]] <- rnorm(K, priors$mu.alpha, sqrt(1/sig2inv.alpha))
  alpha.jk[[k]][k] <- priors$m0
  beta.k[[k]] <- matrix(rmvn(1, priors$mu.beta, priors$Sigma.beta), nrow = q, ncol = 1) # length q
}

################################
### update transition matrix ###
################################

ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)
pi.z = updatePi2(beta = beta.k, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)

# fixed values for new state probabilities 
gampp <- mgamma(nu = nu.df, p = p)
# fixed for "wish", starting value for "ni" because R.mat will change with each iteration as aj.inv changes
# if "ni" then need to update detR.star and log.stuff each iteration
detR.star <- mclapply(1:n, FUN = function(i){
  sapply(1:t.max, FUN = function(t){
    x <- R.mat + priors$lambda*tcrossprod(priors$mu0) + tcrossprod(y[[i]][t,]) - 
      (1/(1+priors$lambda))*tcrossprod(priors$lambda*priors$mu0+y[[i]][t,])
    return(det(x))})
})
log.stuff <- (p/2)*log(priors$lambda/(pi*(priors$lambda+1)))+log(gampp)+(nu.df/2)*log(det(R.mat)); log.stuff

########################
### For MH algorithm ###
########################

# calculate before loop for MH 
ymatrix <- numeric()
for(i in 1:n){
  ymatrix <- rbind(ymatrix, y[[i]])
}
ycomp <- ymatrix[which(rowSums(mism)==0),]

##############################
### MCMC Storage and Setup ###
##############################

z.save <- list()
mu.save <- list()
Sigma.save <- list()
K.save <- rep(NA, niter)
ham <- rep(NA, niter)
mu.mse <- rep(NA, niter)
mu.sse <- rep(NA, niter)
MH.a <- 0 # for MH update a
MH.lam <- 0 # for MH update lams

# missing data sets
imputes <- ceiling(seq.int(nburn, niter, length.out = len.imp))
y.mar.save <- matrix(NA, len.imp, length(which(unlist(mismat)==1)))
y.lod.save <- matrix(NA, len.imp, length(which(unlist(mismat)==2)))
mar.mse <- rep(NA, len.imp)
lod.mse <- rep(NA, len.imp)
mar.sse <- rep(NA, len.imp)
lod.sse <- rep(NA, len.imp)
miss.mse <- rep(NA, len.imp)
mar.bias <- rep(NA, len.imp)
lod.bias <- rep(NA, len.imp)
mar.sum.bias <- rep(NA, len.imp)
lod.sum.bias <- rep(NA, len.imp)
s.imp <- 1

beta.save <- list()

###############
### Sampler ###
###############

for(s in 1:niter){
  
  #####################
  ### initial stuff ### ### done 
  #####################
  
  if(s %% 10 == 0) print(paste("iteration", s))
  #start.time1 <- Sys.time()
  z.prev <- list()
  z.prev <- mclapply(1:n, FUN=function(i) return(z[[i]]))
  
  
  ################
  ### update u ### ### done 
  ################
  
  u <- list()
  for(i in 1:n){
    u[[i]] <- unlist(lapply(1:t.max, FUN = function(t){
      if(t==1){
        return(runif(1, 0, pi.z[[i]][[1]][z[[i]][1]]))
      }else{
        return(runif(1, 0, pi.z[[i]][[t]][z[[i]][t-1], z[[i]][t]]))
      }
    }))
  }
  
  #############################################
  ### update the possible states for each t ###  ### done
  #############################################
  
  state.list <- list()
  # state.list is a list of states that belong to some possible trajectory for each time point
  for(i in 1:n){
    # forward condition: considering the previous possible states, what are the current possible states?
    for.list <- list()
    for.list[[1]] <- which(pi.z[[i]][[1]] >= u[[i]][1])
    for(t in 2:t.max){
      
      for.list[[t]] <- sort(unique(unlist(lapply(intersect(for.list[[t-1]], 1:K),
                                                 FUN = function(k) which(pi.z[[i]][[t]][k,] >= u[[i]][t])))))
    }
    # backward condition
    back.list <- list()
    back.list[[t.max]] <- for.list[[t.max]]
    for(t in (t.max-1):1){
      back.list[[t]] <- sort(unique(unlist(lapply(intersect(back.list[[t+1]], 1:K),
                                                  FUN = function(k) which(pi.z[[i]][[t+1]][,k] >= u[[i]][t+1])))))
    }
    state.list[[i]] <- lapply(1:t.max, FUN = function(t) {
      return(intersect(for.list[[t]], back.list[[t]]))
    })
  }
  
  ###########################
  ### Determine t minus 1 ### 
  ###########################
  
  detzAll <- mclapply(1:n, FUN = function(i){
    detZminus1(i = i, state.list.i = state.list[[i]], pi.z = pi.z[[i]], u.i = u[[i]], t.max = t.max)
  })
  
  ################
  ### Sample Z ### 
  ################
  
  ###################################
  ### make this piece faster next ###
  ###################################
  
  cholSigma <- lapply(1:K, FUN = function(k) chol(Sigma[[k]]))
  K <- length(unique(unlist(z)))
  
  # only the joint NIW option 
  z <- mclapply(1:n, FUN = function(i) {
    updateZ(i=i, y.i=y[[i]], mu=mu, cholSigma=cholSigma, detR.stari=detR.star[[i]],
            nu.df=nu.df, K=K, log.stuff=log.stuff, t.max=t.max, 
            detz = detzAll[[i]], state.list.i = state.list[[i]])
  })

  
  #########################
  ### new state fillers ### 
  #########################
  
  if(any(unlist(z)>K)){
    
    if(algorithm == "MH"){
      lamsNew <- 1/rgamma(3, vj0, rate = deltaj0); lamsNew
      DNew <- diag(lamsNew); DNew
      al.listNew <- list()
      for(j in 2:p){
        al.listNew[[j-1]] <- rnorm(j-1, 0, lamsNew[j])
      }
      alNew <- unlist(al.listNew); alNew # for j = 2 to p
      LNew <- diag(p)
      lowerTriangle(LNew) <- alNew; LNew
      
      SigmaNew <- mhDecomp(LNew, DNew) # cppFunction
      
      #SigmaNew <- solve(LNew)%*%DNew%*%t(solve(LNew)) # Sigma[[k]]
      cholSigmaNew <- chol(SigmaNew)
      muNew <- rmvn(1, priors$mu0, (1/priors$lambda)*SigmaNew)
    }else if(algorithm == "Gibbs"){
      SigmaNew <- invMat(matrix(rWishart(1, df = nu.df, Sigma = invMat(R.mat)),p,p))
      muNew <- rmvn(1, priors$mu0, (1/priors$lambda)*SigmaNew)
    }
    
    # sample starting values if we got a new state
    mu[[K+1]] <- muNew
    if(algorithm == "MH"){
      lams[[K+1]] <- lamsNew
      D[[K+1]] <- DNew
      al[[K+1]] <- alNew
      L[[K+1]] <- LNew
    }
    Sigma[[K+1]] <- SigmaNew
    alpha.0k[K+1] <- 0
    for(k in 1:K){
      alpha.jk[[k]][K+1] <- rnorm(1, priors$mu.alpha, sqrt(1/sig2inv.alpha))
    }
    alpha.jk[[K+1]] <- rnorm(K+1, priors$mu.alpha, sqrt(1/sig2inv.alpha))
    alpha.jk[[K+1]][K+1] <- priors$m0
    beta.k[[K+1]] <- matrix(rmvn(1, priors$mu.beta, priors$Sigma.beta), nrow = q, ncol = 1) # length q
    EnterNew <- TRUE # indicator that we entered into a new state this round 
  }else{
    EnterNew <- FALSE # didn't go into a new state 
  }
  
  #########################################
  ### Relabel the States and Parameters ###
  #########################################
  
  z.new <- recode_Z(unlist(z))
  splits <- seq(1,n*t.max, t.max)
  z.new <- lapply(1:length(splits), FUN = function(i) z.new[splits[i]:(splits[i]+t.max-1)])
  K <- length(sort(unique(unlist(z.new)))) # new K 
  zu <- sort(unique(unlist(z))) # z.unique (old labels for current states to grab from) 
  
  alpha.new <- list()
  beta.new <- list()
  mu.new <- list()
  Sigma.new <- list()
  al.new <- list()
  lams.new <- list()
  L.new <- list()
  D.new <- list()
  
  rj <- 1
  for(k in zu){ # 
    alpha.new[[rj]] <- alpha.jk[[k]][zu]
    beta.new[[rj]] <- beta.k[[k]]
    mu.new[[rj]] <- mu[[k]]
    Sigma.new[[rj]] <- Sigma[[k]]
    if(algorithm == "MH"){
      al.new[[rj]] <- al[[k]]
      lams.new[[rj]] <- lams[[k]]
      D.new[[rj]] <- D[[k]]
      L.new[[rj]] <- L[[k]]
    }
    rj <- rj + 1
  }
  
  alpha.jk <- alpha.new
  beta.k <- beta.new
  z <- z.new
  Sigma <- Sigma.new
  mu <- mu.new
  if(algorithm == "MH"){
    al <- al.new
    lams <- lams.new 
    D <- D.new
    L <- L.new
  }
  
  cholSigma <- lapply(1:K, FUN = function(k) chol(Sigma[[k]]))
  
  ######################
  ### update theta_k ###
  ######################
  
  # first update mu and Sigma 
  for(k in 1:K){
    itimes <- lapply(1:n, FUN = function(i)  which(z[[i]] == k))
    nkk.tilde  <- length(unlist(itimes)) # number in state k 
    y.list <- lapply(1:n, FUN = function(i) matrix(y[[i]][itimes[[i]],], ncol = p))
    yk <- numeric()
    for(i in 1:n){
      yk <- rbind(yk, y.list[[i]])
    }
    ybark <- matrix(apply(yk, 2, mean),p,1)
    nu_nk <- nu.df + nkk.tilde 
    
    if(algorithm == "Gibbs"){
      
      mu_nk <- (priors$lambda*priors$mu0 + nkk.tilde*ybark)/(priors$lambda+nkk.tilde) 
      lambda_nk <- priors$lambda + nkk.tilde 
      if(nkk.tilde == 1){
        M <- R.mat
      }else{
        M <- R.mat + (nkk.tilde-1)*cov(yk) 
      }
      Sigma_nk <- M + (priors$lambda*nkk.tilde)/(nkk.tilde + priors$lambda)*tcrossprod(ybark - priors$mu0) 
      Sigma[[k]] <- invMat(matrix(rWishart(1,df=nu_nk, Sigma=invMat(Sigma_nk))),p,p)
      mu[[k]] <- rmvn(n=1, mu=mu_nk, sigma=chol((1/lambda_nk)*as.matrix(Sigma[[k]], p, p)), isChol = TRUE)  
      
    }else if(algorithm == "MH"){
      # update a
      for(j in 1:length(al[[k]])){
        
        if(resK){
          eta <- rgeom(1, (1/eta.star)) + 1
        }else eta <- 1
        
        if(eta>0){
          for(m in 1:eta){
            al.star <- al[[k]]
            L.star <- L[[k]]
            
            #a.star <- rnorm(1, al[[k]][j], sqrt(tau2)); a.star # proposed value 
            a.star <- rnorm(1, 0, sqrt(tau2)); a.star # proposed value 
            
            al.star[j] <- a.star; al.star
            lowerTriangle(L.star) <- al.star; L.star
            
            SigmaStar <- mhDecomp(L.star, D[[k]]) # cppFunction
            
            #SigmaStar <- solve(L.star)%*%D[[k]]%*%t(solve(L.star)); SigmaStar # function of a.star
            
            # likelihoods
            da.curr <- sum(dmvn(yk, mu = mu[[k]], sigma = cholSigma[[k]], log = TRUE, isChol = TRUE)) +
              dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*Sigma[[k]]), log = TRUE, isChol = TRUE)
            da.star <- sum(dmvn(yk, mu = mu[[k]], sigma = chol(SigmaStar), log = TRUE, isChol = TRUE)) +
              dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*SigmaStar), log = TRUE, isChol = TRUE)
            
            # priors 
            pa.curr <- dnorm(al[[k]][j], 0, sqrt(lams[[k]][which.lams[j]]), log = TRUE)
            pa.star <- dnorm(a.star, 0, sqrt(lams[[k]][which.lams[j]]), log = TRUE)
            
            # proposals
            qa.curr <- dnorm(al[[k]][j], 0, sqrt(tau2), log = TRUE)
            qa.star <- dnorm(a.star, 0, sqrt(tau2), log = TRUE)
            
            mh1 <- pa.star + da.star + qa.curr; mh1
            mh2 <- pa.curr + da.curr + qa.star
            # catch error on da.curr
            
            ar <- mh1-mh2
            
            if(runif(1) < exp(ar)){
              al[[k]][j] <- a.star
              L[[k]] <- L.star # update this too, fxn of a.star
              Sigma[[k]] <- SigmaStar # update this too, fxn of a.star
              MH.a <- MH.a + 1 
            }
          }
        }
      }
      # update lams
      for(j in 1:p){
        D.star <- D[[k]]
        lam.star <- 1/rgamma(1, a.tune, rate = b.tune); lam.star # proposed value 
        D.star[j,j] <- lam.star
        
        SigmaStar <- mhDecomp(L[[k]], D.star) # cppFunction
        
        #SigmaStar <- solve(L[[k]])%*%D.star%*%t(solve(L[[k]]))
        
        # likelihoods
        dlam.curr <- sum(dmvn(yk, mu = mu[[k]], sigma = cholSigma[[k]], log = TRUE, isChol = TRUE)) + 
          dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*Sigma[[k]]), log = TRUE, isChol = TRUE)
        dlam.star <- sum(dmvn(yk, mu = mu[[k]], sigma = chol(SigmaStar), log = TRUE, isChol = TRUE)) +
          dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*SigmaStar), log = TRUE, isChol = TRUE)
        
        # priors
        plam.curr <- dinvgamma(lams[[k]][j], vj0[j]/2, rate = deltaj0[j]/2, log = TRUE)
        plam.star <- dinvgamma(lam.star, vj0[j]/2, rate = deltaj0[j]/2, log = TRUE)
        
        # proposal 
        qlam.curr <- dinvgamma(lams[[k]][j], a.tune, b.tune, log = TRUE)
        qlam.star <- dinvgamma(lam.star, a.tune, b.tune, log = TRUE)
        
        mh1 <- plam.star + dlam.star + qlam.curr; mh1
        mh2 <- plam.curr + dlam.curr + qlam.star; mh2
        ar <- mh1-mh2
        
        if(runif(1) < exp(ar)){
          lams[[k]][j] <- lam.star
          D[[k]] <- D.star # fxn of lam.star
          Sigma[[k]] <- SigmaStar # update this too, fxn of lam.star
          MH.lam <- MH.lam + 1 
        }
      }
      
      # update mu by Gibbs
      mu_nk <- (priors$lambda*priors$mu0 + nkk.tilde*ybark)/(priors$lambda+nkk.tilde) 
      lambda_nk <- priors$lambda + nkk.tilde 
      mu[[k]] <- rmvn(n=1, mu=mu_nk, sigma=chol((1/lambda_nk)*as.matrix(Sigma[[k]], p, p)), isChol = TRUE)  
    } # end if MH 
  } # end sample theta  
  
  
  ######################################################################
  ### if SigmaPrior == "NI": Update aj.inv, detR.star and log.stuff  ###
  ######################################################################
  
  if(SigmaPrior == "non-informative"){
    # then update aj.inv for j = 1 to p 
    shape.aj <- (K*(priors$nu+p-1)+1)/2
    sigmajj <- lapply(1:K, FUN = function(k){
      diag(invMat(Sigma[[k]])) # cppFunction
    }) # diagonal elements of each Sigma.Inv_k
    diags <- numeric()
    for(k in 1:K){
      diags <- rbind(diags, sigmajj[[k]])
    }
    sumdiags <- apply(diags, 2, sum) # sum of diagonals of Sigma.Inv for k = 1 to K 
    rate.aj <- 1/(priors$bj^2) + priors$nu*sumdiags
    aj.inv <- rgamma(p, shape = shape.aj, rate = rate.aj) # these are 1/aj 
    R.mat <- 2*priors$nu*diag(aj.inv)
    
    detR.star <- mclapply(1:n, FUN = function(i){
      sapply(1:t.max, FUN = function(t){
        x <- R.mat + priors$lambda*tcrossprod(priors$mu0) + tcrossprod(y[[i]][t,]) - 
          (1/(1+priors$lambda))*tcrossprod(priors$lambda*priors$mu0+y[[i]][t,])
        return(det(x))
      })})
    
    log.stuff <- (p/2)*log(priors$lambda/(pi*(priors$lambda+1)))+log(gampp)+(nu.df/2)*log(det(R.mat))
  }
  
  ################
  ### update W ### 
  ################
  
  # for each t, w.z gives me a VECTOR based on the previous time point and values up to the current time point
  # so each w.z[[t]] should be a VECTOR of length z_t
  
  w.z <- list()
  for(i in 1:n){
    w.z[[i]] <- list()
    for(t in 1:t.max){
      w.z[[i]][[t]] <- sapply(1:z[[i]][t], FUN = function(l) upWdiff(t=t, l=l, i=i, alpha.0k=alpha.0k, X=X[[i]], 
                                                                     beta.k=beta.k, alpha.jk=alpha.jk, z=z))
    }
  }
  
  #######################
  ### update alpha.0k ### ### done 
  #######################
  
  # number of i's s.t z_i1 >= k for each k 
  itime0 <- sapply(1:K, FUN = function(k) sum(sapply(1:n, FUN = function(i) z[[i]][1] >= k)))
  
  # sum over the i's s.t. z_i1 >= k for each k
  # which.i for each k
  which.i <- sapply(1:K, FUN = function(k) {
    if(any(lapply(1:n, FUN = function(i) z[[i]][1] >= k) == TRUE)){
      return(which(lapply(1:n, FUN = function(i) z[[i]][1] >= k) == TRUE))
    }else{
      return(0)
    }
  })
  # for each k, tells me which i's to sum over, if any
  # sum over w-beta for the i's s.t z_i1 >= k for each k 
  wminusbeta <- sapply(1:K, FUN = function(k) {
    if(any(which.i[[k]] != 0)){
      return(sum(sapply(which.i[[k]], FUN = function(i) wMinusb(i=i, t=1, k=k, w.z=w.z, beta.k=beta.k, X=X[[i]]))))
    }else{
      return(0)
    }
  })
  
  
  # update alpha.0k
  v0k <- 1/(sig2inv.alpha + itime0)
  m0k <- v0k*(priors$mu.alpha*sig2inv.alpha + wminusbeta)
  alpha.0k <- rnorm(K, m0k, sqrt(v0k))
  
  #######################
  ### update alpha.jk ### 
  #######################
  
  alpha.jk <- lapply(1:(K), FUN = function(j){
    unlist(lapply(1:(K), FUN = function(k){
      updateAlphaJKdiff(j=j, k=k, n=n, t.max=t.max, z=z, vinv.alpha=vinv.alpha,
                        sig2inv.alpha = sig2inv.alpha, w.z = w.z, X = X, beta.k = beta.k,
                        m.alpha = m.alpha, mu.alpha = priors$mu.alpha)
    }))
  })
  
  
  if(anyNA(unlist(alpha.jk))) {
    stop("NA's in alpha.jk")
  }
  
  ######################
  ### update m.alpha ###
  ######################
  
  alpha.jj <- list()
  for(k in 1:(K)){
    alpha.jj[[k]] <- alpha.jk[[k]][k]
  }
  sumjj <- sum(unlist(alpha.jj))
  v.star <- 1/(K*vinv.alpha + 1/priors$v0) 
  m.star <- v.star*(sumjj*vinv.alpha + priors$m0/priors$v0)
  m.alpha <- rnorm(1, m.star, sqrt(v.star))
  
  ############################
  ### update sig2inv.alpha ###
  ############################
  
  alpha.jnotk <- list()
  for(k in 1:K){
    alpha.jnotk[[k]] <- alpha.jk[[k]][-k]
  }
  sumAlphajk <- sum((unlist(alpha.jnotk) - priors$mu.alpha)^2)
  sig2inv.alpha <- rgamma(1, priors$a1 + K*(K-1)/2, priors$b1 + .5*sumAlphajk) 
  
  #########################
  ### update vinv.alpha ###
  #########################
  
  sumjj2 <- sum((unlist(alpha.jj) - m.alpha)^2)
  vinv.alpha <- rgamma(1, priors$a2 + K/2, priors$b2 + sumjj2/2) 
  
  #####################
  ### update beta.k ### 
  #####################
  
  beta.k <- mclapply(1:(K), FUN = function(k){
    itimes <- lapply(1:n, FUN = function(i)  which(z[[i]] >= k))
    nk.tilde  <- length(unlist(itimes)) # total number we're looking at the dimension of everything
    if(nk.tilde > 0){
      w.k <- unlist(lapply(1:n, FUN = function(i) {
        if(any(itimes[[i]]>0)){
          sapply(itimes[[i]], FUN = function(t) w.z[[i]][[t]][k])
        }
      }))
      X.k <- numeric()
      for(i in 1:n){
        if(any(itimes[[i]]>0)){
          X.k <- rbind(X.k, X[[i]][itimes[[i]],])
        }
      }
      alpha.k <- list()
      for(i in 1:n){
        alphaTHIS <- numeric()
        if(any(itimes[[i]]>0)){
          for(j in z[[i]][itimes[[i]]-1]){
            alphaTHIS <- c(alphaTHIS, alpha.jk[[j]][k])
          }
          if(1 %in% itimes[[i]]){
            alphaTHIS <- c(alpha.0k[k], alphaTHIS)
          }
        }
        alpha.k[[i]] <- alphaTHIS
      }
      alpha.k <- unlist(alpha.k)
      if(length(alpha.k) != nk.tilde | length(w.k) != nk.tilde | nrow(X.k) != nk.tilde){
        stop("dimension of alpha, w, or X is wrong")
      }
      # update beta.k
      V.k <- chol2inv(chol(priors$SigInv.beta + crossprod(X.k, X.k)))
      m.k <- crossprod(V.k,crossprod(priors$SigInv.beta,priors$mu.beta) + crossprod(X.k,w.k - alpha.k))
      return(matrix(rmvn(n = 1, m.k, V.k), nrow = q))
    }else{ 
      # update from prior
      return(matrix(rmvn(n = 1, mu = priors$mu.beta, sigma = priors$Sigma.beta), nrow = q)) 
    }
  })
  if(anyNA(unlist(beta.k))) {
    stop("NA's in beta.k")
  }
  
  
  ###################
  ### update pi.z ### ### done 
  ###################
  
  ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)
  pi.z = updatePi2(beta = beta.k, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)
  
  #################################
  ### Sample New Missing Values ###
  #################################
  
  # Sample new MAR values conditional on observed data and imputed LOD data ###
  for(i in 1:n){
    if(any(mismat[[i]]==1)){ # MAR = 1 
      nummis <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==1))) # how many missing at each time point 
      for(t in which(nummis>0)){ # only loop through time points with missing data
        whichmis <- which(mismat[[i]][t,]==1) # which ones are missing 
        if(length(whichmis)==p){
          y[[i]][t,] <- rmvn(1, mu[[z[[i]][t]]], chol(Sigma[[z[[i]][t]]]), isChol = TRUE)
        }else{
          y.obs <- y[[i]][t,-whichmis]
          mu.obs <- mu[[z[[i]][t]]][,-whichmis]
          mu.miss <- mu[[z[[i]][t]]][,whichmis]
          Sigma.obs <- matrix(Sigma[[z[[i]][t]]][-whichmis, -whichmis], p-length(whichmis), p-length(whichmis))
          Sigma.miss <- matrix(Sigma[[z[[i]][t]]][whichmis, whichmis], length(whichmis), length(whichmis))
          Sigma.obs.miss <-  matrix(Sigma[[z[[i]][t]]][-whichmis, whichmis], p-length(whichmis), length(whichmis))
          Sigma.miss.obs <- t(Sigma.obs.miss)
          Sigma.mis.obs.inv <- Sigma.miss.obs%*%solve(Sigma.obs)
          mu.mgo <- as.numeric(mu.miss + Sigma.mis.obs.inv%*%(y.obs - mu.obs))
          Sigma.mgo <- Sigma.miss + Sigma.mis.obs.inv%*%Sigma.obs.miss
          y[[i]][t,whichmis] <- rmvn(1, mu.mgo, chol(Sigma.mgo), isChol = TRUE)
        }
      }
    }
  }
  
  # Sample new LOD values conditional on observed data and imputed MAR data ###
  for(i in 1:n){
    if(any(mismat[[i]]==2)){ # LOD = 2
      numlod <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==2))) # how many lod at each time point 
      for(t in which(numlod>0)){ # only loop through time points with missing data
        whichlod <- which(mismat[[i]][t,]==2) # which ones are below lod  
        if(length(whichlod)==p){ 
          y[[i]][t,] <- rtmvn(1, Mean = as.vector(mu[[z[[i]][t]]]), Sigma = Sigma[[z[[i]][t]]], lower = rep(-Inf, p),
                              upper = lod, int = y[[i]][t,], burn = 10, thin = 1)
        }else{
          y.obs <- y[[i]][t,-whichlod]
          mu.obs <- mu[[z[[i]][t]]][,-whichlod]
          mu.miss <- mu[[z[[i]][t]]][,whichlod]
          Sigma.obs <- matrix(Sigma[[z[[i]][t]]][-whichlod, -whichlod], p-length(whichlod), p-length(whichlod))
          Sigma.miss <- matrix(Sigma[[z[[i]][t]]][whichlod, whichlod], length(whichlod), length(whichlod))
          Sigma.obs.miss <-  matrix(Sigma[[z[[i]][t]]][-whichlod, whichlod], p-length(whichlod), length(whichlod))
          Sigma.miss.obs <- t(Sigma.obs.miss)
          Sigma.mis.obs.inv <- Sigma.miss.obs%*%chol2inv(chol(Sigma.obs))
          mu.mgo <- as.numeric(mu.miss + Sigma.mis.obs.inv%*%(y.obs - mu.obs))
          Sigma.mgo <- Sigma.miss + Sigma.mis.obs.inv%*%Sigma.obs.miss
          y[[i]][t,whichlod] <- rtmvn(1, Mean = mu.mgo, Sigma = Sigma.mgo, lower = rep(-Inf, length(whichlod)),
                                      upper = lod[whichlod], int = y[[i]][t, whichlod], burn = 10, thin = 1)
          
        }
      }
    }
  }
  
  #####################
  ### Error Measure ###
  #####################
  
  ## Hamming distance ##
  if(!is.null(unlist(z.true))){
    ham.error <- hamdist(unlist(z.true), unlist(z)) 
    ham[s] <- ham.error/(n*t.max) # proportion of misplaced states
  }else{
    ham <- NULL
  }    
  
  ## MSE for mu ##
  if(!is.null(mu.true)){
    sse <- list()
    for(i in 1:n){
      sse[[i]] <- sapply(1:t.max, FUN = function(t){
        as.numeric(crossprod(unlist(mu[z[[i]][t]]) - mu.true[z.true[[i]][t],]))/p
      })
    }
    mu.sse[s] <- sum(unlist(sse))
    mu.mse[s] <- mean(unlist(sse)) # vector mse for mu, divide by # exposures 
  }else{
    mu.sse <- NULL
    mu.mse <- NULL
  }
  
  
  #####################
  ### Store Results ###
  #####################
  
  z.save[[s]] <- z
  K.save[s] <- K 
  beta.save[[s]] <- beta.k
  mu.save[[s]] <- mu
  
  if(s%in%imputes){
    # imputed values for complete data sets 
    y.mar.save[s.imp,] <- unlist(y)[which(unlist(mismat)==1)] # mar imputations
    y.lod.save[s.imp,] <- unlist(y)[which(unlist(mismat)==2)] # lod imputations
    
    if(!is.null(ycomplete)){
      # MSE
      mar.mse[s.imp] <- mean((unlist(ycomplete)[which(unlist(mismat)==1)] - unlist(y)[which(unlist(mismat)==1)])^2)
      lod.mse[s.imp] <- mean((unlist(ycomplete)[which(unlist(mismat)==2)] - unlist(y)[which(unlist(mismat)==2)])^2)
      
      # SSE
      mar.sse[s.imp] <- sum((unlist(ycomplete)[which(unlist(mismat)==1)] - unlist(y)[which(unlist(mismat)==1)])^2)
      lod.sse[s.imp] <- sum((unlist(ycomplete)[which(unlist(mismat)==2)] - unlist(y)[which(unlist(mismat)==2)])^2)
      
      # mean bias
      mar.bias[s.imp] <- mean((unlist(y)[which(unlist(mismat)==1)] - unlist(ycomplete)[which(unlist(mismat)==1)]))
      lod.bias[s.imp] <- mean((unlist(y)[which(unlist(mismat)==2)] - unlist(ycomplete)[which(unlist(mismat)==2)]))
      
      # sum bias
      mar.sum.bias[s.imp] <- sum((unlist(y)[which(unlist(mismat)==1)] - unlist(ycomplete)[which(unlist(mismat)==1)]))
      lod.sum.bias[s.imp] <- sum((unlist(y)[which(unlist(mismat)==2)] - unlist(ycomplete)[which(unlist(mismat)==2)]))
    }
    
    s.imp <- s.imp+1
  }
  
  #end.time1 <- Sys.time()
  #print(s)
  #if(s %% 10 == 0) print(K)     
  #if(s %% 10 == 0) print(min(unlist(y)))
  #print(mhacc)
  #print(end.time1 - start.time1)
  
  
}
