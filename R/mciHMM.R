#' Fit covariate-dependent PSBP-iHMM to multiple time series 
#' 
#' 
#'
#' @param niter number of total iterations
#' @param nburn number of burn-in iterations
#' @param y list of time series data for each time series 
#' @param X list, matrix of covariates for each time series
#' @param rmlist integer vector identifying repeated time series for the same subject with the same number (e.g. c(1,1,2,2,2,3,3))
#' @param missing logical; if TRUE then the data set y contains missing data, default is FALSE
#' @param lod list of lower limits of detection for p exposures for each time series
#' @param len.imp number of imputations to save. Imputations will be taken at equally spaced iterations between nburn and niter. 
#' @param K.start starting number of hidden states, default is 12
#' @param priors list of priors
#' @param tau2 variance tuning parameter for normal proposal in MH update of lower triangular elements in decomposition of Sigma
#' @param a.tune shape tuning parameter for inverse gamma proposal in MH update of diagonal elements in decomposition of Sigma
#' @param b.tune rate tuning parameter for inverse gamma proposal in MH update of diagonal elements in decomposition of Sigma
#' @param resK logical; if TRUE a resolvent kernel is used in MH update for lower triangular elements in decomposition of Sigma
#' @param eta.star resolvent kernel parameter, must be a real value greater than 1. In the resolvent kernel we take a random draw from the geometric distribution with mean (1-p)/p, eta.star = 1/p.
#' @param z.true list of true hidden states, if known
#' @param mu.true matrix of true exposure means for each true state, if known 
#' @param ycomplete complete data, if available, for evaluating imputations
#' @param holdout list of indicators of missing type in holdout data set, 0 = observed, 1 = MAR, 2 = below LOD, for imputation validation purposes
#'
#' @importFrom parallel mclapply
#' @importFrom stats rnorm runif var rgamma kmeans rWishart cov cov2cor dnorm rgeom pnorm 
#' @importFrom tmvmixnorm rtmvn
#' @importFrom matrixcalc matrix.trace
#' @importFrom mvnfast rmvn dmvn
#' @importFrom invgamma dinvgamma
#' @importFrom gdata lowerTriangle<-
#'
#' @return an object of type "iHMM"
#'
#' @return a list with components
#' \itemize{
#'        \item z.save: list of estimated hidden states for each time series at each iteration
#'        \item K.save: list of estimated number of hidden states for each time series at each iteration
#'        \item ymar: matrix with len.imp rows of imputed values for MAR data
#'        \item ylod: matrix with len.imp rows of imputed values for data below LOD
#'        \item beta.save: list of posterior estimates of beta_k, state-specific regression coefficients in PSBP 
#'        \item gamma.save: list of posterior estimates of gamma_ik, state-specific subject-specific regression coefficients in PSBP
#'        \item mu.save: list of posterior estimates of mu_k, state-specific means 
#'        \item hamming: posterior hamming distance between true and estimated states, if z.true is given
#'        \item mu.mse: mean squared error for estimated state-specific means, if mu.true is given
#'        \item mu.sse: sum of squared errors for estimated state-specific means, if mu.true is given
#'        \item mar.mse: mean squared error of MAR imputations, if ycomplete is given 
#'        \item lod.mse: mean squared error of imputations below LOD, if ycomplete is given 
#'        \item mar.sse: sum of squared errors of MAR imputations, if ycomplete is given 
#'        \item lod.sse: sum of squared errors of imputations below LOD, if ycomplete is given 
#'        \item mar.sum.bias: sum of bias for MAR imputations 
#'        \item lod.sum.bias: sum of bias for below LOD imputations 
#'        \item mar.bias: mean bias for MAR imputations 
#'        \item lod.bias: mean bias for below LOD imputations 
#'        \item mismat: list, each element is a matrix indicating types of missing data for each time series, 0 = complete, 1 = MAR, 2 = below LOD
#'        \item ycomplete: list of complete data
#'        \item MH.arate: MH acceptance rate for lower triangular elements
#'        \item MH.lamrate: MH acceptance rate for diagonal elements 
#' }
#' @export
#'
mciHMM <- function(niter, nburn, y, X, rmlist=NULL, missing = FALSE, 
                   lod = NULL, len.imp = NULL, K.start = 12, priors = NULL,
                   tau2 = NULL, a.tune = NULL, b.tune = NULL, resK = FALSE, eta.star = NULL,
                   z.true = NULL, mu.true = NULL, ycomplete = NULL, holdout = NULL){


  # catch problems with parameter input 
  if(niter <= nburn) stop("niter must be greater than nburn")
  if(nburn < 0) stop("nburn must be greater than or equal to 0")

  if(!is.null(rmlist)){
    if(!is.numeric(rmlist)) stop("rmlist must be a numeric vector")
  }
  
  if(!is.numeric(unlist(y))) stop("y must be numeric")
  if(!is.numeric(unlist(X))) stop("X must be numeric")
  if(!is.numeric(unlist(ycomplete)) & !is.null(ycomplete)) stop("ycomplete must be numeric")

  if(!is.null(holdout) & any(!unique(unlist(holdout)) %in% c(0,1,2)) ){
    stop("holdout data set incorrectly specified, must contain only 0, 1, and 2.")
  }
  
  ### X is a list with a matrix for each i ### 
  if(missing){
    algorithm = "MH"
    SigmaPrior = "wishart"
  }else{
    algorithm = "Gibbs"
    SigmaPrior = "wishart"
  }
  
  if(is.null(rmlist)) beta.sk = NULL
  
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
  }else{
    stop("bad input: y must be a list, matrix, or vector of numeric data")
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
  if(is.null(priors$mu.beta)) priors$mu.beta <- rep(0, q) 
  if(is.null(priors$Sigma.beta)) priors$Sigma.beta <- diag(q) 
  if(!is.null(X)) priors$SigInv.beta <- invMat(priors$Sigma.beta) # cppFunction
  # subject specific beta if repeated measures 
  if(!is.null(rmlist)){
    # beta_sk
    if(is.null(priors$mu.gamma)) priors$mu.gamma <- rep(0, q) 
    if(is.null(priors$Sigma.gamma)) priors$Sigma.gamma <- diag(q) 
    if(!is.null(X)) priors$SigInv.betaS <- invMat(priors$Sigma.gamma) 
   # shrinkage prior on random effects
    if(is.null(priors$a.kappa)) priors$a.kappa <- 1 # shape for kap2inv 
    if(is.null(priors$b.kappa)) priors$b.kappa <- 1 # rate for kap2inv 
  }
  # hyperpriors on alpha
  if(is.null(priors$a1)) priors$a1 <- 1 # shape for sig2inv.alpha
  if(is.null(priors$b1)) priors$b1 <- 1 # rate for sig2inv.alpha
  if(is.null(priors$a2)) priors$a2 <- 1 # shape for vinv.alpha
  if(is.null(priors$b2)) priors$b2 <- 1 # rate for vinv.alpha
  
  # Sigma 
  if(missing) priors$R <- diag(p)
  if(is.null(priors$R)) priors$R <- diag(p) # Sigma_k ~ Inv.Wish(nu, R) hyperparameter for Sigma_k
  if(is.null(priors$nu)) priors$nu <- p+2 # Sigma_k ~ Inv.Wish(nu, R) hyperparameter for Sigma_k, nu > p+1
  nu.df <- priors$nu # 
  R.mat <- priors$R 
  
  if(is.null(priors$lambda)) priors$lambda <- 10 
  
  #############################
  ### Indicate Missing Type ###
  #############################
  
  # indicate missing data: obs = 0, mar = 1, lod = 2
  mismat <- list()
  for(i in 1:n){
    mismat[[i]] <- matrix(sapply(y[[i]], ismissing), ncol = p)
  }
  
  if(is.null(holdout)) holdout = mismat
  
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
      expLod <- exp(lod[[i]])
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
  
  if(is.null(K.start)) K.start = 12
  
  z <- list()
  for(i in 1:n){
    K <- K.start
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
      lams[[k]] <- 1/rgamma(7, vj0, rate = deltaj0)
      D[[k]] <- diag(lams[[k]])
      al.list <- list()
      for(j in 2:p){
        al.list[[j-1]] <- rnorm(j-1, 0, lams[[k]][j])
      }
      al[[k]] <- unlist(al.list) # for j = 2 to p
      which.lams <- unlist(sapply(2:p, FUN = function(j) rep(j,j-1))) # which lams to use for each al 
      L[[k]] <- diag(p)
      lowerTriangle(L[[k]]) <- al[[k]]
      Sigma[[k]] <- mhDecomp(L[[k]], D[[k]]) # cppFunction
      mu[[k]] <- rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
    }else{
      Sigma[[k]] <- chol2inv(chol(matrix(rWishart(1, df = nu.df, Sigma = invMat(R.mat)),p,p)))
      mu[[k]] <- rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
    }
  }
  
  
  # alpha 
  alpha.0k <- rep(0, K) # initial state intercept
  alpha.jk <- list() # state intercepts
  m.alpha <- priors$m0 # mean on alpha.jj
  sig2inv.alpha <- priors$a1/priors$a2 # precision on alpha.jk
  vinv.alpha <- priors$b1/priors$b2 # precision on alpha.jj
  
  # beta
  beta.k <- list() # regression coefficients
  for(k in 1:K){
    alpha.jk[[k]] <- rnorm(K, priors$mu.alpha, sqrt(1/sig2inv.alpha))
    alpha.jk[[k]][k] <- priors$m0
    beta.k[[k]] <- matrix(rmvn(1, priors$mu.beta, priors$Sigma.beta), nrow = q, ncol = 1) # length q
  }
  
  # beta_s random effects 
  if(!is.null(rmlist)){
    kap2inv <- 1
    n.sub = length(unique(rmlist)) # number of unique subjects
    beta.ik = list()
    for(i in 1:n.sub){
      beta.ik[[i]] = list()
      for(k in 1:K){
        beta.ik[[i]][[k]] <- matrix(rmvn(1, priors$mu.beta, (1/kap2inv)*priors$Sigma.beta), nrow = q, ncol = 1) # length q
      }
    }
    beta.sk = list()
    for(i in 1:n){
      beta.sk[[i]] = beta.ik[[rmlist[i]]]
    }
  }

  
  ################################
  ### update transition matrix ###
  ################################
  
  ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)
  
  if(!is.null(rmlist)){
    pi.z = updatePi_rm(beta = beta.k, beta_sk = beta.sk, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)
  }else{
    pi.z = updatePi(beta = beta.k, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)
  }
  
  ################################################
  ### fixed values for new state probabilities ###
  ################################################
  
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
  beta.save <- list()
  beta_rm.save <- list() # check bottom for this 
  K.save <- numeric()
  ham <- numeric()
  mu.mse <- numeric()
  mu.sse <- numeric()
  MH.a <- 0 # for MH update a
  MH.lam <- 0 # for MH update lams
  s.save = 1
  
  # missing data sets
  if(!is.null(len.imp)){
    imputes <- ceiling(seq.int(nburn+1, niter, length.out = len.imp))
    y.mar.save <- matrix(NA, len.imp, length(which(unlist(mismat)==1)))
    y.lod.save <- matrix(NA, len.imp, length(which(unlist(mismat)==2)))
    mar.mse <- numeric()
    lod.mse <- numeric()
    mar.sse <- numeric()
    lod.sse <- numeric()
    miss.mse <- numeric()
    mar.bias <- numeric()
    lod.bias <- numeric()
    mar.sum.bias <- numeric()
    lod.sum.bias <- numeric()
    s.imp <- 1
  }else{
    imputes = 0
    y.mar.save <- NULL
    y.lod.save <- NULL
    mar.mse <- NULL
    lod.mse <- NULL
    mar.sse <- NULL
    lod.sse <- NULL
    miss.mse <- NULL
    mar.bias <- NULL
    lod.bias <- NULL
    mar.sum.bias <- NULL
    lod.sum.bias <- NULL
    s.imp <- NULL
  }

  ###############
  ### Sampler ###
  ###############
  
  #start.time = Sys.time()
  for(s in 1:niter){
    
    if (s%%100==0) print(paste("iteration", s, " number of states =", K))

    z.prev <- list()
    z.prev <- mclapply(1:n, FUN=function(i) return(z[[i]]))
    
    ######################
    ### update theta_k ### 
    ######################
    
    
    cholSigma <- lapply(1:K, FUN = function(k) chol(Sigma[[k]]))
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
        Sigma[[k]] <- chol2inv(chol(matrix(rWishart(1,df=nu_nk, Sigma=invMat(Sigma_nk)),p,p)))
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
              a.star <- rnorm(1, 0, sqrt(tau2)); a.star # proposed value 
              al.star[j] <- a.star; al.star
              lowerTriangle(L.star) <- al.star; L.star
              SigmaStar <- mhDecomp(L.star, D[[k]]) # cppFunction

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
    
    ################
    ### update u ### ### Rcpp
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
    
    ############################
    ### update State List: R ###
    ############################
    
    state.list = lapply(1:n, FUN = function(i){
      return(upStateList_lapply(i, u=u, pi.z = pi.z, K=K, t.max = t.max))
    })
    
    ######################
    ### Sample Z: Rcpp ### 
    ######################
    
    # why is this so slow in Rcpp?
    ## new code with updated handling of NAs in probVec
    z1 <- upZ(stateList = state.list, y = y, mu = mu, Sigma = Sigma, logStuff = log.stuff,
              nudf = nu.df, detRstar = detR.star, piz = pi.z, u = u, tmax = t.max,
              K = K, n = n, d = p)
  
    
    z <- lapply(1:n, FUN = function(i){
      return(as.numeric(z1[[i]]))
    })

    #########################
    ### new state fillers ### 
    #########################
    
    if(any(unlist(z)>K)){
      
      if(algorithm == "MH"){
        lamsNew <- 1/rgamma(7, vj0, rate = deltaj0); lamsNew
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
        SigmaNew <- chol2inv(chol(matrix(rWishart(1, df = nu.df, Sigma = invMat(R.mat)),p,p)))
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
      
      # repeated measures for beta 
      if(!is.null(rmlist)){
        beta.K.new = list()
        for(i in 1:n.sub){
          beta.K.new[[i]] = matrix(rmvn(1, priors$mu.gamma, (1/kap2inv)*priors$Sigma.gamma), nrow = q, ncol = 1)
        }
        for(i in 1:n){
          beta.sk[[i]][[K+1]] = beta.K.new[[rmlist[i]]]
        }
      }
  
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
    
    if(!is.null(rmlist)){
      betaS.new <- list()
      for(i in 1:n){
        betaS.new[[i]] = list()
      }
    }
    
    
    rj <- 1
    for(k in zu){ # 
      alpha.new[[rj]] <- alpha.jk[[k]][zu]
      beta.new[[rj]] <- beta.k[[k]]
      mu.new[[rj]] <- mu[[k]]
      Sigma.new[[rj]] <- Sigma[[k]]
      if(!is.null(rmlist)){
        for(i in 1:n){
          betaS.new[[i]][[rj]] = beta.sk[[i]][[k]]
        }
      }
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
    if(!is.null(rmlist)) beta.sk <- betaS.new
    
    z <- z.new
    Sigma <- Sigma.new
    mu <- mu.new
    if(algorithm == "MH"){
      al <- al.new
      lams <- lams.new 
      D <- D.new
      L <- L.new
    }

    ################
    ### update W ###
    ################
    
    # for each t, w.z gives me a VECTOR based on the previous time point and values up to the current time point
    # so each w.z[[t]] should be a VECTOR of length z_t
    ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)
    if(!is.null(rmlist)){
      w.z = upW_rm(alpha0 = alpha.0k, X = X, beta = beta.k, beta_rm = beta.sk, ajk = ajkmat, z=z, tmax = t.max, n=n)
    }else{
      w.z = upW(alpha0 = alpha.0k, X = X, beta = beta.k, ajk = ajkmat, z=z, tmax = t.max, n=n)
    }

    #######################
    ### update alpha.0k ### 
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
        return(sum(sapply(which.i[[k]], FUN = function(i) wMinusb(i=i, t=1, k=k, w.z=w.z, beta.k=beta.k, beta.sk = beta.sk[[i]], X=X[[i]]))))
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
    
    
    if(is.null(rmlist)){
      alpha.jk = up_ajk(K=K,n=n, tmax = t.max, z = z, vinv_alpha = vinv.alpha, sig2inv_alpha = sig2inv.alpha,
                        w = w.z, X = X, beta_k = beta.k, m_alpha = m.alpha, mu_alpha = priors$mu.alpha)
    }else{
      alpha.jk = up_ajk_rm(K=K,n=n, tmax = t.max, z = z, vinv_alpha = vinv.alpha, sig2inv_alpha = sig2inv.alpha,
                        w = w.z, X = X, beta_k = beta.k, beta_sk = beta.sk, m_alpha = m.alpha, mu_alpha = priors$mu.alpha)
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
        
        if(!is.null(rmlist)){
          Xibetak = numeric()
          for(i in 1:n){
            if(any(itimes[[i]])>0){
              Xibetak = rbind(Xibetak, X[[i]][itimes[[i]],] %*% beta.sk[[i]][[k]])
            }
          }
          V.k <- invMat(priors$SigInv.beta + crossprod(X.k, X.k))
          m.k <- crossprod(V.k,crossprod(priors$SigInv.beta,priors$mu.beta) + crossprod(X.k,w.k - alpha.k - Xibetak))
        }else{
          V.k <- invMat(priors$SigInv.beta + crossprod(X.k, X.k))
          m.k <- crossprod(V.k,crossprod(priors$SigInv.beta,priors$mu.beta) + crossprod(X.k,w.k - alpha.k))
        }
        return(matrix(rmvn(n = 1, m.k, V.k), nrow = q))
      }else{ 
        # update from prior
        return(matrix(rmvn(n = 1, mu = priors$mu.beta, sigma = priors$Sigma.beta), nrow = q)) 
      }
    })
    
    ############################################
    ### update beta.sk for repeated measures ###
    ############################################

    if(!is.null(rmlist)){
      beta.ik = list()
      for(i.sub in 1:n.sub){
        subs = which(rmlist == i.sub) # subjects under consideration 
        num.subs = length(subs)
        
        beta.ik[[i.sub]] <- mclapply(1:K, FUN = function(k){
          itimes <- lapply(subs, FUN = function(i)  which(z[[i]] >= k)) # only subject i.sub
          nk.tilde  <- length(unlist(itimes)) # total number of time points, dimension of everything that follows
          if(nk.tilde > 0){
            w.k <- unlist(lapply(1:num.subs, FUN = function(i) {
              if(any(itimes[[i]]>0)){
                sapply(itimes[[i]], FUN = function(t) w.z[[subs[i]]][[t]][k])
              }
            }))
            X.k <- numeric()
            for(i in 1:num.subs){
              if(any(itimes[[i]]>0)){
                X.k <- rbind(X.k, X[[subs[i]]][itimes[[i]],])
              }
            }
            alpha.k <- list()
            for(i in 1:num.subs){
              alphaTHIS <- numeric()
              if(any(itimes[[i]]>0)){
                for(j in z[[subs[i]]][itimes[[i]]-1]){
                  alphaTHIS <- c(alphaTHIS, alpha.jk[[j]][k])
                }
                if(1 %in% itimes[[i]]){
                  alphaTHIS <- c(alpha.0k[k], alphaTHIS)
                }
              }
              alpha.k[[i]] <- alphaTHIS
            }
            alpha.k <- unlist(alpha.k)

            # update beta.ik
            V.k <- chol2inv(chol(kap2inv*priors$SigInv.betaS + crossprod(X.k, X.k)))
            m.k <- crossprod(V.k,crossprod( (kap2inv*priors$SigInv.betaS) ,priors$mu.gamma) + crossprod(X.k,w.k - alpha.k - crossprod(t(X.k), beta.k[[k]])))
            return(matrix(rmvn(n = 1, m.k, V.k), nrow = q))
          }else{ 
            # update from prior
            return(matrix(rmvn(n = 1, mu = priors$mu.gamma, sigma = (1/kap2inv)*priors$Sigma.gamma), nrow = q)) 
          }
        })
      }

      # relist from 1:n
      beta.sk = list()
      for(i in 1:n){
        beta.sk[[i]] = beta.ik[[rmlist[i]]]
      }
      
      # update kap2inv 
      ssbetaik <- 0 
      for(i in 1:n.sub){
        for(k in 1:K){
          ssbetaik = ssbetaik + t(beta.ik[[i]][[k]] - priors$mu.gamma)%*%priors$SigInv.betaS%*%(beta.ik[[i]][[k]] - priors$mu.gamma)
        }
      }
      a.kap <- priors$a.kappa + n.sub*K/2
      b.kap <- priors$b.kappa + (1/2)*ssbetaik
      kap2inv <- rgamma(1, a.kap, b.kap) 
      
    }
   
    ###################
    ### update pi.z ### 
    ###################
    
    ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)
    
    # why is this so slow in Rcpp?
    if(!is.null(rmlist)){
      pi.z = updatePi_rm(beta = beta.k, beta_sk = beta.sk, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)
    }else{
      pi.z = updatePi(beta = beta.k, X = X, a0 = alpha.0k, ajk = ajkmat, tmax = t.max)
    }

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
                                upper = lod[[i]], int = y[[i]][t,], burn = 10, thin = 1)
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
                                        upper = lod[[i]][whichlod], int = y[[i]][t, whichlod], burn = 10, thin = 1)
            
          }
        }
      }
    }
    
    #####################
    ### Store Results ###
    #####################
    if(s > nburn){
      
      ## Hamming distance ##
      if(!is.null(unlist(z.true))){
        ham.error <- hamdist(unlist(z.true), unlist(z)) 
        ham[s.save] <- ham.error/(n*t.max) # proportion of misplaced states
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
        mu.sse[s.save] <- sum(unlist(sse))
        mu.mse[s.save] <- mean(unlist(sse)) 
      }else{
        mu.sse <- NULL
        mu.mse <- NULL
      }
      
      z.save[[s.save]] <- z
      K.save[s.save] <- K 
      beta.save[[s.save]] <- beta.k
      beta_rm.save[[s.save]] <- beta.sk
      mu.save[[s.save]] <- mu
      
      if(s%in%imputes){
        # imputed values for complete data sets 
        y.mar.save[s.imp,] <- unlist(y)[which(unlist(mismat)==1)] # mar imputations
        y.lod.save[s.imp,] <- unlist(y)[which(unlist(mismat)==2)] # lod imputations
        
        if(!is.null(ycomplete)){
          
          # holdout = mismat if not given 
          # MSE
          mar.mse[s.imp] <- mean((unlist(ycomplete)[which(unlist(holdout)==1)] - unlist(y)[which(unlist(holdout)==1)])^2)
          lod.mse[s.imp] <- mean((unlist(ycomplete)[which(unlist(holdout)==2)] - unlist(y)[which(unlist(holdout)==2)])^2)
          
          # SSE
          mar.sse[s.imp] <- sum((unlist(ycomplete)[which(unlist(holdout)==1)] - unlist(y)[which(unlist(holdout)==1)])^2)
          lod.sse[s.imp] <- sum((unlist(ycomplete)[which(unlist(holdout)==2)] - unlist(y)[which(unlist(holdout)==2)])^2)
          
          # mean bias
          mar.bias[s.imp] <- mean((unlist(y)[which(unlist(holdout)==1)] - unlist(ycomplete)[which(unlist(holdout)==1)]))
          lod.bias[s.imp] <- mean((unlist(y)[which(unlist(holdout)==2)] - unlist(ycomplete)[which(unlist(holdout)==2)]))
          
          # sum bias
          mar.sum.bias[s.imp] <- sum((unlist(y)[which(unlist(holdout)==1)] - unlist(ycomplete)[which(unlist(holdout)==1)]))
          lod.sum.bias[s.imp] <- sum((unlist(y)[which(unlist(holdout)==2)] - unlist(ycomplete)[which(unlist(holdout)==2)]))
        }
        
        s.imp <- s.imp+1
      }
      s.save = s.save + 1
    }
  }

list1 <- list(z.save = z.save, K.save = K.save,
              ymar = y.mar.save, ylod = y.lod.save,
              beta.save = beta.save,
              gamma.save = beta_rm.save,
              mu.save = mu.save,
              hamming = ham, mu.mse = mu.mse, 
              mu.sse = mu.sse,
              mar.mse = mar.mse, lod.mse = lod.mse,
              mar.sse = mar.sse, lod.sse = lod.sse,
              mar.sum.bias = mar.sum.bias, lod.sum.bias = lod.sum.bias,
              mar.bias = mar.bias, lod.bias = lod.bias,
              mismat = mismat, ycomplete = ycomplete,
              MH.arate = MH.a/(length(al)*sum(K.save)),
              MH.lamrate = MH.lam/(p*sum(K.save)))

class(list1) <- "ihmm"
return(list1)


}

