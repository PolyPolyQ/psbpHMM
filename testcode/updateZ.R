
detZminus1 <- function(i,state.list.i,pi.z,u.i,t.max){
  
  return(mclapply(2:t.max, FUN = function(t){
    pass1 <- state.list.i[[t-1]] # possible states for t-1
    pos.zt <- state.list.i[[t]] # possible states for t
    # need loop so we get NULL values where we need them 
    pos.ztminus1 <- list()
    for(k in pos.zt){
      pass2 <- which(pi.z[[t]][, k] >= u.i[t]) # second/backward condition
      pos.ztminus1[[k]] <- intersect(pass1, pass2) # these are the possible values of z_(t-1) given that z_t = k
    }
    return(list(pos.zt = pos.zt, pos.ztminus1 = pos.ztminus1))
  }))
  # pos.zt are possibles states at time t
  # pos.ztminus1 are the specific possible states at t-1 for each k in pos.zt
  # pos.ztminus1 tells me which states to sum over when I calculate prob for z[t] = k
} 


for(t in 2:t.max){
  pass1 = state.list[[i]][[t-1]] # possible states for t-1
  k.prime = state.list[[i]][[t]] # possibles states for t
  
}


updateZ <- function(i,y.i,mu,cholSigma,detR.stari,nu.df,K,log.stuff,t.max, detz, state.list.i){
  
  z.probs <- list()
  # calculate probability for the first time point 
  k.prime <- state.list.i[[1]] # possible states for first time point via slice sampling 
  log.tk.probs <- sapply(intersect(k.prime, 1:K), FUN = function(k) {
    dmvn(y.i[1,], mu[[k]], cholSigma[[k]], log = TRUE, isChol = TRUE) # current states
  })
  if((K+1)%in%k.prime){ # new state
    log.tk.probs <- c(log.tk.probs,log.stuff-((nu.df+1)/2)*log(detR.stari[1]))
  }
  tk.probs <- exp(log.tk.probs-logsum(log.tk.probs))
  z.probs[[1]] <-  data.frame(k.prime, tk.probs)
  
  # this has to be a loop because it depends on previous time
  for(t in 2:t.max){ 
    k.prime <- state.list.i[[t]]; k.prime
    j.prime <- z.probs[[t-1]]$k.prime; j.prime
    j.probs <- z.probs[[t-1]]$tk.probs
    # for each current state k that is possible, sum over the previous possible states and multiply by density of k 
    log.tk.probs <- sapply(k.prime, FUN = function(k){
      log.sumprobs <- log( sum( (j.prime %in% detz[[t-1]]$pos.ztminus1[[k]])*j.probs)) # sum over probs for previous time point 
      if(k <= K){
        loglik_k <- dmvn(y.i[t,], mu[[k]], cholSigma[[k]], log = TRUE, isChol = TRUE) # current state
      }else{
        loglik_k <- log.stuff-((nu.df+1)/2)*log(detR.stari[t])
      }
      return(loglik_k + log.sumprobs)
    })
    tk.probs <- exp(log.tk.probs - logsum(log.tk.probs))
    z.probs[[t]] <- data.frame(k.prime, tk.probs)
  }
  
  upz <- mclapply(1:t.max, FUN = function(t){
    k.prime <- z.probs[[t]]$k.prime
    if(length(k.prime) == 1){
      statez <- k.prime
    }else{
      statez <- sample(k.prime, 1, prob = z.probs[[t]]$tk.probs)
    }
    return(statez) 
  })
  return(unlist(upz))
}
