
# we can write this loop in C++ to replace BOTH detZminus1 and updateZ 

##### parameters that we need for this function are:
# List state.list
# List y 
# List mu 
# List cholSigma
# double log.stuff
# double nu.df 
# List detR.star
# List pi.z 
# List u

##### special functions 
# dmvn (manual)
# logsum (log is in rcpp)
# log (rcpp)
# cbind (join_rows or join_cols)
# sample (rcpp) (vector x, int size, replace = false, probs)
# intersect (rcpp)
# which 
# %in% 



#################################################################################
for(i in 1:n){
  
  
  ### first time point ###
  zprobs <- list()
  kprime <- state.list[[i]][[1]] # possible states for first time point via slice sampling 
  # log likelihood for first state 
  nprime = length(kprime)
  log.tk.probs = numeric(nprime)
  for(idx in 1:nprime){
    k = kprime[idx]
    if(k<=K){
      log.tk.probs[idx] = dmvn(y[[i]][1,], mu[[k]], cholSigma[[k]], log = TRUE, isChol = TRUE) 
    }else{
      log.tk.probs[idx] = log.stuff-((nu.df+1)/2)*log(detR.star[[i]][1])
    }
  }
  
  tk.probs <- exp(log.tk.probs-logsum(log.tk.probs)) # function logsum 
  #sum(tk.probs) # normalized 
  
  zprobs[[1]] <- cbind(kprime, tk.probs)
  
  # sample first state
  if(nprime == 1){
    z[[i]][t] = kprime
  }else{
    z[[i]][t] = sample(kprime, 1, prob = tk.probs)
  }
  ### first time point ### 
  
  ##################################### c++ up to here 
  
  
  # next state probabilities 
  for(t in 2:t.max){ 
    k.prime = state.list[[i]][[t]]; k.prime # current possible states 
    j.prime = zprobs[[t-1]][,1]; j.prime # previous possibles states 
    j.probs = zprobs[[t-1]][,2]; j.probs # previous state probabilities 
    
    # for each current possible state, which are the previous possible states 
    nprime = length(k.prime); nprime
    log.tk.probs = numeric(nprime)
    for(idx in 1:nprime){
      k=k.prime[idx]
      
      # four new lines
      # these are the previous possible states for current possible state k 
      prevK = intersect(which(pi.z[[i]][[t]][,k] >= u[[i]][t]),state.list[[i]][[t-1]]); prevK # intersect and which 
      # if(length(prevK) < length(j.prime)){
      #   print(t)
      #   print(k)
      #   stop()
      # }
      ints = which(j.prime %in% prevK); ints  
      alljprobs = j.probs[ints]; alljprobs
      log.sumprobs = log(sum(alljprobs))
      
      
      # log likelihood 
      if(k<=K){
        loglik_k = dmvn(y[[i]][t,], mu[[k]], cholSigma[[k]], log = TRUE, isChol = TRUE)
      }else{
        loglik_k = log.stuff-((nu.df+1)/2)*log(detR.star[[i]][t])
      }
      
      
      log.tk.probs[idx] = loglik_k + log.sumprobs
    }
    
    tk.probs = exp(log.tk.probs - logsum(log.tk.probs)) # logsum function 
    zprobs[[t]] = cbind(k.prime, tk.probs)
    
    
    # sample next state - same as sample for first time point 
    if(nprime == 1){
      z[[i]][t] = k.prime
    }else{
      z[[i]][t] = sample(k.prime, 1, prob = tk.probs)
    }
  }
}

z
zprobs[[288]]
#################################################################################



