

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
    back.list[[t]] <- sort(unique(unlist(lapply(intersect(back.list[[t+1]], 1:K ),
                                                FUN = function(k) which(pi.z[[i]][[t+1]][,k] >= u[[i]][t+1])))))
  }
  state.list[[i]] <- lapply(1:t.max, FUN = function(t) {
    return(intersect(for.list[[t]], back.list[[t]]))
  })
}



# // reset whichK
# whichKtemp = as<arma::vec>(whichK); 
# whichKtemp.set_size(0);
# whichK = wrap(whichKtemp); 
# 
# idx2 = 0; 
# for(kstar = 0; kstar < K; ++kstar){
#   if(pizitk[kstar] >= uit) {
#     whichK.push_back(kstar+1); // indexing 
#     idx2++; 
#   }
# }


for(t = 1; t < tmax; ++t){
  
  whichKminus1 = as<arma::vec>(forList[t-1]); // vector of previous states 
  whichKprev = wrap(whichKminus1);
  numKminus1 = whichKprev.length();
  uit = ui[t]; 
  
  // reset whichK
  whichKtemp = as<arma::vec>(whichK); 
  whichKtemp.set_size(0);
  whichK = wrap(whichKtemp); 
  
  pizit = as<arma::mat>(pizi[t]); // prob matrix for current time point 
  
  // loop through forList at t-1
  for(idx = 0; idx < numKminus1; ++idx){
    tempK = whichKprev[idx]; // one previous state at a time 
    pizitk = pizit.row(tempK).t(); // get the transition probs for the one previous state 
    for(k = 0; k < (K+1); ++k){
      if(pizitk[k] >= uit){
        whichK.push_back(k); 
      }
    }
  }
  
  // take the unique elements of whichK
  forList[t] = sort_unique(whichK); 
  
}

