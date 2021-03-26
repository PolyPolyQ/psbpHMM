


set.seed(1234)

# Rcpp
w.z <- list()
for(i in 1:n){
  w.z[[i]] <- list()
  for(t in 1:t.max){
    w.z[[i]][[t]] <- sapply(1:z[[i]][t], FUN = function(l) updateW(t=t, l=l, i=i, alpha.0k=alpha.0k, X=X[[i]], 
                                                                   beta.k=beta.k, beta.sk = beta.sk[[i]], 
                                                                   alpha.jk=alpha.jk, z=z))
  }
}

w.rcpp = w.z


set.seed(1234)

# R
w.z <- list()
for(i in 1:n){
  w.z[[i]] <- list()
  for(t in 1:t.max){
    w.z[[i]][[t]] <- sapply(1:z[[i]][t], FUN = function(l) upWsame(t=t, l=l, i=i, alpha.0k=alpha.0k, X=X[[i]], 
                                                                   beta.k=beta.k, alpha.jk=alpha.jk, z=z))
  }
}

w.r = w.z