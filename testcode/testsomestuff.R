compileAttributes()
devtools::build()
devtools::install()
library(psbpHMM)

#################### test rcpp code #######################
i=50
t=8
ajkmat = matrix(unlist(alpha.jk), nrow = K, ncol = K)
updateWrcpp(t=t, i=i, alpha0=alpha.0k, ajk = ajkmat, X=X, z=z, beta=beta.k, beta_sk=beta.sk)
  
############################################################


w.z <- list()
for(i in 1:n){
  w.z[[i]] <- list()
  for(t in 1:t.max){
    w.z[[i]][[t]] <- sapply(1:z[[i]][t], FUN = function(l) updateW(t=t, l=l, i=i, alpha.0k=alpha.0k, X=X[[i]], 
                                                                   beta.k=beta.k, beta.sk = beta.sk[[i]], 
                                                                   alpha.jk=alpha.jk, z=z))
  }
}

# write updateW in C++

x = X[[i]]
beta.ski = beta.sk[[i]]

updateW <- function(t,l,i, alpha.0k, x, beta.k, beta.sk = NULL, alpha.jk, z){
  
  if(!is.null(beta.sk)){
    if(t==1){
      trunc.mean <- as.numeric(alpha.0k[l] + crossprod(x[t,],beta.k[[l]]) + crossprod(x[t,],beta.ski[[l]])) 
    }else{
      trunc.mean <- as.numeric(alpha.jk[[z[[i]][t-1]]][l] + crossprod(X[t,],beta.k[[l]]) + crossprod(X[t,],beta.sk[[l]]))
    }
    return((z[[i]][t]==l)*rtruncnorm(1, a = 0, b = Inf, mean = trunc.mean, sd = 1) +
             (z[[i]][t]>l)*rtruncnorm(1, a = -Inf, b = 0, mean = trunc.mean, sd = 1))
    
  }else{
    if(t==1){
      trunc.mean <- as.numeric(alpha.0k[l] + crossprod(X[t,],beta.k[[l]])) 
    }else{
      trunc.mean <- as.numeric(alpha.jk[[z[[i]][t-1]]][l] + crossprod(X[t,],beta.k[[l]]))
    }
    return((z[[i]][t]==l)*rtruncnorm(1, a = 0, b = Inf, mean = trunc.mean, sd = 1) +
             (z[[i]][t]>l)*rtruncnorm(1, a = -Inf, b = 0, mean = trunc.mean, sd = 1))
  }
}
