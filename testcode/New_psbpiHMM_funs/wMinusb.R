



wMinusb <- function(i, t, k, w.z, beta.k, beta.sk=NULL, X){
  if(!is.null(beta.sk)){
    return(w.z[[i]][[t]][k] - crossprod(X[t,], beta.k[[k]]) - crossprod(X[t,], beta.sk[[k]]))
  }else{
    return(w.z[[i]][[t]][k] - crossprod(X[t,], beta.k[[k]]))
  }
}
