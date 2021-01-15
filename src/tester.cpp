
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

// [[Rcpp::export]]
NumericVector addTwo(NumericVector x) {
  return x + 2;
}

// [[Rcpp::export]]
NumericVector minusTwo(NumericVector x) {
  return x - 2;
}

// [[Rcpp::export]]
arma::mat a1(arma::mat x) {
  return(x) ;
}

// [[Rcpp::export]]
arma::vec a2(arma::vec x) {
  return(x) ;
}

// [[Rcpp::export]]
List a8(int n, int r, double v) {
  arma::mat x1 ;
  x1.print() ;
  x1.reshape(n, r) ;
  x1.fill(v) ;
  arma::mat x2(n, r) ;
  x2.fill(v) ;
  arma::mat x3 = x2 ;
  x3.reshape(n, r) ;
  List ret ;
  ret["x1"] = x1 ;
  ret["x2"] = x2 ;
  ret["x3"] = x3 ;
  return(ret) ;
}

// [[Rcpp::export]]
arma::mat mhDecomp(arma::mat L, arma::mat D) {
  return( inv(L) * D * inv(L).t() );
}

// [[Rcpp::export]]
arma::mat invMat(arma::mat x) {
  return( inv(x));
}

// [[Rcpp::export]]
List updatePi(List beta, arma::mat X, arma::vec a0, arma::mat ajk, int tmax){
  
  //step by step return each piece 
  int K = beta.size();
  arma::vec first(K);
  arma::vec second(K);
  List pimat(tmax); 
  
  for(int t = 0; t < tmax; ++t){
    arma::mat x = X.row(t);
    if(t == 0){
      for(int k = 0; k < K; ++k){
        arma::mat xb = as<arma::rowvec>(beta[k])*x.t();
        arma::mat y = a0[k] + xb; 
        first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
        if (k > 0) second[k] = 1-first[k-1];
      }
      second[0] = 1;
      arma::vec third = cumprod(second);
      arma::vec fourth(K+1);
      double s = 0;
      for(int i = 0; i < K; i++){
        fourth[i] = first[i]*third[i];
        s = s + fourth[i];
      }
      //double last = sum(fourth);
      fourth[K]=1-s; // index starts at 0 
      pimat[t] = fourth;
    }else {
      arma::mat Z;
      for(int j = 0; j < K; ++j){
        for(int k = 0; k < K; ++k){
          arma::mat xb = as<arma::rowvec>(beta[k])*x.t();
          double a = ajk(k,j);
          arma::mat y = a + xb; 
          first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
          if (k > 0) second[k] = 1-first[k-1];
        }
        second[0] = 1;
        arma::vec third = cumprod(second);
        arma::vec fourth(K+1);
        double s = 0;
        for(int i = 0; i < K; i++){
          fourth[i] = first[i]*third[i];
          s = s + fourth[i];
        }
        //double last = sum(fourth);
        fourth[K]=1-s; // index starts a 0 
        // fourth is the jth row
        Z = join_rows(Z, fourth);
        pimat[t] = Z.t();
      }
    }
  }
  return(pimat);
}


// [[Rcpp::export]]
List updatePi2(List beta, List X, arma::vec a0, arma::mat ajk, int tmax){
  
  // X is a List for each i 
  int n = X.size(); // number of people
  int K = beta.size();
  arma::vec first(K);
  arma::vec second(K);
  
  List probMat(n); // each element is a T-length list of matrices 
  
  for(int ind = 0; ind < n; ++ind){
    
    List pimat(tmax);
    
    arma::mat xi = as<arma::mat>(X[ind]);
    for(int t = 0; t < tmax; ++t){
      arma::mat x = xi.row(t);
      if(t == 0){
        for(int k = 0; k < K; ++k){
          arma::mat xb = as<arma::rowvec>(beta[k])*x.t();
          arma::mat y = a0[k] + xb; 
          first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
          if (k > 0) second[k] = 1-first[k-1];
        }
        second[0] = 1;
        arma::vec third = cumprod(second);
        arma::vec fourth(K+1);
        double s = 0;
        for(int i = 0; i < K; i++){
          fourth[i] = first[i]*third[i];
          s = s + fourth[i];
        }
        //double last = sum(fourth);
        fourth[K]=1-s; // index starts at 0 
        
        pimat[t] = fourth;
      }else {
        arma::mat Z;
        for(int j = 0; j < K; ++j){
          for(int k = 0; k < K; ++k){
            arma::mat xb = as<arma::rowvec>(beta[k])*x.t();
            double a = ajk(k,j);
            arma::mat y = a + xb; 
            first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
            if (k > 0) second[k] = 1-first[k-1];
          }
          second[0] = 1;
          arma::vec third = cumprod(second);
          arma::vec fourth(K+1);
          double s = 0;
          for(int i = 0; i < K; i++){
            fourth[i] = first[i]*third[i];
            s = s + fourth[i];
          }
          //double last = sum(fourth);
          fourth[K]=1-s; // index starts a 0 
          // fourth is the jth row
          Z = join_rows(Z, fourth);
          pimat[t] = Z.t();
        }
      }
    }
    probMat[ind] = pimat;
  }
  return(probMat);
}



// [[Rcpp::export]]
double returnPi(){
  double x = M_PI;
  return(x);
}

// [[Rcpp::export]]
double mvndensity(arma::vec y, arma::vec mu, arma::mat Sigma, double d){
  
  // det(Sigma)^(-1/2)
  arma::vec x1; 
  double x2 = det(Sigma); 
  x1 = x2; 
  arma::vec x3 = pow( x1, -0.5 );
  // dmnvorm 
  arma::vec x4 = pow( (2*M_PI), (-d/2))*x3*exp(-0.5 * (y - mu).t() * inv(Sigma) * (y - mu));
  
  double x5 = as_scalar(x4);
  
  
  return(x5);
  
}



// [[Rcpp::export]]
List upZ(List stateList, List y, List mu, List Sigma, double logStuff, double nudf, List detRstar, List piz,
         List u, int tmax, int K, int n, double d){
  
  // initialization 
  List z(1); // set to n 
  List sli; // state list at i 
  arma::mat yi; // y at i 
  arma::vec yit; // y_i at t 
  IntegerVector kprime; 
  int nprime;  
  arma::vec logTKprobs; 
  arma::vec tkProbs;
  double maxlog; 
  double logsum; 
  IntegerVector jprime; 
  arma::vec jprobs; 
  NumericVector probVec; 
  IntegerVector ret; 
  int state; 
  arma::vec zt;
  int idx; 
  int k; 
  arma::mat sigK;
  arma::vec muK;
  arma::vec temp1; 
  double temp2; 
  arma::vec temp3; 
  arma::vec temp4; 
  double temp5; 
  arma::vec detRi; 
  double logR; 
  double newState; 
  double loglik; 
  double uit; 
  
  // for i in 1:n
  int i = 0; 
  sli = stateList[i];
  yi = as<arma::mat>(y[i]);
  zt.set_size(tmax); // change to tmax, put zt in z(i), return z 

  // first time point 
  int t = 0; 
  yit = (yi.row(t)).t(); // redefined for each t 
  kprime = sli[t];
  nprime = kprime.length(); 
  logTKprobs.set_size(nprime); 
  for(idx = 0; idx < nprime; ++idx){
    k = kprime[idx]-1; // because the indexing starts at 0 for mu and Sigma 
    if(k < K){
      sigK = as<arma::mat>(Sigma[k]);
       muK = as<arma::mat>(mu[k]).t();
      temp2 = det(sigK); 
      temp1 = temp2; 
      temp3 = pow( temp1, -0.5 );
      temp4 = pow( (2*M_PI), (-d/2)) * temp3 * exp(-0.5 * (yit - muK).t() * inv(sigK) * (yit - muK));
      temp5 = log(as_scalar(temp4));
      logTKprobs[idx] = temp5;
    }else{
      detRi = as<arma::vec>(detRstar[i]);
      logR = detRi[t]; 
      newState = logStuff - ((nudf + 1)/2)*logR;
      logTKprobs[idx] = newState;
    }
  }
  maxlog = max(logTKprobs);
  logsum = log(sum(exp(logTKprobs - maxlog))) + maxlog;
  tkProbs = exp(logTKprobs - logsum); // becomes jprobs
  // sample t = 0 
  if(nprime == 1){
    zt[t] = kprime[0]; // only one element in kprime
  }else{
    probVec = wrap(tkProbs);
    ret = RcppArmadillo::sample(kprime, 1, 0, probVec);
    state = ret[0]; 
    zt[t] = state;
  }
  
  // the rest of the time points 
  for(t = 1; t < tmax; ++t){
    yit = (yi.row(t)).t(); // check this 
    jprime = kprime; // previous possible states
    jprobs = tkProbs; // previous state probabilities 
    kprime = sli[t]; // update kprime
    nprime = kprime.length(); 
    logTKprobs.set_size(nprime); 
    
    for(idx = 0; idx < nprime; ++idx){
      k = kprime[idx]-1; // because the indexing starts at 0 for mu and Sigma 
      if(k < K){
        sigK = as<arma::mat>(Sigma[k]);
        muK = as<arma::mat>(mu[k]).t();
        temp2 = det(sigK); 
        temp1 = temp2; 
        temp3 = pow( temp1, -0.5 );
        temp4 = pow( (2*M_PI), (-d/2)) * temp3 * exp(-0.5 * (yit - muK).t() * inv(sigK) * (yit - muK));
        loglik = log(as_scalar(temp4)); // loglik 
      }else{
        detRi = as<arma::vec>(detRstar[i]);
        logR = detRi[t]; 
        loglik = logStuff - ((nudf + 1)/2)*logR; // loglik
      }
      
      // calculate logsumprobs
      // these are the previous possible states for current possible state k 
      
      
      // write the following four lines in Rcpp
      
      // 1) prevK = intersect(which(pi.z[[i]][[t]][,k] >= u[[i]][t]),state.list[[i]][[t-1]]) # intersect and which 
      // 2) ints = which(prevK %in% j.prime) # which and %in% "sugar"?
      // 3) alljprobs = j.probs[ints]
      // 4) logsumprobs = log(sum(alljprobs))
      
      
      
      logTKprobs[idx] = loglik;
    }
    maxlog = max(logTKprobs);
    logsum = log(sum(exp(logTKprobs - maxlog))) + maxlog;
    tkProbs = exp(logTKprobs - logsum); // becomes jprobs
    // sample t
    if(nprime == 1){
      zt[t] = kprime[0]; // only one element in kprime
    }else{
      probVec = wrap(tkProbs);
      ret = RcppArmadillo::sample(kprime, 1, 0, probVec);
      state = ret[0]; 
      zt[t] = state;
    }
  }
  
  z[i] = zt; 
  return(z); 
}


// [[Rcpp::export]]
NumericVector csample_num( NumericVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()) {
  NumericVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return(ret);
}

// [[Rcpp::export]]
IntegerVector csample_int( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()) {
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return(ret);
}


// [[Rcpp::export]]
arma::mat a3(NumericMatrix x) {
  arma::mat y = as<arma::mat>(x) ;
  return(y) ;
}

// [[Rcpp::export]]
NumericMatrix a4(arma::mat x) {
  NumericMatrix y = wrap(x) ;
  return(y) ;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
