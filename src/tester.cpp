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
List updatePi(List beta, List X, arma::vec a0, arma::mat ajk, int tmax){
  
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
List updatePi_rm(List beta, List beta_sk, List X, arma::vec a0, arma::mat ajk, int tmax){
  
  // X is a List for each i 
  int n = X.size(); // number of people
  int K = beta.size();
  arma::vec first(K);
  arma::vec second(K);
  arma::mat xi; 
  arma::mat x; 
  arma::mat xb; 
  arma::mat xbs; 
  List bs; 
  arma::mat y; 
  
  List probMat(n); // each element is a T-length list of matrices 
  
  for(int ind = 0; ind < n; ++ind){
    List pimat(tmax);
    xi = as<arma::mat>(X[ind]);
    bs = beta_sk[ind]; 
    
    for(int t = 0; t < tmax; ++t){
      x = xi.row(t);
      if(t == 0){
        for(int k = 0; k < K; ++k){
          xb = as<arma::rowvec>(beta[k])*x.t();
          xbs = as<arma::rowvec>(bs[k])*x.t();
          y = a0[k] + xb + xbs; 
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
            xb = as<arma::rowvec>(beta[k])*x.t();
            xbs = as<arma::rowvec>(bs[k])*x.t();
            double a = ajk(k,j);
            y = a + xb + xbs; 
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
  arma::vec x1; 
  double x2 = det(Sigma); 
  x1 = x2; 
  arma::vec x3 = pow( x1, -0.5 );
  arma::vec x4 = pow( (2*M_PI), (-d/2))*x3*exp(-0.5 * (y - mu).t() * inv(Sigma) * (y - mu));
  double x5 = as_scalar(x4);
  return(x5);
}

// [[Rcpp::export]]
List upZ(List stateList, List y, List mu, List Sigma, double logStuff, double nudf, List detRstar, List piz,
         List u, int tmax, int K, int n, double d){
  
  // initialization 
  List z(n); // set to n 
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
  int s1; 
  int jindex; 
  int s2; 
  double probj;
  List pizi; 
  arma::mat pizit; 
  arma::vec pizitk; 
  arma::vec ui; 
  double uit; 
  IntegerVector slit;
  IntegerVector prevK; 
  IntegerVector whichK(0);
  arma::vec whichKtemp;
  NumericVector allProbs(0);
  arma::vec allProbstemp; 
  double logsumall; 
  int idx2; 
  int kstar; 
  
  for(int i = 0; i < n; ++i){
    sli = stateList[i];
    yi = as<arma::mat>(y[i]);
    zt.set_size(tmax); 
    pizi = piz[i]; 
    ui = as<arma::vec>(u[i]);
    
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
      pizit = as<arma::mat>(pizi[t]); 
      uit = ui[t]; 
      yit = (yi.row(t)).t(); // y_it 
      jprime = kprime; // previous possible states
      jprobs = tkProbs; // previous state probabilities 
      kprime = sli[t]; // update kprime, current possible states
      nprime = kprime.length(); 
      logTKprobs.set_size(nprime); // reset 
      slit = sli[t-1];
      
      // log likelihood for all current possible states kprime 
      for(idx = 0; idx < nprime; ++idx){
        // reset whichK
        whichKtemp = as<arma::vec>(whichK); 
        whichKtemp.set_size(0);
        whichK = wrap(whichKtemp); 
        // reset allProbs 
        allProbstemp = as<arma::vec>(allProbs); 
        allProbstemp.set_size(0); 
        allProbs = wrap(allProbstemp); 
        // for each k we sum over previous possible states for logsumprobs 
        k = kprime[idx]; 
        // because the indexing starts at 0 for mu and Sigma 
        if(k <= K){
          sigK = as<arma::mat>(Sigma[k-1]);
          muK = as<arma::mat>(mu[k-1]).t();
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
        pizitk = pizit.col(k-1); 
        idx2 = 0; 
        for(kstar = 0; kstar < K; ++kstar){
          if(pizitk[kstar] >= uit) {
            whichK.push_back(kstar+1); // indexing 
            idx2++; 
          }
        }
        prevK = intersect(whichK, slit);
        for(s1 = 0; s1 < jprime.length(); s1++){
          jindex = jprime[s1];
          for(s2 = 0; s2 < prevK.length(); s2++){
            if(prevK[s2]==jindex){
              probj = jprobs[s1];
              allProbs.push_back(probj);
            }
          }
        }
        logsumall = log(sum(allProbs)); 
        logTKprobs[idx] = loglik + logsumall;
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
  }
  return(z); 
}

// [[Rcpp::export]]
List upStateList(List piz, List u, int K, int tmax){
  
  List stateList(tmax); 
  List forList(tmax); 
  List backList(tmax);
  List pizi;
  arma::vec ui;
  arma::mat pizit; 
  double uit;
  int k;
  int t; 
  int i; 
  arma::vec whichKtemp; 
  IntegerVector whichK(0);
  int idx; 
  arma::vec whichKminus1; 
  int tempK; 
  arma::vec pizitk;
  IntegerVector whichKprev; 
  arma::vec pizitk0; 
  double uitplus1; 
  arma::vec whichKplus1; 
  IntegerVector whichKnext; 
  arma::mat pizitplus1; 
  arma::vec f1; 
  arma::vec f2; 
  IntegerVector f3; 
  IntegerVector f4; 
  
  
  // loop through individuals 
  i = 0; 
  ui = as<arma::vec>(u[i]);
  pizi = piz[i]; // list
  
  // first time point 
  t = 0; 
  pizitk0 = as<arma::vec>(pizi[t]); 
  uit = ui[t]; 
  whichKtemp = as<arma::vec>(whichK); 
  whichKtemp.set_size(0);
  whichK = wrap(whichKtemp); 
  for(k = 0; k < K+1; ++k){
    if(pizitk0[k] >= uit){
      whichK.push_back(k); 
    }
  }
  forList[t] = whichK; 
  // rest of time points 
  for(t = 1; t < tmax; ++t){
    whichKminus1 = as<arma::vec>(forList[t-1]); // vector of previous states 
    whichKprev = wrap(whichKminus1); // IntegerVector
    uit = ui[t]; 
    whichKtemp = as<arma::vec>(whichK); 
    whichKtemp.set_size(0);
    whichK = wrap(whichKtemp); 
    pizit = as<arma::mat>(pizi[t]); // prob matrix for current time point 
    for(idx = 0; idx < whichKprev.length(); ++idx){
      tempK = whichKprev[idx]; // one previous state at a time 
      if(tempK < K){ // can't go from a new state 
        pizitk = (pizit.row(tempK)).t(); // get the transition probs for the one previous state 
        for(k = 0; k < K+1; ++k){
          if(pizitk[k] >= uit){
            whichK.push_back(k); 
          }
        }
      }
    }
    forList[t] = sort_unique(whichK); 
  }
  
  // last time point 
  backList[tmax-1] = forList[tmax-1]; 
  for(t = tmax-2; t>=0; --t){
    whichKplus1 = as<arma::vec>(backList[t+1]); // vector of previous states 
    whichKnext = wrap(whichKplus1); // IntegerVector
    uitplus1 = ui[t+1]; 
    
    whichKtemp = as<arma::vec>(whichK); 
    whichKtemp.set_size(0);
    whichK = wrap(whichKtemp); 
    
    pizitplus1 = as<arma::mat>(pizi[t+1]); // prob matrix for current time point 
    for(idx = 0; idx < whichKnext.length(); ++idx){
      tempK = whichKnext[idx]; // one previous state at a time 
      pizitk = (pizitplus1.col(tempK)); // get the transition probs for the one previous state 
      for(k = 0; k < K; ++k){
        if(pizitk[k] >= uitplus1){
          whichK.push_back(k); 
        }
      }
    }
    backList[t] = sort_unique(whichK); 
  }
  
  
  // intersect
  for(t = 0; t < tmax; ++t){
    f1 = as<arma::vec>(forList[t]); 
    f2 = as<arma::vec>(backList[t]); 
    f3 = wrap(f1);
    f4 = wrap(f2); 
    // needs to be a vector 
    stateList[t] = sort_unique(intersect(f3, f4))+1;
  }
  
  return(stateList);
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
