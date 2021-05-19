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