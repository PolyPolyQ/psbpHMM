### practice c++ functions



# // [[Rcpp::export]]
# NumericVector fun0(NumericVector a, NumericVector b, int K){
#   NumericVector f1(K);
#   for(int i = 0; i < K; i++){
#     f1[i] = R::pnorm(a[i] + b[i], 0, 1, 1, 0);
#   }
#   return(f1);
# }
# 
# // [[Rcpp::export]]
# List fun3(List beta, arma::vec alpha, arma::vec x1) {
#   int n_list = beta.size();
#   Rcpp::List out(n_list);
#   for (int i = 0; i < n_list; ++i) {
#     out[i] = alpha[i] + as<arma::mat>(beta[i]).t()*x1;
#   }
#   return(out);
# }
# 
# 
# // [[Rcpp::export]]
# arma::vec fun4(List beta, arma::vec alpha, arma::vec x1) {
#   int K = beta.size();
#   arma::vec out(K);
#   for (int i = 0; i < K; ++i) {
#     arma::mat y = alpha[i] + as<arma::mat>(beta[i]).t()*x1;
#     //double z = as_scalar(y);
#     out[i] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
#   }
#   return(out);
# }
# 
# // [[Rcpp::export]]
# arma::vec fun5(List beta, arma::vec alpha, arma::vec x1) {
#   int K = beta.size();
#   arma::vec out(K-1);
#   for (int i = 0; i < K-1; ++i) {
#     arma::mat y = alpha[i] + as<arma::mat>(beta[i]).t()*x1;
#     out[i] = 1-(R::pnorm(as_scalar(y), 0, 1, 1, 0));
#   }
#   return(out);
# }
# 
# // [[Rcpp::export]]
# arma::vec fun6(List beta, arma::vec alpha, arma::vec x1) {
#   int K = beta.size();
#   arma::vec out(K);
#   for (int i = 0; i < K; ++i) {
#     arma::mat y = alpha[i] + as<arma::mat>(beta[i]).t()*x1;
#     out[i] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
#   }
#   
#   // out = first1
#   // z1 = prod1 (almost) 
#   arma::vec z = 1-out;
#   arma::vec z1 = cumprod(z);
#   
#   return(z1);
# }


# // [[Rcpp::export]]
# arma::vec fun7(List beta, arma::vec alpha, arma::vec x1) {
#   int K = beta.size();
#   arma::vec first(K);
#   arma::vec second(K);
#   for (int i = 0; i < K; ++i) {
#     arma::mat y = alpha[i] + as<arma::mat>(beta[i]).t()*x1;
#     first[i] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
#     if (i > 0) second[i] = 1-first[i-1];
#   }
#   second[0] = 1;
#   arma::vec third = cumprod(second);
#   arma::vec fourth(K+1);
#   for(int i = 0; i < K; i++){
#     fourth[i] = first[i]*third[i];
#   }
#   double last = sum(fourth);
#   fourth[K]=1-last;
#   return(fourth);
# }
# 
# 
# // [[Rcpp::export]]
# arma::mat fun8(List beta, arma::mat X, arma::mat ajk, int t){
#   
#   //step by step return each piece 
#   int K = beta.size();
#   arma::mat x = X.row(t);
#   arma::vec first(K);
#   arma::vec second(K);
#   arma::mat Z;
#   
#   for(int j = 0; j < K; ++j){
#     for(int k = 0; k < K; ++k){
#       arma::mat xb = as<arma::rowvec>(beta[k])*x.t();
#       double a = ajk(k,j);
#       arma::mat y = a + xb; 
#       first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
#       if (k > 0) second[k] = 1-first[k-1];
#     }
#     second[0] = 1;
#     arma::vec third = cumprod(second);
#     arma::vec fourth(K+1);
#     double s = 0;
#     for(int i = 0; i < K; i++){
#       fourth[i] = first[i]*third[i];
#       s = s + fourth[i];
#     }
#     //double last = sum(fourth);
#     fourth[K]=1-s; // index starts a 0 
#     // fourth is the jth row
#     
#     Z = join_rows(Z, fourth);
#     
#     
#   }
#   return(Z.t());
# }




