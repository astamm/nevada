// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::mat mean_nvd_impl(const Rcpp::List &z);
