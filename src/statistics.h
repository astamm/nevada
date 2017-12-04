// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double stat_dom_frobenius_impl(const Rcpp::List &x, const Rcpp::List &y, const bool standardize = true);
