// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat aux_nbdsmooth(const arma::mat &A, const unsigned int N);

// [[Rcpp::export]]
arma::mat est_nbdsmooth(const arma::mat &A);
