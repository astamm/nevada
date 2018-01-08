// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::mat repr_adjacency_impl(const unsigned int numberOfVertices, const arma::mat &edgeList, const arma::vec &weights);
