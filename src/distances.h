// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double internal_hamming(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
double internal_frobenius(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
double internal_spectral(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
double internal_root_euclidean(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
arma::mat internal_distance_matrix(const Rcpp::List &z, const std::string distance = "frobenius");
