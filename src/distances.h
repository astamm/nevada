// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double dist_hamming_impl(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
double dist_frobenius_impl(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
double dist_spectral_impl(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
double dist_root_euclidean_impl(const arma::mat &x, const arma::mat &y);

// [[Rcpp::export]]
arma::mat dist_nvd_impl(const Rcpp::List &z, const std::string distance = "frobenius");

// [[Rcpp::export]]
double ipro_frobenius_impl(const arma::mat &x, const arma::mat &y);
