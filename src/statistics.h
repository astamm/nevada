// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double stat_lot_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices);

// [[Rcpp::export]]
double stat_sot_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices);

// [[Rcpp::export]]
double stat_biswas_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices);

// [[Rcpp::export]]
double stat_energy_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices, const unsigned int alpha = 1);

// [[Rcpp::export]]
double stat_t_euclidean_impl(const Rcpp::List &x, const Rcpp::List &y, const bool pooled = true);

// [[Rcpp::export]]
arma::vec stat_edge_count_impl(const arma::mat &E, const arma::vec &indices);
