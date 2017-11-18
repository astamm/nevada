// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double internal_frobenius(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  double squaredDiff = 0.0;
  for (unsigned int i = 0;i < n;i++)
    for (unsigned int j = 0;j < n;j++)
      squaredDiff += (x(i, j) - y(i, j)) * (x(i, j) - y(i, j));

  return std::sqrt(squaredDiff);
}

// [[Rcpp::export]]
double internal_spectral(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  arma::vec xx = arma::eig_sym(x);
  arma::vec yy = arma::eig_sym(y);

  double squaredDiff = 0.0;
  for (unsigned int i = 0;i < n;i++)
    squaredDiff += (xx[i] - yy[i]) * (xx[i] - yy[i]);

  return std::sqrt(squaredDiff);
}
