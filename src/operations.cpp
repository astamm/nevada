#include "operations.h"

arma::mat mean_nvd_impl(const Rcpp::List &z)
{
  unsigned int n = z.size();
  unsigned int vcount = Rcpp::as<arma::mat>(z[0]).n_cols;
  arma::mat out(vcount, vcount, arma::fill::zeros);

  for (unsigned int i = 0;i < n;++i)
  {
    out *= i / (i + 1.0);
    out += Rcpp::as<arma::mat>(z[i]) / (i + 1.0);
  }

  return out;
}
