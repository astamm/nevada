#include "operations.h"
#include "distances.h"

arma::mat mean_nvd_impl(const Rcpp::List &z)
{
  unsigned int n = z.size();
  unsigned int vcount = Rcpp::as<arma::mat>(z[0]).n_cols;
  arma::mat out(vcount, vcount, arma::fill::zeros);

  for (unsigned int i = 0;i < n;++i)
  {
    out *= i / (i + 1.0);
    out += Rcpp::as<arma::mat>(z[i]).replace(arma::datum::nan, 0.0) / (i + 1.0);
  }

  return out;
}

double var_nvd_impl(const Rcpp::List &z, const std::string distance)
{
  unsigned int n = z.size();
  double resVal = 0.0;
  double distanceValue = 0.0;
  arma::mat net1, net2;

  for (unsigned int i = 0;i < n-1;++i)
  {
    net1 = Rcpp::as<arma::mat>(z[i]);

    for (unsigned int j = i+1;j < n;++j)
    {
      net2 = Rcpp::as<arma::mat>(z[j]);

      if (distance == "hamming")
        distanceValue = dist_hamming_impl(net1, net2);
      else if (distance == "frobenius")
        distanceValue = dist_frobenius_impl(net1, net2);
      else if (distance == "spectral")
        distanceValue = dist_spectral_impl(net1, net2);
      else if (distance == "root-euclidean")
        distanceValue = dist_root_euclidean_impl(net1, net2);
      else
        Rcpp::stop("Unavailable distance.\n");

      resVal += distanceValue * distanceValue;
    }
  }

  resVal /= (n * (n - 1.0));
  return resVal;
}
