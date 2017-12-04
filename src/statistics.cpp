#include "statistics.h"
#include "distances.h"

double stat_dom_frobenius_impl(const Rcpp::List &x, const Rcpp::List &y, const bool standardize)
{
  unsigned int n1 = x.size();
  unsigned int n2 = y.size();

  unsigned int numberOfNodes = Rcpp::as<arma::mat>(x[0]).n_rows;

  arma::mat m1(numberOfNodes, numberOfNodes, arma::fill::zeros);
  for (unsigned int i = 0;i < n1;++i)
    m1 += Rcpp::as<arma::mat>(x[i]);
  m1 /= n1;

  arma::mat m2(numberOfNodes, numberOfNodes, arma::fill::zeros);
  for (unsigned int i = 0;i < n2;++i)
    m2 += Rcpp::as<arma::mat>(y[i]);
  m2 /= n2;

  double meanDifference = dist_frobenius_impl(m1, m2);

  if (!standardize)
    return meanDifference;

  double ssd1 = 0;
  for (unsigned int i = 0;i < n1;++i)
  {
    double tmpVal = dist_frobenius_impl(Rcpp::as<arma::mat>(x[i]), m1);
    ssd1 += tmpVal * tmpVal;
  }

  double ssd2 = 0;
  for (unsigned int i = 0;i < n2;++i)
  {
    double tmpVal = dist_frobenius_impl(Rcpp::as<arma::mat>(y[i]), m2);
    ssd2 += tmpVal * tmpVal;
  }

  double pooledVariance = (ssd1 + ssd2) / (n1 + n2 - 2);

  return meanDifference / std::sqrt(pooledVariance);
}
