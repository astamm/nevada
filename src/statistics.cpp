#include "statistics.h"
#include "distances.h"

double stat_t_euclidean_impl(const Rcpp::List &x, const Rcpp::List &y, const bool pooled)
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

  double varianceValue = (pooled) ? (ssd1 + ssd2) / (n1 + n2 - 2.0) * (1.0 / n1 + 1.0 / n2) : ssd1 / (n1 - 1.0) / n1 + ssd2 / (n2 - 1.0) / n2;

  if (varianceValue < std::sqrt(std::numeric_limits<double>::epsilon()))
    return meanDifference;

  return meanDifference / std::sqrt(varianceValue);
}

arma::vec stat_edge_count_impl(const arma::mat &E, const arma::vec &indices)
{
  unsigned int r1 = 0;
  unsigned int r2 = 0;

  for (unsigned int i = 0;i < E.n_rows;++i)
  {
    bool firstIndexIn = false;
    bool secondIndexIn = false;
    for (unsigned int j = 0;j < indices.n_elem;++j)
    {
      if (E(i,0) == indices[j])
        firstIndexIn = true;

      if (E(i,1) == indices[j])
        secondIndexIn = true;

      if (firstIndexIn && secondIndexIn)
      {
        r1++;
        break;
      }
    }

    if (!firstIndexIn && !secondIndexIn)
      r2++;
  }

  arma::vec rValues(2);
  rValues[0] = r1;
  rValues[1] = r2;
  return rValues;
}
