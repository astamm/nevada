#include "distances.h"

double dist_hamming_impl(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  double absDiff = 0.0;
  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = i;j < n;++j)
    {
      double tmpVal = std::abs(x(i, j) - y(i, j));
      absDiff += tmpVal;

      if (i != j)
        absDiff += tmpVal;
    }
  }

  return absDiff / n / (n - 1.0);
}

double dist_frobenius_impl(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  double squaredDiff = 0.0;
  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = i;j < n;++j)
    {
      double tmpVal = (x(i, j) - y(i, j)) * (x(i, j) - y(i, j));
      squaredDiff += tmpVal;

      if (i != j)
        squaredDiff += tmpVal;
    }
  }

  return std::sqrt(squaredDiff);
}

double dist_spectral_impl(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  arma::vec xx = arma::eig_sym(x);
  arma::vec yy = arma::eig_sym(y);

  double squaredDiff = 0.0;
  for (unsigned int i = 0;i < n;i++)
    squaredDiff += (xx[i] - yy[i]) * (xx[i] - yy[i]);

  return std::sqrt(squaredDiff);
}

double dist_root_euclidean_impl(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  arma::vec eigval1, eigval2;
  arma::mat eigvec1, eigvec2;
  arma::eig_sym(eigval1, eigvec1, x);
  arma::eig_sym(eigval2, eigvec2, y);

  for (unsigned int i = 0;i < n;++i)
  {
    eigval1[i] = std::sqrt(std::max(0.0, eigval1[i]));
    eigval2[i] = std::sqrt(std::max(0.0, eigval2[i]));
  }

  double squaredDiff = 0.0;
  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = i;j < n;++j)
    {
      double tmpVal = 0.0;
      for (unsigned int k = 0;k < n;++k)
      {
        double diff = eigval1[k] * eigvec1(i,k) * eigvec1(j,k) - eigval2[k] * eigvec2(i,k) * eigvec2(j,k);
        tmpVal += diff * diff;
      }

      squaredDiff += tmpVal;

      if (i != j)
        squaredDiff += tmpVal;
    }
  }

  return std::sqrt(squaredDiff);
}

Rcpp::NumericVector dist_nvd_impl(const Rcpp::List &z, const std::string distance)
{
  unsigned int n = z.size();
  Rcpp::NumericVector outValue(n * (n - 1) / 2);
  arma::mat net1, net2;

  for (unsigned int i = 0;i < n-1;++i)
  {
    net1 = Rcpp::as<arma::mat>(z[i]);

    for (unsigned int j = i+1;j < n;++j)
    {
      net2 = Rcpp::as<arma::mat>(z[j]);

      double distanceValue = 0.0;
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

      double rowIndex = i + 1;
      double colIndex = j + 1;
      double indexValue = n * (rowIndex - 1) - rowIndex * (rowIndex - 1) / 2 + colIndex - rowIndex -1;
      outValue(indexValue) = distanceValue;
    }
  }

  outValue.attr("class") = "dist";
  outValue.attr("Size") = n;
  outValue.attr("Diag") = false;
  outValue.attr("Upper") = false;
  return outValue;
}

double ipro_frobenius_impl(const arma::mat &x, const arma::mat &y)
{
  return arma::trace(x.t() * y);
}
