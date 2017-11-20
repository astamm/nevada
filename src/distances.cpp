// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double internal_hamming(const arma::mat &x, const arma::mat &y)
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

// [[Rcpp::export]]
double internal_frobenius(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  double squaredDiff = 0.0;
  for (unsigned int i = 0;i < n;++i)
  {
    for (unsigned int j = 0;j < n;++j)
    {
      double tmpVal = (x(i, j) - y(i, j)) * (x(i, j) - y(i, j));
      squaredDiff += tmpVal;

      if (i != j)
        squaredDiff += tmpVal;
    }
  }

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

// [[Rcpp::export]]
double internal_root_euclidean(const arma::mat &x, const arma::mat &y)
{
  unsigned int n = x.n_rows;

  arma::vec eigval1, eigval2;
  arma::mat eigvec1, eigvec2;
  arma::eig_sym(eigval1, eigvec1, x);
  arma::eig_sym(eigval2, eigvec2, y);

  for (unsigned int i = 0;i < n;++i)
  {
    double eigVal = eigval1[i];

    if (eigVal < 0.0)
      eigval1[i] = 0.0;
    else
      eigval1[i] = std::sqrt(eigVal);

    eigVal = eigval2[i];

    if (eigVal < 0.0)
      eigval2[i] = 0.0;
    else
      eigval2[i] = std::sqrt(eigVal);
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

// [[Rcpp::export]]
arma::mat internal_distance_matrix(const Rcpp::List &z, const std::string &distance)
{
  unsigned int n = z.size();
  arma::mat out(n, n, arma::fill::zeros);
  arma::mat net1, net2;

  for (unsigned int i = 0;i < n-1;++i)
  {
    net1 = Rcpp::as<arma::mat>(z[i]);

    for (unsigned int j = i+1;j < n;++j)
    {
      net2 = Rcpp::as<arma::mat>(z[j]);

      double distanceValue = 0.0;
      if (distance == "hamming")
        distanceValue = internal_hamming(net1, net2);
      else if (distance == "frobenius")
        distanceValue = internal_frobenius(net1, net2);
      else if (distance == "spectral")
        distanceValue = internal_spectral(net1, net2);
      else if (distance == "root-euclidean")
        distanceValue = internal_root_euclidean(net1, net2);
      else
      {
        Rcpp::Rcout << "Unavailable distance." << std::endl;
        exit(-1);
      }

      out(i,j) = distanceValue;
      out(j,i) = distanceValue;
    }
  }

  return out;
}
