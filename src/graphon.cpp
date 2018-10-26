#include "graphon.h"

arma::mat aux_nbdsmooth(const arma::mat &A, const unsigned int N)
{
  arma::mat resVal(N, N, arma::fill::zeros);
  arma::mat A2 = (A * A) / N;

  for (unsigned int i = 0;i < N - 1;++i)
  {
    for (unsigned int j = i + 1;j < N;++j)
    {
      double tmpVal = arma::abs(A2.row(i) - A2.row(j)).max();
      resVal(i,j) = tmpVal;
      resVal(j,i) = tmpVal;
    }
  }

  return resVal;
}

arma::mat est_nbdsmooth(const arma::mat &A)
{
  unsigned int N = A.n_rows;
  double h = std::sqrt(std::log(N) / N);
  arma::mat D = aux_nbdsmooth(A, N);
  arma::mat resVal(N, N, arma::fill::zeros);
  arma::rowvec tmpVec;

  for (unsigned int i = 0;i < N;++i)
  {
    // Get h-th quantile of D[i, ]
    tmpVec = arma::sort(D.row(i));
    double quantileValue = tmpVec[std::ceil(N * h)];

    for (unsigned int j = 0;j < N;++j)
      resVal(i, j) = (double)(D(i, j) < quantileValue);
  }

  resVal.each_col() /= arma::sum(resVal, 1);
  resVal = resVal * A;
  resVal = (resVal + resVal.t()) / 2.0;

  return resVal;
}
