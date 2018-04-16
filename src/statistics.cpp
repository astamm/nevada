#include "statistics.h"
#include "distances.h"

double stat_lot_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices)
{
  double xyMean = 0.0, xxMean = 0.0, yyMean = 0.0;
  unsigned int n1 = firstGroupIndices.n_elem;
  unsigned int n2 = secondGroupIndices.n_elem;
  unsigned int counter_xy = 0, counter_xx = 0, counter_yy = 0;
  for (unsigned int i = 0;i < n1;++i)
  {
    for (unsigned int j = 0;j < n2;++j)
    {
      xyMean *= (counter_xy / (counter_xy + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, secondGroupIndices[j] - 1);
      xyMean += distanceValue * distanceValue / (counter_xy + 1.0);
      ++counter_xy;

      if (i == 0)
      {
        for (unsigned int k = j + 1;k < n2;++k)
        {
          yyMean *= (counter_yy / (counter_yy + 1.0));
          double distanceValue = distanceMatrix(secondGroupIndices[j] - 1, secondGroupIndices[k] - 1);
          yyMean += distanceValue * distanceValue / (counter_yy + 1.0);
          ++counter_yy;
        }
      }
    }

    for (unsigned int j = i + 1;j < n1;++j)
    {
      xxMean *= (counter_xx / (counter_xx + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, firstGroupIndices[j] - 1);
      xxMean += distanceValue * distanceValue / (counter_xx + 1.0);
      ++counter_xx;
    }
  }

  double var_x = xxMean / 2.0;
  double var_y = yyMean / 2.0;
  double var_c = var_x / n1 + var_y / n2;
  return (xyMean - var_x - var_y) / var_c;
}

double stat_sot_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices)
{
  double xxMean = 0.0, yyMean = 0.0;
  unsigned int n1 = firstGroupIndices.n_elem;
  unsigned int n2 = secondGroupIndices.n_elem;
  unsigned int counter_xx = 0, counter_yy = 0;
  for (unsigned int i = 0;i < n1;++i)
  {
    for (unsigned int j = 0;j < n2;++j)
    {
      if (i == 0)
      {
        for (unsigned int k = j + 1;k < n2;++k)
        {
          yyMean *= (counter_yy / (counter_yy + 1.0));
          double distanceValue = distanceMatrix(secondGroupIndices[j] - 1, secondGroupIndices[k] - 1);
          yyMean += distanceValue * distanceValue / (counter_yy + 1.0);
          ++counter_yy;
        }
      }
    }

    for (unsigned int j = i + 1;j < n1;++j)
    {
      xxMean *= (counter_xx / (counter_xx + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, firstGroupIndices[j] - 1);
      xxMean += distanceValue * distanceValue / (counter_xx + 1.0);
      ++counter_xx;
    }
  }

  double var_x = xxMean / 2.0;
  double var_y = yyMean / 2.0;
  return std::max(var_x, var_y) / std::min(var_x, var_y);
}

double stat_biswas_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices)
{
  double xyMean = 0.0, xxMean = 0.0, yyMean = 0.0;
  unsigned int n1 = firstGroupIndices.n_elem;
  unsigned int n2 = secondGroupIndices.n_elem;
  unsigned int counter_xy = 0, counter_xx = 0, counter_yy = 0;
  for (unsigned int i = 0;i < n1;++i)
  {
    for (unsigned int j = 0;j < n2;++j)
    {
      xyMean *= (counter_xy / (counter_xy + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, secondGroupIndices[j] - 1);
      xyMean += distanceValue / (counter_xy + 1.0);
      ++counter_xy;

      if (i == 0)
      {
        for (unsigned int k = j + 1;k < n2;++k)
        {
          yyMean *= (counter_yy / (counter_yy + 1.0));
          double distanceValue = distanceMatrix(secondGroupIndices[j] - 1, secondGroupIndices[k] - 1);
          yyMean += distanceValue / (counter_yy + 1.0);
          ++counter_yy;
        }
      }
    }

    for (unsigned int j = i + 1;j < n1;++j)
    {
      xxMean *= (counter_xx / (counter_xx + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, firstGroupIndices[j] - 1);
      xxMean += distanceValue / (counter_xx + 1.0);
      ++counter_xx;
    }
  }

  return (xxMean - xyMean) * (xxMean - xyMean) + (yyMean - xyMean) * (yyMean - xyMean);
}

double stat_energy_impl(const arma::mat &distanceMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices, const unsigned int alpha)
{
  double xyMean = 0.0, xxMean = 0.0, yyMean = 0.0;
  unsigned int n1 = firstGroupIndices.n_elem;
  unsigned int n2 = secondGroupIndices.n_elem;
  unsigned int counter_xy = 0, counter_xx = 0, counter_yy = 0;
  for (unsigned int i = 0;i < n1;++i)
  {
    for (unsigned int j = 0;j < n2;++j)
    {
      xyMean *= (counter_xy / (counter_xy + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, secondGroupIndices[j] - 1);
      if (alpha != 1)
        distanceValue = std::pow(distanceValue, (double)alpha);
      xyMean += distanceValue / (counter_xy + 1.0);
      ++counter_xy;

      if (i == 0)
      {
        for (unsigned int k = 0;k < n2;++k)
        {
          yyMean *= (counter_yy / (counter_yy + 1.0));
          double distanceValue = distanceMatrix(secondGroupIndices[j] - 1, secondGroupIndices[k] - 1);
          if (alpha != 1)
            distanceValue = std::pow(distanceValue, (double)alpha);
          yyMean += distanceValue / (counter_yy + 1.0);
          ++counter_yy;
        }
      }
    }

    for (unsigned int j = 0;j < n1;++j)
    {
      xxMean *= (counter_xx / (counter_xx + 1.0));
      double distanceValue = distanceMatrix(firstGroupIndices[i] - 1, firstGroupIndices[j] - 1);
      if (alpha != 1)
        distanceValue = std::pow(distanceValue, (double)alpha);
      xxMean += distanceValue / (counter_xx + 1.0);
      ++counter_xx;
    }
  }

  return xyMean - (xxMean + yyMean) / 2.0;
}

double stat_cq_impl(const arma::mat &similarityMatrix, const arma::vec &firstGroupIndices, const arma::vec &secondGroupIndices)
{
  double xyMean = 0.0, xxMean = 0.0, yyMean = 0.0;
  unsigned int n1 = firstGroupIndices.n_elem;
  unsigned int n2 = secondGroupIndices.n_elem;
  unsigned int counter_xy = 0, counter_xx = 0, counter_yy = 0;
  for (unsigned int i = 0;i < n1;++i)
  {
    for (unsigned int j = 0;j < n2;++j)
    {
      xyMean *= (counter_xy / (counter_xy + 1.0));
      double similarityValue = similarityMatrix(firstGroupIndices[i] - 1, secondGroupIndices[j] - 1);
      xyMean += similarityValue / (counter_xy + 1.0);
      ++counter_xy;

      if (i == 0)
      {
        for (unsigned int k = 0;k < n2;++k)
        {
          if (j == k)
            continue;

          yyMean *= (counter_yy / (counter_yy + 1.0));
          double similarityValue = similarityMatrix(secondGroupIndices[j] - 1, secondGroupIndices[k] - 1);
          yyMean += similarityValue / (counter_yy + 1.0);
          ++counter_yy;
        }
      }
    }

    for (unsigned int j = 0;j < n1;++j)
    {
      if (i == j)
        continue;

      xxMean *= (counter_xx / (counter_xx + 1.0));
      double similarityValue = similarityMatrix(firstGroupIndices[i] - 1, firstGroupIndices[j] - 1);
      xxMean += similarityValue / (counter_xx + 1.0);
      ++counter_xx;
    }
  }

  return xxMean + yyMean - 2 * xyMean;
}

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
