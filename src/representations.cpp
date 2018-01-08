#include "representations.h"

arma::mat repr_adjacency_impl(const unsigned int numberOfVertices, const arma::mat &edgeList, const arma::vec &weights)
{
  arma::mat adjMatrix(numberOfVertices, numberOfVertices, arma::fill::zeros);
  unsigned int numberOfEdges = edgeList.n_rows;
  for (unsigned int i = 0;i < numberOfEdges;++i)
  {
    unsigned int rowIndex = edgeList(i, 0) - 1;
    unsigned int colIndex = edgeList(i, 1) - 1;
    double w = weights[i];

    adjMatrix(rowIndex, colIndex) = w;

    if (rowIndex != colIndex)
      adjMatrix(colIndex, rowIndex) = w;
  }

  return adjMatrix;
}
