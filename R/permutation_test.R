#' Title Compute the p-value for network-valued data
#'
#' \code{get_pvalue} returns the p-value associated to the chosen test statistic
#' and computed in a permutational framework.
#'
#' All the possibile permutations are performed if their number is smaller than
#' B or if it is smaller than the permutations that has to be performed to
#' achieve tolerance. On contrary, the p-value is estimanted via Monte Carlo
#' simulation.
#'
#'
#' @param x List of adjacency matrix/igraph object of the first sample.
#' @param y List of adjacency matrix/igraph object of the second sample.
#' @param representation The choosen representation for the elements of x and y.
#'   It can be: "adjacency", "laplacian" and "modularity".
#' @param distance The chosen distance for calculate the test statistic. It can
#'   be: "hamming", "frobenius", "spectral", "rooteuclidean".
#' @param statistic The chosen test statistic. It can be: "average" or
#'   "frechet".
#' @param B The number of permutation or the tolerance.
#' @param modality The parameter for fix the number of permutations, either "b" for
#'   the exact number of permutations to be performed or "tol" for the
#'   tolerance.
#' @param alpha The level of significancy.
#'
#' @return \code{get_pvalue} returns a scalar.
#' @export
#'
#' @examples
get_pvalue <- function(x, y, representation, distance, statistic, B, modality, alpha) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  d <- get_paired_distances(x, y, representation, distance)


  #DA FARE: verifica della raggiungibilità di alpha

  T0 <- switch(statistic,
              average <- get_average_statistic(d, 1:n1),
              frechet <- get_frechet_statistic(d, 1:n1))


  all <- FALSE

  Bnew <- switch(modality,
                 b <- B,
                 tol <- (qnorm(alpha/2, lower.tail = FALSE)/tol)^2)

  group1.perm <- t(combn(1:n,n1))
  M <- dim(group1.perm)[1]

  if (Bnew>=M) {
    all <- TRUE
  }


  if (all) {
    Tp <- rep(-1, M)
    for (i in 1:M) {
      p <- group1.perm[i,]
      Tp[i] <- switch(statistic,
                      average <- get_average_statistic(d, p),
                      frechet <- get_frechet_statistic(d, p))
    }
  }

  else {
    Tp <- rep(-1, Bnew)
    for (i in 1:B) {
      p <- sample(1:n, n, replace=FALSE)
      Tp[i] <- switch(statistic,
                      average <- get_average_statistic(d, p[1:n1]),
                      frechet <- get_frechet_statistic(d, p[1:n1]))
    }
  }




  p <- sum(Tp >= T0)/B

  return(p)

}
