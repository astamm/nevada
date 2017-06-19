#' Comparison of Network Distributions
#'
#' This function carries out an hypothesis test where the null hypothesis is
#' that the two populations of networks share the same underlying probabilistic
#' distribution against the alternative hypothesis that the two populations come
#' from different distributions. The test is performed in a non-parametric
#' fashion using a permutational framework in which several statistics can be
#' used, together with several choices of network matrix representations and
#' distances between networks.
#'
#' The exhaustive list of all possible permutations is used if its length is
#' smaller than either the user-defined number of random permutations requested
#' or the number of permutations required for achieving a user-defined p-value
#' resolution. Otherwise, the p-value is estimated via Monte-Carlo simulations.
#'
#' @param x A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects
#'   or matrix representations of underlying networks from a given first
#'   population.
#' @param y A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects
#'   or matrix representations of underlying networks from a given second
#'   population.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"} [default], \code{"laplacian"} and
#'   \code{"modularity"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"} [default], \code{"frobenius"},
#'   \code{"spectral"} and \code{"root-euclidean"}.
#' @param statistic A string specifying the chosen test statistic, among:
#'   \code{"average"} [default] or \code{"frechet"}.
#' @param B The number of permutation or the tolerance (default: \code{1000L}).
#' @param modality A string for determining the number of permutations: either
#'   \code{"b"} for the exact number of permutations to be performed or
#'   \code{"tol"} for the tolerance (default: \code{"b"}).
#' @param alpha The significance level (default: \code{0.05}).
#'
#' @return A \code{\link[base]{list}} with two components: the value of the
#'   statistic for the original two samples and the p-value of the resulting
#'   permutation test.
#' @export
#'
#' @examples
network_test2p <- function(x, y, representation = "adjacency", distance = "hamming",
                           statistic = "average", B = 1000L, modality = "b", alpha = 0.05) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  d <- get_distance_matrix(x, y, representation, distance)

  #DA FARE: verifica della raggiungibilità di alpha

  T0 <- switch(
    statistic,
    average = get_average_statistic(d, 1:n1),
    frechet = get_frechet_statistic(d, 1:n1)
  )

  Bnew <- switch(
    modality,
    b = B,
    tol = (qnorm(alpha / 2, lower.tail = FALSE) / tol)^2
  )

  group1.perm <- t(combn(1:n, n1))
  M <- nrow(group1.perm)

  if (Bnew >= M) {
    Tp <- rep(-1, M)
    for (i in 1:M) {
      p <- group1.perm[i, ]
      Tp[i] <- switch(
        statistic,
        average = get_average_statistic(d, p),
        frechet = get_frechet_statistic(d, p)
      )
    }
  } else {
    Tp <- rep(-1, Bnew)
    for (i in 1:B) {
      p <- sample(1:n, n, replace = FALSE)
      Tp[i] <- switch(
        statistic,
        average = get_average_statistic(d, p[1:n1]),
        frechet = get_frechet_statistic(d, p[1:n1])
      )
    }
  }

  p <- sum(Tp >= T0) / B

  list(statistic = T0, pvalue = p)
}
