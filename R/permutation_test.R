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
#'   If this number is lower than \code{1}, it is intended as a tolerance.
#'   Otherwise, it is intended as the number of required permutations.
#' @param alpha The significance level (default: \code{0.05}).
#'
#' @return A \code{\link[base]{list}} with two components: the value of the
#'   statistic for the original two samples and the p-value of the resulting
#'   permutation test.
#' @export
#'
#' @examples
#' n <- 25L
#'
#' x <- list()
#' y <- list()
#' for (i in 1:10) {
#'   X <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#'   Y <- igraph::barabasi.game(n, m = 3, power = 2, directed = FALSE)
#'   adjX <- get_adjacency(X)
#'   adjY <- get_adjacency(Y)
#'   x[[i]] <- adjX
#'   y[[i]] <- adjY
#' }
#' test1 <- network_test2p(x, y, "adjacency")
#' test1
#'
#' x <- list()
#' y <- list()
#' for (i in 1:10) {
#'   X <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#'   Y <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#'   adjX <- get_adjacency(X)
#'   adjY <- get_adjacency(Y)
#'   x[[i]] <- adjX
#'   y[[i]] <- adjY
#' }
#' test2 <- network_test2p(x, y, "adjacency")
#' test2
network_test2p <- function(x, y, representation = "adjacency", distance = "hamming",
                           statistic = "average", B = 1000L, alpha = 0.05) {
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

  if (B < 1)
    B <- (qnorm(alpha / 2, lower.tail = FALSE) / tol)^2

  if (!pvalue_resolution(n1, n2, B, alpha, TRUE))
    warning("The requested significance level cannot be reached
            given the value for B and the sample sizes.")

  group1.perm <- t(combn(1:n, n1))
  M <- nrow(group1.perm)

  if (B >= M) {
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
    Tp <- rep(-1, B)
    for (i in 1:B) {
      p <- sample(1:n, n, replace = FALSE)
      Tp[i] <- switch(
        statistic,
        average = get_average_statistic(d, p[1:n1]),
        frechet = get_frechet_statistic(d, p[1:n1])
      )
    }
  }

  p <- mean(Tp >= T0)

  list(statistic = T0, pvalue = p)
}
