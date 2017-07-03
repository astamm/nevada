#' Pairwise Distance Matrix Between Two Samples of Networks
#'
#' This function computes the matrix of pairwise distances between all the
#' elements of the two samples put together. The cardinality of the fist sample
#' is denoted by \eqn{n1} and that of the second one is denoted by \eqn{n2}.
#'
#' @param x A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects or matrix
#'   representations of underlying networks from a given first population.
#' @param y A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects or matrix
#'   representations of underlying networks from a given second population.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"} [default], \code{"laplacian"} and
#'   \code{"modularity"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"} [default], \code{"frobenius"},
#'   \code{"spectral"} and \code{"root-euclidean"}.
#'
#' @return A matrix of dimension \eqn{(n1+n2) \times (n1+n2)} containing the
#'   distances between all the elements of the two samples put together.
#' @export
#'
#' @examples
#' n <- 25L
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
#' d <- get_distance_matrix(x, y, "adjacency", "spectral")
#' d
get_distance_matrix <- function(x, y, representation = "adjacency", distance = "hamming") {
  n <- length(x)
  m <- length(y)
  N <- n + m
  x <- purrr::map(x, format_input, representation)
  y <- purrr::map(y, format_input, representation)
  z <- c(x, y)

  dist <- rep(-1, N * (N - 1) / 2)
  k <- 1
  for (i in 1:(N-1)) {
    network1 <- z[[i]]
    for (j in ((i+1):N)) {
      network2 <- z[[j]]
      dist[k] <- switch(
        distance,
        "hamming" = get_hamming_distance(network1, network2, representation),
        "frobenius" = get_frobenius_distance(network1, network2, representation),
        "spectral" = get_spectral_distance(network1, network2, representation),
        "root-euclidean" = get_root_euclidean_distance(network1, network2, representation)
        )
      k <- k + 1
    }
  }

  d <- diag(rep(0, N))
  d[lower.tri(d)] <- dist
  d <- d + t(d)
  d
}
