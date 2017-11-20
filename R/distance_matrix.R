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
  x <- lapply(x, format_input, representation)
  y <- lapply(y, format_input, representation)
  z <- c(x, y)
  internal_distance_matrix(z, distance);
}
