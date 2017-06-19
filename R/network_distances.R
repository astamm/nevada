#' Distances Between Networks
#'
#' This is a collection of functions computing the distance between two
#' networks.
#'
#' Let \eqn{X} be the matrix representation of network \eqn{x} and \eqn{Y} be
#' the matrix representation of network \eqn{y}. The Hamming distance between
#' \eqn{x} and \eqn{y} is given by \deqn{\frac{1}{N(N-1)} \sum_{i,j} |X_{ij} -
#' Y_{ij}|,} where \eqn{N} is the number of vertices in networks \eqn{x} and
#' \eqn{y}. The Frobenius distance between \eqn{x} and \eqn{y} is given by
#' \deqn{\sqrt{\sum_{i,j} (X_{ij} - Y_{ij})^2}.} The spectral distance between
#' \eqn{x} and \eqn{y} is given by \deqn{\sqrt{\sum_i (a_i - b_i)^2},} where
#' \eqn{a} and \eqn{b} of the eigenvalues of \eqn{X} and \eqn{Y}, respectively.
#' This distance gives rise to classes of equivalence. Consider the spectral
#' decomposition of \eqn{X} and \eqn{Y}: \deqn{X=VAV^{-1}} and \deqn{Y =
#' UBU^{-1},} where \eqn{V} and \eqn{U} are the matrices whose colums are the
#' eigenvectors of \eqn{X} and \eqn{Y}, respectively and \eqn{A} and \eqn{B} are
#' the diagonal matrices with elements the eigenvalues of \eqn{X} and \eqn{Y},
#' respectively. The root-Euclidean distance between \eqn{x} and \eqn{y} is
#' given by \deqn{\sqrt{\sum_i (V \sqrt{A} V^{-1} - U \sqrt{B} U^{-1})^2}.}
#'
#' @param x An \code{\link[igraph]{igraph}} object or a matrix representing an
#'   underlying network.
#' @param y An \code{\link[igraph]{igraph}} object or a matrix representing an
#'   underlying network. Should have the same number of vertices as \code{x}.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"} [default], \code{"laplacian"} and
#'   \code{"modularity"}.
#'
#' @return A scalar measuring the distance between the two input networks.
#'
#' @examples
#' g1 <- igraph::erdos.renyi.game(20, 0.1)
#' g2 <- igraph::erdos.renyi.game(20, 0.2)
#' get_hamming_distance(g1, g2, "adjacency")
#' get_frobenius_distance(g1, g2, "adjacency")
#' get_spectral_distance(g1, g2, "laplacian")
#' get_root_euclidean_distance(g1, g2, "laplacian")
#' @name get-distance
NULL

#' @rdname get-distance
#' @export
get_hamming_distance <- function(x, y, representation = "adjacency") {
  l <- format_inputs(x, y, representation)
  x <- l$x
  y <- l$y
  d <- sum(abs(x - y))
  trunc(d * 10^4) / 10^4
}

#' @rdname get-distance
#' @export
get_frobenius_distance <- function(x, y, representation = "adjacency") {
  l <- format_inputs(x, y, representation)
  x <- l$x
  y <- l$y
  d <- sqrt(sum((x - y)^2))
  trunc(d * 10^4) / 10^4
}

#' @rdname get-distance
#' @export
get_spectral_distance <- function(x, y, representation = "adjacency") {
  l <- format_inputs(x, y, representation)
  x <- l$x
  y <- l$y

  dlX <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  dlX[abs(dlX) < 1e-10] <- 0

  dlY <- eigen(y, symmetric = TRUE, only.values = TRUE)$values
  dlY[abs(dlY) < 1e-10] <- 0

  d <- sqrt(sum((dlX - dlY)^2))
  trunc(d * 10^4) / 10^4
}

#' @rdname get-distance
#' @export
get_root_euclidean_distance <- function(x, y, representation = "adjacency") {
  l <- format_inputs(x, y, representation)
  x <- l$x
  y <- l$y

  rX <- eigen(x, symmetric = TRUE)
  vals <- rX$values
  vals[abs(vals) < 1e-10] <- 0
  dlX <- rX$vectors %*% diag(sqrt(vals)) %*% t(rX$vectors)

  rY <- eigen(y, symmetric = TRUE)
  vals <- rY$values
  vals[abs(vals) < 1e-10] <- 0
  dlY <- rY$vectors %*% diag(sqrt(vals)) %*% t(rY$vectors)

  d <- sqrt(sum((dlX - dlY)^2))
  trunc(d * 10^4) / 10^4
}
