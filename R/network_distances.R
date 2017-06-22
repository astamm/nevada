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
#'   among: \code{"adjacency"}, \code{"laplacian"} [default] and
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
get_hamming_distance <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  n <- nrow(x)
  sum(abs(x - y)) / (n * (n - 1))
}

#' @rdname get-distance
#' @export
get_frobenius_distance <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  sqrt(sum((x - y)^2))
}

#' @rdname get-distance
#' @export
get_spectral_distance <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  dlX <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  dlY <- eigen(y, symmetric = TRUE, only.values = TRUE)$values

  sqrt(sum((dlX - dlY)^2))
}

#' @rdname get-distance
#' @export
get_root_euclidean_distance <- function(x, y, representation = "laplacian") {
  if (representation != "lapacian")
    stop("The root-Euclidean distance can only be used with the Laplacian matrix representation.")

  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  rX <- eigen(x, symmetric = TRUE)
  vals <- rX$values
  vals[vals < 0] <- 0
  dlX <- rX$vectors %*% diag(sqrt(vals)) %*% t(rX$vectors)

  rY <- eigen(y, symmetric = TRUE)
  vals <- rY$values
  vals[vals < 0] <- 0
  dlY <- rY$vectors %*% diag(sqrt(vals)) %*% t(rY$vectors)

  sqrt(sum((dlX - dlY)^2))
}
