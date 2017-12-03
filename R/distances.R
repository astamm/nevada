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
#' Root-Euclidean distance can used only with the laplacian matrix
#' representation.
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
#' g1 <- igraph::sample_gnp(20, 0.1)
#' g2 <- igraph::sample_gnp(20, 0.2)
#' dist_hamming(g1, g2, "adjacency")
#' dist_frobenius(g1, g2, "adjacency")
#' dist_spectral(g1, g2, "laplacian")
#' dist_root_euclidean(g1, g2, "laplacian")
#' @name distances
NULL

#' @rdname distances
#' @export
dist_hamming <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  internal_hamming(x, y)
}

#' @rdname distances
#' @export
dist_frobenius <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  internal_frobenius(x, y)
}

#' @rdname distances
#' @export
dist_spectral <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  internal_spectral(x, y)
}

#' @rdname distances
#' @export
dist_root_euclidean <- function(x, y, representation = "laplacian") {
  if (representation != "laplacian")
    stop("The root-Euclidean distance can only be used with the Laplacian matrix representation.")

  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  internal_root_euclidean(x, y)
}

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
#' x <- nvd("smallworld", 10)
#' y <- nvd("pa", 10)
#' dist_nvd(x, y, "adjacency", "spectral")
dist_nvd <- function(x, y = NULL, representation = "adjacency", distance = "hamming") {
  x <- lapply(x, format_input, representation)

  if (!is.null(y)) {
    y <- lapply(y, format_input, representation)
    x <- c(x, y)
  }

  internal_distance_matrix(x, distance)
}