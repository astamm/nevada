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
#' UBU^{-1},} where \eqn{V} and \eqn{U} are the matrices whose columns are the
#' eigenvectors of \eqn{X} and \eqn{Y}, respectively and \eqn{A} and \eqn{B} are
#' the diagonal matrices with elements the eigenvalues of \eqn{X} and \eqn{Y},
#' respectively. The root-Euclidean distance between \eqn{x} and \eqn{y} is
#' given by \deqn{\sqrt{\sum_i (V \sqrt{A} V^{-1} - U \sqrt{B} U^{-1})^2}.}
#' Root-Euclidean distance can used only with the laplacian matrix
#' representation.
#' The match-Frobenius distance is the Frobenius distance considering networks in the Graph Space (alignment approach).
#' Match-Frobenius distance can be used only with adjacency matrix representation.
#'
#' @param x An \code{\link[igraph]{igraph}} object or a matrix representing an
#'   underlying network.
#' @param y An \code{\link[igraph]{igraph}} object or a matrix representing an
#'   underlying network. Should have the same number of vertices as \code{x}.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"}, \code{"laplacian"}, \code{"modularity"} or
#'   \code{"graphon"}. Default is \code{"laplacian"}.
#' @param A string specifying the initialization of the permutation matrix estimate, among: \code{"barycenter"} and \code{"identity"}. Default is
#'   \code{"barycenter"}. It is required for \code{distance == "match-frobenius"}.
#' @param iteration The number of iterations for the Frank-Wolfe algorithm. Default to 20L. It is required for \code{distance == "match-frobenius"}.
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

  dist_hamming_impl(x, y)
}

#' @rdname distances
#' @export
dist_frobenius <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  dist_frobenius_impl(x, y)
}

#' @rdname distances
#' @export
dist_match_frobenius <- function(x, y, representation = "adjacency", start = "barycenter", iteration = 20) {
  if (representation != "adjacency")
    stop("The match-Frobenius distance can only be used with the Adjacency matrix representation.")

  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  num_nodes <- nrow(repr_adjacency(y))

  if (start == "barycenter"){
    start_mat <- matrix(1, num_nodes, num_nodes)/num_nodes
  } else if (start == "identity"){
    start_mat <- diag(rep(1, num_nodes))
  } else stop("Input start matrix is not valid. The initialization of the permutation matrix estimate should be among: barycenter and identity.")


  x <- format_input(x, representation)
  y <- format_input(y, representation)

  distanceValue = -1

  m <- 0 #number of nodes I know a priori that are in correspondence in A and B
  mat1 <- as.matrix(x)
  mat2 <- as.matrix(y)
  perm <- igraph::match_vertices(A = mat1, B = mat2, m = m, start = start_mat, iteration = iteration)
  P <- perm$P
  Pmat2P <- P%*%mat2%*%t(P)

  distanceValue <- dist_frobenius(mat1, Pmat2P, representation = "adjacency")
  distanceValue
}

#' @rdname distances
#' @export
dist_spectral <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  dist_spectral_impl(x, y)
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

  dist_root_euclidean_impl(x, y)
}

#' Inner-Products Between Networks
#'
#' This is a collection of functions computing the inner product between two
#' networks.
#'
#' @param x An \code{\link[igraph]{igraph}} object or a matrix representing an
#'   underlying network.
#' @param y An \code{\link[igraph]{igraph}} object or a matrix representing an
#'   underlying network. Should have the same number of vertices as \code{x}.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"}, \code{"laplacian"}, \code{"modularity"} or
#'   \code{"graphon"}. Default is \code{"laplacian"}.
#'
#' @return A scalar measuring the angle between the two input networks.
#'
#' @examples
#' g1 <- igraph::sample_gnp(20, 0.1)
#' g2 <- igraph::sample_gnp(20, 0.2)
#' ipro_frobenius(g1, g2, "adjacency")
#' @name inner-products
NULL

#' @rdname inner-products
#' @export
ipro_frobenius <- function(x, y, representation = "laplacian") {
  if (!compatible_networks(x, y))
    stop("Input networks are incompatible.")

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  ipro_frobenius_impl(x, y)
}

#' Pairwise Distance Matrix Between Two Samples of Networks
#'
#' This function computes the matrix of pairwise distances between all the
#' elements of the two samples put together. The cardinality of the fist sample
#' is denoted by \eqn{n1} and that of the second one is denoted by \eqn{n2}.
#'
#' @param x A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects
#'   or matrix representations of underlying networks from a given first
#'   population.
#' @param y A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects
#'   or matrix representations of underlying networks from a given second
#'   population.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"}, \code{"laplacian"}, \code{"modularity"} or
#'   \code{"graphon"}. Default is \code{"laplacian"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"}, \code{"frobenius"},
#'   \code{"spectral"}, \code{"root-euclidean"} and \code{"match-frobenius"}. Default is
#'   \code{"frobenius"}.
#' @param start A string specifying the initialization of the permutation matrix estimate, among: \code{"barycenter"} and \code{"identity"}. Default is
#'   \code{"barycenter"}. It is required for \code{distance == "match-frobenius"}.
#' @param iteration The number of iterations for the Frank-Wolfe algorithm. Default to 20L. It is required for \code{distance == "match-frobenius"}.
#'
#' @return A matrix of dimension \eqn{(n1+n2) \times (n1+n2)} containing the
#'   distances between all the elements of the two samples put together.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' k_regular_params <- list(k = 8L)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' y <- nvd(model = "k_regular", n = 10L, model_params = k_regular_params)
#' dist_nvd(x, y, "adjacency", "spectral")
dist_nvd <- function(x,
                     y = NULL,
                     representation = "adjacency",
                     distance = "frobenius",
                     start = "barycenter",
                     iteration = 20L) {
  if(distance == "match-frobenius"){

    x <- append(x, y)
    n <- length(x)
    mat_d <- rep(-1, n * (n - 1) / 2)

    fun <- function(k) {
      j <-floor((2*(n-1)+1-sqrt((2*(n-1)+1)*(2*(n-1)+1)-8*k))/2)
      i <- k-((2*(n-1)-1-j)*j)/2

      net1 = x[[i+2]]
      net2 = x[[j+1]]

      distanceValue <- dist_match_frobenius(net1, net2, representation = representation, start = start, iteration = iteration)
      distanceValue
    }

    mat_d <- sapply(0:((n * (n - 1) / 2)-1), fun)
    attr(mat_d,"class") <- "dist"
    attr(mat_d,"Size") <- n
    attr(mat_d,"Diag") <- FALSE
    attr(mat_d,"Upper") <- FALSE
    mat_d
    } else{
    x <- repr_nvd(x, y, representation = representation)
    dist_nvd_impl(x, distance)
  }
}
