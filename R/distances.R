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
#'
#' @param x An [`igraph::igraph`] object or a matrix representing an underlying
#'   network.
#' @param y An [`igraph::igraph`] object or a matrix representing an underlying
#'   network. Should have the same number of vertices as `x`.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"}, \code{"laplacian"}, \code{"modularity"} or
#'   \code{"graphon"}. Default is \code{"laplacian"}.
#' @param matching_iterations An integer value specifying the maximum number of
#'   runs when looking for the optimal permutation for graph matching. Defaults
#'   to `0L` in which case no matching is done.
#' @param target_matrix A square numeric matrix of size `n` equal to the order
#'   of the graphs specifying a target matrix towards which the initial doubly
#'   stochastic matrix is shrunk each time the graph matching algorithm fails to
#'   provide a good minimum. Defaults to `NULL` in which case the target matrix
#'   is automatically chosen between the identity matrix or the uniform matrix
#'   on the n-simplex.
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
dist_frobenius <- function(x, y,
                           representation = "laplacian",
                           matching_iterations = 0,
                           target_matrix = NULL) {
  if (!compatible_networks(x, y))
    cli::cli_abort("Input networks are incompatible.")

  if (matching_iterations > 0) {
    d0 <- dist_frobenius(x, y, representation = representation)
    Ax <- repr_adjacency(x, validate = FALSE)
    Ay <- repr_adjacency(y, validate = FALSE)

    if (is.null(target_matrix)) {
      n <- igraph::gorder(x)
      target_matrix <- diag(1, n)
      l <- align_graphs(x, y, Ax, Ay, target_matrix = target_matrix, alpha = 1)
      dy <- dist_frobenius(x, l$y, representation = representation)
      dx <- dist_frobenius(l$x, y, representation = representation)
      d1 <- min(dx, dy)
      unif_matrix <- matrix(1, n, n) / n
      l <- align_graphs(x, y, Ax, Ay, target_matrix = unif_matrix, alpha = 1)
      dy <- dist_frobenius(x, l$y, representation = representation)
      dx <- dist_frobenius(l$x, y, representation = representation)
      d2 <- min(dx, dy)
      if (d2 < d1)
        target_matrix <- unif_matrix
      d0 <- min(d0, d1, d2)
    } else {
      l <- align_graphs(x, y, Ax, Ay, target_matrix = target_matrix, alpha = 1)
      dy <- dist_frobenius(x, l$y, representation = representation)
      dx <- dist_frobenius(l$x, y, representation = representation)
      d0 <- min(d0, dx, dy)
    }

    dP <- d0 + 1
    cnt <- 0
    while (dP > d0 && cnt <= matching_iterations) {
      l <- align_graphs(x, y, Ax, Ay, target_matrix = target_matrix, alpha = cnt / matching_iterations)
      dy <- dist_frobenius(x, l$y, representation = representation)
      dx <- dist_frobenius(l$x, y, representation = representation)
      dP <- min(dx, dy)
      cnt <- cnt + 1
    }

    if (cnt > matching_iterations) {
      cli::cli_alert_info("The maximal number of iterations ({matching_iterations}) has been reached.")
      dP <- d0
    }

    return(dP)
  }

  x <- format_input(x, representation)
  y <- format_input(y, representation)

  dist_frobenius_impl(x, y)
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
#' is denoted by \eqn{n_1} and that of the second one is denoted by \eqn{n_2}.
#'
#' @param x A [`base::list`] of [`igraph::igraph`] objects or matrix
#'   representations of underlying networks from a given first population.
#' @param y A [`base::list`] of [`igraph::igraph`] objects or matrix
#'   representations of underlying networks from a given second population.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"}, \code{"laplacian"}, \code{"modularity"} or
#'   \code{"graphon"}. Default is \code{"laplacian"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"}, \code{"frobenius"},
#'   \code{"spectral"} and \code{"root-euclidean"}. Default is
#'   \code{"frobenius"}.
#' @param matching_iterations An integer value specifying the maximum number of
#'   runs when looking for the optimal permutation for graph matching. Defaults
#'   to `0L` in which case no matching is done.
#' @param target_matrix A square numeric matrix of size `n` equal to the order
#'   of the graphs specifying a target matrix towards which the initial doubly
#'   stochastic matrix is shrunk each time the graph matching algorithm fails to
#'   provide a good minimum. Defaults to `NULL` in which case the target matrix
#'   is automatically chosen between the identity matrix or the uniform matrix
#'   on the n-simplex.
#'
#' @return A matrix of dimension \eqn{(n_1+n_2) \times (n_1+n_2)} containing the
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
                     matching_iterations = 0,
                     target_matrix = NULL) {
  if (matching_iterations > 0 && distance == "frobenius") {
    if (!is.null(y)) x <- c(x, y)
    n <- length(x)
    labels <- 1:n
    indices <- linear_index(n)
    D <- furrr::future_map_dbl(indices$k, ~ {
      i <- indices$i[.x]
      j <- indices$j[.x]
      dist_frobenius(
        x[[i]], x[[j]],
        representation = representation,
        matching_iterations = matching_iterations,
        target_matrix = target_matrix
      )
    }, .options = furrr::furrr_options(seed = TRUE))

    attributes(D) <- NULL
    attr(D, "Labels") <- labels
    attr(D, "Size") <- n
    attr(D, "call") <- match.call()
    class(D) <- "dist"
    return(D)
  }

  x <- repr_nvd(x, y, representation = representation)
  dist_nvd_impl(x, distance)
}
