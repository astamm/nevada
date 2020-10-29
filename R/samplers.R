#' Graph samplers using edge distributions
#'
#' A collection of functions to generate random graphs with specified edge
#' distribution.
#'
#' @param n Sample size.
#' @param num_vertices Number of vertices.
#' @param lambda The mean parameter for the Poisson distribution (default: 1).
#' @param rate The rate parameter for the exponential distribution (default: 1).
#' @param size The number of trials for the binomial distribution (default: 1).
#' @param prob The probability of success on each trial for the binomial
#'   distribution (default: 0.5).
#' @param pref_matrix The matrix giving the Bernoulli rates. This is a KxK
#'   matrix, where K is the number of groups. The probability of creating an
#'   edge between vertices from groups i and j is given by element (i,j). For
#'   undirected graphs, this matrix must be symmetric. See
#'   \code{\link[igraph]{sample_sbm}}.
#' @param block_sizes Numeric vector giving the number of vertices in each
#'   group. The sum of the vector must match the number of vertices. See
#'   \code{\link[igraph]{sample_sbm}}.
#'
#' @return A object of class \code{\link{nvd}} containing the sample of graphs.
#' @name samplers
#'
#' @examples
#' nvd <- rexp_network(10, 68)
NULL

#' @rdname samplers
#' @export
rpois_network <- function(n, num_vertices, lambda = 1) {
  as_nvd(replicate(n, {
    A <- diag(0, num_vertices)
    A[upper.tri(A)] <- stats::rpois(
      n = num_vertices * (num_vertices - 1L) / 2L,
      lambda = lambda
    )
    igraph::graph_from_adjacency_matrix(A, mode = "upper", weighted = TRUE)
  }, simplify = FALSE))
}

#' @rdname samplers
#' @export
rexp_network <- function(n, num_vertices, rate = 1) {
  as_nvd(replicate(n, {
    A <- diag(0, num_vertices)
    A[upper.tri(A)] <- stats::rexp(
      n = num_vertices * (num_vertices - 1L) / 2L,
      rate = rate
    )
    igraph::graph_from_adjacency_matrix(A, mode = "upper", weighted = TRUE)
  }, simplify = FALSE))
}

#' @rdname samplers
#' @export
rbinom_network <- function(n, num_vertices, size = 1, prob = 0.5) {
  as_nvd(replicate(n, {
    A <- diag(0, num_vertices)
    A[upper.tri(A)] <- stats::rbinom(
      n = num_vertices * (num_vertices - 1L) / 2L,
      size = size,
      prob = prob
    )
    igraph::graph_from_adjacency_matrix(A, mode = "upper", weighted = TRUE)
  }, simplify = FALSE))
}

#' @rdname samplers
#' @export
rsbm <- function(n, num_vertices, pref_matrix, block_sizes) {
  stopifnot(sum(block_sizes) == num_vertices)
  replicate(n, {
    igraph::sample_sbm(
      n = num_vertices,
      pref.matrix = pref_matrix,
      block.sizes = block_sizes
    )
  }, simplify = FALSE) %>%
    purrr::map(igraph::set_edge_attr, name = "weight", value = 1)
}

sample1_nbs <- function(n, num_vertices, rate = 1) {
  num_vertices <- num_vertices / 2
  as_nvd(
    replicate(n, {
      A11 <- diag(0, num_vertices)
      A11[upper.tri(A11)] <- stats::rexp(
        n = num_vertices * (num_vertices - 1L) / 2L,
        rate = 1
      )
      A22 <- diag(0, num_vertices)
      A22[upper.tri(A22)] <- stats::rexp(
        n = num_vertices * (num_vertices - 1L) / 2L,
        rate = 1
      )
      A12 <- matrix(stats::rexp(
        n = num_vertices * num_vertices,
        rate = rate
      ), nrow = num_vertices, ncol = num_vertices)
      A <- cbind(rbind(A11, A12), rbind(A12, A22))
      igraph::graph_from_adjacency_matrix(A, mode = "upper", weighted = TRUE)
    }, simplify = FALSE)
  )
}

sample2_nbs <- function(n, num_vertices, rate1 = 0.5, rate2 = 3) {
  sample1 <-  sample1_nbs(n, num_vertices, rate1)
  sample2 <-  sample1_nbs(n, num_vertices, rate2)
  sample3 <- rexp_network(n, num_vertices, rate2)
  sample4 <- rexp_network(n, num_vertices, rate1)

  for (i in 1:n) {
    igraph::V(sample1[[i]])$name <- as.character(1:num_vertices)
    igraph::V(sample2[[i]])$name <- as.character(1:num_vertices)
    igraph::V(sample3[[i]])$name <- as.character((num_vertices + 1):(2 * num_vertices))
    igraph::V(sample4[[i]])$name <- as.character((num_vertices + 1):(2 * num_vertices))
  }

  # Disjoint union of graphs
  list(
    x = as_nvd(purrr::map2(sample1, sample3, igraph::`%du%`)),
    y = as_nvd(purrr::map2(sample2, sample4, igraph::`%du%`))
  )
}
