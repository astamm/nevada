#' Two-sample test of the NBS framework
#'
#' @inheritParams test2_global
#' @param threshold Treshold to discard edges with low values of the test
#'   statistic when searching for the connected components (default: \code{0.95}
#'   which means 95th quantile of the null distribution of the test statistic).
#'
#' @return A list containing information about the identified connected
#'   components and corresponding adjusted p-values.
#' @export
#'
#' @examples
#' x <- replicate(10, igraph::sample_gnp(10, 1 / 10), simplify = FALSE)
#' y <- replicate(10, igraph::sample_gnp(10, 1 / 10), simplify = FALSE)
#' test2_nbs(x, y)
test2_nbs <- function(x, y, threshold = 0.95, B = 1000, seed = NULL) {
  set.seed(seed)
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  nv <- igraph::vcount(x[[1]])

  M <- choose(n, n1)
  if (n1 == n2)
    M <- M / 2

  cobs <- .test2_nbs(x, y, threshold)

  # Permutations for finding the null distribution of the maximal component size
  perms <- replicate(B, sample.int(n), simplify = FALSE) %>%
    purrr::map(~c(x, y)[.x]) %>%
    purrr::map(~.test2_nbs(.x[1:n1], .x[(n1+1):n], threshold)) %>%
    purrr::map("csize") %>%
    purrr::map_dbl(max)

  pvals <- cobs$csize %>%
    purrr::map_dbl(~mean(perms >= .x))

  c(cobs, list(pvalues = pvals))
}

.test2_nbs <- function(x, y, threshold) {
  nv <- igraph::vcount(x[[1]])

  x <- x %>%
    purrr::map(repr_adjacency) %>%
    purrr::map(~.x[upper.tri(.x)]) %>%
    purrr::transpose() %>%
    purrr::simplify_all()

  y <- y %>%
    purrr::map(repr_adjacency) %>%
    purrr::map(~.x[upper.tri(.x)]) %>%
    purrr::transpose() %>%
    purrr::simplify_all()

  statval <- purrr::map2_dbl(x, y, ~ (mean(.x) - mean(.y))^2)
  q <- stats::quantile(statval, threshold)
  statval[statval < q] <- 0

  statmat <- matrix(nrow = nv, ncol = nv)
  statmat[upper.tri(statmat)] <- statval

  g <- igraph::graph_from_adjacency_matrix(statmat, mode = "upper", weighted = TRUE, diag = FALSE)
  igraph::components(g)
}

#' Power estimation of the NBS framework
#'
#' @inheritParams test2_nbs
#' @param R Number of simulation runs (default: \code{100L}).
#' @param n_obs Number of observations in the two samples (default: \code{10L}).
#' @param p1 Probability of the presence of an edge in the first sample
#'   (default: \code{0.1}).
#' @param p2 Probability of the presence of an edge in the second sample
#'   (default: \code{0.1}).
#' @param n_nodes Number of vertices in the network structure (default:
#'   \code{10L}).
#' @param alpha Significance level (default: \code{0.05}).
#'
#' @return A matrix of size \code{(n_nodes x n_nodes)} reporting the estimated
#'   probability of rejection.
#' @export
#'
#' @examples
#' power2_nbs(R = 10L)
power2_nbs <- function(
  R = 100,
  n_obs = 10,
  p1 = 1 / 10, p2 = p1, n_nodes = 10,
  threshold = 0.95,
  B = 1000,
  alpha = 0.05,
  seed = NULL
) {
  set.seed(seed)
  sim <- replicate(R, .power2_nbs_gnp(n_obs, p1, p2, n_nodes, threshold, B), simplify = FALSE)
  sim %>%
    purrr::map(.vec2mat) %>%
    purrr::map(`<=`, alpha) %>%
    simplify2array() %>%
    apply(1:2, mean, na.rm = TRUE)
}

.power2_nbs_gnp <- function(n_obs, p1, p2, n_nodes, threshold, B) {
  x <- replicate(n_obs, igraph::sample_gnp(n_nodes, p1), simplify = FALSE)
  y <- replicate(n_obs, igraph::sample_gnp(n_nodes, p2), simplify = FALSE)
  test2_nbs(x, y, threshold, B)
}

.vec2mat <- function(cobs) {
  nv <- sum(cobs$csize)
  pmat <- matrix(NA, nv, nv)
  for (j in 1:length(cobs$csize)) {
    pmat[cobs$membership == j, cobs$membership == j] <- cobs$pvalues[j]
  }
  pmat
}
