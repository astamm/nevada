rpois_network <- function(lambda, n) {
  A <- diag(0, n)
  A[upper.tri(A)] <- rpois(n * (n - 1L) / 2L, lambda)
  A + t(A)
}

get_scenario0_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(5, n), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioA_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(6, n), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioB_dataset <- function(n1, n2 = n1) {
  p1 <- matrix(data = c(0.8, 0.2, 0.2, 0.2), nrow = 2L, ncol = 2L)
  p2 <- matrix(data = c(0.2, 0.2, 0.2, 0.8), nrow = 2L, ncol = 2L)
  x <- replicate(n1, igraph::sample_sbm(n = 25L, pref.matrix = p1, block.sizes = c(12L, 13L)), simplify = FALSE)
  y <- replicate(n2, igraph::sample_sbm(n = 25L, pref.matrix = p2, block.sizes = c(13L, 12L)), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioC_dataset <- function(n1, n2 = n1) {
  x <- replicate(n1, igraph::sample_k_regular(no.of.nodes = 25L, k = 8L), simplify = FALSE)
  y <- replicate(n2, igraph::sample_gnp(n = 25L, p = 1/3), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioD_dataset <- function(n1, n2 = n1) {
  x <- replicate(n1, igraph::sample_smallworld(dim = 1L, size = 25L, nei = 4L, p = 0.15), simplify = FALSE)
  y <- replicate(n2, igraph::sample_pa(n = 25L, power = 2L, m = 4L, directed = FALSE), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioE_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(20, n), simplify = FALSE)
  list(x = x, y = y)
}

perform_single_test <- function(scenario, n_pop, representation, distance, statistic, alpha = 0.05) {
  data <- switch(
    scenario,
    "0" = get_scenario0_dataset(n_pop),
    "A" = get_scenarioA_dataset(n_pop),
    "B" = get_scenarioB_dataset(n_pop),
    "C" = get_scenarioC_dataset(n_pop),
    "D" = get_scenarioD_dataset(n_pop),
    "E" = get_scenarioE_dataset(n_pop)
  )

  test <- twosample_test(
    data$x, data$y,
    representation = representation,
    distance = distance,
    statistic = statistic,
    verbose = FALSE,
    alpha = alpha
  )

  test$pvalue
}

#' Power Simulations for Permutation Tests
#'
#' This function provides a Monte-Carlo estimate of the power of the permutation
#' tests proposed in this package.
#'
#' Currently, six scenarios of pairs of populations are implemented. Scenario 0
#' allows to make sure that all our permutation tests are exact.
#'
#' @param scenario A character specifying the scenario to be simulated (choices
#'   are: "0" [default], "A", "B", "C", "D" or "E").
#' @param n_pop The sample size of each population (the same is used for both
#'   populations, default: 4L)
#' @param representation The chosen network representation as a string (choices
#'   are: "adjacency" [default], "laplacian" or "modularity").
#' @param distance The chosen distance between networks as a string (choices
#'   are: "hamming", "frobenius" [default], "spectral" or "root-euclidean").
#' @param statistic The chosen test statistic as a string (choices are: "mod"
#'   [default], "sdom" or "dom").
#' @param alpha The significance level of the test (default: 5\%).
#' @param R The number of Monte-Carlo runs used to estimate the power (default:
#'   1000L).
#' @param seed An integer specifying the seed to start randomness from (default:
#'   uses clock).
#'
#' @return A vector \code{p} of \code{R} p-values, one for each Monte-Carlo run.
#'   The power can then be estimated as \code{mean(p <= alpha)}.
#' @export
#'
#' @examples
#' alpha <- 0.05
#' p <- twosample_power(alpha = 0.05)
#' mean(p <= alpha)
twosample_power <- function(
  scenario = "0",
  n_pop = 4L,
  representation = "adjacency",
  distance = "frobenius",
  statistic = "mod",
  alpha = 0.05,
  R = 1000L,
  seed = NULL) {
  set.seed(seed)
  replicate(R, perform_single_test(scenario, n_pop, representation, distance, statistic, alpha = alpha))
}
