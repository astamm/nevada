rpois_network <- function(lambda, n) {
  A <- diag(0, n)
  A[upper.tri(A)] <- rpois(n * (n - 1L) / 2L, lambda)
  igraph::graph_from_adjacency_matrix(A + t(A), mode = "undirected")
}

get_scenario0_dataset <- function(n1, n2 = n1) {
  x <- nvd("poisson", 10, lambda = 5)
  y <- nvd("poisson", 10, lambda = 5)
  list(x = x, y = y)
}

get_scenarioA_dataset <- function(n1, n2 = n1) {
  x <- nvd("poisson", 10, lambda = 5)
  y <- nvd("poisson", 10, lambda = 6)
  list(x = x, y = y)
}

get_scenarioB_dataset <- function(n1, n2 = n1) {
  p1 <- matrix(data = c(0.8, 0.2, 0.2, 0.2), nrow = 2L, ncol = 2L)
  p2 <- matrix(data = c(0.2, 0.2, 0.2, 0.8), nrow = 2L, ncol = 2L)
  x <- nvd("sbm", n1, pref.matrix = p1)
  y <- nvd("sbm", n2, pref.matrix = p2)
  list(x = x, y = y)
}

get_scenarioC_dataset <- function(n1, n2 = n1) {
  x <- nvd("k_regular", n1)
  y <- nvd("gnp", n2)
  list(x = x, y = y)
}

get_scenarioD_dataset <- function(n1, n2 = n1) {
  x <- nvd("smallworld", n1)
  y <- nvd("pa", n2)
  list(x = x, y = y)
}

get_scenarioE_dataset <- function(n1, n2 = n1) {
  x <- nvd("poisson", 10, lambda = 5)
  y <- nvd("poisson", 10, lambda = 20)
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

  test <- test_twosample(
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
#' p <- power_twosample(alpha = 0.05)
#' mean(p <= alpha)
power_twosample <- function(
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
