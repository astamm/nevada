rpois_network <- function(lambda, n) {
  A <- diag(0, n)
  A[upper.tri(A)] <- rpois(n * (n - 1L) / 2L, lambda)
  igraph::graph_from_adjacency_matrix(A, mode = "upper", weighted = TRUE)
}

rbinom_network <- function(size, prob, n) {
  A <- diag(0, n)
  A[upper.tri(A)] <- rbinom(n * (n - 1L) / 2L, size, prob)
  igraph::graph_from_adjacency_matrix(A, mode = "upper", weighted = TRUE)
}

# Equal distributions
get_scenario10_dataset <- function(n1, n2 = n1) {
  x <- nvd("binomial", n1, size = 15, prob = 1/3)
  y <- nvd("binomial", n2, size = 15, prob = 1/3)
  list(x = x, y = y)
}

# Distributions with equal means, different variances
get_scenario11_dataset <- function(n1, n2 = n1) {
  x <- nvd("binomial", n1, size = 15, prob = 1/3)
  y <- nvd("binomial", n2, size = 20, prob = 1/4)
  list(x = x, y = y)
}

# Distributions with different means, equal variances
get_scenario12_dataset <- function(n1, n2 = n1) {
  x <- nvd("binomial", n1, size = 15, prob = 31/64)
  y <- nvd("binomial", n2, size = 15, prob = 33/64)
  list(x = x, y = y)
}

# Distributions with different means, different variances
get_scenario13_dataset <- function(n1, n2 = n1) {
  x <- nvd("binomial", n1, size = 15, prob = 1/5)
  y <- nvd("binomial", n2, size = 15, prob = 1/6)
  list(x = x, y = y)
}

get_scenario0_dataset <- function(n1, n2 = n1) {
  x <- nvd("poisson", n1, lambda = 5)
  y <- nvd("poisson", n2, lambda = 5)
  list(x = x, y = y)
}

get_scenarioA_dataset <- function(n1, n2 = n1) {
  x <- nvd("poisson", n1, lambda = 5)
  y <- nvd("poisson", n2, lambda = 6)
  list(x = x, y = y)
}

get_scenarioB_dataset <- function(n1, n2 = n1) {
  p1 <- matrix(data = c(0.8, rep(0.2, 8L)), nrow = 3L)
  p2 <- matrix(data = c(rep(0.2, 8L), 0.8), nrow = 3L)
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
  x <- nvd("poisson", n1, lambda = 5)
  y <- nvd("poisson", n2, lambda = 20)
  list(x = x, y = y)
}

perform_single_test <- function(scenario,
                                n_pop,
                                representation,
                                distance,
                                statistic,
                                alpha = 0.05,
                                test = "exact",
                                k = 5L) {
  data <- switch(
    scenario,
    "0" = get_scenario0_dataset(n_pop),
    "A" = get_scenarioA_dataset(n_pop),
    "B" = get_scenarioB_dataset(n_pop),
    "C" = get_scenarioC_dataset(n_pop),
    "D" = get_scenarioD_dataset(n_pop),
    "E" = get_scenarioE_dataset(n_pop),
    "10" = get_scenario10_dataset(n_pop),
    "11" = get_scenario11_dataset(n_pop),
    "12" = get_scenario12_dataset(n_pop),
    "13" = get_scenario13_dataset(n_pop)
  )

  test_data <- test_twosample(
    data$x,
    data$y,
    representation = representation,
    distance = distance,
    statistic = statistic,
    alpha = alpha,
    test = test,
    k = k,
    verbose = FALSE
  )

  test_data$pvalue
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
#' @inheritParams test_twosample
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
power_twosample <- function(scenario = "0",
                            n_pop = 4L,
                            representation = "adjacency",
                            distance = "frobenius",
                            statistic = "mod",
                            alpha = 0.05,
                            test = "exact",
                            k = 5L,
                            R = 1000L,
                            seed = NULL) {
  set.seed(seed)

  pvalues <- replicate(
    R,
    perform_single_test(
      scenario,
      n_pop,
      representation,
      distance,
      statistic,
      alpha = alpha,
      test = test,
      k = k
    )
  )

  mean(pvalues <= 0.05)
}
