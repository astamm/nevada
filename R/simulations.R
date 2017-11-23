# Utility functions -------------------------------------------------------

rpois_network <- function(lambda, n) {
  A <- diag(0, n)
  A[upper.tri(A)] <- rpois(n * (n - 1L) / 2L, lambda)
  A + t(A)
}

# Scenario 0 --------------------------------------------------------------

get_scenario0_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(5, n), simplify = FALSE)
  list(x = x, y = y)
}

# Scenario A --------------------------------------------------------------

get_scenarioA_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(6, n), simplify = FALSE)
  list(x = x, y = y)
}

# Scenario B --------------------------------------------------------------

get_scenarioB_dataset <- function(n1, n2 = n1) {
  p1 <- matrix(data = c(0.8, 0.2, 0.2, 0.2), nrow = 2L, ncol = 2L)
  p2 <- matrix(data = c(0.2, 0.2, 0.2, 0.8), nrow = 2L, ncol = 2L)
  x <- replicate(n1, igraph::sample_sbm(n = 25L, pref.matrix = p1, block.sizes = c(12L, 13L)), simplify = FALSE)
  y <- replicate(n2, igraph::sample_sbm(n = 25L, pref.matrix = p2, block.sizes = c(13L, 12L)), simplify = FALSE)
  list(x = x, y = y)
}

# Scenario C --------------------------------------------------------------

get_scenarioC_dataset <- function(n1, n2 = n1) {
  x <- replicate(n1, igraph::sample_k_regular(no.of.nodes = 25L, k = 8L), simplify = FALSE)
  y <- replicate(n2, igraph::sample_gnp(n = 25L, p = 1/3), simplify = FALSE)
  list(x = x, y = y)
}

# Scenario D --------------------------------------------------------------

get_scenarioD_dataset <- function(n1, n2 = n1) {
  x <- replicate(n1, igraph::sample_smallworld(dim = 1L, size = 25L, nei = 4L, p = 0.15), simplify = FALSE)
  y <- replicate(n2, igraph::sample_pa(n = 25L, power = 2L, m = 4L, directed = FALSE), simplify = FALSE)
  list(x = x, y = y)
}

# Scenario E --------------------------------------------------------------

get_scenarioE_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(20, n), simplify = FALSE)
  list(x = x, y = y)
}

# Test --------------------------------------------------------------------

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

  test <- network_test2p(
    data$x, data$y,
    representation = representation,
    distance = distance,
    statistic = statistic,
    verbose = FALSE,
    alpha = alpha
  )

  test$pvalue
}

# Power -------------------------------------------------------------------

#' @export
get_power <- function(scenario, n_pop, representation, distance, statistic, alpha = 0.05, R = 10L, seed = NULL) {
  set.seed(seed)
  replicate(R, perform_single_test(scenario, n_pop, representation, distance, statistic, alpha = alpha))
}
