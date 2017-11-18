#######################################
### Simulations in Biometrika paper ###
#######################################

# Scenario 0 --------------------------------------------------------------

rpois_network <- function(lambda, n) {
  A <- diag(0, n)
  A[upper.tri(A)] <- rpois(n * (n - 1L) / 2L, lambda)
  A + t(A)
}

rblocks_network <- function(pois_indices) {
  m <- n * (n - 1L) / 2L
  p <- length(pois_indices)
  if (p > m)
    stop("Too many indices.")

  weights <- numeric(m)
  weights[pois_indices] <- rpois(p, 5)
  weights[-pois_indices] <- rbinom(m - p, 1L, 0.2)

  A <- diag(0, n)
  A[upper.tri(A)] <- weights
  A + t(A)
}

get_scenario0_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(5, n), simplify = FALSE)
  list(x = x, y = y)
}

get_scenario0_datasets <- function(n1, n2 = n1, times = 1L, seed = NULL) {
  set.seed(seed)
  replicate(times, get_scenario0_dataset(n1, n2), simplify = FALSE)
}

# Scenario A --------------------------------------------------------------

get_scenarioA_dataset <- function(n1, n2 = n1) {
  n <- 25L
  x <- replicate(n1, rpois_network(5, n), simplify = FALSE)
  y <- replicate(n2, rpois_network(6, n), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioA_datasets <- function(n1, n2 = n1, times = 1L, seed = NULL) {
  set.seed(seed)
  replicate(times, get_scenarioA_dataset(n1, n2), simplify = FALSE)
}

# Scenario B --------------------------------------------------------------

get_scenarioB_dataset <- function(n1, n2 = n1) {
  n <- 8L
  x <- replicate(n1, rblocks_network(1:6), simplify = FALSE)
  y <- replicate(n2, rblocks_network(c(15, 20:21, 26:28)), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioB_datasets <- function(n1, n2 = n1, times = 1L, seed = NULL) {
  set.seed(seed)
  replicate(times, get_scenarioB_dataset(n1, n2), simplify = FALSE)
}

# Scenario C --------------------------------------------------------------

get_scenarioC_dataset <- function(n1, n2 = n1) {
  x <- replicate(n1, igraph::sample_k_regular(no.of.nodes = 25L, k = 8L), simplify = FALSE)
  y <- replicate(n2, igraph::erdos.renyi.game(n = 25L, p.or.m = 1/3, type = "gnp"), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioC_datasets <- function(n1, n2 = n1, times = 1L, seed = NULL) {
  set.seed(seed)
  replicate(times, get_scenarioC_dataset(n1, n2), simplify = FALSE)
}

# Scenario D --------------------------------------------------------------

get_scenarioD_dataset <- function(n1, n2 = n1) {
  x <- replicate(n1, igraph::watts.strogatz.game(dim = 1L, size = 25L, nei = 4L, p = 0.15), simplify = FALSE)
  y <- replicate(n2, igraph::barabasi.game(n = 25L, power = 2L, m = 4L, directed = FALSE), simplify = FALSE)
  list(x = x, y = y)
}

get_scenarioD_datasets <- function(n1, n2 = n1, times = 1L, seed = NULL) {
  set.seed(seed)
  replicate(times, get_scenarioD_dataset(n1, n2), simplify = FALSE)
}

# Plot simulations --------------------------------------------------------

plot_simulation <- function(df) {
  df %>%
    ggplot2::ggplot(ggplot2::aes(
      x = size,
      y = pvalue_mean,
      color = representation,
      linetype = statistic,
      shape = distance
    )) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_hline(ggplot2::aes(
      yintercept = alpha
    )) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Sample size") +
    ggplot2::ylab("Estimated Power") +
    ggplot2::scale_y_continuous(limits = c(0, 1))
}
