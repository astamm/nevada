#' Comparison of Network Distributions
#'
#' This function carries out an hypothesis test where the null hypothesis is
#' that the two populations of networks share the same underlying probabilistic
#' distribution against the alternative hypothesis that the two populations come
#' from different distributions. The test is performed in a non-parametric
#' fashion using a permutational framework in which several statistics can be
#' used, together with several choices of network matrix representations and
#' distances between networks.
#'
#' The exhaustive list of all possible permutations is used if its length is
#' smaller than either the user-defined number of random permutations requested
#' or the number of permutations required for achieving a user-defined p-value
#' resolution. Otherwise, the p-value is estimated via Monte-Carlo simulations.
#'
#' @param x A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects
#'   or matrix representations of underlying networks from a given first
#'   population.
#' @param y A \code{\link[base]{list}} of \code{\link[igraph]{igraph}} objects
#'   or matrix representations of underlying networks from a given second
#'   population.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"} [default], \code{"laplacian"} and
#'   \code{"modularity"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"} [default], \code{"frobenius"},
#'   \code{"spectral"} and \code{"root-euclidean"}.
#' @param statistic A string specifying the chosen test statistic, among:
#'   \code{"mod"} [default], \code{"dom"} or \code{"sdom"}.
#' @param B The number of permutation or the tolerance (default: \code{1000L}).
#'   If this number is lower than \code{1}, it is intended as a tolerance.
#'   Otherwise, it is intended as the number of required permutations.
#' @param alpha The significance level (default: \code{0.05}).
#' @param test A character string specifying if performing an exact test through
#'   the use of Phipson Smyth p-value or an approximate test through a
#'   Monte-Carlo p-value (default: \code{"exact"}).
#' @param verbose A boolean to switch on/off verbosity (information about
#'   p-value resolution).
#'
#' @return A \code{\link[base]{list}} with three components: the value of the
#'   statistic for the original two samples, the p-value of the resulting
#'   permutation test and a numeric vector storing the values of the permutated
#'   statistics.
#' @export
#'
#' @examples
#' n <- 25L
#'
#' x <- list()
#' y <- list()
#' for (i in 1:10) {
#'   x[[i]] <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#'   y[[i]] <- igraph::barabasi.game(n, m = 3, power = 2, directed = FALSE)
#' }
#' test1 <- network_test2p(x, y, "modularity")
#' test1
#'
#' x <- list()
#' y <- list()
#' for (i in 1:10) {
#'   x[[i]] <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#'   y[[i]] <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#' }
#'
#' test2 <- network_test2p(x, y, "modularity")
#' test2
network_test2p <- function(
  x, y,
  representation = "adjacency", distance = "hamming", statistic = "mod",
  B = 1000L, alpha = 0.05, test = "exact", verbose = TRUE) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  representation <- match.arg(representation, c("adjacency", "laplacian", "modularity"))
  distance <- match.arg(distance, c("hamming", "frobenius", "spectral", "root-euclidean"))
  d <- get_distance_matrix(x, y, representation, distance)

  statistic <- match.arg(statistic, c("mod", "dom", "sdom"))
  T0 <- switch(
    statistic,
    mod = get_mod_statistic(d, 1:n1),
    dom = get_dom_statistic(d, 1:n1),
    sdom = get_sdom_statistic(d, 1:n1)
  )

  if (B < 1)
    B <- (qnorm(alpha / 2, lower.tail = FALSE) / tol)^2

  if (!reachable_significance(n1, n2, B, alpha, verbose))
    warning("The requested significance level cannot be
            reached given the sample sizes.")

  M <- choose(n, n1)
  if (n1 == n2)
    M <- M / 2

  test <- match.arg(test, c("approximate", "exact"))
  if (test == "approximate" & M <= B) {
    B <- M
    group1.perm <- combn(n, n1)[, 1:B]
  } else
    group1.perm <- replicate(B, sample.int(n))[1:n1, ]

  Tp <- switch(
    statistic,
    mod = purrr::map_dbl(1:B, ~ get_mod_statistic(d, group1.perm[, .])),
    dom = purrr::map_dbl(1:B, ~ get_dom_statistic(d, group1.perm[, .])),
    sdom = purrr::map_dbl(1:B, ~ get_sdom_statistic(d, group1.perm[, .]))
  )

  if (test == "approximate")
    p <- mean(Tp >= T0)
  else {
    b <- sum(Tp >= T0)
    p <- phipson_smyth_pvalue(b, B, M)
  }

  list(statistic = T0, pvalue = p, permuted_statistics = Tp)
}
