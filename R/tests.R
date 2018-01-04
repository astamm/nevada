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
#' @param x An \code{\link{nvd}} object listing networks in sample 1.
#' @param y An \code{\link{nvd}} object listing networks in sample 2.
#' @param representation A string specifying the desired type of representation,
#'   among: \code{"adjacency"} [default], \code{"laplacian"} and
#'   \code{"modularity"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"}, \code{"frobenius"} [default],
#'   \code{"spectral"} and \code{"root-euclidean"}.
#' @param statistic A string specifying the chosen test statistic(s), among:
#'   \code{"lot"} [default], \code{"sot"}, \code{"biswas"}, \code{"energy"},
#'   \code{"student"}, \code{"welch"}, \code{"original"}, \code{"generalized"},
#'   \code{"weighted"} or a combination from \code{c("lot", "sot", "biswas",
#'   "energy")}.
#' @param B The number of permutation or the tolerance (default: \code{1000L}).
#'   If this number is lower than \code{1}, it is intended as a tolerance.
#'   Otherwise, it is intended as the number of required permutations.
#' @param alpha The significance level (default: \code{0.05}).
#' @param test A character string specifying if performing an exact test through
#'   the use of Phipson-Smyth estimate of the p-value or an approximate test
#'   through a Monte-Carlo estimate of the p-value (default: \code{"exact"}).
#' @param k An integer specifying the density of the minimum spanning tree used
#'   for the edge count statistics (default: \code{5L}).
#'
#' @return A \code{\link[base]{list}} with three components: the value of the
#'   statistic for the original two samples, the p-value of the resulting
#'   permutation test and a numeric vector storing the values of the permuted
#'   statistics.
#' @export
#'
#' @examples
#' n <- 10L
#'
#' # Two different models for the two populations
#' x <- nvd("smallworld", n)
#' y <- nvd("pa", n)
#' test1 <- test_twosample(x, y, "modularity")
#' test1$pvalue
#'
#' # Same model for the two populations
#' x <- nvd("smallworld", n)
#' y <- nvd("smallworld", n)
#' test2 <- test_twosample(x, y, "modularity")
#' test2$pvalue
test_twosample <- function(x,
                           y,
                           representation = "adjacency",
                           distance = "frobenius",
                           statistic = "lot",
                           B = 1000L,
                           alpha = 0.05,
                           test = "exact",
                           k = 5L) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  representation <- match.arg(
    representation,
    c("adjacency", "laplacian", "modularity", "transitivity")
  )

  distance <- match.arg(
    distance,
    c("hamming", "frobenius", "spectral", "root-euclidean")
  )

  npc <- length(statistic) > 1

  if (!npc) {
    if (statistic %in% c("student", "welch"))
      d <- repr_nvd(x, y, representation = representation)
    else {
      d <- dist_nvd(x, y, representation = representation, distance = distance)
      if (statistic %in% c("original", "generalized", "weighted"))
        d <- edge_count_global_variables(d, n1, k = k)
    }
  } else
    d <- dist_nvd(x, y, representation = representation, distance = distance)

  if (B < 1)
    B <- (qnorm(alpha / 2, lower.tail = FALSE) / tol)^2

  M <- choose(n, n1)
  if (n1 == n2)
    M <- M / 2

  test <- match.arg(test, c("approximate", "exact"))
  if (test == "approximate" & M <= B) {
    B <- M
    group1.perm <- combn(n, n1)[, 1:B]
  } else
    group1.perm <- replicate(B, sample.int(n))[1:n1, ]

  if (!npc)
    Tp <- sapply(0:B, get_permuted_statistic, indices1 = group1.perm, d = d, statistic = statistic)
  else {
    Tp <- statistic %>%
      purrr::map(~ sapply(0:B, get_permuted_statistic, indices1 = group1.perm, d = d, statistic = .)) %>%
      purrr::map(~ sapply(1:(B+1), stats2pvalue, Tp = ., test = test, B = B, M = M)) %>%
      purrr::transpose() %>%
      purrr::simplify_all() %>%
      purrr::map_dbl(combine_pvalues)
  }

  list(
    statistic = Tp[1],
    pvalue = stats2pvalue(1, Tp, test, B, M),
    permuted_statistics = Tp[-1]
  )
}

stats2pvalue <- function(i, Tp, test = "exact", B, M) {
  T0 <- Tp[i]
  Tp <- Tp[-i]
  if (test == "approximate")
    p <- mean(Tp >= T0)
  else {
    b <- sum(Tp >= T0)
    p <- phipson_smyth_pvalue(b, B, M)
  }
}

combine_pvalues <- function(p, method = "tippett") {
  switch (
    method,
    tippett = 1 - min(p),
    fisher = - 2 * sum(log(p))
  )
}
