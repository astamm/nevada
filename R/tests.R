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
#' @param x An \code{\link{nvd}} object listing networks in sample 1.
#' @param y An \code{\link{nvd}} object listing networks in sample 2.
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
                           distance = "hamming",
                           statistic = "mod",
                           B = 1000L,
                           alpha = 0.05,
                           test = "exact",
                           verbose = TRUE) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  representation <-
    match.arg(representation,
              c("adjacency", "laplacian", "modularity", "transitivity"))
  distance <-
    match.arg(distance,
              c("hamming", "frobenius", "spectral", "root-euclidean"))
  statistic <-
    match.arg(statistic,
              c("mod", "dom", "sdom", "dom-frobenius", "sdom-frobenius"))

  if (statistic %in% c("dom-frobenius", "sdom-frobenius"))
    d <- repr_nvd(x, y, representation = representation)
  else
    d <- dist_nvd(x, y, representation = representation, distance = distance)

  T0 <- switch(
    statistic,
    "mod" = stat_mod(d, 1:n1),
    "dom" = stat_dom(d, 1:n1, standardize = FALSE),
    "sdom" = stat_dom(d, 1:n1, standardize = TRUE),
    "dom-frobenius" = stat_dom_frobenius(d, 1:n1, standardize = FALSE),
    "sdom-frobenius" = stat_dom_frobenius(d, 1:n1, standardize = TRUE)
  )

  if (B < 1)
    B <- (qnorm(alpha / 2, lower.tail = FALSE) / tol) ^ 2

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

  Tp <- sapply(1:B, get_permuted_statistic, group1.perm, d, statistic)

  if (test == "approximate")
    p <- mean(Tp >= T0)
  else {
    b <- sum(Tp >= T0)
    p <- phipson_smyth_pvalue(b, B, M)
  }

  list(
    statistic = T0,
    pvalue = p,
    permuted_statistics = Tp
  )
}
