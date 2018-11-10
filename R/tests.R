#' Global Two-Sample Test for Network-Valued Data
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
#' @param test A character string specifying if performing an exact test through
#'   the use of Phipson-Smyth estimate of the p-value or an approximate test
#'   through a Monte-Carlo estimate of the p-value (default: \code{"exact"}).
#' @param k An integer specifying the density of the minimum spanning tree used
#'   for the edge count statistics (default: \code{5L}).
#' @param seed An integer for specifying the seed of the random generator for
#'   result reproducibility (default: \code{NULL}).
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
#' t1 <- test2_global(x, y, "modularity")
#' t1$pvalue
#'
#' # Same model for the two populations
#' x <- nvd("smallworld", n)
#' y <- nvd("smallworld", n)
#' t2 <- test2_global(x, y, "modularity")
#' t2$pvalue
test2_global <- function(x, y,
                         representation = "adjacency",
                         distance = "frobenius",
                         statistic = "lot",
                         B = 1000L,
                         test = "exact",
                         k = 5L,
                         seed = NULL) {

  set.seed(seed)
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

  M <- choose(n, n1)
  if (n1 == n2)
    M <- M / 2

  test <- match.arg(test, c("approximate", "exact"))
  if (test == "approximate" & M <= B) {
    B <- M
    group1.perm <- utils::combn(n, n1)[, 1:B]
  } else
    group1.perm <- replicate(B, sample.int(n))[1:n1, ]

  if (!npc)
    Tp <- sapply(0:B, get_permuted_statistic, indices1 = group1.perm, d = d, statistic = statistic)
  else {
    Tp <- statistic %>%
      purrr::map(~ sapply(0:B, get_permuted_statistic, indices1 = group1.perm, d = d, statistic = .)) %>%
      purrr::map(~ sapply(1:(B+1), stats2pvalue, Tp = ., test = "approximate", B = B, M = M)) %>%
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
  b <- sum(Tp >= T0) - 1
  if (test == "approximate") return(b / B)
  phipson_smyth_pvalue(b, B, M)
}

combine_pvalues <- function(p, method = "tippett") {
  switch (
    method,
    tippett = 1 - min(p),
    fisher = - 2 * sum(log(p))
  )
}

#' Local Two-Sample Test for Network-Valued Data
#'
#' @inheritParams test2_global
#' @param partition Either a list or a vector specifying vertex memberships into
#'   partition elements.
#' @param alpha Significance level for hypothesis testing (default:
#'   \code{0.05}). If set to 1, the function outputs properly adjusted p-values.
#'   If lower than 1, then only p-values lower than alpha are properly adjusted.
#'
#' @return A length-2 list reporting the adjusted p-values of each element of
#'   the partition for the intra- and inter-tests.
#' @export
#'
#' @examples
#' n <- 10
#' p1 <- matrix(
#'   data = c(0.1, 0.4, 0.1, 0.4,
#'            0.4, 0.4, 0.1, 0.4,
#'            0.1, 0.1, 0.4, 0.4,
#'            0.4, 0.4, 0.4, 0.4),
#'   nrow = 4,
#'   ncol = 4,
#'   byrow = TRUE
#' )
#' p2 <- matrix(
#'   data = c(0.1, 0.4, 0.4, 0.4,
#'            0.4, 0.4, 0.4, 0.4,
#'            0.4, 0.4, 0.1, 0.1,
#'            0.4, 0.4, 0.1, 0.4),
#'   nrow = 4,
#'   ncol = 4,
#'   byrow = TRUE
#' )
#' set.seed(1234)
#' sim <- n %>%
#'   purrr::rerun(
#'     x = igraph::sample_sbm(68, p1, c(17,17,17,17)),
#'     y = igraph::sample_sbm(68, p2, c(17,17,17,17))
#'   ) %>%
#'   purrr::transpose() %>%
#'   purrr::map(as_nvd)
#' t0 <- test2_local(sim$x, sim$y, m, statistic = c("lot", "sot"), seed = 1234)
test2_local <- function(x, y, partition,
                        representation = "adjacency",
                        distance = "frobenius",
                        statistic = "lot",
                        B = 1000L,
                        alpha = 0.05,
                        test = "exact",
                        k = 5L,
                        seed = NULL) {

  # Creating sigma-albebra generated by the partition
  sa <- generate_sigma_algebra(partition)

  # Tests on full subgraphs
  p_full <- test2_subgraph(
    x, y, sa,
    subgraph_full,
    representation, distance, statistic, B, test, k, seed
  )

  # Tests on intra subgraphs
  p_intra <- test2_subgraph(
    x, y, sa,
    subgraph_intra,
    representation, distance, statistic, B, test, k, seed
  )

  # Intra-adjusted p-values
  p_full %>%
    dplyr::left_join(p_intra, by = "Elements") %>%
    dplyr::mutate(p_intra = pmax(pvalue.x, pvalue.y)) %>%
    dplyr::select(Elements, p_intra)
}

test2_subgraph <- function(x, y, sa, fun, representation, distance, statistic, B, test, k, seed) {
  x <- x %>%
    purrr::map(function(g) {
      sa %>%
        purrr::modify_depth(2, rlang::as_function(fun), g = g) %>%
        purrr::flatten()
    }) %>%
    purrr::transpose() %>%
    tibble::as_tibble() %>%
    tidyr::gather(Elements, Graphs) %>%
    dplyr::group_by(Elements) %>%
    dplyr::do(x = as_nvd(.$Graphs)) %>%
    dplyr::ungroup()
  y <- y %>%
    purrr::map(function(g) {
      sa %>%
        purrr::modify_depth(2, rlang::as_function(fun), g = g) %>%
        purrr::flatten()
    }) %>%
    purrr::transpose() %>%
    tibble::as_tibble() %>%
    tidyr::gather(Elements, Graphs) %>%
    dplyr::group_by(Elements) %>%
    dplyr::do(y = as_nvd(.$Graphs)) %>%
    dplyr::ungroup()
  p <- x %>%
    dplyr::left_join(y) %>%
    dplyr::mutate(pvalue = purrr::map2(
      x, y,
      test2_global,
      representation = representation,
      distance = distance,
      statistic = statistic,
      B = B,
      test = test,
      k = k,
      seed = seed
    ) %>% purrr::map_dbl("pvalue")) %>%
    dplyr::select(-x, -y) %>%
    dplyr::mutate(Elements = stringr::str_split(Elements, ",")) %>%
    tidyr::unnest() %>%
    dplyr::group_by(Elements) %>%
    dplyr::summarise(pvalue = max(pvalue, na.rm = TRUE))
}
