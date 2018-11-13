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
#' sim <- sample2_sbm(n, 68, p1, c(17, 17, 17, 17), p2, seed = 1234)
#' m <- as.integer(c(rep(1, 17), rep(2, 17), rep(3, 17), rep(4, 17)))
#' t1 <- test2_local(sim$x, sim$y, m, statistic = c("lot", "sot"), seed = 1234, alpha = 1.00)
#' t2 <- test2_local(sim$x, sim$y, m, statistic = c("lot", "sot"), seed = 1234, alpha = 0.05)
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
  partition <- as_vertex_partition(partition)
  E <- names(partition)
  sa <- generate_sigma_algebra(partition)
  psize <- length(sa)

  # Initialize output for intra-adjusted pvalues
  stop_intra <- FALSE
  p_intra <- utils::combn(E, 1, simplify = FALSE) %>%
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    rlang::set_names("E") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(pvalue = 0, truncated = FALSE)

  # Intialize output for inter-adjusted pvalues
  stop_inter <- FALSE
  p_inter <- utils::combn(E, 2, simplify = FALSE) %>%
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    rlang::set_names(c("E1", "E2")) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(pvalue = 0, truncated = FALSE)

  skip_intra <- NULL
  skip_inter <- NULL

  for (i in 1:psize) {
    sas <- sa[[i]]
    compositions <- names(sas)

    for (j in 1:length(sas)) {

      if (stop_intra && stop_inter)
        return(list(intra = p_intra, inter = p_inter))

      if (compositions[j] %in% skip_intra && compositions[j] %in% skip_inter)
        next()

      individuals <- compositions[j] %>%
        stringr::str_split(",") %>%
        purrr::simplify()

      # Tests on full subgraphs
      p <- test2_subgraph(
        x, y, sas[[j]],
        subgraph_full,
        representation, distance, statistic, B, test, k, seed
      )

      # Grab combinations that do not require update anymore
      if (p >= alpha) {
        for (k in 1:length(individuals)) {
          tmp <- individuals %>%
            utils::combn(k, paste0, collapse = ",", simplify = FALSE) %>%
            purrr::simplify()
          skip_intra <- c(skip_intra, tmp) %>% unique()
          skip_inter <- c(skip_inter, tmp) %>% unique()
        }
      }

      # Intra-adjusted p-values from full tests
      if (!stop_intra)
        p_intra <- .update_intra_pvalues(p_intra, individuals, p, alpha)
      stop_intra <- all(p_intra$truncated)

      # Inter-adjusted p-values from full tests
      if (!stop_inter && i < psize)
        p_inter <- .update_inter_pvalues(p_inter, individuals, p, alpha)
      stop_inter <- all(p_inter$truncated)

      if (!stop_intra) {
        if (compositions[j] %in% skip_intra)
          next()

        # Tests on intra subgraphs
        p <- test2_subgraph(
          x, y, sas[[j]],
          subgraph_intra,
          representation, distance, statistic, B, test, k, seed
        )

        # Grab combinations that do not require update anymore
        if (p >= alpha) {
          for (k in 1:length(individuals)) {
            tmp <- individuals %>%
              utils::combn(k, simplify = FALSE) %>%
              purrr::simplify()
            skip_intra <- c(skip_intra, tmp) %>% unique()
          }
        }

        # Intra-adjusted p-values from intra tests
        p_intra <- .update_intra_pvalues(p_intra, individuals, p, alpha)
        stop_intra <- all(p_intra$truncated)
      }

      if (!stop_inter && i < psize) {
        if (compositions[j] %in% skip_inter)
          next()

        # Tests on inter subgraphs
        p <- test2_subgraph(
          x, y, sas[[j]],
          subgraph_inter,
          representation, distance, statistic, B, test, k, seed
        )

        # Grab combinations that do not require update anymore
        if (p >= alpha) {
          for (k in 2:length(individuals)) {
            tmp <- individuals %>%
              utils::combn(k, simplify = FALSE) %>%
              purrr::simplify()
            skip_inter <- c(skip_inter, tmp) %>% unique()
          }
        }

        # Inter-adjusted p-values from inter tests
        p_inter <- .update_inter_pvalues(p_inter, individuals, p, alpha)
        stop_inter <- all(p_inter$truncated)
      }
    }
  }

  list(intra = p_intra, inter = p_inter)
}

.update_intra_pvalues <- function(output, c, p, alpha) {
  output %>%
    dplyr::mutate(
      pvalue = purrr::map2_dbl(E, pvalue, ~ dplyr::if_else(.x %in% c, pmax(.y, p), .y)),
      truncated = pvalue >= alpha
    )
}

.update_inter_pvalues <- function(output, c, p, alpha) {
  output %>%
    dplyr::mutate(
      pvalue = purrr::pmap_dbl(
        list(E1, E2, pvalue),
        ~ dplyr::if_else(all(c(..1, ..2) %in% c), pmax(..3, p), ..3)
      ),
      truncated = pvalue >= alpha
    )
}

test2_subgraph <- function(x, y, subpartition, fun,
                           representation, distance, statistic, B, test, k, seed) {
  x <- x %>%
    purrr::map(rlang::as_function(fun), vids = subpartition) %>%
    as_nvd()
  y <- y %>%
    purrr::map(rlang::as_function(fun), vids = subpartition) %>%
    as_nvd()
  test2_global(
    x, y,
    representation = representation,
    distance = distance,
    statistic = statistic,
    B = B,
    test = test,
    k = k,
    seed = seed
  )$pvalue
}
