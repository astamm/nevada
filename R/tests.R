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
#'   among: \code{"adjacency"}, \code{"laplacian"} and \code{"modularity"}.
#'   Defaults to \code{"adjacency"}.
#' @param distance A string specifying the chosen distance for calculating the
#'   test statistic, among: \code{"hamming"}, \code{"frobenius"},
#'   \code{"spectral"} and \code{"root-euclidean"}. Defaults to
#'   \code{"frobenius"}.
#' @param stats A character vector specifying the chosen test statistic(s),
#'   among: `"original_edge_count"`, `"generalized_edge_count"`,
#'   `"weighted_edge_count"`, `"student_euclidean"`, `"welch_euclidean"` or any
#'   statistics based on inter-point distances available in the **flipr**
#'   package: `"flipr:student_ip"`, `"flipr:fisher_ip"`, `"flipr:bg_ip"`,
#'   `"flipr:energy_ip"`, `"flipr:cq_ip"`. Defaults to `c("flipr:student_ip",
#'   "flipr:fisher_ip")`.
#' @param B The number of permutation or the tolerance. If this number is lower
#'   than \code{1}, it is intended as a tolerance. Otherwise, it is intended as
#'   the number of required permutations. Defaults to `1000L`.
#' @param test A character string specifying the formula to be used to compute
#'   the permutation p-value. Choices are `"estimate"`, `"upper_bound"` and
#'   `"exact"`. Defaults to `"exact"` which provides exact tests.
#' @param k An integer specifying the density of the minimum spanning tree used
#'   for the edge count statistics. Defaults to `5L`.
#' @param seed An integer for specifying the seed of the random generator for
#'   result reproducibility. Defaults to `NULL`.
#'
#' @return A \code{\link[base]{list}} with three components: the value of the
#'   statistic for the original two samples, the p-value of the resulting
#'   permutation test and a numeric vector storing the values of the permuted
#'   statistics.
#' @export
#'
#' @examples
#' n <- 10L
#' gnp_params <- list(p = 1/3)
#' k_regular_params <- list(k = 8L)
#'
#' # Two different models for the two populations
#' x <- nvd(model = "gnp", n = n, model_params = gnp_params)
#' y <- nvd(model = "k_regular", n = n, model_params = k_regular_params)
#' t1 <- test2_global(x, y, representation = "modularity")
#' t1$pvalue
#'
#' # Same model for the two populations
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' y <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' t2 <- test2_global(x, y, representation = "modularity")
#' t2$pvalue
test2_global <- function(x, y,
                         representation = "adjacency",
                         distance = "frobenius",
                         stats = c("flipr:t_ip", "flipr:f_ip"),
                         B = 1000L,
                         test = "exact",
                         k = 5L,
                         seed = NULL) {
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

  use_frechet_stats <- any(grepl("student_euclidean", stats)) ||
    any(grepl("welch_euclidean", stats))
  if (use_frechet_stats &&
      (any(grepl("_ip", stats)) ||
       any(grepl("edge_count", stats))))
    cli::cli_abort("It is not possible to mix statistics based on Frechet means and statistics based on inter-point distances.")

  ecp <- NULL
  if (use_frechet_stats)
    d <- repr_nvd(x, y, representation = representation)
  else {
    d <- dist_nvd(x, y, representation = representation, distance = distance)
    if (any(grepl("edge_count", stats)))
      ecp <- edge_count_global_variables(d, n1, k = k)
  }

  null_spec <- function(y, parameters) {
    return(y)
  }

  stat_functions <- stats %>%
    strsplit(split = ":") %>%
    purrr::map(~ {
      if (length(.x) == 1) {
        s <- paste0("stat_", .x)
        return(rlang::as_function(s))
      }
      s <- paste0("stat_", .x[2])
      getExportedValue(.x[1], s)
    })

  stat_assignments <- list(delta = 1:length(stat_functions))

  if (inherits(d, "dist")) {
    xx <- d
    yy <- as.integer(n1)
  } else {
    xx <- d[1:n1]
    yy <- d[(n1 + 1):(n1 + n2)]
  }

  pf <- flipr::PlausibilityFunction$new(
    null_spec = null_spec,
    stat_functions = stat_functions,
    stat_assignments = stat_assignments,
    xx, yy,
    seed = seed
  )
  pf$set_nperms(B)
  pf$set_pvalue_formula(test)
  pf$set_alternative("right_tail")

  pf$get_value(
    parameters = 0,
    edge_count_prep = ecp,
    keep_null_distribution = TRUE
  )
}

#' Local Two-Sample Test for Network-Valued Data
#'
#' @inheritParams test2_global
#' @param partition Either a list or an integer vector specifying vertex
#'   memberships into partition elements.
#' @param alpha Significance level for hypothesis testing. If set to 1, the
#'   function outputs properly adjusted p-values. If lower than 1, then only
#'   p-values lower than alpha are properly adjusted. Defaults to `0.05`.
#' @param verbose Boolean specifying whether information on intermediate tests
#'   should be printed in the process (default: \code{FALSE}).
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
#' test2_local(sim$x, sim$y, m,
#'             seed = 1234,
#'             alpha = 0.05,
#'             B = 100)
test2_local <- function(x, y, partition,
                        representation = "adjacency",
                        distance = "frobenius",
                        stats = c("flipr:t_ip", "flipr:f_ip"),
                        B = 1000L,
                        alpha = 0.05,
                        test = "exact",
                        k = 5L,
                        seed = NULL,
                        verbose = FALSE) {

  # Creating sigma-algebra generated by the partition
  partition <- as_vertex_partition(partition)
  E <- names(partition)
  sa <- generate_sigma_algebra(partition)
  psize <- length(sa)

  # Initialize output for intra-adjusted pvalues
  stop_intra <- FALSE
  skip_intra <- NULL
  p_intra <- utils::combn(E, 1, simplify = FALSE) %>%
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    rlang::set_names("E") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(pvalue = 0, truncated = FALSE)

  # Intialize output for inter-adjusted pvalues
  stop_inter <- FALSE
  skip_inter <- NULL
  p_inter <- utils::combn(E, 2, simplify = FALSE) %>%
    purrr::transpose() %>%
    purrr::simplify_all() %>%
    rlang::set_names(c("E1", "E2")) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(pvalue = 0, truncated = FALSE)

  for (i in 1:psize) {

    sas <- sa[[i]]
    compositions <- names(sas)

    for (j in 1:length(sas)) {

      if (stop_intra && stop_inter)
        return(list(intra = p_intra, inter = p_inter))

      element_name <- compositions[j]

      update_intra <- !stop_intra && !(element_name %in% skip_intra)
      update_inter <- !stop_inter && i < psize && !(element_name %in% skip_inter)

      if (!update_intra && !update_intra)
        next()

      element_value <- sas[[j]]
      individuals <- element_name %>%
        strsplit(",") %>%
        purrr::simplify()

      # Tests on full subgraphs
      p <- test2_subgraph(
        x, y, element_value,
        subgraph_full,
        representation, distance, stats, B, test, k, seed
      )

      if (verbose) {
        writeLines("- Type of test: FULL")
        writeLines(paste0("Element of the sigma-algebra: ", element_name))
        writeLines(paste0("P-value of the test: ", p))
      }

      # Intra-adjusted p-values from full tests
      if (update_intra)
        p_intra <- .update_intra_pvalues(p_intra, individuals, p, alpha)

      # Inter-adjusted p-values from full tests
      if (update_inter)
        p_inter <- .update_inter_pvalues(p_inter, individuals, p, alpha)

      # Update stopping and skipping conditions
      stop_intra <- all(p_intra$truncated)
      stop_inter <- all(p_inter$truncated)
      if (p >= alpha) {
        skip_intra <- .update_skip_list(skip_intra, individuals)
        skip_inter <- .update_skip_list(skip_inter, individuals)
      }

      update_intra <- !stop_intra && !(element_name %in% skip_intra)

      if (update_intra) {
        # Tests on intra subgraphs
        p <- test2_subgraph(
          x, y, element_value,
          subgraph_intra,
          representation, distance, stats, B, test, k, seed
        )

        if (verbose) {
          writeLines("- Type of test: INTRA")
          writeLines(paste0("Element of the sigma-algebra: ", element_name))
          writeLines(paste0("P-value of the test: ", p))
        }

        # Intra-adjusted p-values from intra tests
        p_intra <- .update_intra_pvalues(p_intra, individuals, p, alpha)

        # Update stopping and skipping conditions
        stop_intra <- all(p_intra$truncated)
        if (p >= alpha)
          skip_intra <- .update_skip_list(skip_intra, individuals)
      }

      update_inter <- !stop_inter && i < psize && !(element_name %in% skip_inter)

      if (update_inter) {
        # Tests on inter subgraphs
        p <- test2_subgraph(
          x, y, element_value,
          subgraph_inter,
          representation, distance, stats, B, test, k, seed
        )

        if (verbose) {
          writeLines("- Type of test: INTER")
          writeLines(paste0("Element of the sigma-algebra: ", element_name))
          writeLines(paste0("P-value of the test: ", p))
        }

        # Inter-adjusted p-values from inter tests
        p_inter <- .update_inter_pvalues(p_inter, individuals, p, alpha)

        # Update stopping and skipping conditions
        stop_inter <- all(p_inter$truncated)
        if (p >= alpha)
          skip_inter <- .update_skip_list(skip_inter, individuals)
      }
    }
  }

  list(intra = p_intra, inter = p_inter)
}

.update_intra_pvalues <- function(output, c, p, alpha) {
  output %>%
    dplyr::mutate(
      pvalue = purrr::map2_dbl(.data$E, .data$pvalue, ~ dplyr::if_else(.x %in% c, pmax(.y, p), .y)),
      truncated = .data$pvalue >= alpha
    )
}

.update_inter_pvalues <- function(output, c, p, alpha) {
  output %>%
    dplyr::mutate(
      pvalue = purrr::pmap_dbl(
        list(.data$E1, .data$E2, .data$pvalue),
        ~ dplyr::if_else(all(c(..1, ..2) %in% c), pmax(..3, p), ..3)
      ),
      truncated = .data$pvalue >= alpha
    )
}

.update_skip_list <- function(skip_list, individuals) {
  for (k in 1:length(individuals)) {
    tmp <- individuals %>%
      utils::combn(k, paste0, collapse = ",", simplify = FALSE) %>%
      purrr::simplify()
    skip_list <- unique(c(skip_list, tmp))
  }
  skip_list
}

test2_subgraph <- function(x, y, subpartition, fun,
                           representation, distance, stats, B, test, k, seed) {
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
    stats = stats,
    B = B,
    test = test,
    k = k,
    seed = seed
  )$pvalue
}
