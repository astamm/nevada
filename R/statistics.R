#' Test Statistics for Network Populations
#'
#' This is a collection of functions that provide statistics for testing
#' equality in distribution between samples of networks.
#'
#' In details, there are three main categories of statistics:
#'
#' - *Euclidean t-Statistics*: both Student `stat_student_euclidean` version for
#' equal variances and Welch `stat_welch_euclidean` version for unequal
#' variances,
#' - *Statistics based on similarity graphs*: 3 types of edge count statistics.
#'
#' @param d Either a matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the
#'   distances between all the elements of the two samples put together (for
#'   distance-based statistics) or the concatenation of the lists of matrix
#'   representations of networks in samples 1 and 2 for Euclidean t-Statistics.
#' @param indices A vector of dimension \eqn{n1} containing the indices of the
#'   elements of the first sample.
#' @param edge_count_prep A list of preprocessed data information used by edge
#'   count statistics and produced by \code{\link{edge_count_global_variables}}.
#' @param ... Extra parameters specific to some statistics.
#'
#' @return A scalar giving the value of the desired test statistic.
#'
#' @examples
#' n1 <- 30L
#' n2 <- 10L
#' x <- nvd("smallworld", n1)
#' y <- nvd("pa", n2)
#' r <- repr_nvd(x, y, representation = "laplacian")
#' stat_student_euclidean(r, 1:n1)
#' stat_welch_euclidean(r, 1:n1)
#' d <- dist_nvd(x, y, representation = "laplacian", distance = "frobenius")
#' ecp <- edge_count_global_variables(d, n1, k = 5L)
#' stat_original_edge_count(d, 1:n1, edge_count_prep = ecp)
#' stat_generalized_edge_count(d, 1:n1, edge_count_prep = ecp)
#' stat_weighted_edge_count(d, 1:n1, edge_count_prep = ecp)
#' @name statistics
NULL

#' @rdname statistics
#' @export
stat_student_euclidean <- function(d, indices, ...) {
  x <- d[indices]
  y <- d[-indices]
  stat_t_euclidean_impl(x, y, pooled = TRUE)
}

#' @rdname statistics
#' @export
stat_welch_euclidean <- function(d, indices, ...) {
  x <- d[indices]
  y <- d[-indices]
  stat_t_euclidean_impl(x, y, pooled = FALSE)
}

stat_edge_count <- function(d, indices, edge_count_prep, edge_count_type = "generalized") {
  E <- edge_count_prep$E
  nE <- edge_count_prep$nE
  n1 <- edge_count_prep$n1
  n2 <- edge_count_prep$n2
  mu0 <- edge_count_prep$mu0
  mu1 <- edge_count_prep$mu1
  mu2 <- edge_count_prep$mu2
  V0 <- edge_count_prep$V0
  V1 <- edge_count_prep$V1
  V2 <- edge_count_prep$V2
  V12 <- edge_count_prep$V12
  Sinv <- edge_count_prep$Sinv

  temp <- stat_edge_count_impl(E, indices)
  R1 <- temp[1]
  R2 <- temp[2]

  switch(edge_count_type,
    original = -(nE - R1 - R2 - mu0) / sqrt(V0),
    generalized = Sinv[1,1] * (R1 - mu1)^2 + Sinv[2, 2] * (R2 - mu2)^2 + 2 * Sinv[1, 2] * (R1 - mu1) * (R2 - mu2),
    weighted = (n2 * (R1 - mu1) + n1 * (R2 - mu2)) / sqrt(n2^2 * V1 + n1^2 * V2 + 2 * n2 * n1 * V12)
  )
}

#' @rdname statistics
#' @export
stat_original_edge_count <- function(d, indices, edge_count_prep, ...) {
  stat_edge_count(d, indices, edge_count_prep, edge_count_type = "original")
}

#' @rdname statistics
#' @export
stat_generalized_edge_count <- function(d, indices, edge_count_prep, ...) {
  stat_edge_count(d, indices, edge_count_prep, edge_count_type = "generalized")
}

#' @rdname statistics
#' @export
stat_weighted_edge_count <- function(d, indices, edge_count_prep, ...) {
  stat_edge_count(d, indices, edge_count_prep, edge_count_type = "weighted")
}

#' Transform distance matrix in edge properties of minimal spanning tree
#'
#' @param d A matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the distances
#'   between all the elements of the two samples put together.
#' @param n1 An integer giving the size of the first sample.
#' @param k An integer specifying the density of the minimal spanning tree to
#'   generate.
#'
#' @return A list of edge properties of the minimal spanning tree.
#' @export
#'
#' @examples
#' n1 <- 30L
#' n2 <- 10L
#' x <- nvd("smallworld", n1)
#' y <- nvd("pa", n2)
#' d <- dist_nvd(x, y, representation = "laplacian", distance = "frobenius")
#' e <- edge_count_global_variables(d, n1, k = 5L)
edge_count_global_variables <- function(d, n1, k = 1L) {
  n <- attr(d, "Size")
  n2 <- n - n1
  k <- min(k, ceiling(n / 4))
  g <- kmst(as.matrix(d), k)
  E <- igraph::as_edgelist(g, names = FALSE)
  Ebynode <- vector("list", n)
  for (i in 1:nrow(E)) {
    Ebynode[[E[i, 1]]] <- c(Ebynode[[E[i, 1]]], E[i, 2])
    Ebynode[[E[i, 2]]] <- c(Ebynode[[E[i, 2]]], E[i, 1])
  }
  nE <- nrow(E)
  nodedeg <- sapply(Ebynode, length)
  nEi <- sum(nodedeg * (nodedeg - 1))
  mu0 <- nE * 2 * n1 * n2 / n / (n - 1)
  mu1 <- nE * n1 * (n1 - 1) / n / (n - 1)
  mu2 <- nE * n2 * (n2 - 1) / n / (n - 1)
  V0 <- nEi * n1 * n2 / n / (n - 1) + (nE * (nE - 1) - nEi) * 4 * n1 * n2 * (n1 - 1) * (n2 - 1) / n / (n - 1) / (n - 2) / (n - 3) + mu0 - mu0^2
  V1 <- nEi * n1 * (n1 - 1) * (n1 - 2) / n / (n - 1) / (n - 2) + (nE * (nE - 1) - nEi) * n1 * (n1 - 1) * (n1 - 2) * (n1 - 3) / n / (n - 1) / (n - 2) / (n - 3) + mu1 - mu1^2
  V2 <- nEi * n2 * (n2 - 1) * (n2 - 2) / n / (n - 1) / (n - 2) + (nE * (nE - 1) - nEi) * n2 * (n2 - 1) * (n2 - 2) * (n2 - 3) / n / (n - 1) / (n - 2) / (n - 3) + mu2 - mu2^2
  V12 <- (nE * (nE - 1) - nEi) * n1 * n2 * (n1 - 1) * (n2 - 1) / n / (n - 1) / (n - 2) / (n - 3) - mu1 * mu2
  S <- matrix(c(V1, V12, V12, V2), nrow = 2)

  list(
    E = E,
    nE = nE,
    n1 = n1,
    n2 = n2,
    mu0 = mu0,
    mu1 = mu1,
    mu2 = mu2,
    V0 = V0,
    V1 = V1,
    V2 = V2,
    V12 = V12,
    Sinv = solve_partial(S)
  )
}
