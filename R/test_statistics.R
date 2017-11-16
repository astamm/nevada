#' Test Statistics for Network Populations
#'
#' This is a collection of functions that provide various test statistics for
#' populations of networks. The \emph{MoD} (mean of the distances) test
#' statistic is given by the mean of all the distances between the \eqn{n1}
#' elements of the fist sample and those (\eqn{n2}) of the second one. The
#' \emph{DoM} (distance of the means) test statistic is given by the distance
#' between the Frechet mean of the first sample (with cardinality \eqn{n1}) and
#' that of the second one (with cardinality \eqn{n2}). The \emph{sDoM}
#' (standardized distance of the means) is the \emph{DoM} statistic standardized
#' by the pooled sample standard deviation. The distance used to compute the
#' value of the test statistic is the same used to find the Fréchet means.
#'
#' @param d A matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the distances
#'   between all the elements of the two samples put together.
#' @param indices A vector of dimension \eqn{n1} containing the indices of the
#'   elements of the first sample.
#'
#' @return A scalar giving the value of the desired test statistic.
#'
#' @examples
#' n <- 25L
#' x <- list()
#' y <- list()
#' for (i in 1:10) {
#'   X <- igraph::watts.strogatz.game(1, n, 3, 0.05)
#'   Y <- igraph::barabasi.game(n, m = 3, power = 2, directed = FALSE)
#'   adjX <- get_adjacency(X)
#'   adjY <- get_adjacency(Y)
#'   x[[i]] <- adjX
#'   y[[i]] <- adjY
#' }
#' d <- get_distance_matrix(x, y, "adjacency", "spectral")
#' stat_mod <- get_mod_statistic(d, 1:10)
#' stat_dom <- get_dom_statistic(d, 1:10)
#' stat_sdom <- get_sdom_statistic(d, 1:10)
#' @name get-statistic
NULL

#' @rdname get-statistic
#' @export
get_mod_statistic <- function(d, indices) {
  mean(d[indices, seq_len(nrow(d))[-indices]])
}

#' @rdname get-statistic
#' @export
get_dom_statistic <- function(d, indices) {
  N <- nrow(d)
  M <- length(indices)
  frechet1 <- rep(NA, M)
  frechet2 <- rep(NA, (N-M))

  for (i in seq_along(indices))
    frechet1[i] <- sum(d[indices[i], indices]^2)
  frechetmean1 <- indices[which.min(frechet1)]

  indices <- seq_len(N)[-indices]
  for (i in seq_along(indices))
    frechet2[i] <- sum(d[indices[i], indices]^2)
  frechetmean2 <- indices[which.min(frechet2)]

  d[frechetmean1, frechetmean2]
}

#' @rdname get-statistic
#' @export
get_sdom_statistic <- function(d, indices) {
  N <- nrow(d)
  M <- length(indices)
  frechet1 <- rep(NA, M)
  frechet2 <- rep(NA, (N-M))

  for (i in seq_along(indices))
    frechet1[i] <- sum(d[indices[i], indices]^2)
  frechetmean1 <- indices[which.min(frechet1)]
  frechetvar1 <- min(frechet1)

  indices <- seq_len(N)[-indices]
  for (i in seq_along(indices))
    frechet2[i] <- sum(d[indices[i], indices]^2)
  frechetmean2 <- indices[which.min(frechet2)]
  frechetvar2 <- min(frechet2)

  pooled_variance <- ((M - 1) * frechetvar1 + (N - M - 1) * frechetvar2) / (N - 2)

  d[frechetmean1, frechetmean2] / sqrt(pooled_variance)
}
