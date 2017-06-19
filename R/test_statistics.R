#' Test Statistics for Network Populations
#'
#' This is a collection of functions that provide various test statistics for
#' populations of networks. The \emph{average} test statistic is given by the
#' mean of all the distances between the \eqn{n1} elements of the fist sample
#' and those (\eqn{n2}) of the second one. The \emph{Fréchet} test statistic is
#' given by the distance between the Frechet mean of the first sample (with
#' cardinality \eqn{n1}) and that of the second one (with cardinality \eqn{n2}).
#' The distance used to compute the value of the test statistic is the same used
#' to find the Fréchet means.
#'
#' @param d A matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the
#'   distances between all the elements of the two samples put together.
#' @param indices A vector of dimension \eqn{n1} containing the indices of the
#'   elements of the first sample.
#'
#' @return A scalar giving the value of the desired test statistic.
#'
#' @examples
#' @name get-statistic
NULL

#' @rdname get-statistic
#' @export
get_average_statistic <- function(d, indices) {
  mean(d[indices, seq_len(nrow(d))[-indices]])
}

#' @rdname get-statistic
#' @export
get_frechet_statistic <- function(d, indices) {
  N <- nrow(d)
  M <- length(indices)
  frechet1 <- rep(NA, M)
  frechet2 <- rep(NA, (N-M))

  t <- 1
  for (i in indices) {
    frechet1[t] <- sum(d[i, indices]^2)
    t <- t+1
  }
  frechetmean1 <- indices[which.min(frechet1)]

  r <- 1
  for (i in seq_len(N)[-indices]) {
    frechet2[t] <- sum(d[i, indices]^2)
    r <- r+1
  }
  frechetmean2 <- seq_len(N)[-indices][which.min(frechet2)]

  statistic <- d[frechetmean1, frechetmean2]
}
