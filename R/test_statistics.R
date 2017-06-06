
#' Title Compute the value of the test statistic 'average'
#'
#' This test statistic is given by the mean of all the distances between the
#' \eqn{n1} elements of the fist sample and those (\eqn{n2}) of the second
#' one.
#'
#' @param d A matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the
#'   distances between all the elements of the two samples put together.
#' @param indices A vector of dimension \eqn{n1} containing the indices of the
#'   elements of the first sample.
#'
#' @return \code{get_average_statistic} returns a scalar.
#' @export
#'
#' @examples
get_average_statistic <- function(d, indices) {
  N <- dim(d)[1]
  statistic <- mean(d[indices, (1:N)[-indices]])
}

#' Title Compute the value of the test statistic 'frechet'
#'
#' This test statistic is given by the distance between the Frechet mean of the
#' first sample (with cardinality \eqn{n1}) and that of the second one (with
#' cardinality \eqn{n2}). The distance used to compute the value of the test
#' statistic is the same used to find the Frechet means.
#'
#' @param d A matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing all the
#'   distances between all the elements of the two samples put together.
#' @param indices A vector of dimension \eqn{n1} containing the indices of the
#'   elements of the first sample.
#'
#' @return \code{get_frechet_statistic} returns a scalar.
#' @export
#'
#' @examples
get_frechet_statistic <- function(d, indices) {
  N <- dim(d)[1]
  M <- length(indices)
  frechet1 <- rep(NA, M)
  frechet2 <- rep(NA, (N-M))

  t <- 1
  for (i in indices) {
    frechet1[t] <- sum(d[i,indices]^2)
    t <- t+1
  }
  frechetmean1 <- indices[which.min(frechet1)]

  r <- 1
  for (i in ((1:N)[-indices])) {
    frechet2[t] <- sum(d[i,indices]^2)
    r <- r+1
  }
  frechetmean2 <- (1:N)[-indices][which.min(frechet2)]


  statistic <- d[frechetmean1, frechetmean2]
}

