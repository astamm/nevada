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
#' @param d Either a matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the
#'   distances between all the elements of the two samples put together or, for
#'   the \code{stat_dom_frobenius} statistic, the concatenation of the lists of
#'   matrix representations of networks in samples 1 and 2.
#' @param indices A vector of dimension \eqn{n1} containing the indices of the
#'   elements of the first sample.
#' @param standardize If \code{TRUE} (default), standardize the test statistic
#'   using the pooled sample variance estimator.
#'
#' @return A scalar giving the value of the desired test statistic.
#'
#' @examples
#' x <- nvd("smallworld", 10)
#' y <- nvd("pa", 10)
#' d <- dist_nvd(x, y, representation = "laplacian", distance = "frobenius")
#' stat_mod(d, 1:10)
#' stat_dom(d, 1:10, FALSE)
#' stat_dom(d, 1:10, TRUE)
#' dm <- repr_nvd(x, y, representation = "laplacian")
#' stat_dom_frobenius(dm, 1:10, FALSE)
#' stat_dom_frobenius(dm, 1:10, TRUE)
#' @name statistics
NULL

#' @rdname statistics
#' @export
stat_mod <- function(d, indices) {
  mean(d[indices, -indices])
}

#' @rdname statistics
#' @export
stat_dom <- function(d, indices, standardize = TRUE) {
  n <- nrow(d)
  n1 <- length(indices)
  n2 <- n - n1

  ssd1_vec <- numeric(n1)
  for (i in seq_along(indices))
    ssd1_vec[i] <- sum(d[indices[i], indices]^2)
  km1 <- indices[which.min(ssd1_vec)]

  ssd2_vec <- numeric(n2)
  indices <- seq_len(n)[-indices]
  for (i in seq_along(indices))
    ssd2_vec[i] <- sum(d[indices[i], indices]^2)
  km2 <- indices[which.min(ssd2_vec)]

  stat <- d[km1, km2]

  if (!standardize)
    return(stat)

  ssd1 <- min(ssd1_vec)
  ssd2 <- min(ssd2_vec)
  pooled_variance <- (ssd1 + ssd2) / (n - 2)
  stat / sqrt(pooled_variance)
}

#' @rdname statistics
#' @export
stat_dom_frobenius <- function(d, indices, standardize = TRUE) {
  x <- d[indices]
  y <- d[-indices]
  stat_dom_frobenius_impl(x, y, standardize)
}
