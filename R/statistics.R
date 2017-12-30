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
#' @param standardize If \code{TRUE} (default), standardize the DoM statistic
#'   using the pooled sample variance estimator.
#' @param type A string specifying the version of the edge count test statistic
#'   to be used. Choices are \code{"original"} [default], \code{"generalized"}
#'   or \code{"weighted"}.
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
stat_lot <- function(d, indices, pooled = TRUE) {
  xy <- d[indices, -indices]^2
  xx <- d[indices, indices]
  xx <- xx[upper.tri(xx)]^2
  yy <- d[-indices, -indices]
  yy <- yy[upper.tri(yy)]^2
  varx <- mean(xx) / 2
  vary <- mean(yy) / 2
  mse <- mean(xy)
  n1 <- length(indices)
  n2 <- nrow(d) - n1
  if (pooled)
    var <- mean(c(xx, yy)) / 2 * (1 / n1 + 1 / n2)
  else
    var <- varx / n1 + vary / n2
  (mse - varx - vary) / var
}

#' @rdname statistics
#' @export
stat_sot <- function(d, indices) {
  xx <- d[indices, indices]
  xx <- xx[upper.tri(xx)]^2
  yy <- d[-indices, -indices]
  yy <- yy[upper.tri(yy)]^2
  varx <- mean(xx) / 2
  vary <- mean(yy) / 2
  max(varx, vary) / min(varx, vary)
}

#' @rdname statistics
#' @export
stat_biswas <- function(d, indices) {
  xy <- d[indices, -indices]
  xx <- d[indices, indices]
  xx <- xx[upper.tri(xx)]
  yy <- d[-indices, -indices]
  yy <- yy[upper.tri(yy)]
  muFF <- mean(xx)
  muGG <- mean(yy)
  muFG <- mean(xy)
  (muFF - muFG)^2 + (muGG - muFG)^2
}

#' @rdname statistics
#' @export
stat_energy <- function(d, indices, alpha = 1) {
  xy <- d[indices, -indices]^alpha
  xx <- d[indices, indices]^alpha
  yy <- d[-indices, -indices]^alpha
  varx <- mean(xx) / 2
  vary <- mean(yy) / 2
  mse <- mean(xy)
  mse - varx - vary
}

#' @rdname statistics
#' @export
stat_mod <- function(d, indices) {
  xy <- d[indices, -indices]
  mean(xy)
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
stat_student_euclidean <- function(d, indices) {
  x <- d[indices]
  y <- d[-indices]
  stat_t_euclidean_impl(x, y, pooled = TRUE)
}

#' @rdname statistics
#' @export
stat_welch_euclidean <- function(d, indices) {
  x <- d[indices]
  y <- d[-indices]
  stat_t_euclidean_impl(x, y, pooled = FALSE)
}

#' @rdname statistics
#' @export
stat_edge_count <- function(d, indices, type = "original") {
  E <- d$E
  nE <- d$nE
  n1 <- d$n1
  n2 <- d$n2
  mu0 <- d$mu0
  mu1 <- d$mu1
  mu2 <- d$mu2
  V0 <- d$V0
  V1 <- d$V1
  V2 <- d$V2
  V12 <- d$V12
  Sinv <- d$Sinv

  temp <- stat_edge_count_impl(E, indices)
  R1 <- temp[1]
  R2 <- temp[2]

  switch(type,
    original = -(nE - R1 - R2 - mu0) / sqrt(V0),
    generalized = Sinv[1,1] * (R1 - mu1)^2 + Sinv[2, 2] * (R2 - mu2)^2 + 2 * Sinv[1, 2] * (R1 - mu1) * (R2 - mu2),
    weighted = (n2 * (R1 - mu1) + n1 * (R2 - mu2)) / sqrt(n2^2 * V1 + n1^2 * V2 + 2 * n2 * n1 * V12)
  )
}
