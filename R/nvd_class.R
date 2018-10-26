#' Network-Valued Data Constructor
#'
#' This is the constructor for objects of class \code{nvd}.
#'
#' @param model A string specifying the model to be used for sampling networks
#'   (current choices are: \code{"sbm"}, \code{"k_regular"}, \code{"gnp"},
#'   \code{"smallworld"} [default], \code{"pa"}, \code{"poisson"} and
#'   \code{"binomial"}).
#' @param n An integer specifying the sample size (default: \code{0L}).
#' @param pref.matrix A matrix giving the Bernoulli rates for the SBM generator
#'   (see \code{\link[igraph]{sample_sbm}} for details). Default is \code{NULL}.
#'   It is required for \code{model == "sbm"}.
#' @param lambda A numeric value specifying the mean value for the Poisson
#'   generator. Default is \code{NULL}. It is required for \code{model ==
#'   "poisson"}.
#' @param size An integer value specifying the number of trials for the binomial
#'   distribution. Default is \code{NULL}. It is required for \code{model ==
#'   "binomial"}.
#' @param prob A numeric value specifying the probability of success of each
#'   trial for the binomial distribution. Default is \code{NULL}. It is required
#'   for \code{model == "binomial"}.
#'
#' @return A \code{nvd} object which is a list of \code{\link[igraph]{igraph}}
#'   objects.
#' @export
#'
#' @examples
#' nvd(n = 10L)
nvd <- function(model = "smallworld",
                n = 0L,
                pref.matrix = NULL,
                lambda = NULL,
                size = NULL,
                prob = NULL) {

  model <- match.arg(
    model,
    c("sbm", "k_regular", "gnp", "smallworld", "pa", "poisson", "binomial")
  )

  if (model == "sbm" & is.null(pref.matrix))
    stop("The pref.matrix argument should be specified to use the SBM generator.")

  if (model == "poisson" & is.null(lambda))
    stop("The lambda argument should be specified to use the Poisson generator.")

  if (model == "binomial" & (is.null(size) | is.null(prob)))
    stop("The size and prob arguments should be specified to use the Binomial generator.")

  obj <- replicate(n, switch(
    model,
    "sbm" = igraph::sample_sbm(n = 25L, pref.matrix = pref.matrix, block.sizes = c(12L, 1L, 12L)),
    "k_regular" = igraph::sample_k_regular(no.of.nodes = 25L, k = 8L),
    "gnp" = igraph::sample_gnp(n = 25L, p = 1/3),
    "smallworld" = igraph::sample_smallworld(dim = 1L, size = 25L, nei = 4L, p = 0.15),
    "pa" = igraph::sample_pa(n = 25L, power = 2L, m = 4L, directed = FALSE),
    "poisson" = rpois_network(lambda = lambda, n = 25L),
    "binomial" = rbinom_network(size = size, prob = prob, n = 25L)
  ), simplify = FALSE)

  as_nvd(obj)
}

#' Coercion to Network-Valued Data Object
#'
#' This function flags a list of \code{\link[igraph]{igraph}} objects as an
#' \code{\link{nvd}} object as defined in this package.
#'
#' @param obj A list of \code{\link[igraph]{igraph}} objects.
#'
#' @return An \code{\link{nvd}} object.
#' @export
#'
#' @examples
#' as_nvd(nvd("smallworld", 10))
as_nvd <- function(obj) {
  if (!is.list(obj))
    stop("Input should be a list.")

  # check that entries are igraph objects
  input_ok <- TRUE
  for (l in obj) {
    if (!igraph::is_igraph(l)) {
      input_ok <- FALSE
      break()
    }
  }

  if (!input_ok)
    stop("List elements should be igraph objects.")

  class(obj) <- c("nvd", "list")
  obj
}

is_nvd <- function(obj) {
  "nvd" %in% class(obj)
}

#' Fréchet Mean of Network-Valued Data
#'
#' This function computes the sample Fréchet mean from an observed sample of
#' network-valued random variables according to a specified matrix
#' representation. It currently only supports the Euclidean geometry i.e. the
#' sample Fréchet mean is obtained as the argmin of the sum of squared Frobenius
#' distances.
#'
#' @param x An \code{\link{nvd}} object.
#' @param representation A string specifying the graph representation to be used
#'   (choices: adjacency [default], laplacian, modularity, graphon).
#' @param ... Other argument to be parsed to the \code{\link[base]{mean}}
#'   function.
#'
#' @return The mean network in the chosen matrix representation assuming
#'   Euclidean geometry for now.
#' @export
#'
#' @examples
#' d <- nvd(n = 10L)
#' mean(d)
mean.nvd <- function(x, representation = "adjacency", weights = NULL, ...) {
  x <- repr_nvd(x, representation = representation)
  if (is.null(weights)) weights <- rep(1, length(x))
  x <- mean_nvd_impl(x, weights)
  switch(
    representation,
    adjacency = as_adjacency(x),
    laplacian = as_laplacian(x),
    modularity = as_modularity(x),
    graphon = as_graphon(x)
  )
}

#' Fréchet Variance of Network-Valued Data Around a Given Network
#'
#' This function computes the Fréchet variance around a specified network from
#' an observed sample of network-valued random variables according to a
#' specified distance. In most cases, the user is willing to compute the sample
#' variance, in which case the Fréchet variance has to be evaluated w.r.t. the
#' sample Fréchet mean. In this case, it is important that the user indicates
#' the same distance as the one (s)he used to separately compute the sample
#' Fréchet mean. This function can also be used as is as the function to be
#' minimized in order to find the Fréchet mean for a given distance.
#'
#' @param x An \code{\link{nvd}} object listing a sample of networks.
#' @param x0 A network already in matrix representation around which to
#'   calculate variance (usually the Fréchet mean but not necessarily). Note
#'   that the chosen matrix representation is extracted from this parameter.
#' @param distance A string specifying the distance to be used. Possible choices
#'   are: hamming, frobenius [default], spectral or root-euclidean. When the
#'   Fréchet mean is used as \code{x0} parameter, the distance should match the
#'   one used to compute the mean. This is not currently checked.
#'
#' @return A positive scalar value evaluating the amount of variability of the
#'   sample around the specified network.
#' @export
#'
#' @examples
#' d <- nvd(n = 10L)
#' m <- mean(d)
#' var_nvd(x = d, x0 = m, distance = "frobenius")
var_nvd <- function(x, x0, distance = "frobenius") {
  if (!is_nvd(x))
    stop("The input x should be of class nvd.")
  if (!is.matrix(x0))
    stop("The input x0 should be of class matrix.")
  representation <- attributes(x0)$representation
  if (representation == "")
    stop("The input x0 matrix should have an attribute named representation.")
  x <- purrr::map(x, format_input, representation = representation)
  ssd <- switch(
    distance,
    hamming = purrr::reduce(x, function(.v, .x, .x0) {
      d <- dist_hamming_impl(.x, .x0)
      .v + d^2
    }, .x0 = x0, .init = 0),
    frobenius = purrr::reduce(x, function(.v, .x, .x0) {
      d <- dist_frobenius_impl(.x, .x0)
      .v + d^2
    }, .x0 = x0, .init = 0),
    spectral = purrr::reduce(x, function(.v, .x, .x0) {
      d <- dist_spectral_impl(.x, .x0)
      .v + d^2
    }, .x0 = x0, .init = 0),
    "root-euclidean" = purrr::reduce(x, function(.v, .x, .x0) {
      d <- dist_root_euclidean_impl(.x, .x0)
      .v + d^2
    }, .x0 = x0, .init = 0)
  )
  ssd / (length(x) - 1)
}

#' Fréchet Variance of Network-Valued Data from Inter-Point Distances
#'
#' This function computes the Fréchet variance using exclusively inter-point
#' distances. As such, it can accommodate any pair of representation and
#' distance.
#'
#' @param x An \code{\link{nvd}} object listing a sample of networks.
#' @param representation A string specifying the graph representation to be used
#'   (choices: adjacency [default], laplacian, modularity, graphon).
#' @param distance A string specifying the distance to be used. Possible choices
#'   are: hamming, frobenius [default], spectral or root-euclidean.
#'
#' @return A positive scalar value evaluating the variance based on inter-point
#'   distances.
#' @export
#'
#' @examples
#' d <- nvd(n = 10L)
#' var2_nvd(x = d, representation = "graphon", distance = "frobenius")
var2_nvd <- function(x, representation = "adjacency", distance = "frobenius") {
  if (!is_nvd(x))
    stop("The input x should be of class nvd.")
  x <- repr_nvd(x, representation = representation)
  var_nvd_impl(x, distance)
}
