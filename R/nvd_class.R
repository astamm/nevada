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

#' Mean of network-valued data
#'
#' @param x An \code{\link{nvd}} object.
#' @param representation A string specifying the graph representation to be used (choices: adjacency [default], laplacian, modularity).
#' @param ... Other argument to be parsed to the \code{\link[base]{mean}} function.
#'
#' @return The mean network in the chosen matrix representation assuming Euclidean geometry for now.
#' @export
#'
#' @examples
#' d <- nvd(n = 10L)
#' mean(d)
mean.nvd <- function(x, representation = "adjacency", ...) {
  x <- switch (representation,
    adjacency = purrr::map(x, repr_adjacency),
    laplacian = purrr::map(x, repr_laplacian),
    modularity = purrr::map(x, repr_modularity)
  )
  x <- mean_nvd_impl(x)
  switch (representation,
    adjacency = as_adjacency(x),
    laplacian = as_laplacian(x),
    modularity = as_modularity(x)
  )
}
