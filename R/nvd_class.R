#' Network-Valued Data Constructor
#'
#' This is the constructor for objects of class \code{nvd}.
#'
#' @param model A string specifying the model to be used for sampling networks
#'   (current choices are: \code{"sbm"}, \code{"k_regular"}, \code{"gnp"},
#'   \code{"smallworld"} [default], \code{"pa"}, \code{"poisson"}).
#' @param n An integer specifying the sample size (default: \code{0L}).
#' @param pref.matrix A matrix giving the Bernoulli rates for the SBM
#'   generator (see \code{\link[igraph]{sample_sbm}} for details).
#' @param lambda A numeric value specifying the mean value for the Poisson
#'   generator.
#'
#' @return A \code{nvd} object which is a list of \code{\link[igraph]{igraph}}
#'   objects.
#' @export
#'
#' @examples
#' nvd(n = 10L)
nvd <- function(model = "smallworld", n = 0L, pref.matrix = NULL, lambda = NULL) {
  model <- match.arg(model, c("sbm", "k_regular", "gnp", "smallworld", "pa", "poisson"))

  if (model == "sbm" & is.null(pref.matrix))
    stop("The pref.matrix argument should be specify to use the SBM generator.")

  if (model == "poisson" & is.null(lambda))
    stop("The lambda argument should be specify to use the Poisson generator.")

  obj <- replicate(n, switch(
    model,
    "sbm" = igraph::sample_sbm(n = 25L, pref.matrix = pref.matrix, block.sizes = c(12L, 1L, 12L)),
    "k_regular" = igraph::sample_k_regular(no.of.nodes = 25L, k = 8L),
    "gnp" = igraph::sample_gnp(n = 25L, p = 1/3),
    "smallworld" = igraph::sample_smallworld(dim = 1L, size = 25L, nei = 4L, p = 0.15),
    "pa" = igraph::sample_pa(n = 25L, power = 2L, m = 4L, directed = FALSE),
    "poisson" = rpois_network(lambda = lambda, n = 25L)
  ), simplify = FALSE)

  as_nvd(obj)
}

as.nvd <- function(obj) {
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

as_nvd <- as.nvd

is.nvd <- function(obj) {
  "nvd" %in% class(obj)
}

is_nvd <- is.nvd
