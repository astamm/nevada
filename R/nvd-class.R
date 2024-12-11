#' Network-Valued Data Constructor
#'
#' This is the constructor for objects of class \code{nvd}.
#'
#' @param sample_size An integer specifying the sample size.
#' @param model A string specifying the model to be used for sampling networks.
#'   All `tidygraph::play_` functions are supported. The model name corresponds
#'   to the name of the function without the `play_` prefix.
#' @param ... Model parameters to be passed to the model function.
#'
#' @return A [`nvd`] object which is a list of [`tidygraph::tbl_graph`] objects.
#' @export
#'
#' @examples
#' params <- list(n_dim = 1L, dim_size = 4L, order = 4L, p_rewire = 0.15)
#' nvd(sample_size = 1L, model = "smallworld", !!!params)
nvd <- function(sample_size, model, ...) {
  model <- rlang::arg_match(
    model,
    c(
      # Graph games from tidygraph
      gsub("play_", "", TIDYGRAPH_PLAY_FUNCTIONS()),
      # User-defined graph games
      "binomial", "exponential", "poisson"
    )
  )

  call_data <- build_play_call(model, ...)

  obj <- replicate(sample_size, {
    rlang::exec(call_data[[1]], !!!call_data[[2]])
  }, simplify = FALSE)

  as_nvd(obj)
}

TIDYGRAPH_PLAY_FUNCTIONS <- function() {
  nms <- names(getNamespace("tidygraph"))
  nms[grepl("^play_*", nms)]
}

#' @importFrom rlang `:=`
build_play_call <- function(model, ...) {
  fn_name <- paste0("play_", model)
  pkg <- if (fn_name %in% TIDYGRAPH_PLAY_FUNCTIONS()) {
    "tidygraph"
  } else {
    "nevada"
  }
  fn_call <-utils::getFromNamespace(fn_name, pkg)
  fn_args <- names(formals(fn_call))
  fn_mandatory_args <- fn_args[purrr::map_lgl(formals(fn_call), \(el) {
    class(el) == "name"
  })]
  fn_optional_args <- fn_args[!fn_args %in% fn_mandatory_args]
  user_args <- rlang::list2(...)
  ok <- all(names(user_args) %in% fn_args)
  if (!ok) {
    cli::cli_abort(c(
      "x" = "The arguments you provided do not match the signature of the {.fn {pkg}::{fn_name}} function.",
      "i" = "You must provide the {.arg {fn_mandatory_args}} argument{?s} and optionally the {.arg {fn_optional_args}} argument{?s}.",
      "i" = "Refer to {.fn {pkg}::{fn_name}} for more information on the model parameters."
    ))
  }
  ok <- all(fn_mandatory_args %in% names(user_args))
  if (!ok) {
    cli::cli_abort(c(
      "x" = "You did not provide all the mandatory arguments to the {.fn {pkg}::{fn_name}} function.",
      "i" = "You must provide the {.arg {fn_mandatory_args}} argument{?s} and optionally the {.arg {fn_optional_args}} argument{?s}.",
      "i" = "See {.fn {pkg}::{fn_name}} for more information on the model parameters."
    ))
  }
  tb <- table(names(user_args))
  dups <- names(tb[tb > 1])
  if (length(dups) > 0) {
    cli::cli_abort(c(
      "x" = "You provided duplicate arguments to the {.fn {pkg}::{fn_name}} function.",
      "i" = "The following arguments are duplicated: {.arg {dups}}."
    ))
  }

  if (length(fn_optional_args) > 0) {
    for (arg in fn_optional_args) {
      if (!arg %in% names(user_args)) {
        user_args <- c(user_args, rlang::list2(
          !!arg := formals(fn_call)[[arg]]
        ))
      }
    }
  }

  cli::cli_alert_info("Calling the {.fn {pkg}::{fn_name}} function with the following arguments:")
  lid <- cli::cli_ul()
  for (i in 1:length(user_args)) {
    argval <- if (is.null(user_args[[i]])) "NULL" else user_args[[i]]
    cli::cli_li("{names(user_args)[i]}: {argval}")
  }
  cli::cli_end(lid)

  list(fn_call, user_args)
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
#' params <- list(n_dim = 1L, dim_size = 4L, order = 4L, p_rewire = 0.15)
#' out <- nvd(sample_size = 1L, model = "smallworld", !!!params)
#' as_nvd(out)
as_nvd <- function(obj) {
  if (!is.list(obj))
    cli::cli_abort("Input should be a list.")

  N <- length(obj)
  if (N == 0)
    cli::cli_abort("Input list should not be empty.")

  # Input must be a list of graph.
  # Each entry can be:
  # - either an igraph object
  # - or a tidygraph::tbl_graph object
  # - or a list itself with entries `nodes`, `edges`, `directed` and `node_key`.

  if (igraph::is_igraph(obj[[1]])) {
    input_ok <- TRUE
    for (i in 1:N) {
      if (!igraph::is_igraph(obj[[i]])) {
        input_ok <- FALSE
        break()
      }
      obj[[i]] <- tidygraph::as_tbl_graph(obj[[i]])
    }

    if (!input_ok)
      cli::cli_abort("All elements in the input list should be {.cls igraph} objects.")
  } else if (tidygraph::is.tbl_graph(obj[[1]])) {
    input_ok <- TRUE
    for (i in 2:N) {
      if (!tidygraph::is.tbl_graph(l[[i]])) {
        input_ok <- FALSE
        break()
      }
    }

    if (!input_ok)
      cli::cli_abort("All elements in the input list should be {.cls tbl_graph} objects.")
  } else if (is.list(obj[[1]])) {
    input_ok <- TRUE
    for (i in 1:N) {
      if (!all(c("nodes", "edges", "directed", "node_key") %in% names(l[[i]]))) {
        input_ok <- FALSE
        break()
      }
      l[[i]] <- tidygraph::tbl_graph(
        nodes = l[[i]]$nodes,
        edges = l[[i]]$edges,
        directed = l[[i]]$directed,
        node_key = l[[i]]$node_key
      )
    }

    if (!input_ok)
      cli::cli_abort("All elements in the input list should be lists with entries `nodes`, `edges`, `directed` and `node_key`.")
  } else {
    cli::cli_abort("All elements in the input list should be either {.cls igraph}, {.cls tbl_graph} or lists with entries `nodes`, `edges`, `directed` and `node_key`.")
  }

  class(obj) <- c("nvd", "list")
  obj
}

is_nvd <- function(obj) {
  "nvd" %in% class(obj)
}

#' Two-Sample Stochastic Block Model Generator
#'
#' This function generates two samples of networks according to the stochastic
#' block model (SBM). This is essentially a wrapper around
#' \code{\link[igraph]{sample_sbm}} which allows to sample a single network from
#' the SBM.
#'
#' @param n Integer scalar giving the sample size.
#' @param nv Integer scalar giving the number of vertices of the generated
#'   networks, common to all networks in both samples.
#' @param p1 The matrix giving the Bernoulli rates for the 1st sample. This is a
#'   KxK matrix, where K is the number of groups. The probability of creating an
#'   edge between vertices from groups i and j is given by element (i,j). For
#'   undirected graphs, this matrix must be symmetric.
#' @param b1 Numeric vector giving the number of vertices in each group for the
#'   first sample. The sum of the vector must match the number of vertices.
#' @param p2 The matrix giving the Bernoulli rates for the 2nd sample (default:
#'   same as 1st sample). This is a KxK matrix, where K is the number of groups.
#'   The probability of creating an edge between vertices from groups i and j is
#'   given by element (i,j). For undirected graphs, this matrix must be
#'   symmetric.
#' @param b2 Numeric vector giving the number of vertices in each group for the
#'   second sample (default: same as 1st sample). The sum of the vector must
#'   match the number of vertices.
#'
#' @return A length-2 list containing the two samples stored as
#'   \code{\link{nvd}} objects.
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
#' sim <- sample2_sbm(n, 68, p1, c(17, 17, 17, 17), p2)
sample2_sbm <- function(n, nv, p1, b1, p2 = p1, b2 = b1) {
  sim <- n |>
    purrr::rerun(
      x = igraph::sample_sbm(nv, p1, b1),
      y = igraph::sample_sbm(nv, p2, b2)
    ) |>
    purrr::transpose() |>
    purrr::map(as_nvd)
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
#' @param weights A numeric vector specifying weights for each observation
#'   (default: equally weighted).
#' @param representation A string specifying the graph representation to be
#'   used. Choices are adjacency, laplacian, modularity, graphon. Default is
#'   adjacency.
#' @param ... Other argument to be parsed to the \code{\link[base]{mean}}
#'   function.
#'
#' @return The mean network in the chosen matrix representation assuming
#'   Euclidean geometry for now.
#' @export
#'
#' @examples
#' params <- list(n = 24L, p = 1/3)
#' x <- nvd(sample_size = 1L, model = "gnp", !!!params)
#' mean(x)
mean.nvd <- function(x, weights = rep(1, length(x)), representation = "adjacency", ...) {
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
#' @param weights A numeric vector specifying weights for each observation
#'   (default: equally weighted).
#' @param distance A string specifying the distance to be used. Possible choices
#'   are: hamming, frobenius, spectral or root-euclidean. Default is frobenius.
#'   When the Fréchet mean is used as \code{x0} parameter, the distance should
#'   match the one used to compute the mean. This is not currently checked.
#'
#' @return A positive scalar value evaluating the amount of variability of the
#'   sample around the specified network.
#' @export
#'
#' @examples
#' params <- list(n = 24L, p = 1/3)
#' x <- nvd(sample_size = 1L, model = "gnp", !!!params)
#' m <- mean(x)
#' var_nvd(x = x, x0 = m, distance = "frobenius")
var_nvd <- function(x, x0, weights = rep(1, length(x)), distance = "frobenius") {
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
    hamming = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_hamming_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0),
    frobenius = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_frobenius_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0),
    spectral = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_spectral_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0),
    "root-euclidean" = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_root_euclidean_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0)
  )
  ssd / sum(weights)
}

#' Fréchet Variance of Network-Valued Data from Inter-Point Distances
#'
#' This function computes the Fréchet variance using exclusively inter-point
#' distances. As such, it can accommodate any pair of representation and
#' distance.
#'
#' @param x An \code{\link{nvd}} object listing a sample of networks.
#' @param representation A string specifying the graph representation to be
#'   used. Choices are adjacency, laplacian, modularity, graphon. Default is
#'   adjacency.
#' @param distance A string specifying the distance to be used. Possible choices
#'   are: hamming, frobenius, spectral or root-euclidean. Default is frobenius.
#'
#' @return A positive scalar value evaluating the variance based on inter-point
#'   distances.
#' @export
#'
#' @examples
#' params <- list(n = 24L, p = 1/3)
#' x <- nvd(sample_size = 1L, model = "gnp", !!!params)
#' var2_nvd(x = x, representation = "graphon", distance = "frobenius")
var2_nvd <- function(x, representation = "adjacency", distance = "frobenius") {
  if (!is_nvd(x))
    stop("The input x should be of class nvd.")
  x <- repr_nvd(x, representation = representation)
  var_nvd_impl(x, distance)
}
