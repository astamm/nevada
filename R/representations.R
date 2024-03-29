#' Network Representation Functions
#'
#' This is a collection of functions that convert a graph stored as an
#' \code{\link[igraph]{igraph}} object into a desired matrix representation
#' among adjacency matrix, graph laplacian, modularity matrix or graphon (edge
#' probability matrix).
#'
#' @param network An \code{\link[igraph]{igraph}} object.
#' @param validate A boolean specifying whether the function should check the
#'   class of its input (default: \code{TRUE}).
#'
#' @return A numeric square matrix giving the desired network representation
#'   recorded in the object's class.
#'
#' @examples
#' g <- igraph::sample_smallworld(1, 25, 3, 0.05)
#' repr_adjacency(g)
#' repr_laplacian(g)
#' repr_modularity(g)
#' repr_graphon(g)
#' @name representations
NULL

#' @rdname representations
#' @export
repr_adjacency <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  if ("weight" %in% igraph::edge_attr_names(network))
    repr <- repr_adjacency_impl(
      igraph::gorder(network),
      igraph::as_edgelist(network, names = FALSE),
      igraph::E(network)$weight
    )
  else
    repr <- igraph::as_adjacency_matrix(network, type = "both", sparse = FALSE)
  as_adjacency(repr)
}

#' @rdname representations
#' @export
repr_laplacian <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  repr <- igraph::laplacian_matrix(network, sparse = FALSE)
  as_laplacian(repr)
}

#' @rdname representations
#' @export
repr_modularity <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }

  repr <- igraph::modularity_matrix(network)
  repr[is.nan(repr)] <- 0
  as_modularity(repr)
}

#' @rdname representations
#' @export
repr_graphon <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  A <- repr_adjacency(network, FALSE)
  A <- matrix(as.numeric(A > 0), igraph::gorder(network))
  repr <- est_nbdsmooth(A)
  as_graphon(repr)
}

as_adjacency <- function(x) {
  l <- attributes(x)
  l$representation <- "adjacency"
  attributes(x) <- l
  x
}

as_laplacian <- function(x) {
  l <- attributes(x)
  l$representation <- "laplacian"
  attributes(x) <- l
  x
}

as_modularity <- function(x) {
  l <- attributes(x)
  l$representation <- "modularity"
  attributes(x) <- l
  x
}

as_graphon <- function(x) {
  l <- attributes(x)
  l$representation <- "graphon"
  attributes(x) <- l
  x
}

is_adjacency <- function(x) {
  "adjacency" == attributes(x)$representation
}

is_laplacian <- function(x) {
  "laplacian" == attributes(x)$representation
}

is_modularity <- function(x) {
  "modularity" == attributes(x)$representation
}

is_graphon <- function(x) {
  "graphon" == attributes(x)$representation
}

as_transitivity <- function(x) {
  x <- as.matrix(x)
  l <- attributes(x)
  l$representation <- "transitivity"
  attributes(x) <- l
  x
}

repr_transitivity <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  repr <- igraph::transitivity(network, type = "undirected")
  as_transitivity(repr)
}

#' Network-Valued to Matrix-Valued Data
#'
#' @param x An \code{\link{nvd}} object.
#' @param y An \code{\link{nvd}} object. If \code{NULL} (default), it is not
#'   taken into account.
#' @param representation A string specifying the requested matrix
#'   representation. Choices are: \code{"adjacency"}, \code{"laplacian"} or
#'   \code{"modularity"}. Default is \code{"adjacency"}.
#'
#' @return A list of matrices.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' xm <- repr_nvd(x)
repr_nvd <- function(x, y = NULL, representation = "adjacency") {
  x <- lapply(x, format_input, representation)
  if (!is.null(y)) {
    y <- lapply(y, format_input, representation)
    x <- c(x, y)
  }
  x
}
