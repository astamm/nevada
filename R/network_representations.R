#' Network Representation Functions
#'
#' This is a collection of functions that convert a graph stored as an
#' \code{\link[igraph]{igraph}} object into a desired matrix representation
#' among adjacency matrix, graph laplacian or modularity matrix.
#'
#' @param network An \code{\link[igraph]{igraph}} object.
#' @param validate A boolean specifying whether the function should check the
#'   class of its input (default: \code{TRUE}).
#'
#' @return A numeric square matrix giving the desired network representation
#'   recorded in the object's class.
#'
#' @examples
#' X <- igraph::watts.strogatz.game(1, 25, 3, 0.05)
#' adjX <- get_adjacency(X)
#' lapX <- get_laplacian(X)
#' modX <- get_modularity(X)
#' @name get-representation
NULL

#' @rdname get-representation
#' @export
get_adjacency <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  repr <- igraph::as_adjacency_matrix(network, type = "both", sparse = FALSE)
  as_adjacency(repr)
}

#' @rdname get-representation
#' @export
get_laplacian <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  repr <- igraph::laplacian_matrix(network, sparse = FALSE)
  as_laplacian(repr)
}

#' @rdname get-representation
#' @export
get_modularity <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  repr <- igraph::modularity_matrix(network, rep(1, igraph::vcount(network)))
  as_modularity(repr)
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

is_adjacency <- function(x) {
  "adjacency" == attributes(x)$representation
}

is_laplacian <- function(x) {
  "laplacian" == attributes(x)$representation
}

is_modularity <- function(x) {
  "modularity" == attributes(x)$representation
}

as_transitivity <- function(x) {
  x <- as.matrix(x)
  l <- attributes(x)
  l$representation <- "transitivity"
  attributes(x) <- l
  x
}

get_transitivity <- function(network, validate = TRUE) {
  if (validate) {
    if (!igraph::is_igraph(network))
      stop("Input network should be of class igraph.")
  }
  repr <- igraph::transitivity(network, type = "undirected")
  as_transitivity(repr)
}
