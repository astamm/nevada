#' Full, intra and inter subgraph generators
#'
#' This is a collection of functions for extracting full, intra and inter
#' subgraphs of a graph given a list of vertex subsets.
#'
#' @param g An \code{\link[igraph]{igraph}} object.
#' @param vids A list of integer vectors identifying vertex subsets.
#'
#' @return An \code{\link[igraph]{igraph}} object storing a subgraph of type
#'   full, intra or inter.
#'
#' @examples
#' g <- igraph::make_ring(10)
#' g_full  <- subgraph_full (g, list(1:3, 4:5, 8:10))
#' g_intra <- subgraph_intra(g, list(1:3, 4:5, 8:10))
#' g_inter <- subgraph_inter(g, list(1:3, 4:5, 8:10))
#' g_all   <- subgraph_all  (g, list(1:3, 4:5, 8:10))
#'
#' @name subgraphs
NULL

#' @export
#' @rdname subgraphs
subgraph_full <- function(g, vids) {
  if (!igraph::is_named(g)) igraph::V(g)$name <- seq_len(igraph::vcount(g))
  igraph::induced_subgraph(g, unlist(vids)) %>%
    igraph::set_graph_attr("atoms", names(vids))
}

#' @export
#' @rdname subgraphs
subgraph_intra <- function(g, vids) {
  if (!igraph::is_named(g)) igraph::V(g)$name <- seq_len(igraph::vcount(g))
  # Handle graph attributes
  ga_names <- igraph::graph_attr_names(g)
  ga_values <- ga_names %>%
    purrr::map(~ igraph::graph_attr(g, .x)) %>%
    rlang::set_names(ga_names)
  vids %>%
    purrr::map(~ igraph::induced_subgraph(g, .x)) %>%
    purrr::map(delete_graph_attributes) %>%
    purrr::reduce(igraph::disjoint_union) %>%
    add_graph_attributes(ga_values) %>%
    igraph::set_graph_attr("atoms", names(vids))
}

#' @export
#' @rdname subgraphs
subgraph_inter <- function(g, vids) {
  if (!igraph::is_named(g)) igraph::V(g)$name <- seq_len(igraph::vcount(g))
  vids %>%
    purrr::map(~ igraph::induced_subgraph(g, .x)) %>%
    purrr::reduce(igraph::disjoint_union) %>%
    igraph::difference(big = subgraph_full(g, vids)) %>%
    igraph::set_graph_attr("atoms", names(vids))
}

#' @export
#' @rdname subgraphs
subgraph_all <- function(g, vids) {
  if (!igraph::is_named(g)) igraph::V(g)$name <- seq_len(igraph::vcount(g))
  # Handle graph attributes
  ga_names <- igraph::graph_attr_names(g)
  ga_values <- ga_names %>%
    purrr::map(~ igraph::graph_attr(g, .x)) %>%
    rlang::set_names(ga_names)
  g_full <- subgraph_full(g, vids)
  g_intra <- vids %>%
    purrr::map(~ igraph::induced_subgraph(g, .x)) %>%
    purrr::map(delete_graph_attributes) %>%
    purrr::reduce(igraph::disjoint_union) %>%
    add_graph_attributes(ga_values) %>%
    igraph::set_graph_attr("atoms", names(vids))
  g_inter <- igraph::difference(g_full, g_intra)
  list(full = g_full, intra = g_intra, inter = g_inter)
}

delete_graph_attributes <- function(g) {
  attr_names <- igraph::graph_attr_names(g)
  for (name in attr_names) g <- igraph::delete_graph_attr(g, name)
  g
}

add_graph_attributes <- function(g, attribute_list) {
  attr_names <- names(attribute_list)
  for (i in seq_along(attr_names))
    g <- igraph::set_graph_attr(g, attr_names[i], attribute_list[[i]])
  g
}
