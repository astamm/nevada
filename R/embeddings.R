#' Graph Sample Embeddings
#'
#' A collection of functions to embed a sample of graphs into suitable spaces
#' for further statistical analysis.
#'
#' @param obj An object of class [`nvd`] containing the sample of graphs.
#'
#' @return An object of class [`nvd`] containing the sample of graphs in graph space.
#'
#' @export
#' @examples
#' x <- nvd(sample_size = 5L, model = "gnp", n = 24L, p = 1/3)
#' x <- push_to_graph_space(x)
push_to_graph_space <- function(obj) {
  if (!is_nvd(obj))
    cli::cli_abort("Input should be an object of class {.cls nvd}.")

  # Here we handle graphs of different orders by adding as many vertices as
  # necessary to match the order of the largest graph (a.k.a pushing graphs into
  # graph space structure).
  num_vertices <- purrr::map_int(obj, igraph::gorder)
  cli::cli_alert_info("Graphs have the following number of vertices: {num_vertices}")
  max_num_vertices <- max(num_vertices)
  n_chars <- nchar(max_num_vertices)
  obj <- purrr::map2(obj, num_vertices, \(g, n) {
    g <- dplyr::mutate(
      g,
      name = paste0("O", formatC(1:dplyr::n(), width = n_chars, flag = "0"))
    )
    n_extra <- max_num_vertices - n
    if (n_extra == 0)
      return(g)
    extra_nodes <- data.frame(
      name = paste0("V", formatC(1:n_extra, width = n_chars, flag = "0"))
    )
    tidygraph::bind_nodes(g, extra_nodes)
  })

  obj
}
