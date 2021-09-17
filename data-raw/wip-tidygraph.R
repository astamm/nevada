library(nevada)
library(ggraph)
library(tidygraph)
library(graphlayouts)

as_tbl_graph.nvd <- function(x, directed = TRUE, ...) {
  nm <- names(x)
  if (is.null(nm)) nm <- seq_len(length(x))
  x |>
    purrr::map(igraph::as_edgelist, names = FALSE) |>
    purrr::map2(nm, cbind) |>
    purrr::map(`colnames<-`, c("from", "to", "id")) |>
    purrr::map(tibble::as_tibble) |>
    purrr::reduce(rbind) |>
    as_tbl_graph(directed = directed)
}

x <- nvd("pa", 6)
xx <- purrr::map(x, igraph::as_edgelist, names = FALSE) |>
  purrr::imap(~ cbind(.x, .y)) |>
  purrr::map(`colnames<-`, c("from", "to", "id")) |>
  purrr::map(tibble::as_tibble) |>
  bind_rows() |>
  as_tbl_graph(directed = FALSE) |>
  mutate(Popularity = centrality_degree(mode = 'in'))

ggraph(xx, layout = 'kk') +
  geom_edge_fan(aes(alpha = stat(index)), show.legend = FALSE) +
  geom_node_point(aes(size = Popularity)) +
  facet_edges(~id, nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

xxx <- as_tbl_graph(x)
ggraph(xxx, layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  # geom_node_point(aes(size = Popularity)) +
  facet_edges(vars(id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
