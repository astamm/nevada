library(nevada)
library(ggraph)
library(tidygraph)
library(graphlayouts)

as_tbl_graph.nvd <- function(x, directed = TRUE, memberships = NULL, ...) {
  nm <- names(x)
  if (is.null(nm)) {
    if (is.null(memberships))
      nm <- seq_len(length(x))
    else {
      nm <- numeric(length(x))
      groups <- unique(memberships)
      for (.g in groups)
        nm[memberships == .g] <- seq_len(sum(memberships == .g))
    }
  }
  print(nm)
  x <- x |>
    purrr::map(igraph::as_edgelist, names = FALSE) |>
    purrr::map2(nm, cbind)
  if (is.null(memberships)) {
    x <- purrr::map(x, `colnames<-`, c("from", "to", "id"))
  } else {
    x <- x |>
      purrr::map2(memberships, cbind) |>
      purrr::map(`colnames<-`, c("from", "to", "id", "membership"))
  }
  x |>
    purrr::map(tibble::as_tibble) |>
    purrr::reduce(rbind) |>
    as_tbl_graph(directed = directed)
}

x <- nvd("pa", 6, model_params = list(power = 1, m = 1, directed = FALSE))
xx <- purrr::map(x, igraph::as_edgelist, names = FALSE) |>
  purrr::imap(~ cbind(.x, .y)) |>
  purrr::map(`colnames<-`, c("from", "to", "id")) |>
  purrr::map(tibble::as_tibble) |>
  bind_rows() |>
  as_tbl_graph(directed = FALSE) |>
  mutate(Popularity = centrality_degree(mode = 'in'))

ggraph(xx, layout = 'kk') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  geom_node_point(aes(size = Popularity)) +
  facet_edges(~id, nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

xxx <- as_tbl_graph(x)
ggraph(xxx, layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  # geom_node_point(aes(size = Popularity)) +
  facet_edges(vars(id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

library(nevada)
library(tidygraph)
library(ggraph)
n <- 10
p1 <- matrix(
  data = c(0.1, 0.4, 0.1, 0.4,
           0.4, 0.4, 0.1, 0.4,
           0.1, 0.1, 0.4, 0.4,
           0.4, 0.4, 0.4, 0.4),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)
p2 <- matrix(
  data = c(0.1, 0.4, 0.4, 0.4,
           0.4, 0.4, 0.4, 0.4,
           0.4, 0.4, 0.1, 0.1,
           0.4, 0.4, 0.1, 0.4),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)
sim <- sample2_sbm(n, 68, p1, c(17, 17, 17, 17), p2, seed = 1234)
m <- as.integer(c(rep(1, 17), rep(2, 17), rep(3, 17), rep(4, 17)))

ggraph(as_tbl_graph(sim[[1]]), layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  facet_edges(vars(id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

ggraph(as_tbl_graph(sim[[2]]), layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  facet_edges(vars(id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

res <- test2_local(sim$x, sim$y, m,
                   seed = 1234,
                   # alpha = 0.05,
                   B = 100)
alpha <- 0.05
edge_width_min <- 1 - max(res$inter$pvalue)
edge_width_max <- 1 - min(res$inter$pvalue)
node_width_min <- 1 - max(res$intra$pvalue)
node_width_max <- 1 - min(res$intra$pvalue)
g <- tbl_graph(
  nodes = res$intra |>
    mutate(
      signif = pvalue <= alpha,
      weight = 1 - pvalue,
      weight = (weight - node_width_min) / (node_width_max - node_width_min)
    ),
  edges = res$inter |>
    rename(from = E1, to = E2) |>
    mutate(
      signif = pvalue <= alpha,
      weight = 1 - pvalue,
      weight = (weight - edge_width_min) / (edge_width_max - edge_width_min)
    ),
  directed = FALSE,
  node_key = "E"
)

ggraph(g, layout = 'stress') +
  geom_edge_link(aes(width = weight, edge_colour = signif), show.legend = FALSE) +
  geom_node_point(aes(size = weight, fill = signif), shape = 21) +
  geom_node_text(aes(label = E), repel = TRUE) +
  scale_size(range = c(0, 5)) +
  scale_edge_width(range = c(0, 1)) +
  theme_graph()




n <- 5
p1 <- matrix(
  data = c(0.5, 0.5, 0.1, 0.1,
           0.5, 0.5, 0.1, 0.1,
           0.1, 0.1, 0.5, 0.5,
           0.1, 0.1, 0.5, 0.5),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)
p2 <- matrix(
  data = c(0.1, 0.1, 0.1, 0.1,
           0.1, 0.5, 0.5, 0.1,
           0.1, 0.5, 0.5, 0.1,
           0.1, 0.1, 0.1, 0.1),
  nrow = 4,
  ncol = 4,
  byrow = TRUE
)
sim <- sample2_sbm(n, 20, p1, rep(5, 4), p2, seed = 1234)
m <- as.integer(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)))

ggg <- as_tbl_graph(
  as_nvd(c(sim[[1]], sim[[2]])),
  memberships = c(rep(1, 5), rep(2, 5))
)
ggraph(ggg, layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  facet_edges(vars(membership, id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

sim <- sample2_sbm(30, 20, p1, rep(5, 4), p2, seed = 1234)
res <- test2_local(sim$x, sim$y, m,
                   seed = 1234,
                   # alpha = 0.05,
                   B = 100)
res


ggraph(as_tbl_graph(sim[[1]]), layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  facet_edges(vars(id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

ggraph(as_tbl_graph(sim[[2]]), layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  facet_edges(vars(id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
