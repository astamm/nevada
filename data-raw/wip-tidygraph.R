library(nevada)
library(ggraph)
library(tidygraph)
library(graphlayouts)

as_tbl_graph.nvd <- function(x, ...) {
  x |>
    map(as_tbl_graph) |>
    imap(~ mutate(activate(.x, "edges"), id = .y)) |>
    # imap(~ mutate(activate(.x, "nodes"), id = .y)) |>
    reduce(graph_join)
}

sample_islands <- function(n, n_islands, size_islands, p_within, m_between) {
  l <- replicate(n, {
    g <- igraph::sample_islands(n_islands, size_islands, p_within, m_between)
    g <- igraph::simplify(g)
    igraph::V(g)$grp <- paste0("I", rep(1:n_islands, each = size_islands))
    g
  }, simplify = FALSE)
  as_nvd(l, matched = TRUE)
}

# Only intra changes ------------------------------------------------------

# Data generation
x <- sample_islands(
  n = 5,
  n_islands = 3,
  size_islands = 10,
  p_within = 0.1,
  m_between = 5
)
y <- sample_islands(
  n = 5,
  n_islands = 3,
  size_islands = 10,
  p_within = 0.9,
  m_between = 5
)

# Viz graphs
z <- bind_nvd(x, y)
g <- as_tbl_graph(z)
bb <- layout_as_backbone(y[[1]], keep = 0.4)
igraph::E(g)$col <- FALSE
igraph::E(g)$col[bb$backbone] <- TRUE
ggraph(g, layout = "manual", x = bb$xy[, 1], y = bb$xy[, 2]) +
  geom_edge_link0(colour = rgb(0, 0, 0, 0.5), width = 0.1) +
  geom_node_point(aes(col = grp)) +
  scale_color_brewer(palette = "Set1") +
  facet_edges(vars(id), nrow = 2) +
  theme_graph() +
  theme(legend.position = "none")

# Test
m <- rep(1:3, each = 10)
res <- test2_local(x, y, m,
                   seed = 1234,
                   alpha = 0.05,
                   B = 1000)
res
alpha <- 0.05
width_min <- 1 - max(res$inter$pvalue, res$intra$pvalue)
width_max <- 1 - min(res$inter$pvalue, res$intra$pvalue)
g2 <- tbl_graph(
  nodes = res$intra |>
    mutate(
      signif = factor(pvalue <= alpha, levels = c(TRUE, FALSE)),
      weight = 1 - pvalue,
      weight = (weight - width_min) / (width_max - width_min)
    ),
  edges = res$inter |>
    rename(from = E1, to = E2) |>
    mutate(
      signif = factor(pvalue <= alpha, levels = c(TRUE, FALSE)),
      weight = 1 - pvalue,
      weight = (weight - width_min) / (width_max - width_min)
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

ggraph(g, layout = 'stress') +
  geom_edge_link(aes(colour = signif), show.legend = FALSE) +
  geom_node_point(aes(fill = signif), shape = 21) +
  geom_node_text(aes(label = E), repel = TRUE) +
  # scale_size(range = c(0, 5)) +
  # scale_edge_width(range = c(0, 1)) +
  theme_graph()
