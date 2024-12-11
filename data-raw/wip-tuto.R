library(nevada)
# library(igraph)
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

withr::with_seed(1234, {
  sim <- sample2_sbm(n, 20, p1, rep(5, 4), p2, seed = 1234)
})
m <- as.integer(c(rep(1, 5), rep(2, 5), rep(3, 5), rep(4, 5)))

g <- as_tbl_graph(
  as_nvd(c(sim[[1]], sim[[2]])),
  memberships = c(rep(1, 5), rep(2, 5))
)
ggraph(g, layout = 'stress') +
  geom_edge_fan(aes(alpha = after_stat(index)), show.legend = FALSE) +
  facet_edges(vars(membership, id), nrow = 2) +
  theme_graph(foreground = 'steelblue', fg_text_colour = 'white')

withr::with_seed(1234, {
  sim_all <- sample2_sbm(50, 20, p1, rep(5, 4), p2, seed = 1234)
})
res <- test2_local(sim_all$x, sim_all$y, m, seed = 1234, B = 100)
res
