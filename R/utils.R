format_inputs <- function(x, y, representation = "adjacency") {
  if (igraph::is_igraph(x))
    x <- switch(representation,
                adjacency = get_adjacency(x),
                laplacian = get_laplacian(x),
                modularity = get_modularity(x)
    )

  if (igraph::is_igraph(y))
    y <- switch(representation,
                adjacency = get_adjacency(y),
                laplacian = get_laplacian(y),
                modularity = get_modularity(y)
    )

  if (representation != attributes(x)$representation || representation != attributes(y)$representation)
    stop("The input networks are not in the desired representation.")

  if (nrow(x) != nrow(y))
    stop("The input networks do not have the same number of vertices.")

  list(x = x, y = y)
}
