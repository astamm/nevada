format_input <- function(x, representation = "adjacency") {
  if (igraph::is_igraph(x))
    x <- switch(
      representation,
      adjacency = get_adjacency(x),
      laplacian = get_laplacian(x),
      modularity = get_modularity(x)
    )
  else {
    if ("representation" %in% names(attributes(x))) {
      if (attributes(x)$representation != "adjacency" ||
          (attributes(x)$representation == "adjacency" && representation == "adjacency"))
        return(x)
    } else {
      if (any(x < 0))
        stop("All entries of an adjacency matrix should be non-negative.")
      if (any(x != t(x)))
        stop("The input adjacency matrix should be symmetric.")
      if (any(diag(x) != 0))
        stop("The input adjacency matrix should have a value of 0 on the diagonal.")

      x <- switch(
        representation,
        adjacency = get_adjacency(igraph::graph_from_adjacency_matrix(x, mode = "undirected")),
        laplacian = get_laplacian(igraph::graph_from_adjacency_matrix(x, mode = "undirected")),
        modularity = get_modularity(igraph::graph_from_adjacency_matrix(x, mode = "undirected"))
      )
    }
  }

  if (!("representation" %in% names(attributes(x))))
    stop("The input matrix representation should have a representation
         attribute that specifies which network representation has
         been used.")

  if (representation != attributes(x)$representation)
    stop("The input network is not in the desired representation.")

  if (nrow(x) != ncol(x))
    stop("A matrix representation of a network should be a square matrix.")

  x
}

compatible_networks <- function(x, y) {
  compatible <- TRUE

  if (igraph::is_igraph(x)) {
    nx <- igraph::vcount(x)
    ny <- igraph::vcount(y)
  } else {
    nx <- nrow(x)
    ny <- nrow(y)
  }

  if (nx != ny) {
    message("The input networks do not have the same number of vertices.")
    compatible <- FALSE
  }

  compatible
}

reachable_significance <- function(nx, ny, B, alpha = 0.05, verbose = FALSE) {
  n_comb <- choose(nx + ny, nx)
  p_min <- 1 / min(B, n_comb)

  if (verbose) {
    writeLines(paste(" - P-value resolution:", p_min))
    if (B >= n_comb) { # Case of exact test
      writeLines(" - Computing exact p-value.")
      writeLines(paste(" - P-value will never drop below", p_min))
    } else { # Case of approximate test
      writeLines(paste(" - Computing approximate p-value using", B,
                       "random permutations."))
      writeLines(paste(" - P-value will not drop below", 1 / n_comb,
                       "on average over repeated Monte-Carlo estimates."))
    }
  }

  alpha * n_comb >= 1
}

