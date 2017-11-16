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

rperm <- function(m, size=2) { # Obtain m unique permutations of 1:size

  # Function to obtain a new permutation.
  newperm <- function() {
    count <- 0                # Protects against infinite loops
    repeat {
      # Generate a permutation and check against previous ones.
      p <- sample(1:size)
      hash.p <- paste(p, collapse="")
      if (is.null(cache[[hash.p]])) break

      # Prepare to try again.
      count <- count+1
      if (count > 1000) {   # 1000 is arbitrary; adjust to taste
        p <- NA           # NA indicates a new permutation wasn't found
        hash.p <- ""
        break
      }
    }
    cache[[hash.p]] <<- TRUE  # Update the list of permutations found
    p                         # Return this (new) permutation
  }

  # Obtain m unique permutations.
  cache <- list()
  replicate(m, newperm())
} # Returns a `size` by `m` matrix; each column is a permutation of 1:size.

rperm2 <- function(m, size=2) { # Obtain m unique permutations of 1:size
  max.failures <- 10

  # Function to index into the upper-level cache.
  prefix <- function(p, k) {    # p is a permutation, k is the prefix size
    sum((p[1:k] - 1) * (size ^ ((1:k)-1))) + 1
  } # Returns a value from 1 through size^k

  # Function to obtain a new permutation.
  newperm <- function() {
    # References cache, k.head, and failures in parent context.
    # Modifies cache and failures.

    count <- 0                # Protects against infinite loops
    repeat {
      # Generate a permutation and check against previous ones.
      p <- sample(1:size)
      k <- prefix(p, k.head)
      ip <- cache[[k]]
      hash.p <- paste(tail(p,-k.head), collapse="")
      if (is.null(ip[[hash.p]])) break

      # Prepare to try again.
      n.failures <<- n.failures + 1
      count <- count+1
      if (count > max.failures) {
        p <- NA           # NA indicates a new permutation wasn't found
        hash.p <- ""
        break
      }
    }
    if (count <= max.failures) {
      ip[[hash.p]] <- TRUE      # Update the list of permutations found
      cache[[k]] <<- ip
    }
    p                         # Return this (new) permutation
  }

  # Initialize the cache.
  k.head <- min(size-1, max(1, floor(log(m / log(m)) / log(size))))
  cache <- as.list(1:(size^k.head))
  for (i in 1:(size^k.head)) cache[[i]] <- list()

  # Count failures (for benchmarking and error checking).
  n.failures <- 0

  # Obtain (up to) m unique permutations.
  s <- replicate(m, newperm())
  s[is.na(s)] <- NULL
  list(failures=n.failures, sample=matrix(unlist(s), ncol=size))
} # Returns an m by size matrix; each row is a permutation of 1:size.
