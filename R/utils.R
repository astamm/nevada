format_input <- function(x, representation = "adjacency") {
  if (igraph::is_igraph(x))
    x <- switch(
      representation,
      adjacency = repr_adjacency(x),
      laplacian = repr_laplacian(x),
      modularity = repr_modularity(x),
      transitivity = repr_transitivity(x)
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
        adjacency = repr_adjacency(igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE)),
        laplacian = repr_laplacian(igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE)),
        modularity = repr_modularity(igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE)),
        transitivity = repr_transitivity(igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE))
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

phipson_smyth_pvalue <- function(b, B, M) {
  if (M <= 10000) {
    pt <- (1:M) / M
    return(mean(pbinom(q = b, size = B, prob = pt)))
  }

  corr <- integrate(pbinom, 0, 0.5 / M, q = b, size = B)$value
  (b + 1) / (B + 1) - corr
}

get_permuted_statistic <- function(i, indices1, d, statistic) {
  if (i == 0)
    indices <- seq_len(nrow(indices1))
  else
    indices <- indices1[, i]

  switch(
    statistic,
    plot = stat_lot(d, indices, pooled = TRUE),
    lot = stat_lot(d, indices, pooled = FALSE),
    sot = stat_sot(d, indices),
    biswas = stat_biswas(d, indices),
    energy = stat_energy(d, indices),
    student = stat_student_euclidean(d, indices),
    welch = stat_welch_euclidean(d, indices),
    original = stat_edge_count(d, indices, type = "original"),
    generalized = stat_edge_count(d, indices, type = "generalized"),
    weighted = stat_edge_count(d, indices, type = "weighted")
  )
}

capitalize <- function(x) {
  gsub("(?<=\\b)([a-z])", "\\U\\1", tolower(x), perl = TRUE)
}

edge_count_global_variables <- function(d, n1, k = 1L) {
  k <- min(k, ceiling(nrow(d) / 4))
  g <- kmst(d, k)
  E <- igraph::as_edgelist(g)
  n <- nrow(d)
  n2 <- n - n1
  Ebynode <- vector("list", n)
  for (i in 1:nrow(E)) {
    Ebynode[[E[i, 1]]] <- c(Ebynode[[E[i, 1]]], E[i, 2])
    Ebynode[[E[i, 2]]] <- c(Ebynode[[E[i, 2]]], E[i, 1])
  }
  nE <- nrow(E)
  nodedeg <- sapply(Ebynode, length)
  nEi <- sum(nodedeg * (nodedeg - 1))
  mu0 <- nE * 2 * n1 * n2 / n / (n - 1)
  mu1 <- nE * n1 * (n1 - 1) / n / (n - 1)
  mu2 <- nE * n2 * (n2 - 1) / n / (n - 1)
  V0 <- nEi * n1 * n2 / n / (n - 1) + (nE * (nE - 1) - nEi) * 4 * n1 * n2 * (n1 - 1) * (n2 - 1) / n / (n - 1) / (n - 2) / (n - 3) + mu0 - mu0^2
  V1 <- nEi * n1 * (n1 - 1) * (n1 - 2) / n / (n - 1) / (n - 2) + (nE * (nE - 1) - nEi) * n1 * (n1 - 1) * (n1 - 2) * (n1 - 3) / n / (n - 1) / (n - 2) / (n - 3) + mu1 - mu1^2
  V2 <- nEi * n2 * (n2 - 1) * (n2 - 2) / n / (n - 1) / (n - 2) + (nE * (nE - 1) - nEi) * n2 * (n2 - 1) * (n2 - 2) * (n2 - 3) / n / (n - 1) / (n - 2) / (n - 3) + mu2 - mu2^2
  V12 <- (nE * (nE - 1) - nEi) * n1 * n2 * (n1 - 1) * (n2 - 1) / n / (n - 1) / (n - 2) / (n - 3) - mu1 * mu2
  S <- matrix(c(V1, V12, V12, V2), nrow = 2)

  list(
    E = E,
    nE = nE,
    n1 = n1,
    n2 = n2,
    mu0 = mu0,
    mu1 = mu1,
    mu2 = mu2,
    V0 = V0,
    V1 = V1,
    V2 = V2,
    V12 = V12,
    Sinv = solve_partial(S)
  )
}

kmst <- function(d, k = 1L) {
  for (i in 1:k) {
    g <- igraph::graph_from_adjacency_matrix(d, mode = "undirected", weighted = TRUE)
    m <- igraph::mst(g, algorithm = "prim")
    if (i == 1)
      res <- m
    else
      res <- igraph::union(res, m)
    e <- igraph::as_edgelist(m)
    for (j in 1:nrow(e)) {
      d[e[j, 1], e[j, 2]] <- .Machine$double.xmax
      d[e[j, 2], e[j, 1]] <- .Machine$double.xmax
    }
  }
  res
}

solve_partial <- function(M) {
  eig <- eigen(M, symmetric = TRUE)
  Minv <- matrix(0, nrow = nrow(M), ncol = ncol(M))
  for (i in seq_along(eig$values)) {
    l <- eig$values[i]
    if (l < .Machine$double.eps)
      next()
    v <- eig$vectors[, i]
    Minv <- Minv + (v %*% t(v)) / l
  }
  Minv
}

xi <- function(th) {
  2 + th^2 - pi / 8 * ( (2 + th^2) * besselI(th^2 / 4, 0, expon.scaled = TRUE) + th^2 * besselI(th^2 / 4, 1, expon.scaled = TRUE) )^2
}

fixed_point <- function(th, r) {
  sqrt(xi(th) * (1 + r^2) - 2)
}

rice_params <- function(x) {
  mu <- mean(x)
  sig <- sd(x)
  r <- mu / sig
  lb <- sqrt(2 / xi(0) - 1)
  if (r <= lb)
    th <- 0
  else {
    th <- r - lb
    th_new <- 0
    iter <- 0
    while (abs(th_new - th) > 1.0e-6) {
      iter <- iter + 1
      th_new <- th
      th <- fixed_point(th, r)
    }
  }
  sigma <- sig / sqrt(xi(th))
  nu <- sqrt( mu^2 + (xi(th) - 2) * sig^2 )
  list(nu = nu, sigma = sigma)
}
