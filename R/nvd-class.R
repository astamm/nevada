#' Network-Valued Data Constructor
#'
#' This is the constructor for objects of class \code{nvd}.
#'
#' @param model A string specifying the model to be used for sampling networks
#'   (current choices are: `"sbm"`, `"k_regular"`, `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"`). Defaults to `"smallworld"`.
#' @param n An integer specifying the sample size. Defaults to `1L`.
#' @param num_vertices An integer specifying the order of the graphs to be
#'   generated (i.e. the number of nodes). Defaults to `25L`.
#' @param model_params A named list setting the parameters of the model you are
#'   considering. Defaults to `NULL`.
#' @param seed An integer specifying the random generator seed. Defaults to
#'   `1234`.
#' @param rand_num_vertices A boolean specifying whether the number of vertices of the networks in the sample should be random. In particular N_1,...,N_sample_size iid Poisson(num_vertices). It is compatible with `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"` models. Defaults to `FALSE`.
#' @param  rand_num_vertices A string specifying whether the number of vertices of the networks in the first sample should be random. Choices are `"poisson"`, `"uniform"`. In particular, if `"poisson"` N_1,...,N_sample_size iid Poisson(num_vertices), if `"uniform"` N_1,...,N_sample_size iid Uniform(floor(num_vertices*(1-rate)), ceiling(num_vertices*(1+rate))). It is compatible with `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"` models. Defaults to `NULL`.
#' @param rate A rate. See `rand_num_vertices` for further information. Defaults to `0.25`.
#'
#' @return A \code{nvd} object which is a list of \code{\link[igraph]{igraph}}
#'   objects.
#' @export
#'
#' @examples
#' smallworld_params <- list(dim = 1L, nei = 4L, p = 0.15)
#' nvd(model_params = smallworld_params)
nvd <- function(model = "smallworld",
                n = 1L,
                num_vertices = 25L,
                model_params = NULL,
                seed = 1234,
                rand_num_vertices = NULL,
                rate = 0.25) {
  if (!is.null(seed))
    withr::local_seed(seed)

  model <- match.arg(
    model,
    c("sbm", "k_regular", "gnp", "smallworld", "pa", "poisson", "binomial")
  )

  if (is.null(rand_num_vertices)) {
    if (!rlang::is_named2(model_params))
      cli::cli_abort("The {.code model_params} list should be named after the parameters of the model you are considering but is not.")

    if (model == "poisson") {
      if (!all(c("lambda") %in% names(model_params)))
        cli::cli_abort("The {.code model_params} list should contain the field {.field lambda} to use the Poisson generator.")
      return(rlang::eval_tidy(rlang::quo(
        rpois_network(n = n, num_vertices = num_vertices, !!!model_params)
      )))
    }

    if (model == "binomial") {
      if (!all(c("size", "prob") %in% names(model_params)))
        cli::cli_abort("The {.code model_params} list should contain the fields {.field size} and {.field prob} to use the binomial generator.")
      return(rlang::eval_tidy(rlang::quo(
        rbinom_network(n = n, num_vertices = num_vertices, !!!model_params)
      )))
    }

    if (model == "sbm" && !all(c("pref.matrix", "block.sizes") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the fields {.field pref.matrix} and {.field block.sizes} to use the SBM generator.")

    if (model == "k_regular" && !all(c("k") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the field {.field k} to use the k-regular model generator.")

    if (model == "gnp" && !all(c("p") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the field {.field p} to use the GNP model generator.")

    if (model == "smallworld" && !all(c("dim", "nei", "p") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the fields {.field dim}, {.field nei} and {.field p} to use the Watts-Strogatz small-world model generator.")

    if (model == "pa" && !all(c("power", "m", "directed") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the fields {.field power}, {.field m} and {.field directed} to use the Barabasi-Albert scale-free model generator.")

    obj <- replicate(n, {
      rlang::eval_tidy(rlang::quo(switch(
        model,
        "sbm" = igraph::sample_sbm(
          n = num_vertices,
          !!!model_params
        ),
        "k_regular" = igraph::sample_k_regular(
          no.of.nodes = num_vertices,
          !!!model_params
        ),
        "gnp" = igraph::sample_gnp(
          n = num_vertices,
          !!!model_params
        ),
        "smallworld" = igraph::sample_smallworld(
          size = num_vertices,
          !!!model_params
        ),
        "pa" = igraph::sample_pa(
          n = num_vertices,
          !!!model_params
        )
      )))
    }, simplify = FALSE)
  } else {

    if (!rlang::is_named2(model_params))
      cli::cli_abort("The {.code model_params} list should be named after the parameters of the model you are considering but is not.")

    if (model == "sbm" || model == "k_regular") {
      cli::cli_abort("It is not possible to generate a sample with random number of vertices from this model.")
    }

    if (model == "poisson") {
      if (!all(c("lambda") %in% names(model_params)))
        cli::cli_abort("The {.code model_params} list should contain the field {.field lambda} to use the Poisson generator.")
    }

    if (model == "binomial") {
      if (!all(c("size", "prob") %in% names(model_params)))
        cli::cli_abort("The {.code model_params} list should contain the fields {.field size} and {.field prob} to use the binomial generator.")
    }

    if (model == "gnp" && !all(c("p") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the field {.field p} to use the GNP model generator.")

    if (model == "smallworld" && !all(c("dim", "nei", "p") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the fields {.field dim}, {.field nei} and {.field p} to use the Watts-Strogatz small-world model generator.")

    if (model == "pa" && !all(c("power", "m", "directed") %in% names(model_params)))
      cli::cli_abort("The {.code model_params} list should contain the fields {.field power}, {.field m} and {.field directed} to use the Barabasi-Albert scale-free model generator.")
    if (rand_num_vertices == "poisson") {
      obj <- replicate(n, {
        rlang::eval_tidy(rlang::quo(switch(
          model,
          "binomial" = rbinom_network2(
            num_vertices = stats::rpois(n = 1, lambda = num_vertices),
            !!!model_params),
          "poisson" = rpois_network2(
            num_vertices = stats::rpois(n = 1, lambda = num_vertices),
            !!!model_params),
          "gnp" = igraph::sample_gnp(
            n = stats::rpois(n = 1, lambda = num_vertices),
            !!!model_params
          ),
          "smallworld" = igraph::sample_smallworld(
            size = stats::rpois(n = 1, lambda = num_vertices),
            !!!model_params
          ),
          "pa" = igraph::sample_pa(
            n = stats::rpois(n = 1, lambda = num_vertices),
            !!!model_params
          )
        )))
      }, simplify = FALSE)
    } else if (rand_num_vertices == "uniform") {
      obj <- replicate(n, {
        rlang::eval_tidy(rlang::quo(switch(
          model,
          "binomial" = rbinom_network2(
            num_vertices = sample(floor(num_vertices*(1-rate)):ceiling(num_vertices*(1+rate)), size = 1),
            !!!model_params),
          "poisson" = rpois_network2(
            num_vertices = sample(floor(num_vertices*(1-rate)):ceiling(num_vertices*(1+rate)), size = 1),
            !!!model_params),
          "gnp" = igraph::sample_gnp(
            n = sample(floor(num_vertices*(1-rate)):ceiling(num_vertices*(1+rate)), size = 1),
            !!!model_params
          ),
          "smallworld" = igraph::sample_smallworld(
            size = sample(floor(num_vertices*(1-rate)):ceiling(num_vertices*(1+rate)), size = 1),
            !!!model_params
          ),
          "pa" = igraph::sample_pa(
            n = sample(floor(num_vertices*(1-rate)):ceiling(num_vertices*(1+rate)), size = 1),
            !!!model_params
          )
        )))
      }, simplify = FALSE)
    } else {
      cli::cli_abort("rand_num_vertices is not correctly initialized.")
    }

  }

  as_nvd(obj)
}

#' Coercion to Network-Valued Data Object
#'
#' This function flags a list of \code{\link[igraph]{igraph}} objects as an
#' \code{\link{nvd}} object as defined in this package.
#'
#' @param obj A list of \code{\link[igraph]{igraph}} objects.
#'
#' @return An \code{\link{nvd}} object.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' as_nvd(nvd(model = "gnp", n = 10L, model_params = gnp_params))
as_nvd <- function(obj) {
  if (!is.list(obj))
    cli::cli_abort("Input should be a list.")

  # check that entries are igraph objects
  input_ok <- TRUE
  for (l in obj) {
    if (!igraph::is_igraph(l)) {
      input_ok <- FALSE
      break()
    }
  }

  if (!input_ok)
    cli::cli_abort("List elements should be igraph objects.")

  class(obj) <- c("nvd", "list")
  obj
}

is_nvd <- function(obj) {
  "nvd" %in% class(obj)
}

#' Two-Sample Stochastic Block Model Generator
#'
#' This function generates two samples of networks according to the stochastic
#' block model (SBM). This is essentially a wrapper around
#' \code{\link[igraph]{sample_sbm}} which allows to sample a single network from
#' the SBM.
#'
#' @param n Integer scalar giving the sample size.
#' @param nv Integer scalar giving the number of vertices of the generated
#'   networks, common to all networks in both samples.
#' @param p1 The matrix giving the Bernoulli rates for the 1st sample. This is a
#'   KxK matrix, where K is the number of groups. The probability of creating an
#'   edge between vertices from groups i and j is given by element (i,j). For
#'   undirected graphs, this matrix must be symmetric.
#' @param b1 Numeric vector giving the number of vertices in each group for the
#'   first sample. The sum of the vector must match the number of vertices.
#' @param p2 The matrix giving the Bernoulli rates for the 2nd sample (default:
#'   same as 1st sample). This is a KxK matrix, where K is the number of groups.
#'   The probability of creating an edge between vertices from groups i and j is
#'   given by element (i,j). For undirected graphs, this matrix must be
#'   symmetric.
#' @param b2 Numeric vector giving the number of vertices in each group for the
#'   second sample (default: same as 1st sample). The sum of the vector must
#'   match the number of vertices.
#' @param seed The seed for the random number generator (default: \code{NULL}).
#'
#' @return A length-2 list containing the two samples stored as
#'   \code{\link{nvd}} objects.
#' @export
#'
#' @examples
#' n <- 10
#' p1 <- matrix(
#'   data = c(0.1, 0.4, 0.1, 0.4,
#'            0.4, 0.4, 0.1, 0.4,
#'            0.1, 0.1, 0.4, 0.4,
#'            0.4, 0.4, 0.4, 0.4),
#'   nrow = 4,
#'   ncol = 4,
#'   byrow = TRUE
#' )
#' p2 <- matrix(
#'   data = c(0.1, 0.4, 0.4, 0.4,
#'            0.4, 0.4, 0.4, 0.4,
#'            0.4, 0.4, 0.1, 0.1,
#'            0.4, 0.4, 0.1, 0.4),
#'   nrow = 4,
#'   ncol = 4,
#'   byrow = TRUE
#' )
#' sim <- sample2_sbm(n, 68, p1, c(17, 17, 17, 17), p2, seed = 1234)
sample2_sbm <- function(n, nv, p1, b1, p2 = p1, b2 = b1, seed = NULL) {
  withr::local_seed(seed)
  sim <- n %>%
    purrr::rerun(
      x = igraph::sample_sbm(nv, p1, b1),
      y = igraph::sample_sbm(nv, p2, b2)
    ) %>%
    purrr::transpose() %>%
    purrr::map(as_nvd)
}

#' Fréchet Mean of Network-Valued Data
#'
#' This function computes the sample Fréchet mean from an observed sample of
#' network-valued random variables according to a specified matrix
#' representation. It currently only supports the Euclidean geometry i.e. the
#' sample Fréchet mean is obtained as the argmin of the sum of squared Frobenius
#' distances.
#'
#' @param x An \code{\link{nvd}} object.
#' @param weights A numeric vector specifying weights for each observation
#'   (default: equally weighted).
#' @param representation A string specifying the graph representation to be
#'   used. Choices are adjacency, laplacian, modularity, graphon. Default is
#'   adjacency.
#' @param aac A boolean specifying whether the function should perform the Align All and Compute algorithm (AAC). If it is \code{TRUE}, the alignment is performed via FAQ (Frank-Wolfe algorithm), with barycenter as initialization of the permutation matrix and 20 iterations. Defaults to `FALSE`.
#' @param tol Tolerance for AAC. Default to `0.001`.
#' @param max_iteration Maximum number of iteration for AAC. Default to `200`.
#' @param seed An integer specifying the random generator seed for AAC. Defaults to
#'   `1234`.
#' @param ... Other argument to be parsed to the \code{\link[base]{mean}}
#'   function.
#'
#' @return The mean network in the chosen matrix representation assuming
#'   Euclidean geometry for now.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' mean(x)
mean.nvd <- function(x, weights = rep(1, length(x)), representation = "adjacency", aac = FALSE, tol = 0.001, max_iteration = 200, seed = 1234, ...) {
  if(!aac){
    x <- repr_nvd(x, representation = representation)
    if (is.null(weights)) weights <- rep(1, length(x))
    x <- mean_nvd_impl(x, weights)
    switch(
      representation,
      adjacency = as_adjacency(x),
      laplacian = as_laplacian(x),
      modularity = as_modularity(x),
      graphon = as_graphon(x)
    )
  } else {
    if (!is.null(seed))
      withr::local_seed(seed)
    x <- repr_nvd(x, representation = representation)
    if (is.null(weights)) weights <- rep(1, length(x))

    first_id <- sample(1:length(x), 1)
    m <- x[[first_id]]
    randP <- randRotation::randpermut(nrow(m))
    m <- randP%*%m%*%t(randP)

    s <- tol + 1
    k <- 1
    while (s > tol & k <= max_iteration) {
      cat("Start iteration ", k, "\n")
      x_aligned <- align_networks(m, x)
      fm <- mean_nvd_impl(x_aligned, weights)
      s <- dist_frobenius_impl(m, fm)
      m <- fm
      k <- k+1
    }

    if(k >= max_iteration)
      print("Max number of iteration reached")

    switch(
      representation,
      adjacency = as_adjacency(m),
      laplacian = as_laplacian(m),
      modularity = as_modularity(m),
      graphon = as_graphon(m)
    )
  }
}

#' Fréchet Variance of Network-Valued Data Around a Given Network
#'
#' This function computes the Fréchet variance around a specified network from
#' an observed sample of network-valued random variables according to a
#' specified distance. In most cases, the user is willing to compute the sample
#' variance, in which case the Fréchet variance has to be evaluated w.r.t. the
#' sample Fréchet mean. In this case, it is important that the user indicates
#' the same distance as the one (s)he used to separately compute the sample
#' Fréchet mean. This function can also be used as is as the function to be
#' minimized in order to find the Fréchet mean for a given distance.
#'
#' @param x An \code{\link{nvd}} object listing a sample of networks.
#' @param x0 A network already in matrix representation around which to
#'   calculate variance (usually the Fréchet mean but not necessarily). Note
#'   that the chosen matrix representation is extracted from this parameter.
#' @param weights A numeric vector specifying weights for each observation
#'   (default: equally weighted).
#' @param distance A string specifying the distance to be used. Possible choices
#'   are: hamming, frobenius, spectral or root-euclidean. Default is frobenius.
#'   When the Fréchet mean is used as \code{x0} parameter, the distance should
#'   match the one used to compute the mean. This is not currently checked.
#'
#' @return A positive scalar value evaluating the amount of variability of the
#'   sample around the specified network.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' m <- mean(x)
#' var_nvd(x = x, x0 = m, distance = "frobenius")
var_nvd <- function(x, x0, weights = rep(1, length(x)), distance = "frobenius") {
  if (!is_nvd(x))
    stop("The input x should be of class nvd.")
  if (!is.matrix(x0))
    stop("The input x0 should be of class matrix.")
  representation <- attributes(x0)$representation
  if (representation == "")
    stop("The input x0 matrix should have an attribute named representation.")
  x <- purrr::map(x, format_input, representation = representation)
  ssd <- switch(
    distance,
    hamming = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_hamming_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0),
    frobenius = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_frobenius_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0),
    spectral = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_spectral_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0),
    "root-euclidean" = purrr::reduce2(x, weights, function(.v, .x, .y, .x0) {
      d <- dist_root_euclidean_impl(.x, .x0)
      .v + .y * d^2
    }, .x0 = x0, .init = 0)
  )
  ssd / sum(weights)
}

#' Fréchet Variance of Network-Valued Data from Inter-Point Distances
#'
#' This function computes the Fréchet variance using exclusively inter-point
#' distances. As such, it can accommodate any pair of representation and
#' distance.
#'
#' @param x An \code{\link{nvd}} object listing a sample of networks.
#' @param representation A string specifying the graph representation to be
#'   used. Choices are adjacency, laplacian, modularity, graphon. Default is
#'   adjacency.
#' @param distance A string specifying the distance to be used. Possible choices
#'   are: hamming, frobenius, spectral or root-euclidean. Default is frobenius.
#'
#' @return A positive scalar value evaluating the variance based on inter-point
#'   distances.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' var2_nvd(x = x, representation = "graphon", distance = "frobenius")
var2_nvd <- function(x, representation = "adjacency", distance = "frobenius") {
  if (!is_nvd(x))
    stop("The input x should be of class nvd.")
  x <- repr_nvd(x, representation = representation)
  var_nvd_impl(x, distance)
}

#' From adjacency matrices to \code{nvd} object
#'
#' This function creates an \code{nvd} object from a list of adjacency matrices.
#'
#' @param adjmat A list of adjacency matrices.
#' @param directed A boolean specifying whether the function should consider the input matrices as directed or undirected. Default to `FALSE`.
#' @param weighted A boolean specifying whether the function should consider the input matrices as weighted or unweighted Default to `FALSE`.
#'
#' @return An \code{nvd} object.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' x_adj <- repr_nvd(x, representation = "adjacency")
#' nvd_from_adjacency_matrix(x_adj, directed = FALSE, weighted = FALSE)
nvd_from_adjacency_matrix <- function(adjmat, directed = FALSE, weighted = FALSE){
  sample_size <- length(adjmat)
  y <- lapply(1:sample_size, function(x) NULL)

  if (!directed) {
    mode <- "undirected"
  } else{
    mode <- "directed"
  }
  if (!weighted) {
    weight <- NULL
  } else{
    weight <- TRUE
  }

  for (i in 1:sample_size) {
    attr(adjmat[[i]], "representation") <- "adjacency"
    y[[i]] <- igraph::graph_from_adjacency_matrix(adjmat[[i]], weighted = weight, mode = mode)
  }
  as_nvd(y)
}

#' Summary of number of vertices in a sample of networks
#'
#' This function summarizes information related to the number of vertices in an \code{nvd} object.
#'
#' @param x An \code{nvd} object.
#'
#' @return A \code{\link[base]{list}} with 5 components: a vector with the number of vertices of each network, the maximum, the minimum, the mean and the variance of the number of vertices of the networks.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' nodes <- nvd_num_vertices(x)
#' nodes$num_vertices
#' nodes$max
nvd_num_vertices <- function(x){
  sample_size <- length(x)
  vec <- rep(-1, sample_size)
  x_adj <- repr_nvd(x, representation = "adjacency")
  for (i in 1:sample_size) {
    vec[i] <- nrow(x_adj[[i]])
  }
  list("num_vertices" = vec, "max" = max(vec), "min" = min(vec), "mean" = mean(vec), "var" = var(vec))
}

#' Adding null nodes to networks
#'
#' This function adds null nodes to networks stored in an \code{nvd} object.
#'
#' @param x An \code{nvd} object.
#' @param num_vertices An integer specifying the desired number of vertices of the networks. If necessary, null nodes are added in each network to reach this number.  If `NULL`, the maximum number of vertices in the sample is considered. Default to `NULL`.
#' @param directed A boolean specifying whether the function should consider the input networks as directed or undirected. Default to `FALSE`.
#' @param weighted A boolean specifying whether the function should consider the input networks as weighted or unweighted Default to `FALSE`.
#'
#' @return An \code{nvd} object.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' x <- nvd(model = "gnp", n = 10L, num_vertices = 10L, model_params = gnp_params)
#' add_null_nodes(x, num_vertices = 12L, directed = FALSE, weighted = FALSE)
add_null_nodes <- function(x, num_vertices = NULL, directed = FALSE, weighted = FALSE){
  max_num_vertices <- nvd_num_vertices(x)$max

  if (is.null(num_vertices)) {
    num_vertices <- max_num_vertices
  } else if (num_vertices < max_num_vertices){
    stop("The input number of vertices is smaller than the maximum number of vertices in the sample.")
  }

  sample_size <- length(x)

  zero_mat <- matrix(0, num_vertices, num_vertices)
  attr(zero_mat, "representation") <- "adjacency"
  x_null <- lapply(1:sample_size, function(x) zero_mat)

  x_adj <- repr_nvd(x, representation = "adjacency")
  for (k in 1:sample_size) {
    num_vertices_old <- nrow(x_adj[[k]])
    if(num_vertices_old != 0){
      for (l in 1:num_vertices_old) {
        for (g in 1:num_vertices_old) {
          x_null[[k]][l,g] <- x_adj[[k]][l,g]
        }
      }
    }
  }

  nvd_from_adjacency_matrix(x_null, directed = directed, weighted = weighted)
}
