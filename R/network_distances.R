#' Title Compute the Hamming distance between two networks
#'
#' If \eqn{X} is the matrix representation for the graph \eqn{x} and \eqn{Y} is
#' that of the graph \eqn{y}, the Hamming distance between \eqn{x} and \eqn{y}
#' is given by \deqn{1/N(N-1)sum_{i,j} |X_ij - Y_ij|} where \eqn{N} is the
#' number of vertices of \eqn{x} and \eqn{y}.
#'
#' @param x Adjacency matrix/igraph object with N vertices.
#' @param y Adjacency matrix/igraph object with N vertices.
#' @param representation The choosen matrix representation for the elements of x
#'   and y. It can be: "adjacency", "laplacian" and "modularity".
#'
#' @return \code{get_hamming_distance} returns a scalar.
#' @export
#'
#' @examples g1 <- erdos.renyi.game(20, 0.1)
#' g2 <- erdos.renyi.game(20, 0.2)
#' get_hamming_distance(g1, g2, "adjacency")
get_hamming_distance <- function(x, y, representation) {
  X <- switch(representation,
              adjacency <- get_adjacency(x),
              laplacian <- get_laplacian(x),
              modularity <- get_modularity(x))
  Y <- switch(representation,
              adjacency <- get_adjacency(Y),
              laplacian <- get_laplacian(Y),
              modularity <- get_modularity(Y))

  d <- sum(abs(X-Y))
  return(trunc(d*10^4)/10^4)

}


#' Title Compute the Frobenius distance between two networks
#'
#' If \eqn{X} is the matrix representation for the graph \eqn{x} and \eqn{Y}
#' is that of the graph \eqn{y}, the Frobenius distance between \eqn{x} and
#' \eqn{y} is given by \deqn{sqrt(sum_{i,j} (X_ij - Y_ij)^2).}
#'
#' @param x Adjacency matrix/igraph object with N vertices.
#' @param y Adjacency matrix/igraph object with N vertices.
#' @param representation The choosen matrix representation for the elements of x
#'   and y. It can be: "adjacency", "laplacian" and "modularity".
#'
#' @return \code{get_frobenius_distance} returns a scalar.
#' @export
#'
#' @examples g1 <- erdos.renyi.game(20, 0.1)
#' g2 <- erdos.renyi.game(20, 0.2)
#' get_frobenius_distance(g1, g2, "adjacency")
get_frobenius_distance <- function(x, y, representation) {
  X <- switch(representation,
              adjacency <- get_adjacency(x),
              laplacian <- get_laplacian(x),
              modularity <- get_modularity(x))
  Y <- switch(representation,
              adjacency <- get_adjacency(Y),
              laplacian <- get_laplacian(Y),
              modularity <- get_modularity(Y))

  d <- sqrt(sum((X-Y)^2))
  return(trunc(d*10^4)/10^4)

}


#' Title  Compute the Spectral distance between two networks
#'
#' If \eqn{X} is the matrix representation for the graph \eqn{x} and \eqn{Y}
#' is that of the graph \eqn{y}, consider the vectors \eqn{a} and \eqn{b} of
#' the eigenvalues of \eqn{X} and \eqn{Y}, respectively. The Spectral distance
#' between \eqn{x} and \eqn{y} is given by \deqn{sqrt(sum_{i} (a_i - b_i)^2).}
#' This distance gives rise to classes of equivalence.
#'
#' @param x Adjacency matrix/igraph object with N vertices.
#' @param y Adjacency matrix/igraph object with N vertices.
#' @param representation The choosen matrix representation for the elements of x
#'   and y. It can be: "adjacency", "laplacian" and "modularity".
#'
#' @return \code{get_spectral_distance} returns a scalar.
#' @export
#'
#' @examples g1 <- erdos.renyi.game(20, 0.1)
#' g2 <- erdos.renyi.game(20, 0.2)
#' get_spectral_distance(g1, g2, "laplacian")
get_spectral_distance <- function(x, y, representation) {
  X <- switch(representation,
              adjacency <- get_adjacency(x),
              laplacian <- get_laplacian(x),
              modularity <- get_modularity(x))
  Y <- switch(representation,
              adjacency <- get_adjacency(Y),
              laplacian <- get_laplacian(Y),
              modularity <- get_modularity(Y))

  rX <- eigen(X, symmetric=TRUE);
  for (t in 1:n) {
    if (abs(rX$values[t])<1e-10) {
      rX$values[t] <- 0;
    }
  }
  dlX <- rX$values

  rY <- eigen(Y, symmetric=TRUE);
  for (t in 1:n) {
    if (abs(rY$values[t])<1e-10) {
      rY$values[t] <- 0;
    }
  }
  dlY <- rY$values

  d <- sqrt(sum((dlX-dlY)^2))
  return(trunc(d*10^4)/10^4)

}


#' Title Compute the Root Euclidean distance between two networks
#'
#' If \eqn{X} is the matrix representation for the graph \eqn{x} and \eqn{Y} is
#' that of the graph \eqn{y}, consider the spectral decomposition of \eqn{X} and
#' \eqn{Y}: \deqn{X=VAV^(-1)} \deqn{X=UBU^(-1)} where \eqn{V} and \eqn{U} are
#' the matrix whose colums are the eigenvectors of \eqn{X} and \eqn{Y},
#' respectively; \eqn{A} and \eqn{B} are the diagonal matrix with elements the
#' eigenvalues of \eqn{X} and \eqn{Y}, respectively. The Root Euclidean distance
#' between \eqn{x} and \eqn{y} is given by \deqn{sqrt(sum_{i} (Vsqrt(A)V^(-1) -
#' Usqrt(B)U^(-1))^2)} where the square root of a matrix is the matrix with the
#' square rotted elements.
#'
#' @param x Adjacency matrix/igraph object with N vertices.
#' @param y Adjacency matrix/igraph object with N vertices.
#' @param representation The choosen matrix representation for the elements of x
#'   and y. It can be: "adjacency", "laplacian" and "modularity".
#'
#' @return \code{get_spectral_distance} returns a scalar.
#' @export
#'
#' @examples g1 <- erdos.renyi.game(20, 0.1)
#' g2 <- erdos.renyi.game(20, 0.2)
#' get_rooteuclidean_distance(g1, g2, "laplacian")
get_rooteuclidean_distance <- function(x, y, representation) {
  X <- switch(representation,
              adjacency <- get_adjacency(x),
              laplacian <- get_laplacian(x),
              modularity <- get_modularity(x))
  Y <- switch(representation,
              adjacency <- get_adjacency(Y),
              laplacian <- get_laplacian(Y),
              modularity <- get_modularity(Y))

  rX <- eigen(X, symmetric=TRUE)
  for (t in 1:n) {
    if (abs(rX$values[t])<1e-10) {
      rX$values[t] <- 0
    }
  }
  dlX <- rX$vectors%*%sqrt(diag(rX$values))%*%solve(rX$vectors)

  rY <- eigen(Y, symmetric=TRUE)
  for (t in 1:n) {
    if (abs(rY$values[t])<1e-10) {
      rY$values[t] <- 0;
    }
  }
  dlY <- rY$vectors%*%sqrt(diag(rY$values))%*%solve(rY$vectors)

  d <- sqrt(sum((dlX-dlY)^2))
  return(trunc(d*10^4)/10^4)

}
