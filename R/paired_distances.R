#' Title Compute all the paired distances between networks of two samples
#'
#' \code{get_paired_distances} computes the paired distances between all the
#' elements of the two samples put together. The cardinality of the fist sample
#' is indicated with \eqn{n1} and that of the second one is indicated with
#' \eqn{n2}.
#'
#' @param x List of adjacency matrix/igraph object of the first sample.
#' @param y List of adjacency matrix/igraph object of the second sample.
#' @param representation The choosen representation for the elements of x and y.
#'   It can be: "adjacency", "laplacian" and "modularity".
#' @param distance The chosen distance for calculate the test statistic. It can
#'   be: "hamming", "frobenius", "spectral", "rooteuclidean".
#'
#' @return A matrix of dimension \eqn{(n1+n2)x(n1+n2)} containing the distances
#'   between all the elements of the two samples put together.
#' @export
#'
#' @examples x <- NULL
#' y <- NULL
#' for (i in 1:10) {
#' X <- watts.strogatz.game(1, 25, 3, 0.05)
#' Y <- barabasi.game(n, m=3, power=2, directed = FALSE)
#' adjX <- get_adjacency(X)
#' adjY <- get_adjacency(Y)
#' if (i==1) {
#' x <- list(adjX)
#' y <- list(adjY)
#' }
#' else {
#' x[[i]] <- adjX
#' y[[i]] <- adjY
#' }
#' }
#' get_paired_distances(x,y,"laplacian","spectral")
#'


get_paired_distances <- function(x, y, representation, distance) {
  n <- length(x)
  m <- length(y)
  z <- c(x, y)
  N <- n + m
  d <- matrix(data = rep(0, N*N), nrow = N, ncol = N)


  #DUE PROBLEMI: 1.non uso la funzione map del pacchetto purrr 2.ricalcolo la rappresentazione della stessa matrice più volte

  dist=rep(-1,N*(N-1)/2);
  k <- 1
  for (i in 1:(N-1)) {
    for (j in ((i+1):N)) {
      dist[k] <- switch(distance,
                        hamming <- get_hamming_distance(z[[i]], z[[j]], representation),
                        frobenius <- get_frobenius_distance(z[[i]], z[[j]], representation),
                        spectral <- get_spectral_distance(z[[i]], z[[j]], representation),
                        rooteuclidean <- get_rooteuclidean_distance(z[[i]], z[[j]], representation))
      k <- k+1
    }
  }

  d[lower.tri(d)]=dist
  d=t(d)
  d=d+t(d)


  return(d)

}
