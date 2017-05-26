as_adjacency <- function(matrix) {
  class(matrix) <- c(class(matrix), "adjacency")
}

as_laplacian <- function(matrix) {
  class(matrix) <- c(class(matrix), "laplacian")
}

as_modularity <- function(matrix) {
  class(matrix) <- c(class(matrix), "modularity")
}

is_adjacency <- function(matrix) {
  "adjacency" %in% class(matrix)
}

is_laplacian <- function(matrix) {
  "laplacian" %in% class(matrix)
}

is_modularity <- function(matrix) {
  "modularity" %in% class(matrix)
}
