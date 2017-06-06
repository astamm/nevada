
get_adjacency <- function(network){
  adjacent <- get.adjacency(network, type = "both")
}

get_laplacian <- function(network){
  laplacian <- laplacian_matrix(network)
}

get_modularity <- function(network){
  modularity <- modularity_matrix(network, rep(1,n)) #NUMERO DI VERTICI?
}

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
