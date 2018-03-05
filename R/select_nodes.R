select_nodes_tot <- function(x, zone, vertex) {
  vertici <- vertex[[zone[1]]]
  r <- 2
  for (l in zone[-1]) {
    vertici <- c(vertici, vertex[[zone[r]]])
    r <- r+1
  }

  xtot <- x[sort(vertici), sort(vertici)]
  return(xtot)
}

select_nodes_int <- function(x, zone, vertex) {
  vertici <- vertex[[zone[1]]]
  r <- 2
  for (l in zone[-1]) {
    vertici <- c(vertici, vertex[[zone[r]]])
    r <- r+1
  }

  ZERO <- NULL
  r <- 1
  for (k in zone) {
    for (h in zone[-r]) {
      zeros <- rbind(matrix(unlist(expand.grid(vertex[[k]], vertex[[h]])), nrow = length(vertex[[k]])*length(vertex[[h]]), ncol = 2), matrix(unlist(expand.grid(vertex[[h]], vertex[[k]])), nrow = length(vertex[[k]])*length(vertex[[h]]), ncol = 2))
      ZERO <- rbind(ZERO, zeros)
    }
    r <- r+1
  }

  xint <- x
  xint[ZERO] <- 0
  xint <- xint[sort(vertici), sort(vertici)]

  return(xint)
}

select_nodes_ext <- function(x, zone, vertex) {
  xtot <- select_nodes_tot(x, zone, vertex)
  xint <- select_nodes_int(x, zone, vertex)
  xext <- xtot - xint

  return(xext)
}
