subgroups <- function(Y, dimens) {
  Z <- apply(Y, MARGIN = 2, function(x) utils::combn(x, dimens))
  W <- matrix(data = -1, nrow = dimens, ncol = dim(Z)[2] * dim(Z)[1] / dimens)
  k <- 1
  for (i in seq(1, dim(Z)[1], by = dimens)) {
    W[, k:(k + dim(Z)[2] - 1)] <- Z[i:(i+dimens-1), ]
    k <- k + dim(Z)[2]
  }
  W[, !duplicated(W, MARGIN = 2), drop = FALSE]
}

subgroups_orig <- function(Y, dimens) {

  if (dimens == 1) {
    Z <- apply(as.matrix(Y), MARGIN = 2, function(x) utils::combn(x, dimens))
    W <- matrix(data = -1, nrow = dimens, ncol = (dim(Z)[2])*((dim(Z)[1])/dimens))
    k <- 1
    for (i in 1:2) {
      for (j in 1:dim(Z)[2]) {
        W[,k] <- Z[i, j]
        k <- k+1
      }
    }
    W <- as.matrix(W[,!duplicated.matrix(W, MARGIN = 2)])
  }

  if (dimens>1) {
    Z <- apply(as.matrix(Y), MARGIN = 2, function(x) utils::combn(x, dimens))
    W <- matrix(data = -1, nrow = dimens, ncol = (dim(Z)[2])*((dim(Z)[1])/dimens))
    k <- 1
    for (i in seq(1, dim(Z)[1]-1, by = dimens)) {
      for (j in 1:dim(Z)[2]) {
        W[,k] <- Z[i:(i+(dimens-1)), j]
        k <- k+1
      }
    }
    W <- as.matrix(W[,!duplicated.matrix(W, MARGIN = 2)])
  }

  return(W)
}
