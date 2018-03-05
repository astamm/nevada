possibly_significant_areas_int <- function(x, y, Y, V, dimens, p_int) {

  m <- as.numeric(dim(Y)[2])
  pvalue <- rep(-1, m)

  for(j in 1:m) {

    xx <- sapply(x, select_nodes_int, zone = Y[,j], vertex = V, simplify = FALSE)
    yy <- sapply(y, select_nodes_int, zone = Y[,j], vertex = V, simplify = FALSE)

    pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"))[[2]]
    p_int[1,Y[,j]] <- apply(t(as.matrix(p_int[1,Y[,j]])), 2, function(x) max(x, pvalue[j]))
  }

  Z <- NULL
  t <- 1
  for (j in 1:m) {
    if (pvalue[j] < 0.05) {
      Z[[t]] <- combn(Y[,j], dimens)
      t <- t + 1
    }
  }

  possibly_significant_areas <- t(unique.matrix(t(do.call("cbind", Z))))

  results <- list(possibly_significant_areas, p_int)
  return(results)
}

possibly_significant_areas_ext <- function(x, y, Y, V, dimens, p_ext) {

  m <- as.numeric(dim(Y)[2])
  pvalue <- rep(-1, m)

  for(j in 1:m) {

    xx <- sapply(x, select_nodes_ext, zone = Y[,j], vertex = V, simplify = FALSE)
    yy <- sapply(y, select_nodes_ext, zone = Y[,j], vertex = V, simplify = FALSE)

    pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"))[[2]]
    inter <- combn(Y[,j],2)
    mm <- dim(inter)[2]
    for (k in 1:mm) {
      p_ext[inter[1,k], inter[2,k]] <- max(p_ext[inter[1,k], inter[2,k]], pvalue[j])
    }
  }

  Z <- NULL
  t <- 1
  for (j in 1:m) {
    if (pvalue[j] < 0.05) {
      Z[[t]] <- combn(Y[,j], dimens)
      t <- t + 1
    }

  }
  possibly_significant_areas <- t(unique.matrix(t(do.call("cbind", Z))))

  results <- list(possibly_significant_areas, p_ext)
  return(results)
}
