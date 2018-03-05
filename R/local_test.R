network_localtest2p <- function(x, y, location) {
  p_global <- test_twosample(x, y, statistic = c("sot", "lot"))[2]

  areas <- sort(unique(location))
  N <- length(areas)

  V <- NULL #lista in cui l'elemento i contiene i vertici che stanno nella zona i

  for (i in areas) {
    V[[i]] <- which(location==i)
  }


  if (p_global < 0.05) {

    xINT <- sapply(x, select_nodes_int, zone = 1:N, vertex = V, simplify = FALSE)
    yINT <- sapply(y, select_nodes_int, zone = 1:N, vertex = V, simplify = FALSE)

    p_global_int <- test_twosample(xINT,yINT, statistic = c("sot", "lot"))[2]


    xEXT <- sapply(x, select_nodes_ext, zone = 1:N, vertex = V, simplify = FALSE)
    yEXT <- sapply(y, select_nodes_ext, zone = 1:N, vertex = V, simplify = FALSE)

    p_global_ext <- test_twosample(xEXT,yEXT, statistic = c("sot", "lot"))[2]

  } else {
    p_global_int <- 1
    p_global_ext <- 1
    p_int <- matrix(data = rep(1,N), nrow=1, ncol = N)
    p_ext <- matrix(data = rep(1 ,N*N), nrow = N, ncol = N)
  }





  ## INTRA LOCAL TESTS

  p_int <- matrix(data = rep(-1,N), nrow=1, ncol = N)

  if (p_global_int < 0.05) {

    results_int <- possibly_significant_areas_int(x, y, combn(1:N, (N-1)), V, N-2, p_int)
    areas_to_be_tested_int <- results_int[[1]]
    p_int <- results_int[[2]]

    if (N>3){
      for(i in (N-2):2) {
        results_int <- possibly_significant_areas_int(x, y, areas_to_be_tested_int, V, i-1, p_int)
        areas_to_be_tested_int <- results_int[[1]]
        p_int <- results_int[[2]]
      }


      m <- as.numeric(dim(areas_to_be_tested_int)[2])
      pvalue <- rep(-1, m)

      for(j in 1:m) {

        xx <- sapply(x, select_nodes_int, zone = areas_to_be_tested_int[,j], vertex = V, simplify = FALSE)
        yy <- sapply(y, select_nodes_int, zone = areas_to_be_tested_int[,j], vertex = V, simplify = FALSE)

        pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"))[[2]]
        p_int[1,areas_to_be_tested_int[,j]] <- apply(t(as.matrix(p_int[1,areas_to_be_tested_int[,j]])), 2, function(x) max(x, pvalue[j]))
      }
    }

  }




  ## INTER LOCAL TESTS

  p_ext <- matrix(data = rep(-1 ,N*N), nrow = N, ncol = N)

  if (p_global_ext < 0.05) {

    results_ext <- possibly_significant_areas_ext(x, y, combn(1:N, (N-1)), V, N-2, p_ext)
    areas_to_be_tested_ext <- results_ext[[1]]
    p_ext <- results_ext[[2]]

    if (N>4){
      for(i in (N-2):3) {
        results_ext <- possibly_significant_areas_ext(x, y, areas_to_be_tested_ext, V, i-1, p_ext)
        areas_to_be_tested_ext <- results_ext[[1]]
        p_ext <- results_ext[[2]]
      }
    }


    m <- as.numeric(dim(areas_to_be_tested_ext)[2])
    pvalue <- rep(-1, m)

    for(j in 1:m) {

      xx <- sapply(x, select_nodes_ext, zone = areas_to_be_tested_ext[,j], vertex = V, simplify = FALSE)
      yy <- sapply(y, select_nodes_ext, zone = areas_to_be_tested_ext[,j], vertex = V, simplify = FALSE)

      pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"))[[2]]
      inter <- combn(areas_to_be_tested_ext[,j],2)
      m <- dim(inter)[2]
      for (k in 1:m) {
        p_ext[inter[1,k], inter[2,k]] <- max(p_ext[inter[1,k], inter[2,k]], pvalue[j])
      }
    }

  }

  return(list(p_global_int, p_int, p_global_ext, p_ext))

}
