network_localtest2p <- function(x, y, location, representation = "adjacency", distance = "frobenius", alpha = 0.05, B = 1000) {
  set.seed(1234)
  p_global <- test_twosample(x, y, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[2]

  areas <- sort(unique(location))
  N <- length(areas)

  V <- list() #lista in cui l'elemento i contiene i vertici che stanno nella zona i

  for (i in areas) {
    V[[i]] <- which(location==i)
  }

  dimens1 <- which(data.frame(table(location))[,2] == 1)


  if (N == dim(x[[1]])[1]) {
    if (p_global < alpha) {

      p_single_node <- matrix(data = rep(as.numeric(p_global), N*N), nrow = N, ncol = N)

      if (N == 2) {
        pvalue <- rep(-1, N)
        Y <- combn(c(1,2),1)
        m <- as.numeric(dim(Y)[2])

        for(j in 1:m) {

          xx <- sapply(x, select_nodes_tot, zone = Y[,j], vertex = V, simplify = FALSE)
          yy <- sapply(y, select_nodes_tot, zone = Y[,j], vertex = V, simplify = FALSE)

          pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
          p_single_node[Y[,j],Y[,j]] <- max(as.numeric(p_single_node[Y[,j], Y[,j]]), pvalue[j])
        }
      }
      else {

        results_tot <- possibly_significant_areas_tot(x, y, combn(1:N, (N-1)), V, N-2, p_single_node, representation = representation, distance = distance, alpha = alpha, N = N, B = B)


        if(results_tot[[1]]){
          p_single_node <- results_tot[[2]]
        }
        else {
          areas_to_be_tested_tot <- as.matrix(results_tot[[2]])
          p_single_node <- results_tot[[3]]

          if (N>3) {

            for (i in (N-2):2) {
              results_tot <- possibly_significant_areas_tot(x, y, areas_to_be_tested_tot, V, i-1, p_single_node, representation = representation, distance = distance, alpha = alpha, N = N, B = B)
              if (results_tot[[1]]) {
                p_single_node <- results_tot[[2]]
                break
              }
              else {
                areas_to_be_tested_tot <- results_tot[[2]]
                p_single_node <- results_tot[[3]]
              }
            }
          }
        }
      }

    }
    else
      p_single_node <- matrix(data = rep(1,N), nrow = N, ncol = N)
    return(p_single_node)
  }

  if (N != dim(x[[1]])[1]) {


  if (p_global < alpha) {

    xINT <- sapply(x, select_nodes_int, zone = 1:N, vertex = V, simplify = FALSE)
    yINT <- sapply(y, select_nodes_int, zone = 1:N, vertex = V, simplify = FALSE)

    set.seed(1234)
    p_global_int <- test_twosample(xINT,yINT, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[2]


    xEXT <- sapply(x, select_nodes_ext, zone = 1:N, vertex = V, simplify = FALSE)
    yEXT <- sapply(y, select_nodes_ext, zone = 1:N, vertex = V, simplify = FALSE)

    set.seed(1234)
    p_global_ext <- test_twosample(xEXT,yEXT, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[2]

  } else {
    p_global_int <- 1
    p_global_ext <- 1
    p_int <- matrix(data = rep(1,N), nrow=1, ncol = N)
    p_ext <- matrix(data = rep(1 ,N*N), nrow = N, ncol = N)
  }


  exit_tot <- FALSE
  exit_int <- FALSE
  exit_ext <- FALSE



  ## TOTAL LOCAL TESTS

  p_tot <- matrix(data = rep(as.numeric(p_global), N*N), nrow = N, ncol = N)

  if (p_global < alpha) {

    if (N == 2) {
      pvalue <- rep(-1, N)
      Y <- combn(c(1,2),1)
      m <- as.numeric(dim(Y)[2])

      for(j in 1:m) {

        xx <- sapply(x, select_nodes_tot, zone = Y[,j], vertex = V, simplify = FALSE)
        yy <- sapply(y, select_nodes_tot, zone = Y[,j], vertex = V, simplify = FALSE)

        set.seed(1234)
        pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
        p_tot[Y[,j],Y[,j]] <- max(as.numeric(p_tot[Y[,j], Y[,j]]), pvalue[j])
      }
    }
    else {

        results_tot <- possibly_significant_areas_tot(x, y, combn(1:N, (N-1)), V, N-2, p_tot, representation = representation, distance = distance, alpha = alpha, N = N, B = B)


    if(results_tot[[1]]){
      p_tot <- results_tot[[2]]
    }
    else {
      areas_to_be_tested_tot <- as.matrix(results_tot[[2]])
      p_tot <- results_tot[[3]]

      if (N>3) {

    for (i in (N-2):2) {
      results_tot <- possibly_significant_areas_tot(x, y, areas_to_be_tested_tot, V, i-1, p_tot, representation = representation, distance = distance, alpha = alpha, N = N, B = B)
      if (results_tot[[1]]) {
        p_tot <- results_tot[[2]]
        break
      }
      else {
        areas_to_be_tested_tot <- results_tot[[2]]
        p_tot <- results_tot[[3]]
      }
    }
      }
    }


    if (results_tot[[1]])
      p_tot <- p_tot
    else {
      no_test <- which(match(areas_to_be_tested_tot, dimens1, nomatch = 0) != 0)
      if (length(no_test) > 0) {
        areas_to_be_tested_tot <- areas_to_be_tested_tot[-no_test]
      }
      m <- length(areas_to_be_tested_tot)

      if (m > 0) {
        pvalue <- rep(-1, m)

        for(j in 1:m) {

          xx <- sapply(x, select_nodes_tot, zone = areas_to_be_tested_tot[j], vertex = V, simplify = FALSE)
          yy <- sapply(y, select_nodes_tot, zone = areas_to_be_tested_tot[j], vertex = V, simplify = FALSE)

          set.seed(1234)
          pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
          p_tot[areas_to_be_tested_tot[j],areas_to_be_tested_tot[j]] <- max(as.numeric(p_tot[areas_to_be_tested_tot[j], areas_to_be_tested_tot[j]]), pvalue[j])
        }
      }

    }
    }

  }


  ## INTRA LOCAL TESTS

  p_int <- matrix(data = rep(as.numeric(p_global_int),N), nrow=1, ncol = N)

  if (p_global_int < alpha) {

    if (N == 2) {
      pvalue <- rep(-1, N)
      Y <- combn(c(1,2),1)
      m <- as.numeric(dim(Y)[2])
      for(j in 1:m) {

        xx <- sapply(x, select_nodes_int, zone = Y[,j], vertex = V, simplify = FALSE)
        yy <- sapply(y, select_nodes_int, zone = Y[,j], vertex = V, simplify = FALSE)

        set.seed(1234)
        pvalue[j] <- as.numeric(test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]])
        p_int[1,j] <- max(as.numeric(p_int[1,j]), pvalue[j])
      }
      p_int <- as.vector(p_int)
      p_int <- pmax(diag(p_tot), p_int)
    }
    else {
      results_int <- possibly_significant_areas_int(x, y, combn(1:N, (N-1)), V, N-2, p_int, dimens1, representation = representation, distance = distance, alpha = alpha, N = N, B = B)

    if(results_int[[1]]){
      p_int <- results_int[[2]]
      p_int <- pmax(p_int, diag(p_tot))
    }
    else {
      areas_to_be_tested_int <- as.matrix(results_int[[2]])
      p_int <- results_int[[3]]

      if (N>3){
        for(i in (N-2):2) {
          results_int <- possibly_significant_areas_int(x, y, areas_to_be_tested_int, V, i-1, p_int, dimens1, representation = representation, distance = distance, alpha = alpha, N = N, B = B)
          if (results_int[[1]]) {
            p_int <- results_int[[2]]
            break
          }
          else {
            areas_to_be_tested_int <- results_int[[2]]
            p_int <- results_int[[3]]
          }
        }
      }
    }

    if (results_int[[1]])
      p_int <- pmax(p_int, diag(p_tot))
    else {
      no_test <- which(match(areas_to_be_tested_int, dimens1, nomatch = 0) != 0)
      if (length(no_test) > 0) {
        areas_to_be_tested_int <- areas_to_be_tested_int[-no_test]
      }
      m <- length(areas_to_be_tested_int)

      if (m > 0) {
        pvalue <- rep(-1, m)

        for(j in 1:m) {

          xx <- sapply(x, select_nodes_int, zone = areas_to_be_tested_int[j], vertex = V, simplify = FALSE)
          yy <- sapply(y, select_nodes_int, zone = areas_to_be_tested_int[j], vertex = V, simplify = FALSE)

          set.seed(1234)
          pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
          p_int[areas_to_be_tested_int[j]] <- max(p_int[areas_to_be_tested_int[j]], pvalue[j])
          #p_int[1,areas_to_be_tested_int[j]] <- apply(t(as.matrix(p_int[1,areas_to_be_tested_int[j]])), 2, function(x) max(as.numeric(x), pvalue[j]))
        }
      }

    }
    # p_int[1,] <- apply(t(as.matrix(p_int[1,])), 2, function(x) max(as.numeric(x), pvalue))
      p_int <- pmax(diag(p_tot), p_int)
    }

  }




  ## INTER LOCAL TESTS

  p_ext <- matrix(data = rep(as.numeric(p_global_ext), N*N), nrow = N, ncol = N)

  if (N>2) {

  if (p_global_ext < alpha) {

    if (N == 3) {
      pvalue <- rep(-1, N)
      Y <- combn(1:3, 2)
      m <- as.numeric(dim(Y)[2])
      for(j in 1:m) {

        xx <- sapply(x, select_nodes_ext, zone = Y[,j], vertex = V, simplify = FALSE)
        yy <- sapply(y, select_nodes_ext, zone = Y[,j], vertex = V, simplify = FALSE)

        set.seed(1234)
        pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
        inter <- combn(Y[,j],2)
        m <- dim(inter)[2]
        for (k in 1:m) {
          p_ext[inter[1,k], inter[2,k]] <- max(as.numeric(p_ext[inter[1,k], inter[2,k]]), as.numeric(pvalue[j]))
        }
      }
      p_ext <- pmax(p_ext, p_tot)
    }
    else {
       results_ext <- possibly_significant_areas_ext(x, y, combn(1:N, (N-1)), V, N-2, p_ext, representation = representation, distance = distance, alpha = alpha, N = N, B = B)

    if(results_ext[[1]]){
      p_ext <- results_ext[[2]]
      p_ext <- pmax(p_ext, p_tot)
    }
    else {
      areas_to_be_tested_ext <- as.matrix(results_ext[[2]])
      p_ext <- results_ext[[3]]


    if (N>4){
      for(i in (N-2):3) {
        results_ext <- possibly_significant_areas_ext(x, y, areas_to_be_tested_ext, V, i-1, p_ext, representation = representation, distance = distance, alpha = alpha, N = N, B = B)
        if (results_ext[[1]]) {
          p_ext <- results_ext[[2]]
          break
        }
        else {
          areas_to_be_tested_ext <- results_ext[[2]]
          p_ext <- results_ext[[3]]
        }

      }
    }
  }


    if (results_ext[[1]])
      p_ext <- pmax(p_ext, p_tot)
    else {
      m <- as.numeric(dim(as.matrix(areas_to_be_tested_ext))[2])

      if (m > 0) {
        pvalue <- rep(-1, m)

        for(j in 1:m) {

          xx <- sapply(x, select_nodes_ext, zone = as.matrix(areas_to_be_tested_ext)[,j], vertex = V, simplify = FALSE)
          yy <- sapply(y, select_nodes_ext, zone = as.matrix(areas_to_be_tested_ext)[,j], vertex = V, simplify = FALSE)

          set.seed(1234)
          pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
          inter <- combn(sort(as.matrix(areas_to_be_tested_ext)[,j]),2)
          mm <- dim(inter)[2]
          for (k in 1:mm) {
            p_ext[inter[1,k], inter[2,k]] <- max(as.numeric(p_ext[inter[1,k], inter[2,k]]), as.numeric(pvalue[j]))
          }
        }
      }

      p_ext <- pmax(p_ext, p_tot)
    }


    }

  }
  }

  p_ext[lower.tri(p_ext, diag = TRUE)] <- NA

  return(list(p_global_int, p_int, p_global_ext, p_ext))

}
}
