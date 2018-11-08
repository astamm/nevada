possibly_significant_areas_tot <-
  function(x,
           y,
           Y,
           vertex_set,
           dimens,
           p_tot,
           representation = representation,
           distance = distance,
           statistic = statistic,
           alpha = alpha,
           B = B,
           seed = seed) {
    N <- length(vertex_set)
    m <- length(Y)
    pvalue <- rep(-1, m)

    for (j in 1:m) {
      xx <- lapply(x,
                   select_nodes_tot,
                   subsets = Y[[j]],
                   vertex_set = V)
      yy <-
        lapply(y,
               select_nodes_tot,
               subsets = Y[[j]],
               vertex_set = V)

      pvalue[j] <- test2_global(
        xx,
        yy,
        representation = representation,
        distance = distance,
        statistic = statistic,
        alpha = alpha,
        B = B,
        seed = seed
      )$pvalue

      for (k in 1:(dimens + 1))
        for (l in k:(dimens + 1))
          p_tot[Y[[j]][k], Y[[j]][l]] <-
        max(p_tot[Y[[j]][k], Y[[j]][l]], pvalue[j])
    }

    Z <- unique(unlist(Y[pvalue < alpha]))

    if (length(Z) == 0) {
      exit_tot <- TRUE
      return(list(exit_tot, p_tot))
    }

    # Y is now a list, not a matrix anymore
    exit_tot <- FALSE
    nonSign <- Y[pvalue >= alpha]
    if (length(nonSign) > 0) {
      if (dimens > 1) {
        previously_excluded <-
          which(prodlim::row.match(data.frame(t(combn(
            1:N, (dimens + 1)
          ))), data.frame(t(Y)), nomatch = 0) == 0)
        if (length(previously_excluded) > 0) {
          to_be_excluded <-
            subgroups(combn(1:N, (dimens + 1))[, previously_excluded], dimens)
          eliminate1 <-
            which(prodlim::row.match(data.frame(t(
              subgroups(Y[, which(pvalue < alpha)], dimens)
            )), data.frame(t(
              to_be_excluded
            )), nomatch = 0) > 0)
          eliminate2 <-
            which(prodlim::row.match(data.frame(t(
              subgroups(Y[, which(pvalue < alpha)], dimens)
            )), data.frame(t(
              subgroups(Y[, which(pvalue >= alpha)], dimens)
            )), nomatch = 0) > 0)
          eliminate <- unique(sort(c(eliminate1, eliminate2)))
        }
        else
          eliminate <-
            which(prodlim::row.match(data.frame(t(
              subgroups(Y[, which(pvalue < alpha)], dimens)
            )), data.frame(t(
              subgroups(Y[, which(pvalue >= alpha)], dimens)
            )), nomatch = 0) > 0)
      }

      if (dimens == 1) {
        previously_excluded <-
          which(prodlim::row.match(data.frame(t(combn(
            1:N, (dimens + 1)
          ))), data.frame(t(Y)), nomatch = 0) == 0)
        if (length(previously_excluded) > 0) {
          to_be_excluded <-
            subgroups(combn(1:N, (dimens + 1))[, previously_excluded], dimens)
          eliminate1 <-
            which(prodlim::row.match(
              data.frame(subgroups(Y[, which(pvalue < alpha)], dimens)),
              data.frame(to_be_excluded),
              nomatch = 0
            ) > 0)
          eliminate2 <-
            which(prodlim::row.match(
              data.frame(subgroups(Y[, which(pvalue < alpha)], dimens)),
              data.frame(subgroups(Y[, which(pvalue >= alpha)], dimens)),
              nomatch = 0
            ) > 0)
          eliminate <- unique(sort(c(eliminate1, eliminate2)))
        }
        else
          eliminate <-
            which(prodlim::row.match(data.frame(t(
              subgroups(Y[, which(pvalue < alpha)], dimens)
            )), data.frame(t(
              subgroups(Y[, which(pvalue >= alpha)], dimens)
            )), nomatch = 0) > 0)
      }

      if (length(eliminate) == dim(subgroups(Y[, which(pvalue < alpha)], dimens))[2]) {
        exit_tot <- TRUE
        return(list(exit_tot, p_tot))
      }
      if (length(eliminate) != dim(subgroups(Y[, which(pvalue < alpha)], dimens))[2]) {
        if (dimens > 1)
          c <-
            subgroups(Y[, which(pvalue < alpha)], dimens)[, -eliminate]
        if (dimens == 1)
          c <-
            t(subgroups(Y[, which(pvalue < alpha)], dimens))[, -eliminate]
      }

    }

    if (dim(nonSign)[2] == 0) {
      previously_excluded <-
        which(prodlim::row.match(data.frame(t(combn(
          1:N, (dimens + 1)
        ))), data.frame(t(Y)), nomatch = 0) == 0)
      if (length(previously_excluded) > 0) {
        to_be_excluded <-
          subgroups(combn(1:N, (dimens + 1))[, previously_excluded], dimens)
        eliminate <-
          which(prodlim::row.match(data.frame(t(
            subgroups(Y[, which(pvalue < alpha)], dimens)
          )), data.frame(t(
            to_be_excluded
          )), nomatch = 0) > 0)

        if (length(eliminate) == dim(subgroups(Y[, which(pvalue < alpha)], dimens))[2]) {
          exit_tot <- TRUE
          return(list(exit_tot, p_tot))
        }
        if (length(eliminate) != dim(subgroups(Y[, which(pvalue < alpha)], dimens))[2]) {
          c <-
            as.matrix(subgroups(Y[, which(pvalue < alpha)], dimens)[, -eliminate])
        }
      }
      else
        c <- subgroups(Y[, which(pvalue < alpha)], dimens)
    }

    possibly_significant_areas <- c
    results <- list(exit_tot, possibly_significant_areas, p_tot)
    results
  }

possibly_significant_areas_int <- function(x, y, Y, V, dimens, p_int, dimens1, representation = representation, distance = distance, alpha = alpha, N = N, B = B) {

  m <- as.numeric(dim(Y)[2])
  pvalue <- rep(-1, m)

  for(j in 1:m) {

    xx <- sapply(x, select_nodes_int, zone = Y[,j], vertex = V, simplify = FALSE)
    yy <- sapply(y, select_nodes_int, zone = Y[,j], vertex = V, simplify = FALSE)

    set.seed(1234)
    pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
    p_int[1,Y[,j]] <- apply(t(as.matrix(p_int[1,Y[,j]])), 2, function(x) max(as.numeric(x), pvalue[j]))
  }

  #Z <- NULL
  #t <- 1
  #for (j in 1:m) {
   # if (pvalue[j] < 0.05) {
    #  Z[[t]] <- combn(Y[,j], dimens)
     # t <- t + 1
    #}
  #}

  #possibly_significant_areas <- t(unique.matrix(t(do.call("cbind", Z))))

  Z <- unique(as.vector(Y[,which(pvalue<alpha)]))

  if (length(Z) == 0) {
    exit_int <- TRUE
    return(list(exit_int, p_int))
  }
  else {
    exit_int <- FALSE
    nonSign <- as.matrix(Y[,which(pvalue>=alpha)])
    if (dim(nonSign)[2] > 0) {
      if (dimens > 1) {
        previously_excluded <- which(prodlim::row.match(data.frame(t(combn(1:N, (dimens+1)))), data.frame(t(Y)), nomatch = 0)==0)
        if (length(previously_excluded)>0) {
          to_be_excluded <- subgroups(combn(1:N, (dimens+1))[, previously_excluded], dimens)
          eliminate1 <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(to_be_excluded)), nomatch = 0)>0)
          eliminate2 <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(subgroups(Y[,which(pvalue>=alpha)], dimens))), nomatch = 0)>0)
          eliminate <- unique(sort(c(eliminate1, eliminate2)))
        }
        else
          eliminate <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(subgroups(Y[,which(pvalue>=alpha)], dimens))), nomatch = 0)>0)
      }
      #eliminate <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(subgroups(Y[,which(pvalue>=alpha)], dimens))), nomatch = 0)>0)
      if (dimens == 1) {
        previously_excluded <- which(prodlim::row.match(data.frame(t(combn(1:N, (dimens+1)))), data.frame(t(Y)), nomatch = 0)==0)
        if (length(previously_excluded)>0) {
          to_be_excluded <- subgroups(combn(1:N, (dimens+1))[, previously_excluded], dimens)
          eliminate1 <- which(prodlim::row.match(data.frame(subgroups(Y[,which(pvalue<alpha)], dimens)), data.frame(to_be_excluded), nomatch = 0)>0)
          eliminate2 <- which(prodlim::row.match(data.frame(subgroups(Y[,which(pvalue<alpha)], dimens)), data.frame(subgroups(Y[,which(pvalue>=alpha)], dimens)), nomatch = 0)>0)
          eliminate <- unique(sort(c(eliminate1, eliminate2)))
        }
        else
          eliminate <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(subgroups(Y[,which(pvalue>=alpha)], dimens))), nomatch = 0)>0)
      }
      if (length(eliminate) == dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
        exit_int <- TRUE
        return(list(exit_int, p_int))
      }
      if (length(eliminate) != dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
        if (dimens > 1)
          c <- subgroups(Y[,which(pvalue < alpha)], dimens)[,-eliminate]
        if (dimens == 1)
          c <- t(subgroups(Y[,which(pvalue < alpha)], dimens))[,-eliminate]
      }

    }

    if (dim(nonSign)[2] == 0) {
      previously_excluded <- which(prodlim::row.match(data.frame(t(combn(1:N, (dimens+1)))), data.frame(t(Y)), nomatch = 0)==0)
      if (length(previously_excluded)>0) {
        to_be_excluded <- subgroups(combn(1:N, (dimens+1))[, previously_excluded], dimens)
        eliminate <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(to_be_excluded)), nomatch = 0)>0)

        if (length(eliminate) == dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
          exit_int <- TRUE
          return(list(exit_int, p_int))
        }
        if (length(eliminate) != dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
          c <- as.matrix(subgroups(Y[,which(pvalue < alpha)], dimens)[,-eliminate])
        }
      }
      else
        c <- subgroups(Y[,which(pvalue<alpha)], dimens)
    }

    if (dimens <= length(dimens1)) {
      no_int <- combn(dimens1, dimens)
      eliminate_int <- which(prodlim::row.match(data.frame(t(c)), data.frame(t(no_int)), nomatch = 0)>0)
      if (length(eliminate_int)>0) {
        if (length(eliminate_int) == dim(c)[2]) {
          exit_int <- TRUE
          return(list(exit_int, p_int))
        }
        if (length(eliminate_int) != dim(c)[2]) {
          c <- c[,-eliminate_int]
        }
      }


    }

    possibly_significant_areas <- c
    #possibly_significant_areas <- combn(Z, dimens)
    results <- list(exit_int, possibly_significant_areas, p_int)
    return(results)
  }
}

possibly_significant_areas_ext <- function(x, y, Y, V, dimens, p_ext, representation = representation, distance = distance, alpha = alpha, N = N, B = B) {

  m <- as.numeric(dim(Y)[2])
  pvalue <- rep(-1, m)

  for(j in 1:m) {

    xx <- sapply(x, select_nodes_ext, zone = Y[,j], vertex = V, simplify = FALSE)
    yy <- sapply(y, select_nodes_ext, zone = Y[,j], vertex = V, simplify = FALSE)

    set.seed(1234)
    pvalue[j] <- test_twosample(xx, yy, statistic = c("sot", "lot"), representation = representation, distance = distance, alpha = alpha, B = B)[[2]]
    inter <- combn(sort(Y[,j]),2)
    mm <- dim(inter)[2]
    for (k in 1:mm) {
      p_ext[inter[1,k], inter[2,k]] <- max(as.numeric(p_ext[inter[1,k], inter[2,k]]), as.numeric(pvalue[j]))
    }
  }

  #possibly_significant_areas <- t(unique.matrix(t(do.call("cbind", Z))))

  Z <- unique(as.vector(Y[,which(pvalue<alpha)]))

  if (length(Z) == 0) {
    exit_ext <- TRUE
    return(list(exit_ext, p_ext))
  }
  else
  {
    exit_ext <- FALSE
    nonSign <- as.matrix(Y[,which(pvalue>=alpha)])
    if (dim(nonSign)[2] > 0) {
      #c <- cbind(combn(Z,dimens),combn(Y[,which(pvalue>=0.05)],dimens))
      #c <- c[,!duplicated(c, MARGIN = 2)]
      #aggiungere ciclio?

      previously_excluded <- which(prodlim::row.match(data.frame(t(combn(1:N, (dimens+1)))), data.frame(t(Y)), nomatch = 0)==0)
      if (length(previously_excluded)>0) {
        to_be_excluded <- subgroups(combn(1:N, (dimens+1))[, previously_excluded], dimens)
        eliminate1 <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(to_be_excluded)), nomatch = 0)>0)
        eliminate2 <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(subgroups(Y[,which(pvalue>=alpha)], dimens))), nomatch = 0)>0)
        eliminate <- unique(sort(c(eliminate1, eliminate2)))
      }
      else
        eliminate <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(subgroups(Y[,which(pvalue>=alpha)], dimens))), nomatch = 0)>0)



      if (length(eliminate) == dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
        exit_ext <- TRUE
        return(list(exit_ext, p_ext))
      }
      if (length(eliminate) != dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
        c <- subgroups(Y[,which(pvalue < alpha)], dimens)[,-eliminate]
      }

    }

    if (dim(nonSign)[2] == 0) {
      previously_excluded <- which(prodlim::row.match(data.frame(t(combn(1:N, (dimens+1)))), data.frame(t(Y)), nomatch = 0)==0)
      if (length(previously_excluded)>0) {
        to_be_excluded <- subgroups(combn(1:N, (dimens+1))[, previously_excluded], dimens)
        eliminate <- which(prodlim::row.match(data.frame(t(subgroups(Y[,which(pvalue<alpha)], dimens))), data.frame(t(to_be_excluded)), nomatch = 0)>0)

        if (length(eliminate) == dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
          exit_ext <- TRUE
          return(list(exit_ext, p_ext))
        }
        if (length(eliminate) != dim(subgroups(Y[,which(pvalue<alpha)], dimens))[2]) {
          c <- as.matrix(subgroups(Y[,which(pvalue < alpha)], dimens)[,-eliminate])
        }
      }
      else
        c <- subgroups(Y[,which(pvalue<alpha)], dimens)
    }


    possibly_significant_areas <- c
    #possibly_significant_areas <- combn(Z, dimens)
    results <- list(exit_ext, possibly_significant_areas, p_ext)
    return(results)
  }


}
