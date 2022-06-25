#' Power Simulations for Permutation Tests
#'
#' This function provides a Monte-Carlo estimate of the power of the permutation
#' tests proposed in this package.
#'
#' Currently, six scenarios of pairs of populations are implemented. Scenario 0
#' allows to make sure that all our permutation tests are exact.
#'
#' Match-Frobenius distance is the Frobenius distance considering networks in the Graph Space.
#' Match-Frobenius distance is computed via the graph matching algorithm with indefinite relaxation (via Frank-Wolfe), using \code{\link[iGraphMatch]{iGraphMatch}} \code{\link[iGraphMatch]{gm}} function.
#' Match-Frobenius distance can be used only with adjacency matrix representation.
#'
#' To perform it in a parallel fashion, use \code{future::plan(future::multisession)} before the call.
#'
#' @param model1 A string specifying the model to be used for generating the
#'   first sample. Choices are `"sbm"`, `"k_regular"`, `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"`. Defaults to `"gnp"`.
#' @param model2 A string specifying the model to be used for generating the
#'   second sample. Choices are `"sbm"`, `"k_regular"`, `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"`. Defaults to `"k_regular"`.
#' @param n1 The size of the first sample. Defaults to `20L`.
#' @param n2 The size of the second sample. Defaults to `20L`.
#' @param num_vertices1 The number of nodes in the generated graphs of the first sample. Defaults to
#'   `25L`.
#' @param num_vertices2 The number of nodes in the generated graphs of the second sample. Defaults to
#'   `25L`.
#' @param model1_params A named list setting the parameters of the first chosen
#'   model. Defaults to `list(p = 1/3)`.
#' @param model2_params A named list setting the parameters of the second chosen
#'   model. Defaults to `list(k = 8L)`.
#' @inheritParams test2_global
#' @param alpha Significance level for hypothesis testing. Defaults to `0.05`.
#' @param R Number of Monte-Carlo trials used to estimate the power. Defaults to
#'   `1000L`.
#' @param seed An integer specifying the random generator seed. Defaults to
#'   `1234`.
#' @param  rand_num_vertices1 A string specifying whether the number of vertices of the networks in the first sample should be random. Choices are `"poisson"`, `"uniform"`. In particular, if `"poisson"` N_1,...,N_sample_size iid Poisson(num_vertices1), if `"uniform"` N_1,...,N_sample_size iid Uniform(floor(num_vertices1*(1-rate)), ceiling(num_vertices1*(1+rate))). It is compatible with `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"` models. Defaults to `NULL`.
#' @param  rand_num_vertices2 A string specifying whether the number of vertices of the networks in the first sample should be random. Choices are `"poisson"`, `"uniform"`. In particular, if `"poisson"` N_1,...,N_sample_size iid Poisson(num_vertices2), if `"uniform"` N_1,...,N_sample_size iid Uniform(floor(num_vertices2*(1-rate)), ceiling(num_vertices2*(1+rate))). It is compatible with `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"` models. Defaults to `NULL`.
#' @param rate1 A rate. See `rand_num_vertices1` for further information. Defaults to `0.25`.
#' @param rate2 A rate. See `rand_num_vertices2` for further information. Defaults to `0.25`.
#'
#' @return A numeric value estimating the power of the test.
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' k_regular_params <- list(k = 8L)
#' power2(
#'   model1_params = gnp_params,
#'   model2_params = k_regular_params,
#'   R = 10,
#'   B = 100,
#'   seed = 1234
#' )
power2 <- function(model1 = "gnp", model2 = "k_regular",
                   n1 = 20L, n2 = 20L,
                   num_vertices1 = 25L,
                   num_vertices2 = 25L,
                   model1_params = NULL, model2_params = NULL,
                   representation = "adjacency",
                   distance = "frobenius",
                   stats = c("flipr:t_ip", "flipr:f_ip"),
                   B = 1000L,
                   start = "barycenter",
                   iteration = 20L,
                   alpha = 0.05,
                   test = "exact",
                   k = 5L,
                   R = 1000L,
                   seed = 1234,
                   rand_num_vertices1 = NULL,
                   rand_num_vertices2 = NULL,
                   rate1 = 0.25,
                   rate2 = 0.25) {
  if (!is.null(seed))
    withr::local_seed(seed)


  if (num_vertices1 != num_vertices2 | !is.null(rand_num_vertices1) | !is.null(rand_num_vertices2)){

    pvalues_function <- function(xs) {

      test_function <- function(x){
        x = nvd(
          model = model1,
          n = n1,
          num_vertices = num_vertices1,
          model_params = model1_params,
          seed = NULL,
          rand_num_vertices = rand_num_vertices1,
          rate = rate1
        )
        y = nvd(
          model = model2,
          n = n2,
          num_vertices = num_vertices2,
          model_params = model2_params,
          seed = NULL,
          rand_num_vertices = rand_num_vertices2,
          rate = rate2
        )
        maximum_num_vertices <- max(c(nvd_num_vertices(x)$max, nvd_num_vertices(y)$max))
        x <- add_null_nodes(x = x, num_vertices = maximum_num_vertices, directed = FALSE, weighted = ifelse(model1 == "binomial" | model1 == "poisson", TRUE, FALSE))
        y <- add_null_nodes(x = y, num_vertices = maximum_num_vertices, directed = FALSE, weighted = ifelse(model2 == "binomial" | model2 == "poisson", TRUE, FALSE))

        test2_global(
          x = x,
          y = y,
          representation = representation,
          distance = distance,
          stats = stats,
          B = B,
          test = test,
          k = k,
          seed = 1234,
          start = start,
          iteration = iteration
        )$pvalue
      }
      y <- furrr::future_map_dbl(1:xs, test_function, .options = furrr::furrr_options(seed = TRUE))
    }
  } else {
    pvalues_function <- function(xs) {
      test_function <- function(x){
        x = nvd(
          model = model1,
          n = n1,
          num_vertices = num_vertices1,
          model_params = model1_params,
          seed = NULL,
          rand_num_vertices = rand_num_vertices1
        )
        y = nvd(
          model = model2,
          n = n2,
          num_vertices = num_vertices2,
          model_params = model2_params,
          seed = NULL,
          rand_num_vertices = rand_num_vertices2
        )

        test2_global(
          x = x,
          y = y,
          representation = representation,
          distance = distance,
          stats = stats,
          B = B,
          test = test,
          k = k,
          seed = 1234,
          start = start,
          iteration = iteration
        )$pvalue
      }
      y <- furrr::future_map_dbl(1:xs, test_function, .options = furrr::furrr_options(seed = TRUE))
    }
  }

  pvalues <- pvalues_function(R)

  mean(pvalues <= alpha)
}
