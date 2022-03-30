#' Power Simulations for Permutation Tests
#'
#' This function provides a Monte-Carlo estimate of the power of the permutation
#' tests proposed in this package.
#'
#' Currently, six scenarios of pairs of populations are implemented. Scenario 0
#' allows to make sure that all our permutation tests are exact.
#'
#' @param model1 A string specifying the model to be used for generating the
#'   first sample. Choices are `"sbm"`, `"k_regular"`, `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"`. Defaults to `"gnp"`.
#' @param model2 A string specifying the model to be used for generating the
#'   second sample. Choices are `"sbm"`, `"k_regular"`, `"gnp"`, `"smallworld"`,
#'   `"pa"`, `"poisson"` and `"binomial"`. Defaults to `"k_regular"`.
#' @param n1 The size of the first sample. Defaults to `20L`.
#' @param n2 The size of the second sample. Defaults to `20L`.
#' @param num_vertices The number of nodes in the generated graphs. Defaults to
#'   `25L`.
#' @param model1_params A named list setting the parameters of the first chosen
#'   model. Defaults to `list(p = 1/3)`.
#' @param model2_params A named list setting the parameters of the second chosen
#'   model. Defaults to `list(k = 8L)`.
#' @inheritParams test2_global
#' @param alpha Significance level for hypothesis testing. Defaults to `0.05`.
#' @param R Number of Monte-Carlo trials used to estimate the power. Defaults to
#'   `1000L`.
#' @param seed An integer specifying the seed to start randomness from. Defaults
#'   to using clock time.
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
                   num_vertices = 25L,
                   model1_params = NULL, model2_params = NULL,
                   representation = "adjacency",
                   distance = "frobenius",
                   stats = c("flipr:t_ip", "flipr:f_ip"),
                   B = 1000L,
                   alpha = 0.05,
                   test = "exact",
                   k = 5L,
                   R = 1000L,
                   seed = NULL) {

  withr::local_seed(seed)

  pvalues <- replicate(
    R,
    test2_global(
      x = nvd(
        model = model1,
        n = n1,
        num_vertices = num_vertices,
        model_params = model1_params
      ),
      y = nvd(
        model = model2,
        n = n2,
        num_vertices = num_vertices,
        model_params = model2_params
      ),
      representation = representation,
      distance = distance,
      stats = stats,
      B = B,
      test = test,
      k = k
    )$pvalue
  )

  mean(pvalues <= alpha)
}
