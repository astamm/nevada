#' Power Simulations for Permutation Tests
#'
#' This function provides a Monte-Carlo estimate of the power of the permutation
#' tests proposed in this package.
#'
#' Currently, six scenarios of pairs of populations are implemented. Scenario 0
#' allows to make sure that all our permutation tests are exact.
#'
#' @param sample_size1 An integer specifying the size of the first sample.
#' @param model1 A string specifying the model to be used for sampling networks
#'   in the first sample. All `tidygraph::play_` functions are supported. The
#'   model name corresponds to the name of the function without the `play_`
#'   prefix.
#' @param params1 A list specifying the parameters to be passed to the model
#'  function that will generate the first sample.
#' @param sample_size2 An integer specifying the size of the second sample.
#' @param model2 A string specifying the model to be used for sampling networks
#'  in the second sample. All `tidygraph::play_` functions are supported. The
#'  model name corresponds to the name of the function without the `play_`
#'  prefix.
#' @param params2 A list specifying the parameters to be passed to the model
#' @inheritParams test2_global
#' @param alpha Significance level for hypothesis testing. Defaults to `0.05`.
#' @param R Number of Monte-Carlo trials used to estimate the power. Defaults to
#'   `1000L`.
#'
#' @return A numeric value estimating the power of the test.
#' @export
#'
#' @examples
#' gnp_params <- list(n = 24L, p = 1/3)
#' degree_params <- list(out_degree = rep(2, 24L), method = "configuration")
#' power2(
#'   sample_size1 = 10L, model1 = "gnp", params1 = gnp_params,
#'   sample_size2 = 10L, model2 = "degree", params2 = degree_params,
#'   R = 10L,
#'   B = 100L
#' )
power2 <- function(sample_size1, model1, params1,
                   sample_size2, model2, params2,
                   representation = "adjacency",
                   distance = "frobenius",
                   stats = c("flipr:t_ip", "flipr:f_ip"),
                   B = 1000L,
                   alpha = 0.05,
                   test = "exact",
                   k = 5L,
                   R = 1000L) {
  pvalues <- replicate(
    R,
    test2_global(
      x = nvd(sample_size = sample_size1, model = model1, !!!params1),
      y = nvd(sample_size = sample_size2, model = model2, !!!params2),
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
