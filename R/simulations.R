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
#' @param seed An integer specifying the random generator seed. Defaults to
#'   `1234`.
#' @param parallel A boolean specifying whether the function should be performed in parallel. Defaults to `FALSE`.
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
                   start = "barycenter",
                   iteration = 20L,
                   alpha = 0.05,
                   test = "exact",
                   k = 5L,
                   R = 1000L,
                   seed = 1234,
                   parallel = FALSE) {
  if (!is.null(seed))
    withr::local_seed(seed)

  if(!parallel){

  pvalues <- replicate(
    R,
    test2_global(
      x = nvd(
        model = model1,
        n = n1,
        num_vertices = num_vertices,
        model_params = model1_params,
        seed = NULL
      ),
      y = nvd(
        model = model2,
        n = n2,
        num_vertices = num_vertices,
        model_params = model2_params,
        seed = NULL
      ),
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
  )
  }else{
    pvalues=numeric(R)
    n.cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n.cores)
    fun <- function(dummy){
      pvalues <- nevada::test2_global(
        x = nevada::nvd(
          model = model1,
          n = n1,
          num_vertices = num_vertices,
          model_params = model1_params,
          seed = NULL
        ),
        y = nevada::nvd(
          model = model2,
          n = n2,
          num_vertices = num_vertices,
          model_params = model2_params,
          seed = NULL
        ),
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
    parallel::clusterExport(cl=cl,list('model1','model2','n1','n2','num_vertices','model1_params','model2_params','representation','distance','stats','B','start','iteration','test','k'), envir=environment())
    #pvalues <- parallel::parLapply(pvalues, fun, cl=cl)
    pvalues <- pbapply::pbsapply(pvalues, fun, cl=cl) # with a progress bar
    parallel::stopCluster(cl)
  }

  mean(pvalues <= alpha)
}
