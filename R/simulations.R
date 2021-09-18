#' Power Simulations for Permutation Tests
#'
#' This function provides a Monte-Carlo estimate of the power of the permutation
#' tests proposed in this package.
#'
#' Currently, six scenarios of pairs of populations are implemented. Scenario 0
#' allows to make sure that all our permutation tests are exact.
#'
#' @param model1 A string specifying the model to be used for sample 1 among
#'   \code{"sbm"}, \code{"k_regular"}, \code{"gnp"}, \code{"smallworld"},
#'   \code{"pa"}, \code{"poisson"} and \code{"binomial"}. Default is
#'   \code{"gnp"}.
#' @param model2 A string specifying the model to be used for sample 2 among
#'   \code{"sbm"}, \code{"k_regular"}, \code{"gnp"}, \code{"smallworld"},
#'   \code{"pa"}, \code{"poisson"} and \code{"binomial"}. Default is
#'   \code{"k_regular"}.
#' @param n1 The size of sample 1 (default: 20L).
#' @param n2 The size of sample 2 (default: 20L).
#' @param size1 The number of trials for the binomial distribution of sample 1
#'   (default: \code{NULL}). Must be specified if \code{model == "binomial"}.
#' @param prob1 The probability of success of each trial for the binomial
#'   distribution of sample 1 (default: \code{NULL}). Must be specified if
#'   \code{model == "binomial"}.
#' @param size2 The number of trials for the binomial distribution of sample 2
#'   (default: \code{NULL}). Must be specified if \code{model == "binomial"}.
#' @param prob2 The probability of success of each trial for the binomial
#'   distribution of sample 2 (default: \code{NULL}). Must be specified if
#'   \code{model == "binomial"}.
#' @param lambda1 The mean of the Poisson distribution of sample 1 (default:
#'   \code{NULL}). Must be specified if \code{model == "poisson"}.
#' @param lambda2 The mean of the Poisson distribution of sample 2 (default:
#'   \code{NULL}). Must be specified if \code{model == "poisson"}.
#' @param pref.matrix1 A matrix giving the Bernoulli rates for the SBM generator
#'   of sample 1 (default: \code{NULL}). Must be specified if \code{model ==
#'   "sbm"}.
#' @param pref.matrix2 A matrix giving the Bernoulli rates for the SBM generator
#'   of sample 2 (default: \code{NULL}). Must be specified if \code{model ==
#'   "sbm"}.
#' @inheritParams test2_global
#' @param alpha Significance level for hypothesis testing (default:
#'   \code{0.05}).
#' @param R The number of Monte-Carlo runs used to estimate the power (default:
#'   1000L).
#' @param seed An integer specifying the seed to start randomness from (default:
#'   uses clock).
#'
#' @return A numeric value estimating the power of the test.
#' @export
#'
#' @examples
#' power2(R = 10, B = 100, seed = 1234)
power2 <- function(model1 = "gnp", model2 = "k_regular",
                   n1 = 20L, n2 = 20L,
                   size1 = NULL, prob1 = NULL, size2 = NULL, prob2 = NULL,
                   lambda1 = NULL, lambda2 = NULL,
                   pref.matrix1 = NULL, pref.matrix2 = NULL,
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
      x = nvd(model1, n1, pref.matrix1, lambda1, size1, prob1),
      y = nvd(model2, n2, pref.matrix2, lambda2, size2, prob2),
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
