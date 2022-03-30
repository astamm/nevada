#' MDS Visualization of Network Distributions
#'
#' This function generates 2-dimensional plots of samples of networks via
#' multi-dimensional scaling using all representations and distances included in
#' the package.
#'
#' @param x A \code{\link{nvd}} object.
#' @param y A \code{\link{nvd}} object.
#' @param ... Extra arguments to be passed to the plot function.
#'
#' @return Invisibly returns a \code{\link[ggplot2]{ggplot}} object. In
#'   particular, the data set computed to generate the plot can be retrieved via
#'   `$data`. This is a \code{\link[tibble]{tibble}} containing the following
#'   variables:
#'
#' - `V1`: the x-coordinate of each observation in the plane,
#' - `V2`: the y-coordinate of each observation in the plane,
#' - `Label`: the sample membership of each observation,
#' - `Representation`: the type of matrix representation used to manipulate each
#' observation,
#' - `Distance`: the distance used to measure how far each observation is from
#' the others.
#'
#' @export
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' k_regular_params <- list(k = 8L)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' y <- nvd(model = "k_regular", n = 10L, model_params = k_regular_params)
#' plot(x, y)
plot.nvd <- function(x, y, ...) {
  rchoices <- c("adjacency", "laplacian", "modularity")
  dchoices <- c("hamming", "frobenius", "spectral", "root-euclidean")
  tidyr::crossing(Representation = rchoices, Distance = dchoices) %>%
    dplyr::mutate(
      mds = purrr::map2(
        .x = .data$Representation,
        .y = .data$Distance,
        .f = dist_nvd,
        x = x,
        y = y
      ) %>%
        purrr::map(stats::cmdscale) %>%
        purrr::map(`colnames<-`, c("V1", "V2")) %>%
        purrr::map(tibble::as_tibble) %>%
        purrr::map(
          .f = dplyr::mutate,
          Label = c(rep("1", length(x)), rep("2", length(y)))
        )
    ) %>%
    tidyr::unnest(cols = .data$mds) %>%
    dplyr::mutate(
      Representation = .data$Representation %>%
        forcats::as_factor() %>%
        forcats::fct_relabel(capitalize),
      Distance = .data$Distance %>%
        forcats::as_factor() %>%
        forcats::fct_relabel(capitalize),
      Label = forcats::as_factor(.data$Label)
    ) %>%
    ggplot(aes(x = .data$V1, y = .data$V2, color = .data$Label)) +
    geom_point() +
    theme_bw() +
    facet_wrap(
      facets = vars(.data$Representation, .data$Distance),
      scales = "free",
      labeller = "label_both"
    ) +
    theme(legend.position = "none") +
    labs(
      x = "First Principal Coordinate",
      y = "Second Principal Coordinate"
    )
}
