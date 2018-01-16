#' MDS Visualization of Network Distributions
#'
#' This function generates 2-dimensional plots of samples of networks via
#' multi-dimensional scaling using all represenations and distances included in
#' the package.
#'
#' @param x A \code{\link{nvd}} object.
#' @param y A \code{\link{nvd}} object.
#' @param ... Extra arguments to be passed to the plot function.
#'
#' @return Invisibly returns the dataset computed to generate the plot.
#' @importFrom dplyr %>%
#' @export
#'
#' @examples
#' x <- nvd("smallworld", 10)
#' y <- nvd("pa", 10)
#' plot(x, y)
plot.nvd <- function(x, y, ...) {
  rchoices <- c("adjacency", "laplacian", "modularity")
  dchoices <- c("hamming", "frobenius", "spectral", "root-euclidean")
  tidyr::crossing(Representation = rchoices, Distance = dchoices) %>%
    dplyr::mutate(
      dist_matrix = purrr::map2(Representation, Distance, dist_nvd, x = x, y = y),
      mds = purrr::map(dist_matrix, stats::cmdscale) %>%
        purrr::map(tibble::as_tibble) %>%
        purrr::map(dplyr::mutate, Label = c(rep("1", length(x)), rep("2", length(y))))
    ) %>%
    dplyr::select(-dist_matrix) %>%
    tidyr::unnest() %>%
    dplyr::mutate(
      Representation = Representation %>%
        forcats::as_factor() %>%
        forcats::fct_relabel(capitalize),
      Distance = Distance %>%
        forcats::as_factor() %>%
        forcats::fct_relabel(capitalize),
      Label = forcats::as_factor(Label)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = V1, y = V2, color = Label)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~ Representation + Distance, scales = "free", labeller = "label_both") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("First Principal Coordinate") +
    ggplot2::ylab("Second Principal Coordinate")
}
