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
#' @return Invisibly returns the dataset computed to generate the plot.
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
      mds = purrr::map2(
        .x = Representation,
        .y = Distance,
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
      Representation = Representation %>%
        forcats::as_factor() %>%
        forcats::fct_relabel(capitalize),
      Distance = Distance %>%
        forcats::as_factor() %>%
        forcats::fct_relabel(capitalize),
      Label = forcats::as_factor(Label)
    ) %>%
    ggplot(aes(x = V1, y = V2, color = Label)) +
    geom_point() +
    theme_bw() +
    facet_wrap(
      facets = ~ Representation + Distance,
      scales = "free",
      labeller = "label_both"
    ) +
    theme(legend.position = "none") +
    labs(
      x = "First Principal Coordinate",
      y = "Second Principal Coordinate"
    )
}
