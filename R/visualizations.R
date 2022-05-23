.project_data <- function(x, method = "mds") {
  if (sum(x) < sqrt(.Machine$double.eps)) {
    n <- attr(x, "Size")
    return(tibble(V1 = rep(0, n), V2 = rep(0, n)))
  }
  x <- stats::cmdscale(x, add = TRUE)$points
  colnames(x) <-  c("V1", "V2")
  tibble::as_tibble(x)
}

nvd_data <- function(x, y, representation = NULL, distance = NULL) {
  if (is.null(representation)) {
    rchoices <- c("adjacency", "laplacian", "modularity")
  } else {
    rchoices <- representation
  }

  if (is.null(distance)) {
    dchoices <- c("hamming", "frobenius", "spectral", "root-euclidean")
  } else {
    dchoices <- distance
  }

  tidyr::crossing(Representation = rchoices, Distance = dchoices) %>%
    dplyr::mutate(
      mds = purrr::map2(
        .x = .data$Representation,
        .y = .data$Distance,
        .f = dist_nvd,
        x = x, y = y
      ) %>%
        purrr::map(.project_data) %>%
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
    )
}

#' MDS Visualization of Network Distributions
#'
#' This function generates 2-dimensional plots of samples of networks via
#' multi-dimensional scaling using all representations and distances included in
#' the package.
#'
#' @param x A \code{\link{nvd}} object.
#' @param y A \code{\link{nvd}} object.
#' @param representation A character vector specifying the chosen representation(s),
#'   among: `"adjacency"`, `"laplacian"`,
#'   `"modularity"`. If `NULL`, each representation is chosen. Defaults to `NULL`.
#' @param distance A character vector specifying the chosen distance(s),
#'   among: `"frobenius"`, `"hamming"`,
#'   `"spectral"`, `"root-euclidean"` and `"match-frobenius"`. If `NULL`, the first four distances are chosen. Defaults to `NULL`.
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
#' @name nvd-plot
#'
#' @examples
#' gnp_params <- list(p = 1/3)
#' k_regular_params <- list(k = 8L)
#' x <- nvd(model = "gnp", n = 10L, model_params = gnp_params)
#' y <- nvd(model = "k_regular", n = 10L, model_params = k_regular_params)
#' ggplot2::autoplot(x, y)
#' plot(x, y)
NULL

#' @export
#' @rdname nvd-plot
#' @importFrom ggplot2 autoplot
autoplot.nvd <- function(x, y, representation = NULL, distance = NULL, ...) {
  nvd_data(x, y, representation = representation, distance = distance) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data$V1,
      y = .data$V2,
      color = .data$Label
    )) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$Representation, .data$Distance),
      scales = "free",
      labeller = "label_both"
    ) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      x = "First Principal Coordinate",
      y = "Second Principal Coordinate"
    )
}

#' @export
#' @rdname nvd-plot
#' @importFrom graphics plot
plot.nvd <- function(x, y, representation = NULL, distance = NULL, ...) {
  print(autoplot(x, y, representation = representation, distance = distance, ...))
}
