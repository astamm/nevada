.project_data <- function(x, method = "mds") {
  if (sum(x) < sqrt(.Machine$double.eps)) {
    n <- attr(x, "Size")
    return(tibble(V1 = rep(0, n), V2 = rep(0, n)))
  }
  x <- switch (method,
    mds = stats::cmdscale(x, add = TRUE)$points,
    tsne = suppressMessages(tsne::tsne(x)),
    umap = {
      custom_settings <- umap::umap.defaults
      custom_settings$input <- "dist"
      umap::umap(as.matrix(x), config = custom_settings)$layout
    }
  )
  colnames(x) <-  c("V1", "V2")
  tibble::as_tibble(x)
}

nvd_data <- function(x, memberships, method = "mds") {
  rchoices <- c("adjacency", "laplacian", "modularity")
  dchoices <- c("hamming", "frobenius", "spectral", "root-euclidean")
  tidyr::crossing(Representation = rchoices, Distance = dchoices) |>
    dplyr::mutate(
      mds = purrr::map2(.data$Representation, .data$Distance, \(.repr, .dist) {
        dist_nvd(x, representation = .repr, distance = .dist)
      }) |>
        purrr::map(.project_data, method = method) |>
        purrr::map(dplyr::mutate, Label = memberships)
    ) |>
    tidyr::unnest(cols = .data$mds) |>
    dplyr::mutate(
      Representation = .data$Representation |>
        forcats::as_factor() |>
        forcats::fct_relabel(capitalize),
      Distance = .data$Distance |>
        forcats::as_factor() |>
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
#' @param object,x A list containing two samples of network-valued data stored
#'   as objects of class [`nvd`].
#' @param memberships An integer vector specifying the membership of each
#'   network to a specific sample. Defaults to `rep(1, length(nvd))` which
#'   assumes that all networks in the input [`nvd`] object belong to a single
#'   group.
#' @param method A string specifying which dimensionality reduction method to
#'   use for projecting the samples into the cartesian plane. Choices are
#'   `"mds"`, `"tsne"` or `"umap"`. Defaults to `"mds"`.
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
#' gnp_params <- list(n = 24L, p = 1/3)
#' degree_params <- list(out_degree = rep(2, 24L), method = "configuration")
#' x <- nvd(sample_size = 10L, model = "gnp", !!!gnp_params)
#' y <- nvd(sample_size = 10L, model = "degree", !!!degree_params)
#' mb <- c(rep(1, length(x)), rep(2, length(y)))
#' z <- as_nvd(c(x, y))
#' ggplot2::autoplot(z, memberships = mb)
#' plot(z, memberships = mb)
NULL

#' @export
#' @rdname nvd-plot
#' @importFrom ggplot2 autoplot
autoplot.nvd <- function(object,
                         memberships = rep(1, length(object)),
                         method = "mds",
                         ...) {
  nvd_data(x = object, memberships = memberships, method = method) |>
    ggplot2::ggplot(ggplot2::aes(
      x = .data$V1,
      y = .data$V2,
      color = .data$Label
    )) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$Representation),
      cols = ggplot2::vars(.data$Distance),
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
plot.nvd <- function(x, method = "mds", ...) {
  print(autoplot(x, method = method, ...))
}
