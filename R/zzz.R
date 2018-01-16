.onLoad <- function(libname = find.package("grattan"), pkgname = "grattan"){

  # CRAN Note avoidance
  if(getRversion() >= "2.15.1")
    utils::globalVariables(
      # for plot.nvd
      c("Representation", "Distance", "dist_matrix", "Label", "V1", "V2")
    )
  invisible()

}

.onUnload <- function (libpath) {
  library.dynam.unload("nevada", libpath)
}
