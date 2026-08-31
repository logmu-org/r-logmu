# Extracted from test-aev-chevron.R:316

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
mixed_aev <- function() {
  aev <- create_aev(
    A = c(1100, 2000, 1900,   0,  400, 0, NaN),
    E = c(1000, 1000, 1000,  50, 1000, 0, NaN),
    V = c(2500, 2000, 90000, 125,  400, 0, NaN)
  )
  names(aev) <- c(
    "ordinary", "off the top", "bar reaches in", "no deaths",
    "off the bottom", "nobody", "unknown"
  )
  aev
}
chevron_layer <- function(p) {
  index <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) ==
                   "GeomAevChevron")
  p$layers[[index[[1]]]]
}
chevron_data <- function(p) chevron_layer(p)$data
marker_data <- function(p) {
  index <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) ==
                   "GeomAevInterval")
  p$layers[[index[[1]]]]$data
}
chevron_children <- function(p, width = 6, height = 4) {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(height, "inches")
  ))

  index <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) ==
                   "GeomAevChevron")
  grob <- ggplot2::layer_grob(p, index[[1]])[[1]]
  grid::makeContent(grob)$children
}

# test -------------------------------------------------------------------------
render <- function(p) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
