# Extracted from test-aev-residual.R:298

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
mixed_aev <- function() {
  aev <- create_aev(
    A = c(1100,  970,    0, 1500),
    E = c(1000, 1000,   50, 1000),
    V = c(2500, 1387,  125,  100)
  )
  names(aev) <- c("ok", "small", "zero A", "loud")
  aev
}
layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
layer_index <- function(p, geom) which(layer_geoms(p) == geom)[[1]]
residual_layer <- function(p) p$layers[[layer_index(p, "GeomAevResidual")]]
residual_children <- function(p, width = 4) {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(aev_residual$height_px / 96, "inches")
  ))

  # `[[2]]`: `layer_grob()` returns one grob per panel and this layer lives in
  # the second. Asking for the first gives a `zeroGrob` and no children at all.
  grob <- ggplot2::layer_grob(p, layer_index(p, "GeomAevResidual"))[[2]]
  grid::makeContent(grob)$children
}
residual_grob <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::layer_grob(p, layer_index(p, "GeomAevResidual"))[[2]]
}
panel_params <- function(p, panel) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplot_build(p)$layout$panel_params[[panel]]
}

# test -------------------------------------------------------------------------
just_under <- create_aev(A = 1000, E = 1000, V = 1600)
just_over <- create_aev(A = 1030, E = 1000, V = 1600)
expect_lt(abs(just_under$deviance_residual), 0.25)
expect_gt(abs(just_over$deviance_residual), 0.25)
expect_s3_class(residual_children(autoplot(just_under))[[1]], "null")
