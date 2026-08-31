# Extracted from test-aev-plot.R:167

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
prototype_aev <- function() {
  ratio <- c(1.218, 1.075, 1.001, 0.970)
  confidence <- c(0.058, 0.015, 0.101, 0.073)
  E <- rep(1000, 4)
  aev <- create_aev(A = ratio * E, E = E, V = (confidence / 1.959963984540 * E)^2)
  names(aev) <- c("[65, 75)", "[75, 85)", "[85, 95)", "[1930, 1940)")
  aev
}
awkward_aev <- function() {
  aev <- create_aev_unchecked(
    A = c(1100, 0, NaN,   0, 3000),
    E = c(1000, 0, NaN,  50, 1000),
    V = c(2500, 0, NaN, 125, 2000)
  )
  names(aev) <- c("ok", "empty", "missing", "zero A", "too high")
  aev
}
panel_params <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplot_build(p)$layout$panel_params[[1]]
}
layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
layer_index <- function(p, geom) which(layer_geoms(p) == geom)[[1]]
marker_layer <- function(p) p$layers[[layer_index(p, "GeomAevInterval")]]
interval_grob <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::layer_grob(p, layer_index(p, "GeomAevInterval"))[[1]]
}
measured <- function(p, width = 4, height = 3) {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(height, "inches")
  ))

  # Not `interval_grob()`: that opens its own device, and the viewport just
  # pushed belongs to this one.
  grob <- ggplot2::layer_grob(p, layer_index(p, "GeomAevInterval"))[[1]]
  children <- grid::makeContent(grob)$children

  in_pixels <- function(length) {
    grid::convertWidth(length, "inches", valueOnly = TRUE) * 96
  }

  list(
    grob = grob,
    children = children,
    radius_npc = grid::convertHeight(
      grid::unit(aev_geometry$marker_radius_px / 96, "inches"), "npc",
      valueOnly = TRUE
    ),
    radius_px = in_pixels(children[[1]]$r),
    cap_px = in_pixels(children[[4]]$x1[1] - children[[4]]$x0[1])
  )
}
as_npc <- function(values) (values - aev_log_ratio_limits[1]) / diff(aev_log_ratio_limits)

# test -------------------------------------------------------------------------
aev <- create_aev(A = 1500, E = 1000, V = 100)
expect_identical(as.character(aev_status(aev)), "ok")
