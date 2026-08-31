# Extracted from test-aev-legend.R:199

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
legend_aev <- function() {
  aev <- create_aev(A = c(1100, 970, 1010), E = rep(1000, 3), V = c(2500, 1387, 900))
  names(aev) <- c("one", "two", "three")
  aev
}
layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
built <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplotGrob(p)
}
key_grob <- function(entry) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(aev_legend$key_width_px / 96, "inches"),
    height = grid::unit(aev_legend$key_height_px / 96, "inches")
  ))
  aev_draw_key(data.frame(shape = entry), NULL, NULL)
}

# test -------------------------------------------------------------------------
expect_gt(aev_legend$key_width_px, aev_geometry$cap_width_px)
