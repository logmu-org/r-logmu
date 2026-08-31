# Extracted from test-aev-mark-scale.R:437

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
two_records <- function() {
  create_aev(A = c(1100, 970), E = c(1000, 1000), V = c(2500, 1600))
}
marker_widths <- function(size = NULL) {

  plot <- autoplot(two_records())

  if (!is.null(size)) {
    plot <- plot + ggplot2::theme(text = ggplot2::element_text(size = size))
  }

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(plot))
  marker <- Filter(function(g) inherits(g, "aev_interval"), grobs_in(forced))[[1]]

  vapply(marker$children, function(child) child$gp$lwd %||% NA_real_, numeric(1))
}
residual_depth <- function(size = NULL) {

  plot <- autoplot(two_records())

  if (!is.null(size)) {
    plot <- plot + ggplot2::theme(text = ggplot2::element_text(size = size))
  }

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  table <- ggplot2::ggplotGrob(plot)
  lowest <- max(aev_panel_rows(table))

  grid::convertHeight(table$heights[lowest], "inches", valueOnly = TRUE) * 96
}

# test -------------------------------------------------------------------------
expect_true(isTRUE(ggplot2::calc_element(aev_signature, theme_aev())))
