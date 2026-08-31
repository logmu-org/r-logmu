# Extracted from test-aev-mark-scale.R:469

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
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)
aev <- create_aev(A = c(1100, 970), E = rep(1000, 2), V = c(2500, 1600))
quiet <- list(
    "the chart as drawn" = autoplot(aev),
    "text resized" = autoplot(aev) + ggplot2::theme(text = ggplot2::element_text(size = 18)),
    "a caption added" = autoplot(aev) + ggplot2::labs(caption = "Overdispersion 2.0"),
    "the legend hidden" = autoplot(aev) + ggplot2::theme(legend.position = "none"),
    "the theme put back" = autoplot(aev) + ggplot2::theme_minimal() + theme_aev(),
    "resized properly" = autoplot(aev) + theme_aev(base_size = 13)
  )
