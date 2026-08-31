# Extracted from test-aev-wrap.R:111

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
with_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}
axis_text <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(p)
  cell <- which(table$layout$name == "axis-b-1-1")
  aev_text_within(table$grobs[[cell]])
}
drawn_labels <- function(p, width = 8, height = 7) {

  grDevices::pdf(NULL, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(p))

  wrappers <- Filter(
    function(grob) inherits(grob, "forcedgrob") && inherits(grob, "aev_axis_labels"),
    grobs_in(forced)
  )

  if (length(wrappers) == 0L) return(NULL)

  as.character(aev_text_within(wrappers[[1]])$label)
}

# test -------------------------------------------------------------------------
expect_identical(aev_label_lines(character(0)), 1L)
