# Extracted from test-aev-notes.R:33

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
noted_aev <- function() {
  aev <- create_aev(A = c(1100, 970), E = c(1000, 1000), V = c(2500, 1600))
  names(aev) <- c("a", "b")
  aev
}
caption_of <- function(plot) plot$labels$caption
drawn_caption <- function(plot) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(plot)
  cell <- which(table$layout$name == "caption")
  if (length(cell) == 0L) return(character(0))
  labels_in(table$grobs[[cell]])
}

# test -------------------------------------------------------------------------
expect_null(aev_notes_text(NULL))
