# Extracted from test-aev-group-bracket.R:79

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
grouped_aev <- function() {
  aev <- create_aev(
    A = c(1100, 970, 1050, 880, 1200),
    E = rep(1000, 5),
    V = rep(2500, 5)
  )
  names(aev) <- c("a", "b", "c", "d", "e")
  group_names(aev) <- c("one", "one", "two", "two", "two")
  aev
}
bracket_cells <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(p)
  table$layout$name[startsWith(table$layout$name, "aev-group-")]
}
bracket_grob <- function(p, which = 1L) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(p)
  table$grobs[[which(table$layout$name == paste0("aev-group-", which))]]
}
bracket_children <- function(p, which = 1L, width = 2) {

  grob <- bracket_grob(p, which)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(0.2, "inches")
  ))

  grid::makeContent(grob)$children
}

# test -------------------------------------------------------------------------
ends <- bracket_children(autoplot(grouped_aev()))[[1]]
