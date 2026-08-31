# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# THE CHART READS ggplot2'S OWN LAYOUT BY NAME, and those names are ggplot2's to
# change. Tim ruled on 2026-08-26 that every structural lookup must fail LOUDLY:
# each one used to miss quietly, and a quiet miss draws a chart that looks right
# and is not.
#
# Every test here renames a cell to simulate the ggplot2 release that moves it,
# then checks the chart refuses rather than carries on. Renaming is the closest
# thing to that release we can write down.

guarded_aev <- function() {
  aev <- create_aev(A = c(1100, 970, 1050), E = rep(1000, 3), V = c(2500, 1600, 900))
  names(aev) <- c("a", "b", "c")
  aev
}

grouped_guarded_aev <- function() {
  aev <- guarded_aev()
  group_names(aev) <- c("one", "one", "two")
  aev
}

built <- function(aev = guarded_aev(), ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplotGrob(autoplot(aev, ...))
}

# The release that renames a cell, as far as a test can stage it.
renaming <- function(table, from, to) {
  table$layout$name <- sub(paste0("^", from), to, table$layout$name)
  table
}

#### The helper itself ####

test_that("the guard passes what it is given straight back", {

  expect_identical(aev_require_layout(c(3L, 7L), "something"), invisible(c(3L, 7L)))
})

test_that("the guard names what was missing and which ggplot2 is installed", {

  expect_error(aev_require_layout(integer(0), "the widget"), "the widget")
  expect_error(aev_require_layout(integer(0), "the widget"), "could not find")
  expect_error(
    aev_require_layout(integer(0), "the widget"),
    as.character(utils::packageVersion("ggplot2")),
    fixed = TRUE
  )
})

test_that("a row number that came back NA counts as missing", {

  expect_error(aev_require_layout(NA_integer_, "a row"), "could not find")
})

test_that("a theme full of legitimate NAs is not missing", {

  # `is.atomic` guards the NA test for exactly this: the chart's own theme has
  # `fill = NA` on the legend box, and testing a list for NAs would refuse every
  # chart there is.
  expect_no_error(aev_require_layout(theme_aev(), "the theme"))
  # Read the property directly: a ggplot2 4.0 element is an S7 object, so
  # `unlist()` does not reach into it and `anyNA()` warns about the attempt.
  expect_true(is.na(theme_aev()$legend.background$fill))
})

#### Each structural lookup ####

test_that("a renamed panel cell is an error, not a wrong chart", {

  expect_error(
    aev_panel_rows(renaming(built(), "panel", "pane1")),
    "any panel"
  )
})

test_that("a missing lower x axis is an error, not record names printed twice", {

  # Quietly, this left the axis under the residual panel in place, so the record
  # names appeared under both panels.
  expect_error(
    aev_drop_lowest_x_axis(renaming(built(), "axis-b", "axis-B")),
    "x axis below the residual panel"
  )
})

test_that("a missing panel title is an error, not a chart without one", {

  expect_error(
    aev_require_panel_titles(renaming(built(), "strip-l", "strip-L")),
    "A/E panel's title"
  )
})

test_that("missing residual axis numbers are an error, not two scales the same size", {

  expect_error(
    aev_shrink_lower_panel(renaming(built(), "axis-l", "axis-L")),
    "residual panel's axis numbers"
  )
})

test_that("group strips with nowhere to go are an error, not brackets left off", {

  # The strips are found but the row to put them in is not, which is the case
  # that used to `return(table)` and drop the brackets silently.
  expect_error(
    aev_move_group_strips(renaming(built(grouped_guarded_aev()), "axis-b", "axis-B")),
    "row for the group names"
  )
})

#### And the chart still draws when nothing is renamed ####

test_that("none of the guards fires on an ordinary chart", {

  for (aev in list(guarded_aev(), grouped_guarded_aev())) {
    for (residuals in c(TRUE, FALSE)) {
      expect_no_error(built(aev, residuals = residuals))
    }
  }
})
