# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `facet_aev()` does two things a theme cannot: it puts the shared x axis
# BETWEEN the two panels rather than under both, and it makes the residual
# panel's numbers and title smaller than the A/E panel's.
#
# BOTH LEAN ON GGPLOT2 INTERNALS -- that `axes = "all_x"` draws an x axis under
# every panel, and that the assembled table names its cells `panel`, `axis-b`,
# `axis-l` and `strip-l`. This file exists so that a ggplot2 change breaks the
# suite instead of quietly producing a chart with its labels in the wrong place.
# Every assertion below is on the assembled table, which is the thing that would
# change.

facet_aev_aev <- function() {
  aev <- create_aev(A = c(1100, 970, 1010), E = rep(1000, 3), V = c(2500, 1387, 900))
  names(aev) <- c("[65, 75)", "[75, 85)", "[85, 95)")
  aev
}

# Five records in three groups of unequal size, which is what the gaps and the
# proportional widths are for.
#
# THE GROUPS ARE IN AN ORDER ALPHABETICAL SORTING WOULD CHANGE, deliberately:
# period, age, cohort, as sample B has them. A fixture of Age then Cohort could
# not tell record order from sorted order, and a recheck showed that it did
# not -- replacing the ordering with a plain `factor()` left the suite green.
grouped_aev <- function() {
  aev <- create_aev(
    A = c(1100, 970, 1010, 1040, 990),
    E = rep(1000, 5),
    V = c(2500, 1387, 900, 1600, 1200)
  )
  names(aev) <- c(
    "[2010, 2015)", "[65, 75)", "[75, 85)", "[85, 95)", "[1940, 1960)"
  )
  group_names(aev) <- c("Period", "Age", "Age", "Age", "Cohort")
  aev
}

# The table ggplot2 assembles, which is what the facet edits.
assembled <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplotGrob(p)
}

rows_of <- function(table, prefix) {
  table$layout$t[startsWith(table$layout$name, prefix)]
}

cell <- function(table, prefix, row) {
  table$grobs[[which(startsWith(table$layout$name, prefix) & table$layout$t == row)[[1]]]]
}

#### It is still a plot ####

test_that("the facet leaves an ordinary extensible ggplot", {

  # The reason this is a facet subclass and not `ggplotGrob()` surgery after
  # the fact: surgery returns a `gtable`, which a user cannot add to.
  p <- autoplot(facet_aev_aev())

  expect_s3_class(p, "ggplot")
  expect_s3_class(p + ggplot2::labs(subtitle = "extended"), "ggplot")
  expect_s3_class(p + ggplot2::theme(plot.title = ggplot2::element_text(size = 20)), "ggplot")
})

test_that("the facet is a FacetGrid that draws an axis under every panel", {

  facet <- facet_aev(groups = FALSE)

  expect_s3_class(facet, "FacetGrid")
  expect_true(facet$params$draw_axes$x)
  expect_false(facet$params$draw_axes$y)
  expect_true(facet$params$free$y)
})

test_that("groups add a second facet dimension, and only when there are any", {

  grouped <- facet_aev(groups = TRUE)

  # The columns are what draw the gaps, and freeing the x space is what keeps a
  # record the same width whichever group it is in.
  expect_length(grouped$params$cols, 1L)
  expect_true(grouped$params$free$x)
  expect_true(grouped$params$space_free$x)

  # Without groups there are no columns and so no gaps.
  expect_length(facet_aev(groups = FALSE)$params$cols, 0L)
})

#### Where the x axis ends up ####

test_that("the category labels sit between the two panels", {

  # Slide 2's layout stack, and what both rendered samples show.
  table <- assembled(autoplot(facet_aev_aev()))
  panels <- sort(rows_of(table, "panel"))
  axes <- rows_of(table, "axis-b")

  expect_length(panels, 2L)
  expect_true(any(axes > panels[1] & axes < panels[2]))
})

test_that("nothing is left under the residual panel", {

  table <- assembled(autoplot(facet_aev_aev()))
  panels <- sort(rows_of(table, "panel"))
  below <- which(startsWith(table$layout$name, "axis-b") & table$layout$t > max(panels))

  # The cell is still in the table; it is emptied and its row closed up.
  for (i in below) {
    expect_s3_class(table$grobs[[i]], "zeroGrob")
    expect_equal(as.numeric(table$heights[table$layout$t[i]]), 0)
  }
})

test_that("the axis that stays is the one with the record names on it", {

  # Emptying the wrong one of the two would leave a chart with no labels at
  # all, which this catches and the row positions alone would not.
  table <- assembled(autoplot(facet_aev_aev()))
  panels <- sort(rows_of(table, "panel"))
  between <- table$layout$t[startsWith(table$layout$name, "axis-b")]
  between <- between[between > panels[1] & between < panels[2]][[1]]

  expect_true(all(names(facet_aev_aev()) %in% labels_in(cell(table, "axis-b", between))))
})

test_that("a chart with no residual panel keeps its axis at the bottom", {

  table <- assembled(autoplot(facet_aev_aev(), residuals = FALSE))
  panels <- rows_of(table, "panel")
  axes <- rows_of(table, "axis-b")

  expect_length(panels, 1L)
  expect_true(all(axes > panels))

  # And nothing is emptied, because there is nothing below to empty.
  for (i in which(startsWith(table$layout$name, "axis-b"))) {
    expect_false(inherits(table$grobs[[i]], "zeroGrob"))
  }
})

#### How big the lower panel's text is ####

test_that("the residual panel's numbers and title are smaller than the A/E panel's", {

  table <- assembled(autoplot(facet_aev_aev()))
  panels <- sort(rows_of(table, "panel"))

  upper_scale <- cell(table, "axis-l", panels[1])
  lower_scale <- cell(table, "axis-l", panels[2])
  upper_title <- cell(table, "strip-l", panels[1])
  lower_title <- cell(table, "strip-l", panels[2])

  # The A/E panel's numbers are left alone. ITS TITLE IS ENLARGED, from
  # 2026-08-25: Tim asked for the upper title to grow rather than the lower one
  # to shrink, so the residual panel's title stays where he approved it.
  expect_null(upper_scale$gp$cex)
  # NO LONGER SCALED AT ALL. Until 2026-08-29 the theme sized the panel titles
  # like axis numbers and this scaled the upper one up afterwards; the theme now
  # carries `YAxisTitleFontSize` directly, so the grob is left as built.
  expect_null(upper_title$gp$cex)

  expect_equal(lower_scale$gp$cex, aev_residual$scale_ratio)

  # THE TITLE IS SHRUNK TOO, from 2026-08-29. Both panels named themselves in
  # the same size until then, because nobody had found `DevYAxisTitleRatio` in
  # `AEChart.cs` -- the design shrinks the residual panel's title to 0.8 of the
  # A/E panel's, which is the distinction Tim had already asked for by eye.
  expect_equal(lower_title$gp$cex, aev_residual$title_ratio)
})

test_that("the residual panel differs from the A/E panel in two ways, not one", {

  # THE TWO RATIOS GO OPPOSITE WAYS, which is the thing to hold on to.
  #
  # The NUMBERS are the same size as the A/E panel's, which is what the design
  # says. This ratio has moved twice: 0.7 from 2026-08-21 put them at 5.04 pt,
  # under everything else on the chart; 10/9 on 2026-08-29 crowded seven numbers
  # into a fixed 84 px strip. Parity is where it settled, and a change here
  # should be a decision rather than a drift.
  expect_equal(aev_residual$scale_ratio, 1)

  # The TITLE is smaller: `DevYAxisTitleRatio` in `AEChart.cs`.
  expect_equal(aev_residual$title_ratio, 0.8)
  expect_lt(aev_residual$title_ratio, 1)
})

test_that("shrinking keeps whatever else the grob's settings said", {

  # `cex` is merged in rather than replacing the settings, and it multiplies
  # the size each piece of text already has rather than overriding it.
  grob <- grid::textGrob("x", gp = grid::gpar(col = "red"))
  scaled <- aev_scale_text(grob, 0.5)

  expect_identical(scaled$gp$col, "red")
  expect_equal(scaled$gp$cex, 0.5)
  expect_s3_class(scaled$gp, "gpar")

  # And it works on a grob with no settings at all.
  bare <- aev_scale_text(grid::textGrob("x"), 0.5)
  expect_equal(bare$gp$cex, 0.5)
  expect_s3_class(bare$gp, "gpar")
})

#### The helpers on their own ####

test_that("neither edit touches a single-panel table", {

  table <- assembled(autoplot(facet_aev_aev(), residuals = FALSE))

  expect_identical(aev_drop_lowest_x_axis(table), table)
  expect_identical(aev_shrink_lower_panel(table), table)
})

test_that("the panel rows are found in order, top to bottom", {

  table <- assembled(autoplot(facet_aev_aev()))
  rows <- aev_panel_rows(table)

  expect_length(rows, 2L)
  expect_false(is.unsorted(rows))
})

#### Groups ####

test_that("each group is its own column, in record order", {

  table <- assembled(autoplot(grouped_aev()))
  columns <- sort(unique(table$layout$l[startsWith(table$layout$name, "panel")]))

  expect_length(columns, 3L)

  # Record order, which sorting would not give.
  strips <- which(startsWith(table$layout$name, "aev-group"))
  strips <- strips[order(table$layout$l[strips])]
  expect_identical(
    unlist(lapply(table$grobs[strips], labels_in)),
    c("Period", "Age", "Cohort")
  )
})

test_that("a record is the same width whichever group it is in", {

  # One record, then three, then one: the columns should be one to three to one.
  table <- assembled(autoplot(grouped_aev()))
  columns <- sort(unique(table$layout$l[startsWith(table$layout$name, "panel")]))
  widths <- as.numeric(table$widths[columns])

  expect_identical(as.character(grid::unitType(table$widths[columns])), rep("null", 3))
  expect_equal(widths / widths[1], c(1, 3, 1))
})

test_that("the gap between groups is the design's, and runs through both panels", {

  # A column spans the whole grid, so one gap serves both charts -- which is
  # the prototype's rule that a gap "applies to both A/E and residuals".
  table <- assembled(autoplot(grouped_aev()))
  columns <- sort(unique(table$layout$l[startsWith(table$layout$name, "panel")]))
  for (column in columns[-length(columns)]) {
    expect_equal(pixels(table$widths[column + 1L]), aev_geometry$group_gap_px)
  }
})

test_that("the group names sit under the record names, not at the bottom", {

  # Slide 2's stack: plot area, x scale, group names, residual panel. The facet
  # puts column strips at the very bottom and they are moved.
  table <- assembled(autoplot(grouped_aev()))
  panels <- sort(rows_of(table, "panel"))
  axis <- aev_interior_axis_row(table)
  moved <- rows_of(table, "aev-group")

  # Counted first: `all()` of nothing is TRUE, so without this the test passed
  # happily when the strips had not been moved at all.
  expect_length(moved, 3L)
  expect_true(all(moved > axis))
  expect_true(all(moved < max(panels)))
})

test_that("the row the strips came from is emptied and closed up", {

  table <- assembled(autoplot(grouped_aev()))

  for (i in which(startsWith(table$layout$name, "strip-b"))) {
    expect_s3_class(table$grobs[[i]], "zeroGrob")
    expect_equal(as.numeric(table$heights[table$layout$t[i]]), 0)
  }
})

test_that("an unlabelled aev gets no group row and no gaps", {

  table <- assembled(autoplot(facet_aev_aev()))

  expect_length(rows_of(table, "strip-b"), 0L)
  expect_length(rows_of(table, "aev-group"), 0L)
  expect_length(unique(table$layout$l[startsWith(table$layout$name, "panel")]), 1L)
})

test_that("moving the strips is a no-op when there are none", {

  table <- assembled(autoplot(facet_aev_aev()))

  expect_identical(aev_move_group_strips(table), table)
})
