# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The bracket round each group name: a vertical rule at each end of the group's
# span and a horizontal rule across the middle, broken either side of the name.
#
# The rule is `XAxisGroups.Draw` in `AEChart.cs`, and `Test AE chart.svg` is a
# rendering of it. Both were read; the C# is what the arithmetic comes from.

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

# Resolved inside a viewport of a stated size, so the npc lengths mean
# something a test can check.
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

#### Where the brackets are ####

test_that("each group gets a bracket and an ungrouped chart gets none", {

  expect_identical(bracket_cells(autoplot(grouped_aev())), c("aev-group-1", "aev-group-2"))

  plain <- create_aev(A = 1100, E = 1000, V = 2500)
  expect_length(bracket_cells(autoplot(plain)), 0L)
})

test_that("the bracket keeps the strip rather than redrawing the name", {

  # The name stays in whatever the theme says, which is the reason the strip is
  # carried as a child instead of the text being drawn again here.
  children <- bracket_children(autoplot(grouped_aev()))

  expect_length(children, 3L)
  expect_identical(labels_in(children[[3]]), "one")
})

#### The shape ####

test_that("a vertical rule sits at each end of the group's span", {

  ends <- bracket_children(autoplot(grouped_aev()))[[1]]

  expect_s3_class(ends, "segments")
  expect_equal(as.numeric(ends$x0), c(0, 1))
  expect_equal(as.numeric(ends$x1), c(0, 1))
  # The full depth of the name band.
  expect_equal(as.numeric(ends$y0), c(0, 0))
  expect_equal(as.numeric(ends$y1), c(1, 1))
})

test_that("the horizontal rule is halfway down and broken about the middle", {

  across <- bracket_children(autoplot(grouped_aev()))[[2]]

  expect_s3_class(across, "segments")
  expect_equal(as.numeric(across$y0), c(0.5, 0.5))
  expect_equal(as.numeric(across$y1), c(0.5, 0.5))

  left <- as.numeric(across$x1[1]) - as.numeric(across$x0[1])
  right <- as.numeric(across$x1[2]) - as.numeric(across$x0[2])

  # The two segments share what is left equally, and each starts at an end.
  expect_equal(left, right)
  expect_equal(as.numeric(across$x0[1]), 0)
  expect_equal(as.numeric(across$x1[2]), 1)
})

test_that("the break is the name plus a margin at each end", {

  width <- 2
  children <- bracket_children(autoplot(grouped_aev()), width = width)
  across <- children[[2]]

  gap <- as.numeric(across$x0[2]) - as.numeric(across$x1[1])

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"), height = grid::unit(0.2, "inches")
  ))

  text <- aev_text_within(children[[3]])
  label <- grid::convertWidth(grid::grobWidth(text), "npc", valueOnly = TRUE)
  margin <- grid::convertWidth(
    inches(aev_geometry$group_gap_px), "npc", valueOnly = TRUE
  )

  expect_equal(gap, label + 2 * margin)
})

test_that("a name with no room for a rule gets none, and keeps its verticals", {

  # The verticals carry the grouping on their own. The C# drops the horizontal
  # the same way, on `spaceForLine > 0`.
  # 0.2 inches is 19 px, which the name and its two 6 px margins overrun.
  children <- bracket_children(autoplot(grouped_aev()), width = 0.2)

  expect_s3_class(children[[1]], "segments")
  expect_s3_class(children[[2]], "null")
})

test_that("the bracket is drawn in the grid line's colour and weight", {

  children <- bracket_children(autoplot(grouped_aev()))

  for (index in 1:2) {
    expect_identical(children[[index]]$gp$col, aev_palette$grid_line)
    expect_identical(children[[index]]$gp$lwd, aev_geometry$grid_width_px)
    expect_identical(children[[index]]$gp$lineend, "butt")
  }
})

#### That it reaches the page ####

test_that("the bracket resolves when the chart is drawn", {

  # Same trap as the chevron: a layer-level test finds `makeContent` by scope,
  # and `grid` finds it only if the package registers it.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(autoplot(grouped_aev())))

  brackets <- Filter(
    function(grob) inherits(grob, "forcedgrob") && inherits(grob, "aev_group_bracket"),
    grobs_in(forced)
  )

  expect_length(brackets, 2L)

  # NOT `length(children) > 0`. A bracket holds the strip it wraps as a child
  # from the moment it is constructed, so that is true whether or not
  # `makeContent` ever ran -- it asserted the constructor, not the drawing.
  # What exists only after the content is made is the rule ink.
  ink <- unlist(lapply(brackets, function(bracket) {
    lapply(grobs_in(bracket), function(grob) {
      if (inherits(grob, "segments")) grob$gp$col else NULL
    })
  }))

  expect_gte(sum(ink == aev_palette$grid_line), 2L)
})

#### Runs, not names ####

# THE DEFECT THIS COVERS, found by Tim on 2026-08-25. The chart used to be split
# on the group NAME, which gathered every record sharing a name into one place
# however far apart they were. A name that occurred twice left two stretches of
# the axis overlapping: the record names came out twice and records were drawn
# under other records' names, which read as records having gone missing.
#
# `XAxisGroups.Draw` closes a run when the NEXT name differs, and that is the
# rule: a group is a run of adjacent records, not every record sharing a name.

repeated_aev <- function() {
  aev <- create_aev(A = rep(1000, 6), E = rep(1000, 6), V = rep(2500, 6))
  names(aev) <- letters[1:6]
  group_names(aev) <- c("one", "one", "two", "two", "one", "one")
  aev
}

test_that("a run ends where the group name changes", {

  expect_identical(as.integer(aev_runs(c("a", "a", "b", "b", "a"))), c(1L, 1L, 2L, 2L, 3L))
  expect_identical(levels(aev_runs(c("a", "a", "b"))), c("1", "2"))
  expect_length(aev_runs(character(0)), 0L)
})

test_that("records with no group name are in the same run as each other", {

  # `NA != NA` is NA and cannot start a run, so the comparison is made by hand.
  expect_identical(as.integer(aev_runs(c(NA, NA))), c(1L, 1L))
  expect_identical(as.integer(aev_runs(c(NA, "a", NA))), c(1L, 2L, 3L))
})

test_that("a group name that occurs twice gets a bracket each time", {

  expect_identical(bracket_cells(autoplot(repeated_aev())), paste0("aev-group-", 1:3))

  labels <- vapply(
    1:3,
    function(k) labels_in(bracket_grob(autoplot(repeated_aev()), k)),
    character(1)
  )

  expect_identical(labels, c("one", "two", "one"))
})

test_that("records stay in the order they were given", {

  # Three runs of two, so the columns must cover the axis in order and without
  # overlapping. Splitting on the name gave two columns covering 0.5 to 6.5 and
  # 2.5 to 4.5, one inside the other.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  built <- ggplot2::ggplot_build(autoplot(repeated_aev()))
  ranges <- lapply(built$layout$panel_params, function(panel) panel$x.range)

  expect_equal(ranges[[1]], c(0.5, 2.5))
  expect_equal(ranges[[2]], c(2.5, 4.5))
  expect_equal(ranges[[3]], c(4.5, 6.5))
})

test_that("each record name is written once", {

  # The symptom Tim saw: every column drew every break falling inside it, so a
  # name appeared under a column that did not hold its record.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  table <- ggplot2::ggplotGrob(autoplot(repeated_aev()))
  axes <- which(startsWith(table$layout$name, "axis-b"))
  written <- unlist(lapply(table$grobs[axes], labels_in))

  expect_setequal(written, letters[1:6])
  expect_length(written, 6L)
})


#### The chart drawn without its residual panel ####

test_that("groups are bracketed even when the residual panel is off", {

  # FOUND BY REVIEW, 2026-08-26. The band was placed by looking for the axis row
  # BETWEEN the panels, and with one panel there is none -- so the strips were
  # abandoned wherever the facet had left them, with no bracket and no gap. The
  # documentation promises the bracket unconditionally.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (residuals in c(TRUE, FALSE)) {

    table <- ggplot2::ggplotGrob(autoplot(grouped_aev(), residuals = residuals))
    brackets <- which(startsWith(table$layout$name, "aev-group-"))

    expect_length(brackets, 2L)
    expect_true(all(vapply(
      table$grobs[brackets], inherits, logical(1), "aev_group_bracket"
    )))
  }
})
