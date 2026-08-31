# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The off-scale chevron: the mark against the top or bottom of the A/E panel
# saying that a record went that way, drawn INSTEAD of the record.
#
# THE TEST IS THE INNER INTERVAL -- A/E give or take one standard deviation --
# and not the point and not the 95% arms. Tim's rule of 2026-08-25. A record
# whose inner interval misses the panel is replaced by a chevron even where its
# outer arms would have shown, because a bar with no marker and no cap says less
# than the chevron does. `aev_status()` asks the same question the same way, so
# the two agree by construction and there is a test below that says so.

# One of each case the rule has to tell apart.
#
# * ordinary       A/E 110%, drawn normally.
# * off the top    A/E 200% with a tight interval -- nothing reaches back in.
# * bar reaches in A/E 190% but sqrt(V)/E of 0.3, so the inner interval comes
#                  down to 141% and the record is drawn, clipped at the edge.
# * no deaths      A of zero, so log A/E is -Inf whatever the range.
# * off the bottom A/E 40% with a tight interval.
# * nobody         A, E and V all zero.
# * unknown        all missing.
mixed_aev <- function() {
  aev <- create_aev(
    A = c(1100, 2000, 1900,   0,  400, 0, NaN),
    E = c(1000, 1000, 1000,  50, 1000, 0, NaN),
    V = c(2500, 2000, 90000, 125,  400, 0, NaN)
  )
  names(aev) <- c(
    "ordinary", "off the top", "bar reaches in", "no deaths",
    "off the bottom", "nobody", "unknown"
  )
  aev
}

chevron_layer <- function(p) {
  index <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) ==
                   "GeomAevChevron")
  p$layers[[index[[1]]]]
}

chevron_data <- function(p) chevron_layer(p)$data

marker_data <- function(p) {
  index <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) ==
                   "GeomAevInterval")
  p$layers[[index[[1]]]]$data
}

# The resolved shapes, measured inside a viewport of a known size so that the
# pixel rules come out in the units the design states them in.
chevron_children <- function(p, width = 6, height = 4) {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(height, "inches")
  ))

  index <- which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) ==
                   "GeomAevChevron")
  grob <- ggplot2::layer_grob(p, index[[1]])[[1]]
  grid::makeContent(grob)$children
}

#### Which records are marked ####

test_that("only records whose inner interval misses the panel are marked", {

  marked <- aev_off_scale(aev_plot_frame(mixed_aev()), c(-0.5, 0.5))

  expect_identical(marked$name, c("off the top", "no deaths", "off the bottom"))
  expect_identical(marked$top, c(TRUE, FALSE, FALSE))
})

test_that("a record whose interval reaches into the panel is not marked", {

  # A/E of 190% with sqrt(V)/E of 0.3: the bottom of the interval is at
  # log(1.9) - 0.3, which is inside the default range.
  reaching <- create_aev(A = 1900, E = 1000, V = 90000)
  expect_lt(log(1.9) - 0.3, 0.5)

  expect_identical(nrow(aev_off_scale(aev_plot_frame(reaching), c(-0.5, 0.5))), 0L)

  # And it is drawn, marker and all, even though its centre is off the top.
  # What shows is the lower part of its interval, clipped at the edge.
  expect_identical(as.character(aev_status(reaching)), "ok")
})

test_that("the test is one standard deviation, not the 95% arms", {

  # Tim's rule, 2026-08-25. This record sits at log A/E of 0.80 with a standard
  # deviation of 0.2, so its inner interval starts at 0.60 and is off a scale
  # that stops at 0.50 -- but its outer arm reaches back to 0.41 and would be
  # partly visible. The chevron replaces it.
  outer <- create_aev(A = 2226, E = 1000, V = 40000)

  expect_gt(log(2.226) - 0.2, 0.5)
  expect_lt(log(2.226) - aev_z_975 * 0.2, 0.5)

  marked <- aev_off_scale(aev_plot_frame(outer), c(-0.5, 0.5))
  expect_identical(nrow(marked), 1L)
  expect_true(marked$top)
})

test_that("a marked record is not drawn as well", {

  # The chevron REPLACES the marker. Both would be a record counted twice.
  outer <- create_aev(A = 2226, E = 1000, V = 40000)
  names(outer) <- "replaced"

  expect_identical(as.character(aev_status(outer)), "off_scale")
  expect_identical(nrow(marker_data(autoplot(outer))), 0L)
  expect_identical(nrow(chevron_data(autoplot(outer))), 1L)
})

test_that("status and the chevron never disagree", {

  # One rule, asked in two places: whether a record is off the scale, and which
  # way it went. A record marked by one and not the other is a chart with a
  # category missing or a category counted twice.
  # LOG LIMITS, not a ratio. They were a ratio until 2026-08-26, when the
  # round trip through exp() and log() turned out to lose a unit in the last
  # place on some ranges and put this test and the chevron's on opposite sides
  # of it.
  frame <- aev_plot_frame(mixed_aev(), c(-0.5, 0.5))
  marked <- aev_off_scale(frame, c(-0.5, 0.5))

  expect_setequal(marked$name, frame$name[frame$status == "off_scale"])
})

test_that("a record with no deaths is marked at the bottom like any other", {

  # No special case anywhere: -Inf plus a finite reach is still -Inf.
  zero <- create_aev(A = 0, E = 50, V = 125)
  marked <- aev_off_scale(aev_plot_frame(zero), c(-0.5, 0.5))

  expect_identical(nrow(marked), 1L)
  expect_false(marked$top)
})

test_that("empty and missing records are not marked", {

  # Both read NaN in every calculated property, so every comparison gives NA.
  # A chevron says which way a record went and neither of these went anywhere.
  nothing <- create_aev(A = c(0, NaN), E = c(0, NaN), V = c(0, NaN))

  expect_identical(nrow(aev_off_scale(aev_plot_frame(nothing), c(-0.5, 0.5))), 0L)
})

test_that("the axis decides what is off it", {

  aev <- create_aev(A = 2000, E = 1000, V = 2000)

  expect_identical(nrow(aev_off_scale(aev_plot_frame(aev), c(-0.5, 0.5))), 1L)
  # log(2) is 0.693, and the interval reaches down to 0.606.
  expect_identical(nrow(aev_off_scale(aev_plot_frame(aev), c(-0.6, 0.6))), 1L)
  expect_identical(nrow(aev_off_scale(aev_plot_frame(aev), c(-0.7, 0.7))), 0L)
})

#### The layer ####

test_that("the chart carries a chevron layer even when nothing is off scale", {

  # The layer is always there, so a plot's structure does not depend on its
  # data. It draws nothing when it has no rows.
  ordinary <- create_aev(A = 1100, E = 1000, V = 2500)

  expect_identical(nrow(chevron_data(autoplot(ordinary))), 0L)
  expect_identical(nrow(chevron_data(autoplot(mixed_aev()))), 3L)
})

test_that("chevrons are drawn in the A/E panel only", {

  data <- chevron_data(autoplot(mixed_aev()))

  expect_identical(levels(droplevels(data$panel)), aev_panels[1])
})

test_that("a layer with no rows draws nothing", {

  ordinary <- create_aev(A = 1100, E = 1000, V = 2500)
  children <- chevron_children(autoplot(ordinary))

  expect_length(children, 0L)
})

#### What grid is given ####

test_that("each chevron is one three-point polyline", {

  children <- chevron_children(autoplot(mixed_aev()))

  expect_length(children, 1L)
  expect_s3_class(children[[1]], "polyline")
  expect_length(children[[1]]$x, 9L)
  expect_identical(children[[1]]$id, rep(1:3, each = 3))
})

test_that("the chevron is the design's colour and weight, with butt ends", {

  style <- chevron_children(autoplot(mixed_aev()))[[1]]$gp

  expect_identical(style$col, aev_palette$off_scale)
  expect_identical(style$lwd, aev_chevron$line_width_px)
  expect_identical(style$lineend, "butt")
  # A mitred join is what keeps the apex sharp. Two butt-ended segments would
  # leave a notch in it.
  expect_identical(style$linejoin, "mitre")
})

test_that("the chevron is as wide and as tall as the design says", {

  height <- 4
  child <- chevron_children(autoplot(mixed_aev()), height = height)[[1]]

  x <- as.numeric(child$x)
  y <- as.numeric(child$y)

  # First chevron: three points at the left arm, the apex and the right arm.
  expect_equal(
    (x[3] - x[1]) * pixels(grid::unit(6, "inches")),
    aev_chevron$width_px
  )
  expect_equal(
    abs(y[2] - y[1]) * height * 96,
    aev_chevron$height_px
  )
  # The apex is on the axis of the arms.
  expect_equal(x[2], mean(c(x[1], x[3])))
})

test_that("the apex clears the panel edge by the inset and half the mitre", {

  height <- 4
  child <- chevron_children(autoplot(mixed_aev()), height = height)[[1]]

  y <- as.numeric(child$y)
  clearance <- aev_chevron$inset_px + aev_chevron$line_width_px / 2

  # The first chevron is the one at the top, so its apex is that far below 1.
  expect_equal((1 - y[2]) * height * 96, clearance)
  # The other two are at the bottom, the same distance above 0.
  expect_equal(y[5] * height * 96, clearance)
  expect_equal(y[8] * height * 96, clearance)
})

test_that("a chevron at the top points up and one at the bottom points down", {

  child <- chevron_children(autoplot(mixed_aev()))[[1]]
  y <- as.numeric(child$y)

  # The apex is the middle point of each three. Above its arms at the top of
  # the panel, below them at the bottom.
  expect_gt(y[2], y[1])
  expect_lt(y[5], y[4])
  expect_lt(y[8], y[7])
})

test_that("the chevron survives being drawn, not only being built", {

  # THE ONE THAT MATTERS, and the one every other test here misses.
  #
  # Asking a layer for its grob and calling `grid::makeContent()` on it finds
  # `makeContent.aev_chevron` by ordinary scoping, registered or not. `grid`
  # dispatches from its own namespace when it draws, and finds nothing unless
  # the package registers the method -- so a missing `S3method()` line leaves
  # every test above passing and the chart blank. It happened on 2026-08-25.
  #
  # `forceGrob` resolves the tree the way drawing does, which is what makes
  # this test see what the others cannot.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(autoplot(mixed_aev())))
  drawn <- unlist(lapply(grobs_in(forced), function(grob) grob$gp$col))

  expect_true(aev_palette$off_scale %in% drawn)
})

test_that("every grob the chart resolves at draw time actually resolves", {

  # The general form of the bug above, covering all three custom geoms rather
  # than only the chevron. A `makeContent` grob that grid cannot dispatch to
  # comes out of the forcing pass with no children and draws nothing at all.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(autoplot(mixed_aev())))

  # `forceGrob` keeps the unresolved original alongside the resolved copy, so
  # the resolved ones are those that also carry its own class.
  ours <- Filter(
    function(grob) {
      inherits(grob, "forcedgrob") &&
        inherits(grob, c("aev_interval", "aev_residual_marker", "aev_chevron"))
    },
    grobs_in(forced)
  )

  expect_length(ours, 3L)
  expect_true(all(vapply(ours, function(grob) length(grob$children) > 0L, logical(1))))
})

test_that("the whole chart draws without complaint", {

  render <- function(p) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }

  expect_no_warning(render(autoplot(mixed_aev())))
  expect_no_warning(render(autoplot(mixed_aev(), residuals = FALSE)))
})


#### Records the two tests used to disagree about ####

test_that("a record with exposure of zero is marked, not dropped", {

  # FOUND BY REVIEW, 2026-08-26. `aev_status()` and the chevron asked the same
  # question through opposite polarities of `%in% TRUE`, so a row whose
  # comparisons are all NA was off scale to one and not off scale to the other.
  # It was drawn in neither panel and marked by nothing: it simply vanished,
  # keeping only its slot on the x axis.
  #
  # Not reachable through `create_aev()`, which validates. Reachable through a
  # computed `aev`, which is not, and through `exp()` underflow.
  odd <- create_aev_unchecked(A = 3, E = 0, V = 0)

  expect_identical(as.character(aev_status(odd)), "off_scale")

  marks <- aev_off_scale(aev_plot_frame(odd, c(-0.5, 0.5)), c(-0.5, 0.5))

  expect_identical(nrow(marks), 1L)
  expect_true(marks$top)
})

test_that("every off-scale record with a direction is marked, and nothing else is", {

  # The invariant, stated properly rather than checked against one fixture.
  # A record is marked exactly when it is off scale AND its A/E says which way,
  # which is every off-scale record whose A/E is not itself NaN.
  aev <- create_aev_unchecked(
    A = c(1100, 2000,   0,   3,   0, NaN, 400),
    E = c(1000, 1000,  50,   0,   0, NaN, 1000),
    V = c(2500,    1, 125,   0,   5, NaN,  400)
  )

  frame <- aev_plot_frame(aev, c(-0.5, 0.5))
  marks <- aev_off_scale(frame, c(-0.5, 0.5))

  off <- frame$status == "off_scale"
  directed <- off & is.finite(frame$log_ratio | 0) | (off & is.infinite(frame$log_ratio))

  expect_identical(sum(marks$top) + sum(!marks$top), nrow(marks))
  expect_identical(nrow(marks), sum(directed))

  # And no record is marked twice, or marked while drawable.
  expect_false(any(marks$status == "ok"))
  expect_identical(anyDuplicated(marks$position), 0L)
})

test_that("a limit that does not survive exp and log still decides one way", {

  # `log(exp(0.35))` is a unit in the last place short of 0.35. Status used to
  # be judged against `exp(log_range)` and re-logged while the chevron used the
  # raw range, so a record sitting exactly on such a limit could be off scale
  # to one and on scale to the other. Both now read the same number.
  reach <- 0.05
  aev <- create_aev(A = 1000 * exp(0.35 + reach), E = 1000, V = (reach * 1000)^2)

  frame <- aev_plot_frame(aev, c(-0.35, 0.35))
  marks <- aev_off_scale(frame, c(-0.35, 0.35))

  drawn <- as.character(frame$status) == "ok"
  expect_identical(sum(drawn) + nrow(marks), 1L)
})
