# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The main A/E panel. Nothing here compares pixels as images -- the assertions
# are on the object **ggplot2** builds and on the shapes it hands `grid`, which
# is where every decision the prototype fixes actually lives.
#
# THREE FAILURES THIS FILE EXISTS TO CATCH, all found by rendering and looking:
#
# * A record that cannot be drawn losing its place on the x axis, so that five
#   records of which one was drawable came out as a chart of one.
# * Scale limits censoring whole shapes rather than clipping them, which
#   silently deleted the caps of the first and last records.
# * The thin line running through the circle instead of stopping at its edge.
#
# None was visible to any test of the numbers, and all three are pinned below.

# Slide 8's own table, turned back into a triple: E is 1000, A is the stated
# A/E, and V is whatever makes 1.96 * sqrt(V) / E the stated 95% confidence.
#
# The second record matters beyond being a third row: its interval is narrower
# than the marker, which is the case the thin line has to be left out of.
prototype_aev <- function() {
  ratio <- c(1.218, 1.075, 1.001, 0.970)
  confidence <- c(0.058, 0.015, 0.101, 0.073)
  E <- rep(1000, 4)
  aev <- create_aev(A = ratio * E, E = E, V = (confidence / 1.959963984540 * E)^2)
  names(aev) <- c("[65, 75)", "[75, 85)", "[85, 95)", "[1930, 1940)")
  aev
}

# NO GROUP NAMES ON PURPOSE. Groups become a second facet dimension, so a
# grouped aev is a grid of columns and `layer_grob(p, i)[[1]]` is the first
# column rather than the whole row. Everything in this file is about the marker
# and the axes, so it uses the single-column layout and `test-aev-facet.R` owns
# the grouped one.

# One of each kind of record that cannot be drawn, plus one that can.
awkward_aev <- function() {
  aev <- create_aev_unchecked(
    A = c(1100, 0, NaN,   0, 3000),
    E = c(1000, 0, NaN,  50, 1000),
    V = c(2500, 0, NaN, 125, 2000)
  )
  names(aev) <- c("ok", "empty", "missing", "zero A", "too high")
  aev
}

# ON A NULL DEVICE for the same reason `interval_grob()` is: building a plot
# needs text metrics, and without a device R opens the default and leaves an
# `Rplots.pdf` behind.
panel_params <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplot_build(p)$layout$panel_params[[1]]
}

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))

# LAYERS ARE FOUND BY GEOM, NOT BY NUMBER. The residual panel added three more
# and moved every index, which is a silly way for a test to break.
layer_index <- function(p, geom) which(layer_geoms(p) == geom)[[1]]

marker_layer <- function(p) p$layers[[layer_index(p, "GeomAevInterval")]]

# The marker layer's grob, before `grid` has resolved anything. Its fields are
# the data-space positions in npc; the pixel sizes are not in it yet.
#
# ON A NULL DEVICE, because building a grob needs text metrics and so needs a
# device. Without one R opens the default, which writes an `Rplots.pdf` into
# `tests/testthat` and leaves it behind for `git status` to find.
interval_grob <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::layer_grob(p, layer_index(p, "GeomAevInterval"))[[1]]
}

# And after resolving, which needs a viewport of a known size, because that is
# the whole point: the pixel lengths only mean something once there is one.
#
# EVERY MEASUREMENT HAS TO HAPPEN INSIDE THAT VIEWPORT. Converting an npc length
# to inches after leaving it silently uses the device instead, and the first
# version of these tests did exactly that -- reporting a cap of 31.5 px, which
# is 18 scaled by the ratio of the 7-inch device to the 4-inch viewport. So the
# helper returns the derived numbers rather than letting a test compute them.
measured <- function(p, width = 4, height = 3) {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(height, "inches")
  ))

  # Not `interval_grob()`: that opens its own device, and the viewport just
  # pushed belongs to this one.
  grob <- ggplot2::layer_grob(p, layer_index(p, "GeomAevInterval"))[[1]]
  children <- grid::makeContent(grob)$children

  in_pixels <- function(length) {
    grid::convertWidth(length, "inches", valueOnly = TRUE) * 96
  }

  list(
    grob = grob,
    children = children,
    radius_npc = grid::convertHeight(
      grid::unit(aev_geometry$marker_radius_px / 96, "inches"), "npc",
      valueOnly = TRUE
    ),
    radius_px = in_pixels(children[[1]]$r),
    cap_px = in_pixels(children[[4]]$x1[1] - children[[4]]$x0[1])
  )
}

# npc of a value in the transformed (log A/E) space, worked out here rather
# than asked of the code under test.
as_npc <- function(values) (values - aev_log_ratio_limits[1]) / diff(aev_log_ratio_limits)

#### Status ####

test_that("status tells the four cases apart", {

  status <- aev_status(awkward_aev())

  expect_identical(
    as.character(status),
    c("ok", "empty", "missing", "off_scale", "off_scale")
  )
  expect_identical(levels(status), c("ok", "empty", "missing", "off_scale"))
})

test_that("status reads the triple, which is the only thing that can tell empty from missing", {

  # Both are NaN in every calculated property. If this ever starts agreeing,
  # the chart has lost the distinction.
  empty <- create_aev(A = 0, E = 0, V = 0)
  missing <- create_aev_unchecked(A = NaN, E = NaN, V = NaN)

  expect_true(is.nan(empty$A_over_E) && is.nan(missing$A_over_E))
  expect_identical(as.character(aev_status(empty)), "empty")
  expect_identical(as.character(aev_status(missing)), "missing")
})

test_that("missing beats empty beats off scale", {

  # A missing record can satisfy more than one test at once, so the order the
  # labels are applied in decides the answer. Missing is the strongest claim.
  expect_identical(as.character(aev_status(create_aev_unchecked(A = NaN, E = 0, V = 0))), "missing")
  expect_identical(as.character(aev_status(create_aev(A = 0, E = 0, V = 0))), "empty")
})

test_that("a zero-death record with exposure is off scale, not empty", {

  # Finite in the residual panel, infinite in this one.
  expect_identical(as.character(aev_status(create_aev(A = 0, E = 50, V = 125))), "off_scale")
})

test_that("status is judged against the limits it is given", {

  aev <- create_aev(A = 1500, E = 1000, V = 100)

  expect_identical(as.character(aev_status(aev)), "ok")
  expect_identical(as.character(aev_status(aev, limits = log(c(0.9, 1.1)))), "off_scale")
})

test_that("the plot frame numbers records that carry no label", {

  frame <- aev_plot_frame(create_aev(A = c(1, 2), E = c(1, 2), V = c(1, 2)))

  expect_identical(frame$position, 1:2)
  expect_identical(frame$name, c("1", "2"))
})

#### The chart as a whole ####

test_that("autoplot returns a ggplot", {

  p <- autoplot(prototype_aev())

  expect_s3_class(p, "ggplot")
  expect_identical(class(p + ggplot2::labs(subtitle = "extended"))[1], class(p)[1])
})

test_that("the chart is the layers the two panels need, in order", {

  # The whole marker is one grob rather than four layers, because the design
  # mixes data lengths with pixel lengths and only `grid` can hold both. Same
  # for the residual marker below it.
  # The two `GeomBlank`s draw nothing; they fix how far each panel reaches. The
  # `GeomPoint` draws nothing either; it exists so that a legend exists.
  # The chevron layer is always present and usually has no rows: a plot's
  # structure does not depend on whether anything happens to be off the scale.
  expect_identical(
    unname(layer_geoms(autoplot(prototype_aev()))),
    c("GeomRect", "GeomBlank", "GeomBlank", "GeomPoint", "GeomHline",
      "GeomHline", "GeomAevInterval", "GeomAevChevron", "GeomHline",
      "GeomAevResidual")
  )

  # Without the residual panel, only the A/E half of that.
  expect_identical(
    unname(layer_geoms(autoplot(prototype_aev(), residuals = FALSE))),
    c("GeomBlank", "GeomBlank", "GeomPoint", "GeomHline", "GeomHline",
      "GeomAevInterval", "GeomAevChevron")
  )
})

test_that("the 100% line is thick and coloured, not another grid line", {

  layers <- autoplot(prototype_aev())$layers
  reference <- layers[[which(vapply(
    layers,
    function(l) identical(l$aes_params$colour, aev_palette$reference_line),
    logical(1)
  ))[[1]]]]

  expect_identical(reference$aes_params$colour, aev_palette$reference_line)
  expect_equal(reference$aes_params$linewidth, line_width(1.5))

  # Three times the grid line, which is the point of it. It was five until
  # 2026-08-25, when the line was lightened because it was swallowing the dark
  # middle of any interval that crossed it, and the width moved with the colour.
  expect_equal(aev_geometry$reference_width_px / aev_geometry$grid_width_px, 3)
})

test_that("only drawable records reach the marker layer", {

  data <- marker_layer(autoplot(awkward_aev()))$data

  expect_identical(nrow(data), 1L)
  expect_identical(data$name, "ok")
})

#### Where the interval ends ####

test_that("the marker grob holds the five positions the design needs", {

  aev <- prototype_aev()
  centre <- log(aev$A_over_E)
  sigma <- aev$log_A_over_E_stddev
  grob <- interval_grob(autoplot(aev))

  expect_equal(grob$y, as_npc(centre))
  expect_equal(grob$inner_lower, as_npc(centre - sigma))
  expect_equal(grob$inner_upper, as_npc(centre + sigma))
  expect_equal(grob$outer_lower, as_npc(centre - aev_z_975 * sigma))
  expect_equal(grob$outer_upper, as_npc(centre + aev_z_975 * sigma))
})

test_that("the interval is the one logmu already computes", {

  # `log_A_over_E_95pc` is the accessor for this quantity, and the chart must
  # not have its own idea of what it is.
  aev <- prototype_aev()
  grob <- interval_grob(autoplot(aev))
  reach <- (grob$outer_upper - grob$y) * diff(aev_log_ratio_limits)

  expect_equal(reach, aev$log_A_over_E_95pc)
})

test_that("the interval is symmetric about the marker in log space", {

  grob <- interval_grob(autoplot(prototype_aev()))

  expect_equal(grob$inner_upper - grob$y, grob$y - grob$inner_lower)
  expect_equal(grob$outer_upper - grob$y, grob$y - grob$outer_lower)
})

#### What grid is actually given ####

test_that("the pieces are drawn in slide 7's order", {

  children <- measured(autoplot(prototype_aev()))$children

  expect_identical(
    unname(vapply(children, function(k) class(k)[1], character(1))),
    c("circle", "segments", "segments", "segments")
  )
})

test_that("the circle is five pixels, open and transparent", {

  m <- measured(autoplot(prototype_aev()))
  circle <- m$children[[1]]

  expect_equal(m$radius_px, 5)
  expect_true(is.na(circle$gp$fill))
  expect_identical(circle$gp$col, aev_palette$marker_weak)
  expect_equal(circle$gp$lwd, 0.75)
})

test_that("the thin line stops at the circle instead of running through it", {

  # THE DEFECT THIS REPLACES: the thin line was drawn straight from one end of
  # the interval to the other, crossing the marker. It was hidden for a while
  # by filling the circle with the plot background, and reappeared the moment
  # the centres were made transparent.
  m <- measured(autoplot(prototype_aev()))
  grob <- m$grob
  thin <- m$children[[2]]

  radius <- m$radius_npc
  drawn <- abs(grob$inner_upper - grob$y) > radius

  # Every arm starts exactly one radius from the centre, never at the centre.
  expect_equal(
    as.numeric(thin$y0),
    c(grob$y[drawn] + radius, grob$y[drawn] - radius)
  )
  expect_equal(
    as.numeric(thin$y1),
    c(grob$inner_upper[drawn], grob$inner_lower[drawn])
  )
})

test_that("an interval narrower than the marker gets no thin line at all", {

  # The second prototype record is that case: 1.5% either side of the point
  # against a five-pixel radius. Three records of four keep their arms.
  thin <- measured(autoplot(prototype_aev()))$children[[2]]

  expect_identical(length(as.numeric(thin$y0)), 6L)

  # And with nothing else on the chart there is no thin grob to draw.
  tight <- create_aev(A = 1000, E = 1000, V = 1e-6)
  expect_s3_class(measured(autoplot(tight))$children[[2]], "null")
})

test_that("the thick arms are drawn whether or not they intrude", {

  # The design says so explicitly, and it is what distinguishes them from the
  # thin line. Two arms per record, always.
  thick <- measured(autoplot(prototype_aev()))$children[[3]]

  expect_identical(length(as.numeric(thick$y0)), 8L)
  expect_identical(thick$gp$col, aev_palette$marker_medium)
})

test_that("the cap is eighteen pixels wide however many records there are", {

  # WHY THE CUSTOM GEOM EXISTS. An ordinary layer can only be given a width in
  # x units, so this was `0.03 * records` -- a fraction fitted to the two
  # sources rather than the width they both actually state.
  cap_width <- function(records, width) {
    aev <- create_aev(A = rep(1000, records), E = rep(1000, records), V = rep(2500, records))
    measured(autoplot(aev), width = width)$cap_px
  }

  expect_equal(cap_width(4, width = 4), 18)
  expect_equal(cap_width(16, width = 4), 18)

  # And however wide the chart is drawn, which is the other half of the claim.
  expect_equal(cap_width(4, width = 9), 18)
})

test_that("the three weights use the three greys and the stated pixel widths", {

  children <- measured(autoplot(prototype_aev()))$children

  expect_identical(children[[2]]$gp$col, aev_palette$marker_weak)
  expect_identical(children[[3]]$gp$col, aev_palette$marker_medium)
  expect_identical(children[[4]]$gp$col, aev_palette$marker_cap)

  # `grid`'s `lwd` is already in ninety-sixths of an inch, so the design's
  # pixel widths go straight in with no conversion.
  expect_equal(children[[2]]$gp$lwd, aev_geometry$thin_width_px)
  expect_equal(children[[3]]$gp$lwd, aev_geometry$thick_width_px)
  expect_equal(children[[4]]$gp$lwd, aev_geometry$cap_height_px)
})

test_that("the lines have square ends, as rectangles do", {

  # `grid` rounds line ends by default, which gave the arms and caps rounded
  # corners the design does not have -- `AEChart.cs` draws them as rectangles.
  # The default also lengthened them: a round end adds half a line width at
  # each end, so an 18 px cap 3 px thick was covering 21.
  children <- measured(autoplot(prototype_aev()))$children

  for (piece in children[2:4]) {
    expect_identical(piece$gp$lineend, "butt")
  }
})

test_that("the arm and cap weights are the ones the C# and the SVG carry", {

  # `AEChart.cs` and `Test AE chart.svg` both say 2.5 and 4, and both are back.
  #
  # They were lightened to 1.75 and 3 on 2026-08-20 and restored on 2026-08-25:
  # that lightening was calibrated against the panel rather than against one
  # record's width, so on a chart with fewer records than the prototype's
  # fourteen it left the marker too faint. The note in `aev_plot.R` has the
  # arithmetic. Anything lightening them again should fail here and read it.
  expect_equal(aev_geometry$thick_width_px, 2.5)
  expect_equal(aev_geometry$cap_height_px, 4)

  # The pieces still get heavier in the order they are read.
  expect_true(aev_geometry$thin_width_px < aev_geometry$thick_width_px)
  expect_true(aev_geometry$thick_width_px < aev_geometry$cap_height_px)
})

test_that("a linewidth for ggplot2 is not the same number as a width for grid", {

  # Measured by rendering to SVG and reading the numbers back: a `linewidth` of
  # 1 draws `.pt` pixels wide. Only the reference line and the theme need this;
  # everything inside the marker goes to `grid` in pixels as it stands.
  expect_equal(line_width(72.27 / 25.4), 1)
  expect_equal(line_width(0.75) * (72.27 / 25.4), 0.75)
})

#### Axes ####

test_that("the y axis is even in log A/E, not in A/E", {

  params <- panel_params(autoplot(prototype_aev()))

  # Slide 3's table: log A/E from minus 50% to plus 50% in tens. The breaks are
  # held in transformed space, so exponentiating recovers them, and it is the
  # LOGS that are evenly spaced.
  expect_equal(params$y$get_breaks(), log(aev_ratio_breaks))
  expect_identical(
    params$y$get_labels(),
    c("61%", "67%", "74%", "82%", "90%", "100%",
      "111%", "122%", "135%", "149%", "165%")
  )
})

test_that("the y axis runs tick to tick with nothing added at the ends", {

  params <- panel_params(autoplot(prototype_aev()))

  expect_equal(params$y.range, c(-0.5, 0.5))
})

test_that("the x axis leaves half a record at each end", {

  params <- panel_params(autoplot(prototype_aev()))

  expect_equal(params$x.range, c(0.5, 4.5))
})

test_that("every record keeps its place and its label, drawable or not", {

  # The failure this replaces: the axis followed the data, an undrawable record
  # contributed none, and four of five categories vanished from the chart.
  params <- panel_params(autoplot(awkward_aev()))

  expect_equal(params$x$get_breaks(), 1:5)
  expect_identical(params$x$get_labels(), c("ok", "empty", "missing", "zero A", "too high"))
  expect_equal(params$x.range, c(0.5, 5.5))
})

test_that("nothing is censored, so drawing emits no warning", {

  # Scale limits would drop any interval reaching past the end of the axis.
  # Coordinate limits clip instead.
  #
  # IT HAS TO BE DRAWN, NOT MERELY BUILT. `ggplot_build()` sets the censored
  # values to NA and says nothing; the "removed N rows" warning is raised when
  # the layer is rendered. A version of this test that only built the plot
  # watched nothing at all -- putting the limits back on the scales left it
  # green.
  render <- function(aev) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(autoplot(aev))
  }

  expect_no_warning(render(awkward_aev()))
  expect_no_warning(render(prototype_aev()))

  # A record whose marker is on the chart but whose interval runs off both
  # ends of it. Its A/E of 150% is inside the axis; its cap is at 400%.
  wide <- create_aev(A = 1500, E = 1000, V = 250000)
  expect_identical(as.character(aev_status(wide)), "ok")
  expect_no_warning(render(wide))
})

#### Odds and ends ####

test_that("an aev with no records builds and draws", {

  p <- autoplot(prototype_aev()[integer(0)], residuals = FALSE)

  expect_s3_class(p, "ggplot")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(ggplot2::ggplot_build(p))

  # ggplot2 drops an empty layer before the geom is ever asked to draw it.
  expect_s3_class(interval_grob(p), "zeroGrob")
})

test_that("plot draws and returns its input invisibly", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  aev <- prototype_aev()
  result <- withVisible(plot(aev))

  expect_false(result$visible)
  expect_identical(result$value, aev)
})

test_that("autoplot refuses anything that is not an aev", {

  expect_error(autoplot.aev(data.frame(A = 1)), "must be an `aev`")
})

test_that("the title is the caller's and there is none by default", {

  expect_null(autoplot(prototype_aev())$labels$title)
  expect_identical(autoplot(prototype_aev(), title = "Scheme A")$labels$title, "Scheme A")
})
