# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The deviance residual panel: the strip below the A/E chart that carries the
# actual test, where the A/E panel carries magnitude and precision.
#
# THE PANEL IS NOT A SECOND VIEW OF THE SAME RECORDS. A cell with no deaths has
# no position on a log scale and a perfectly good residual, so it appears here
# and not above. That difference is the reason both panels exist.
#
# Nothing here compares images. The assertions are on the layers, on the fixed
# panel geometry, and on the shapes handed to `grid`.

# Two records that can be drawn in both panels, one with no deaths at all
# (absent above, present below), and one whose residual is far off this panel's
# scale (drawn but clipped, with its value still printed).
mixed_aev <- function() {
  aev <- create_aev(
    A = c(1100,  970,    0, 1500),
    E = c(1000, 1000,   50, 1000),
    V = c(2500, 1387,  125,  100)
  )
  names(aev) <- c("ok", "small", "zero A", "loud")
  aev
}

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))

layer_index <- function(p, geom) which(layer_geoms(p) == geom)[[1]]

residual_layer <- function(p) p$layers[[layer_index(p, "GeomAevResidual")]]

# As with the marker, everything pixel-valued is measured inside a viewport of
# a known size. Here that size is the one the theme actually gives the panel,
# so the derived rules come out in the units the design states them in.
residual_children <- function(p, width = 4) {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(width, "inches"),
    height = grid::unit(aev_residual$height_px / 96, "inches")
  ))

  # `[[2]]`: `layer_grob()` returns one grob per panel and this layer lives in
  # the second. Asking for the first gives a `zeroGrob` and no children at all.
  grob <- ggplot2::layer_grob(p, layer_index(p, "GeomAevResidual"))[[2]]
  grid::makeContent(grob)$children
}

residual_grob <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::layer_grob(p, layer_index(p, "GeomAevResidual"))[[2]]
}

panel_params <- function(p, panel) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplot_build(p)$layout$panel_params[[panel]]
}

#### The two panels ####

test_that("there are two panels by default and one on request", {

  panel_count <- function(p) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    length(ggplot2::ggplot_build(p)$layout$panel_params)
  }

  expect_identical(panel_count(autoplot(mixed_aev())), 2L)
  expect_identical(panel_count(autoplot(mixed_aev(), residuals = FALSE)), 1L)
})

test_that("the panel names are the y-axis titles, A/E first", {

  # They are facet strips rather than axis titles, because a plot has one
  # `axis.title.y` and two panels need two.
  expect_identical(aev_panels[1], "A/E\n(log\u00A0scale)")
  # `~` rather than U+2248, which is not in latin1 and so cannot be drawn on the
  # device `R CMD check` uses. Tim's ruling of 2026-08-29, and it matches the
  # legend, which has always written "(~68%)". U+00B1 stays: latin1 holds it.
  expect_identical(aev_panels[2], "Deviance residual\n(\u00B12 ~ 95%)")
  expect_false(grepl("\u2248", aev_panels[2], fixed = TRUE))

  # The qualifier is on its own line in both, and the residual title is the
  # only place on the chart that ties the two panels' statistics together.
  expect_true(all(grepl("\n", aev_panels, fixed = TRUE)))
  expect_true(grepl("95%", aev_panels[2], fixed = TRUE))

  levels <- levels(residual_layer(autoplot(mixed_aev()))$data$panel)
  expect_identical(levels, aev_panels)
})

test_that("the residual panel is a fixed depth and the A/E panel takes the rest", {

  # `MeasureVariableSize` in `AEChart.cs` returns a constant, so the strip is
  # 84 px whatever the chart and does not grow with it.
  heights <- aev_panel_heights(TRUE)

  expect_length(heights, 2L)
  expect_identical(as.character(grid::unitType(heights)), c("null", "inches"))
  expect_equal(as.numeric(heights)[2] * 96, 84)

  expect_identical(as.character(grid::unitType(aev_panel_heights(FALSE))), "null")
})

test_that("each panel gets its own range from its own data", {

  # `scales = \"free_y\"` calls the limits function once per panel. Without that
  # both panels would share one range and one of them would be unreadable.
  p <- autoplot(mixed_aev())

  expect_equal(panel_params(p, 1)$y.range, aev_log_ratio_limits)
  expect_equal(panel_params(p, 2)$y.range, c(-3.5, 3.5))
})

test_that("the limits dispatch on the range, not on equality with it", {

  # The range reaching this function may carry the scale's expansion. An
  # equality test was tried first and fell through to the A/E limits when it
  # did, which put percentage labels on the residual axis.
  expect_equal(aev_panel_limits(c(-0.3, 0.2)), aev_log_ratio_limits)
  expect_equal(aev_panel_limits(c(-1.38, 2.44)), c(-3.5, 3.5))
  expect_equal(aev_panel_limits(c(-3.85, 3.85)), c(-3.5, 3.5))
})

test_that("one scale labels the two panels differently", {

  p <- autoplot(mixed_aev())

  expect_equal(panel_params(p, 1)$y$get_breaks(), log(aev_ratio_breaks))
  expect_identical(panel_params(p, 1)$y$get_labels()[6], "100%")

  expect_equal(panel_params(p, 2)$y$get_breaks(), -3:3)
  # PLOTMATH, not strings: the sign is drawn from the Symbol font so that a true
  # minus survives a latin1 device. `deparse` shows the construction -- a quoted
  # number so the formatting is kept, with the sign as an operator outside it.
  labels <- panel_params(p, 2)$y$get_labels()

  # A LIST HERE, not an expression vector: `aev_panel_labels()` returns an
  # expression and ggplot2 stores it unpacked into a list of language objects.
  # The content is what matters, and `deparse` reads it either way.
  expect_type(labels, "list")
  expect_identical(
    vapply(labels, deparse, character(1)),
    c("-\"3\"", "-\"2\"", "-\"1\"", "\"0\"", "+\"1\"", "+\"2\"", "+\"3\"")
  )
})

test_that("positive residuals are at the top", {

  # `AEChart.cs` places its markers that way but writes its scale labels the
  # other way up, and the two rendered samples disagree with each other because
  # of it. The markers are right: sample A has +3 at the top.
  params <- panel_params(autoplot(mixed_aev()), 2)

  expect_identical(deparse(params$y$get_labels()[[1]]), "-\"3\"")
  expect_identical(deparse(params$y$get_labels()[[7]]), "+\"3\"")
})

test_that("numbers are written with a minus sign, and zero with no sign", {

  # Slide 1: never a hyphen. The character is U+2212.
  # THE DIGITS ARE QUOTED AND THE SIGN IS NOT. `parse(text = "-2.60")` renders
  # "2.6" -- plotmath reads the digits as a number and drops the trailing zero,
  # silently reformatting every residual. Quoting them keeps the formatting and
  # leaves the sign as an operator, which is what draws it from Symbol.
  expect_identical(
    vapply(aev_signed_math(c(-1.5, 0, 2.25), 2), deparse, character(1)),
    c("-\"1.50\"", "\"0.00\"", "+\"2.25\"")
  )
  expect_identical(
    vapply(aev_signed_math(-3:3, 0), deparse, character(1)),
    c("-\"3\"", "-\"2\"", "-\"1\"", "\"0\"", "+\"1\"", "+\"2\"", "+\"3\"")
  )
})

#### Which records appear ####

test_that("the residual panel shows a record the A/E panel cannot", {

  # A cell with no deaths: log(A/E) is infinite, so it has no marker above, but
  # its residual is finite and it is drawn below.
  p <- autoplot(mixed_aev())
  above <- p$layers[[layer_index(p, "GeomAevInterval")]]$data
  below <- residual_layer(p)$data

  expect_false("zero A" %in% above$name)
  expect_true("zero A" %in% below$name)
})

test_that("records with no residual at all are left out", {

  # Empty and missing cells read NaN, and NaN is not a residual.
  aev <- create_aev_unchecked(A = c(1100, 0, NaN), E = c(1000, 0, NaN), V = c(2500, 0, NaN))
  names(aev) <- c("ok", "empty", "missing")

  expect_identical(residual_layer(autoplot(aev))$data$name, "ok")
})

#### The ramp ####

test_that("weight and colour follow the size of the residual", {

  # `AEChart.cs`: half a pixel below 0.5 and doubling at 1.5. BOTH ENDS OF THE
  # TOP SEGMENT ARE TIM'S, on 2026-08-25: the width it reaches went from the
  # C#'s 2 to 4, so that a cell fitting this badly cannot be missed, and the
  # residual it reaches it at went from 2.5 to 3, which is where the outer tint
  # band ends. Flat after that.
  expect_equal(aev_residual_width(c(0, 0.25, 0.5)), rep(0.5, 3))
  expect_equal(aev_residual_width(1.5), 1)
  expect_equal(aev_residual_width(3), 4)
  expect_equal(aev_residual_width(c(6, -6)), c(4, 4))

  # Halfway between the breakpoints, halfway between the widths.
  expect_equal(aev_residual_width(1), 0.75)
  expect_equal(aev_residual_width(2.25), 2.5)
})

test_that("the ramp is symmetric in the sign of the residual", {

  expect_equal(aev_residual_width(-1.7), aev_residual_width(1.7))
  expect_identical(aev_residual_colour(-1.7), aev_residual_colour(1.7))
})

test_that("the colour ramp runs weak to medium to strong", {

  # The endpoints are the palette's own greys; only the interpolation between
  # them is computed, and it is done in CIELAB to match `LerpLab`.
  expect_identical(toupper(aev_residual_colour(0)), aev_palette$marker_weak)
  expect_identical(toupper(aev_residual_colour(1.5)), aev_palette$marker_medium)
  expect_identical(toupper(aev_residual_colour(3)), aev_palette$marker_strong)
  expect_identical(toupper(aev_residual_colour(9)), aev_palette$marker_strong)
})

test_that("a bigger residual is never a lighter one", {

  greys <- vapply(
    aev_residual_colour(seq(0, 4, by = 0.1)),
    function(hex) mean(grDevices::col2rgb(hex)),
    numeric(1)
  )

  expect_false(is.unsorted(rev(greys)))
})

#### What grid is given ####

test_that("the residual grob carries the marker, the zero line and the residual", {

  aev <- mixed_aev()
  grob <- residual_grob(autoplot(aev))
  z <- aev$deviance_residual[is.finite(aev$deviance_residual)]

  expect_identical(grob$z, z)
  expect_equal(grob$zero, rep(0.5, length(z)))
  expect_equal(grob$y, (z + 3.5) / 7)
})

test_that("the panel draws a drop line, a circle and a printed value", {

  children <- residual_children(autoplot(mixed_aev()))

  expect_identical(
    unname(vapply(children, function(k) class(k)[1], character(1))),
    c("segments", "circle", "text")
  )
})

test_that("the marker is three pixels and the drop line stops at its edge", {

  children <- residual_children(autoplot(mixed_aev()))
  circle <- children[[2]]

  expect_equal(pixels(circle$r), 3)
  expect_identical(circle$gp$col, aev_palette$marker_medium)
  expect_true(is.na(circle$gp$fill))
  expect_equal(circle$gp$lwd, 1)
})

test_that("a residual too small to clear its own marker gets no drop line", {

  # At the panel's fixed depth one unit of residual is 12 px and the marker is
  # 3 px across, so the line appears from a quarter of a unit upwards.
  just_under <- create_aev(A = 1000, E = 1000, V = 1600)
  just_over <- create_aev(A = 1030, E = 1000, V = 1600)

  expect_lt(abs(just_under$deviance_residual), 0.25)
  expect_gt(abs(just_over$deviance_residual), 0.25)

  expect_s3_class(residual_children(autoplot(just_under))[[1]], "null")
  expect_s3_class(residual_children(autoplot(just_over))[[1]], "segments")
})

test_that("the printed value sits across the zero line from its marker", {

  # Which is what makes the slide's claim that the values always fit true: the
  # marker occupies one half of the panel and the number the other.
  children <- residual_children(autoplot(mixed_aev()))
  grob <- residual_grob(autoplot(mixed_aev()))
  text <- children[[3]]

  above <- grob$z >= 0
  text_y <- as.numeric(text$y)

  expect_true(all(text_y[above] < 0.5))
  expect_true(all(text_y[!above] > 0.5))

  # And it is anchored so it grows away from the line, not across it.
  expect_equal(text$vjust, ifelse(above, 1, 0))
})

test_that("the printed value is the residual to two places", {

  aev <- mixed_aev()
  text <- residual_children(autoplot(aev))[[3]]
  z <- aev$deviance_residual[is.finite(aev$deviance_residual)]

  expect_identical(
    vapply(text$label, deparse, character(1)),
    vapply(aev_signed_math(z, 2), deparse, character(1))
  )
})

test_that("a residual off the end of the panel is clipped, not dropped", {

  # Sample B has one at +5.51: its line runs off the top, its marker is not
  # visible, and its value is still printed. All three records stay in the grob.
  loud <- create_aev(A = 1500, E = 1000, V = 100)
  expect_gt(loud$deviance_residual, aev_residual$limit)

  children <- residual_children(autoplot(loud))
  expect_identical(
    deparse(children[[3]]$label[[1]]),
    deparse(aev_signed_math(loud$deviance_residual, 2)[[1]])
  )

  render <- function(p) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }
  expect_no_warning(render(autoplot(loud)))
  expect_no_warning(render(autoplot(mixed_aev())))
})

#### The tint bands ####

test_that("the bands are the four the design has, mirrored about the line", {

  bands <- aev_residual_bands(factor(aev_panels[2], levels = aev_panels))

  expect_equal(bands$lower, c(1, 2, -2, -3))
  expect_equal(bands$upper, c(2, 3, -1, -2))

  # The further band is the darker one, on both sides.
  expect_equal(bands$alpha, c(0.1, 0.2, 0.1, 0.2))
  expect_true(all(bands$alpha[c(2, 4)] > bands$alpha[c(1, 3)]))
})

test_that("nothing is tinted between minus one and plus one", {

  bands <- aev_residual_bands(factor(aev_panels[2], levels = aev_panels))

  expect_false(any(bands$lower < 1 & bands$upper > -1))
})

test_that("the bands are drawn behind everything else in the panel", {

  # They are a background, and the first layer is the one at the back.
  expect_identical(layer_geoms(autoplot(mixed_aev()))[[1]], "GeomRect")
})
