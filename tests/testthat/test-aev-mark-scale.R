# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Marks sized against the chart's type rather than against the device.
#
# WHY THESE TESTS MEASURE THE FORCED GROB. Everything else about the marker can
# be checked by asking a layer for its grob, but the scale is not in the layer:
# the facet reads it off the finished theme and stamps it on while the table is
# assembled. So a test that stops short of assembling the table sees a scale of
# one whatever the theme says, which is exactly the failure it is meant to
# catch.
#
# `lwd` is what they assert on, because it is a plain number that means the same
# thing in any viewport. The lengths are in npc by the time `grid` has them, and
# npc outside the panel means nothing.

two_records <- function() {
  create_aev(A = c(1100, 970), E = c(1000, 1000), V = c(2500, 1600))
}

# Every line width inside the resolved marker, in the order the design draws
# them: circle, thin, thick, caps.
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

#### The scale itself ####

test_that("no theme and a broken theme both mean no scaling", {

  expect_equal(aev_mark_scale(NULL), 1)
  expect_equal(aev_mark_scale(ggplot2::theme(text = ggplot2::element_text(size = 0))), 1)
})

test_that("the scale is the theme's text size against the reference", {

  expect_equal(aev_reference_text_size, 9)
  expect_equal(aev_mark_scale(theme_aev()), 1)
  expect_equal(aev_mark_scale(theme_aev(base_size = 18)), 2)
  expect_equal(aev_mark_scale(theme_aev(base_size = 4.5)), 0.5)
})

test_that("scaling names the fields it touches and leaves the rest alone", {

  # `breakpoints` are residuals, not lengths. Scaling this list wholesale would
  # quietly move where the ramp saturates when someone enlarged the type.
  scaled <- aev_scaled(aev_residual, 2, c("marker_radius_px", "text_px"))

  expect_equal(scaled$marker_radius_px, aev_residual$marker_radius_px * 2)
  expect_equal(scaled$text_px, aev_residual$text_px * 2)
  expect_identical(scaled$breakpoints, aev_residual$breakpoints)
  expect_identical(scaled$limit, aev_residual$limit)
})

test_that("a scale of one is the identity, so adopting this changed nothing", {

  # THE NO-OP PROPERTY, which is the whole reason this was safe to adopt: at the
  # reference size every constant has to come back untouched, not merely close.
  expect_identical(aev_scaled(aev_geometry, 1, names(aev_geometry)), aev_geometry)
})

#### Reaching the marks ####

test_that("the default chart draws the marker at its stated pixel widths", {

  widths <- marker_widths()

  expect_equal(unname(widths[2]), aev_geometry$thin_width_px)
  expect_equal(unname(widths[3]), aev_geometry$thick_width_px)
  expect_equal(unname(widths[4]), aev_geometry$cap_height_px)
})

test_that("larger type gives a heavier marker, in proportion", {

  # The one that fails if the facet stops stamping the scale, or if the marks
  # stop reading it. Nothing else in the suite would notice.
  doubled <- marker_widths(size = aev_reference_text_size * 2)

  expect_equal(unname(doubled[2]), aev_geometry$thin_width_px * 2)
  expect_equal(unname(doubled[3]), aev_geometry$thick_width_px * 2)
  expect_equal(unname(doubled[4]), aev_geometry$cap_height_px * 2)
})

test_that("smaller type gives a lighter marker", {

  halved <- marker_widths(size = aev_reference_text_size / 2)

  expect_equal(unname(halved[3]), aev_geometry$thick_width_px / 2)
})

test_that("the residual panel scales with the same factor", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot <- autoplot(two_records()) +
    ggplot2::theme(text = ggplot2::element_text(size = aev_reference_text_size * 2))

  forced <- grid::forceGrob(ggplot2::ggplotGrob(plot))
  marker <- Filter(
    function(g) inherits(g, "aev_residual_marker"), grobs_in(forced)
  )[[1]]

  # The printed value is type and scales too, or a doubled chart would carry
  # residuals in the original 8 px against everything else at twice the size.
  label <- Filter(function(g) inherits(g, "text"), marker$children)[[1]]

  expect_equal(
    label$gp$fontsize, aev_residual$text_px * pt_per_px * 2
  )
})

#### The argument list the scale arrives through ####

test_that("a positional call still finds the theme", {

  # ggplot2 calls `draw_panels` with no argument names at all, so this is the
  # only thing standing between the scale and silently always being one.
  theme <- theme_aev(base_size = 18)
  found <- aev_facet_arguments(
    list(), data.frame(), list(), list(), list(), NULL, list(), theme, list()
  )

  expect_identical(found$theme, theme)
  expect_equal(aev_mark_scale(found$theme), 2)
})

test_that("the argument names are read from ggplot2 rather than written down", {

  # If a future ggplot2 reorders `draw_panels`, this follows it. The test is
  # that the names come from that function and not from a list in our source.
  inner <- environment(ggplot2::FacetGrid$draw_panels)$f
  expected <- setdiff(names(formals(inner %||% ggplot2::FacetGrid$draw_panels)), "self")

  expect_true("theme" %in% expected)
  expect_identical(
    names(aev_facet_arguments(1, 2, 3)),
    expected[1:3]
  )
})

#### The layout follows too ####

# The depth of the residual strip, in pixels, as the table actually carries it.
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

test_that("the residual panel is 84 px at the reference size", {
  expect_equal(residual_depth(), aev_residual$height_px)
})

test_that("the residual panel's depth follows the type", {

  # WITHOUT THIS THE STRIP STAYS 84 PX while its contents grow, so the residual
  # axis tightens up as the A/E panel above it opens out. That is what a fixed
  # depth looked like, and it is why the heights are set by the facet rather
  # than by `theme(panel.heights = )`, which ggplot2 applies afterwards and
  # would overwrite.
  expect_equal(
    residual_depth(size = aev_reference_text_size * 2),
    aev_residual$height_px * 2
  )
  expect_equal(
    residual_depth(size = aev_reference_text_size * 1.5),
    aev_residual$height_px * 1.5
  )
})

test_that("a chart with no residual panel has one stretchy row", {

  # CONVERTING A UNIT NEEDS A DEVICE even when the unit is absolute. Without
  # one R opens the default and leaves an `Rplots.pdf` in `tests/testthat` --
  # which is exactly what happened here, and it reached a commit.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_identical(as.character(grid::unitType(aev_panel_heights(FALSE, 2))), "null")
  expect_equal(
    grid::convertHeight(aev_panel_heights(TRUE, 2)[2], "inches", valueOnly = TRUE) * 96,
    aev_residual$height_px * 2
  )
})

#### The legend and the bracket ####

# The legend is the one part the facet never sees: it is added to the table
# after `draw_panels` has returned. So the key is sized against its own box,
# which `theme_aev()` sizes against the type.

test_that("the key glyph is sized against the box it is drawn in", {

  reference <- aev_legend$key_width_px / 96 * 25.4

  expect_equal(aev_key_scale(c(reference, 1)), 1)
  expect_equal(aev_key_scale(c(reference * 2, 1)), 2)

  # A box of no size, or none at all, must not divide by zero or scale by NA.
  expect_equal(aev_key_scale(numeric(0)), 1)
  expect_equal(aev_key_scale(c(0, 1)), 1)
  expect_equal(aev_key_scale(c(NA_real_, 1)), 1)
})

test_that("a doubled key box gives a doubled glyph", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  single <- aev_key_confidence(1)
  double <- aev_key_confidence(2)

  expect_equal(double$children[[1]]$gp$lwd, single$children[[1]]$gp$lwd * 2)
  expect_equal(double$children[[2]]$gp$lwd, single$children[[2]]$gp$lwd * 2)

  marker <- aev_key_marker(2)
  expect_equal(marker$children[[2]]$gp$lwd, aev_geometry$thin_width_px * 2)

  residual <- aev_key_residual(2)
  expect_equal(residual$children[[2]]$gp$lwd, aev_residual$marker_width_px * 2)
})

test_that("the key box follows the theme's type size", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  in_pixels <- function(unit) {
    grid::convertWidth(unit, "inches", valueOnly = TRUE) * 96
  }

  expect_equal(
    in_pixels(theme_aev()$legend.key.width), aev_legend$key_width_px
  )
  expect_equal(
    in_pixels(theme_aev(base_size = aev_reference_text_size * 2)$legend.key.width),
    aev_legend$key_width_px * 2
  )
})

# The bracket is drawn by the facet, so it gets the scale directly.

test_that("the group bracket's rules scale with the type", {

  rule_width <- function(scale) {

    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      width = grid::unit(2, "inches"), height = grid::unit(0.3, "inches")
    ))

    content <- grid::makeContent(
      aev_group_bracket(grid::textGrob("group"), scale)
    )
    content$children[[1]]$gp$lwd
  }

  expect_equal(rule_width(1), aev_geometry$grid_width_px)
  expect_equal(rule_width(2), aev_geometry$grid_width_px * 2)
})

test_that("a chart at larger type carries the scale into its brackets", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  aev <- create_aev(A = c(1100, 970), E = c(1000, 1000), V = c(2500, 1600))
  group_names(aev) <- c("one", "two")

  table <- ggplot2::ggplotGrob(
    autoplot(aev) +
      ggplot2::theme(text = ggplot2::element_text(size = aev_reference_text_size * 2))
  )

  brackets <- table$grobs[startsWith(table$layout$name, "aev-group-")]

  expect_length(brackets, 2L)
  expect_true(all(vapply(brackets, function(g) g$scale, numeric(1)) == 2))
})

test_that("the key function passes its own box size through to the glyph", {

  # THE WIRING, not the arithmetic. Testing `aev_key_confidence(2)` directly
  # leaves `aev_draw_key()` free to drop the scale on the floor and every
  # assertion above still passes -- checked by breaking exactly that.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  reference <- aev_legend$key_width_px / 96 * 25.4
  confidence <- data.frame(shape = 2L)

  single <- aev_draw_key(confidence, list(), c(reference, 1))
  double <- aev_draw_key(confidence, list(), c(reference * 2, 1))

  expect_equal(double$children[[1]]$gp$lwd, single$children[[1]]$gp$lwd * 2)
  expect_equal(double$children[[2]]$gp$lwd, single$children[[2]]$gp$lwd * 2)
})

#### The exported way to resize ####

test_that("adding the theme back larger scales the marks and the legend key", {

  # THE DIFFERENCE BETWEEN THE TWO ROUTES, which the documentation promises.
  # Both move the marks, because the facet reads the finished theme's text size.
  # Only `theme_aev()` moves the legend key, because the key is sized by
  # `legend.key.width` and nothing else sets it.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  aev <- create_aev(A = c(1100, 970), E = c(1000, 1000), V = c(2500, 1600))
  larger <- aev_reference_text_size * 2

  arm <- function(plot) {
    forced <- grid::forceGrob(ggplot2::ggplotGrob(plot))
    marker <- Filter(function(g) inherits(g, "aev_interval"), grobs_in(forced))[[1]]
    marker$children[[3]]$gp$lwd
  }

  key <- function(plot) {
    grid::convertWidth(
      ggplot2::calc_element("legend.key.width", plot$theme), "inches",
      valueOnly = TRUE
    ) * 96
  }

  plain <- autoplot(aev)
  by_text <- plain + ggplot2::theme(text = ggplot2::element_text(size = larger))
  by_theme <- plain + theme_aev(base_size = larger)

  expect_equal(arm(plain), aev_geometry$thick_width_px)
  expect_equal(arm(by_text), aev_geometry$thick_width_px * 2)
  expect_equal(arm(by_theme), aev_geometry$thick_width_px * 2)

  expect_equal(key(plain), aev_legend$key_width_px)
  expect_equal(key(by_text), aev_legend$key_width_px)
  expect_equal(key(by_theme), aev_legend$key_width_px * 2)
})

test_that("the theme is still a theme that can be extended", {

  # Exporting it makes it part of the interface, so `+` has to keep working on
  # the result rather than only on the chart.
  extended <- theme_aev() + ggplot2::theme(plot.title = ggplot2::element_blank())

  expect_s3_class(extended, "theme")
  expect_s3_class(
    autoplot(create_aev(A = 1100, E = 1000, V = 2500)) + theme_aev(base_size = 11),
    "ggplot"
  )
})

#### What a standard theme does to the chart ####

test_that("a standard theme keeps the structure and loses the styling", {

  # The documented behaviour, measured rather than asserted. Nothing errors and
  # nothing warns -- a ggplot cannot stop its theme being replaced -- so the
  # docs promise exactly this and the promise is checked here.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  aev <- create_aev(A = c(1100, 970, 1050), E = rep(1000, 3), V = c(2500, 1600, 900))
  group_names(aev) <- c("one", "one", "two")

  standard <- autoplot(aev) + ggplot2::theme_minimal()

  # It warns, and that is the subject of its own tests below. Here it is noise.
  table <- suppressWarnings(ggplot2::ggplotGrob(standard))

  # Structure survives: two panels, both brackets, the residual strip.
  expect_length(aev_panel_rows(table), 2L)
  expect_length(which(startsWith(table$layout$name, "aev-group-")), 2L)

  # Styling does not.
  expect_false(inherits(
    ggplot2::calc_element("panel.grid.major", standard$theme), "element_blank"
  ))
  expect_equal(ggplot2::calc_element("strip.text.y.left", standard$theme)$angle, 90)
  expect_identical(ggplot2::calc_element("strip.placement", standard$theme), "inside")

  # And the type it brings takes the marks with it.
  expect_equal(aev_mark_scale(standard$theme), 11 / aev_reference_text_size)
})

test_that("adding theme_aev after a standard theme puts the styling back", {

  restored <- autoplot(create_aev(A = 1100, E = 1000, V = 2500)) +
    ggplot2::theme_minimal() +
    theme_aev()

  expect_true(inherits(
    ggplot2::calc_element("panel.grid.major", restored$theme), "element_blank"
  ))
  expect_equal(ggplot2::calc_element("strip.text.y.left", restored$theme)$angle, 0)
  expect_identical(ggplot2::calc_element("strip.placement", restored$theme), "outside")
  expect_equal(aev_mark_scale(restored$theme), 1)
})


#### The chart notices when its theme has been replaced ####

test_that("the signature survives a merge and is lost to a replacement", {

  # The whole mechanism in four lines. `+ theme(...)` MERGES, which is how
  # anybody customises a chart and must never be flagged; a complete theme
  # REPLACES, which is the thing worth saying something about.
  expect_true(isTRUE(ggplot2::calc_element(aev_signature, theme_aev())))
  expect_true(isTRUE(ggplot2::calc_element(
    aev_signature,
    theme_aev() + ggplot2::theme(text = ggplot2::element_text(size = 18))
  )))
  expect_null(ggplot2::calc_element(
    aev_signature, theme_aev() + ggplot2::theme_minimal()
  ))
  expect_null(ggplot2::calc_element(aev_signature, ggplot2::theme_minimal()))
})

test_that("the element is registered, which is what makes the signature legible", {

  # Registered in `.onLoad`. Without it `calc_element()` errors on an unknown
  # element and the check would decide it cannot tell.
  expect_true(aev_signature %in% names(ggplot2::get_element_tree()))
})

test_that("a replaced theme warns and an extended one does not", {

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

  for (label in names(quiet)) {
    expect_no_warning(ggplot2::ggplotGrob(quiet[[label]]))
  }

  for (standard in list(ggplot2::theme_minimal(), ggplot2::theme_bw(), ggplot2::theme_grey())) {
    expect_warning(
      ggplot2::ggplotGrob(autoplot(aev) + standard),
      "theme has been replaced"
    )
  }
})

test_that("the warning names the fix", {

  # A warning that does not say what to do is a warning people learn to ignore.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  aev <- create_aev(A = 1100, E = 1000, V = 2500)

  expect_warning(
    ggplot2::ggplotGrob(autoplot(aev) + ggplot2::theme_minimal()),
    "theme_aev()", fixed = TRUE
  )
})

test_that("a signature that cannot be read means cannot tell, not replaced", {

  # If the signature can ever not be read -- a load that half failed, or a
  # `reset_theme_settings()` that dropped the registration -- the chart should
  # say nothing rather than warn about a theme that is perfectly fine. Anything
  # `calc_element()` refuses has to read as "cannot tell".
  #
  # NOT an empty theme: the element is registered with a default of FALSE, so
  # `calc_element()` answers FALSE for a theme that simply lacks it, which is
  # the replaced case and not this one.
  expect_false(aev_theme_replaced("not a theme at all"))
  expect_false(aev_theme_replaced(NULL))
})
