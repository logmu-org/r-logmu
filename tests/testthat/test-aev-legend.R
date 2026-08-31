# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The legend takes the marker apart. Three entries name its three pieces and
# the fourth names the panel below, and each key is drawn at the size that
# piece is drawn on the chart, so the key is a sample rather than a diagram.
#
# THE WHOLE THING HANGS OFF ONE LAYER. ggplot2 draws every layer's key into
# every entry, so four layers would give four entries each showing all four
# glyphs stacked up. A single layer whose `draw_key` reads which entry it has
# been handed is the way round it. That is what most of this file is about.

legend_aev <- function() {
  aev <- create_aev(A = c(1100, 970, 1010), E = rep(1000, 3), V = c(2500, 1387, 900))
  names(aev) <- c("one", "two", "three")
  aev
}

layer_geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))

built <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplotGrob(p)
}

# Every key glyph, resolved at a known size so the pixel geometry means
# something.
key_grob <- function(entry) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(aev_legend$key_width_px / 96, "inches"),
    height = grid::unit(aev_legend$key_height_px / 96, "inches")
  ))
  aev_draw_key(data.frame(shape = entry), NULL, NULL)
}

#### What the entries say ####

test_that("the four entries are the design's, in the design's order", {

  # A LIST, because the sigma entry is plotmath and the rest are strings. Sigma
  # cannot be drawn from a string on a latin1 device -- see `aev_signed_math()`.
  entries <- aev_legend_entries(TRUE)

  expect_type(entries, "list")
  expect_identical(entries[[1]], "log A/E")
  expect_identical(entries[[2]], "95% confidence")
  expect_identical(entries[[4]], "Deviance residual")

  expect_true(is.expression(entries[[3]]))
  expect_identical(
    deparse(entries[[3]][[1]]),
    "paste(pm, sigma, \" (~68%) confidence\")"
  )
})

test_that("the residual entry goes when the residual panel does", {

  expect_length(aev_legend_entries(FALSE), 3L)
  expect_false("Deviance residual" %in% aev_legend_entries(FALSE))
})

test_that("the gloss on the residual is not repeated here", {

  # It lives in the panel title, against the thing it explains. Having it in
  # both would be saying the same thing twice on one chart.
  expect_false(any(grepl("95%", aev_legend_entries(TRUE)[4], fixed = TRUE)))
  expect_true(grepl("95%", aev_panels[2], fixed = TRUE))
})

#### One layer, four glyphs ####

test_that("the legend comes from a single layer that draws nothing", {

  p <- autoplot(legend_aev())
  points <- which(layer_geoms(p) == "GeomPoint")

  expect_length(points, 1L)

  # Its data is entirely missing, which is why nothing appears in the panel.
  data <- p$layers[[points]]$data
  expect_true(all(is.na(data$x)))
  expect_true(all(is.na(data$y)))
  expect_identical(nrow(data), 4L)
})

test_that("the key is chosen by the entry, not by the layer", {

  # Four layers would each draw into all four keys. This is the alternative,
  # and it only works because `draw_key` can see which entry it has.
  glyphs <- lapply(1:4, key_grob)
  kinds <- lapply(glyphs, function(g) {
    sort(unname(vapply(g$children, function(k) class(k)[1], character(1))))
  })

  expect_identical(kinds[[1]], c("circle", "segments"))
  expect_identical(kinds[[2]], c("segments", "segments"))
  expect_identical(kinds[[3]], c("segments", "segments"))
  expect_identical(kinds[[4]], c("circle", "segments"))

  # And they are not all the same glyph.
  expect_false(identical(glyphs[[2]], glyphs[[3]]))
})

test_that("each glyph is the piece of the marker it names", {

  in_pixels <- function(length) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::convertWidth(length, "inches", valueOnly = TRUE) * 96
  }

  # The point: a five-pixel circle in the weak grey, as on the chart.
  marker <- key_grob(1)
  expect_equal(in_pixels(marker$children[[2]]$r), aev_geometry$marker_radius_px)
  expect_identical(marker$children[[2]]$gp$col, aev_palette$marker_weak)

  # The 95% cap: eighteen pixels of the strong grey, over an arm of the medium.
  confidence <- key_grob(2)
  expect_identical(confidence$children[[1]]$gp$col, aev_palette$marker_medium)
  expect_identical(confidence$children[[2]]$gp$col, aev_palette$marker_cap)
  expect_equal(
    in_pixels(confidence$children[[2]]$x1[1] - confidence$children[[2]]$x0[1]),
    aev_geometry$cap_width_px
  )

  # The sigma entry: the heavier arm over the lighter line, which is the
  # transition it names.
  sigma <- key_grob(3)
  expect_identical(sigma$children[[1]]$gp$col, aev_palette$marker_medium)
  expect_identical(sigma$children[[2]]$gp$col, aev_palette$marker_weak)
  expect_gt(sigma$children[[1]]$gp$lwd, sigma$children[[2]]$gp$lwd)

  # The residual: the smaller circle, at the residual panel's own size.
  residual <- key_grob(4)
  expect_equal(in_pixels(residual$children[[2]]$r), aev_residual$marker_radius_px)
  expect_identical(residual$children[[2]]$gp$col, aev_palette$marker_medium)
})

test_that("the two circles differ, as they do on the chart", {

  expect_gt(aev_geometry$marker_radius_px, aev_residual$marker_radius_px)
})

#### Where it sits ####

test_that("the legend is above the panels and below the title", {

  table <- built(autoplot(legend_aev(), title = "A title"))
  row_of <- function(prefix) table$layout$t[startsWith(table$layout$name, prefix)]

  legend <- row_of("guide-box-top")
  panels <- row_of("panel")

  expect_length(legend, 1L)
  expect_true(all(legend < panels))
  expect_true(all(row_of("title") < legend))
})

test_that("the legend is left-aligned with the y axis, not with the plot edge", {

  # "panel" rather than "plot": Tim wants its left edge on the y axis, where
  # "plot" would put it out at the far left with the title.
  theme <- theme_aev()

  expect_identical(theme$legend.position, "top")
  expect_identical(theme$legend.justification, "left")
  expect_identical(theme$legend.location, "panel")
})

test_that("the key is small and faded, in a light box", {

  theme <- theme_aev()

  expect_identical(theme$legend.text$colour, aev_palette$marker_medium)
  expect_s3_class(theme$legend.text$size, "rel")
  expect_lt(unclass(theme$legend.text$size), 1)

  # The box the panels no longer have. Tim asked for both changes together.
  expect_s3_class(theme$legend.background, "element_rect")
  expect_identical(theme$legend.background$colour, aev_palette$grid_line)
  expect_true(is.na(theme$legend.background$fill))
  expect_s3_class(theme$panel.border, "element_blank")

  expect_s3_class(theme$legend.key, "element_blank")
  expect_s3_class(theme$legend.title, "element_blank")
})

test_that("the key is wide enough for the widest glyph", {

  # The cap is the widest thing in it, and a key narrower than the cap would
  # crop the entry it is meant to be a sample of.
  expect_gt(aev_legend$key_width_px, aev_geometry$cap_width_px)
  expect_gt(aev_legend$key_height_px, 2 * aev_geometry$marker_radius_px)
})

#### It does not disturb the chart ####

test_that("the legend layer contributes nothing to either axis", {

  # Its data is all missing, so the ranges are the ones the skeleton sets.
  with_legend <- ggplot2::ggplot_build(autoplot(legend_aev()))$layout$panel_params

  expect_equal(with_legend[[1]]$y.range, aev_log_ratio_limits)
  expect_equal(with_legend[[1]]$x.range, c(0.5, 3.5))
})

test_that("the legend layer joins the runs rather than reordering them", {

  # THE REGRESSION THIS REPLACES: the layer passed its group as a character,
  # which made the combined facet variable a character, and the columns came
  # out alphabetically instead of in record order. The chart splits on the run
  # a record is in rather than on its group name, so that is what the layer has
  # to carry, and it has to carry it as the factor.
  aev <- create_aev(A = rep(1000, 3), E = rep(1000, 3), V = rep(900, 3))
  names(aev) <- c("p", "a", "c")
  group_names(aev) <- c("Period", "Age", "Cohort")

  p <- autoplot(aev)
  legend_layer <- p$layers[[which(layer_geoms(p) == "GeomPoint")]]

  expect_s3_class(legend_layer$data$run, "factor")
  expect_identical(levels(legend_layer$data$run), c("1", "2", "3"))
})

test_that("drawing the chart with a legend emits no warning", {

  render <- function(p) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    print(p)
  }

  # `na.rm = TRUE` on the legend layer: without it the all-missing data draws a
  # "removed 4 rows" warning on every chart.
  expect_no_warning(render(autoplot(legend_aev())))
  expect_no_warning(render(autoplot(legend_aev(), residuals = FALSE)))
})
