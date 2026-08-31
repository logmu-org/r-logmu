# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `log_range` and `log_step` are the only way to change the A/E axis, and they
# exist because the obvious alternative is a trap.
#
# ADDING `+ scale_y_continuous(...)` REPLACES THE WHOLE SCALE, taking with it
# the per-panel limits function and the `oob` passthrough that the two panels
# depend on. Measured before these arguments existed: the residual panel's range
# went from [-3.5, 3.5] to [-4.16, 10.30] and its labels became raw logarithms.
# Silently wrong, with no error. The last test in this file holds that line.

scale_aev <- function() {
  aev <- create_aev(A = c(1218, 1075, 1001), E = rep(1000, 3), V = c(876, 59, 2656))
  names(aev) <- c("a", "b", "c")
  aev
}

params <- function(p, panel = 1L) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplot_build(p)$layout$panel_params[[panel]]
}

#### The breaks that come out ####

test_that("the default axis is unchanged by the arguments existing", {

  expect_equal(aev_log_breaks(c(-0.5, 0.5), 0.1), seq(-0.5, 0.5, by = 0.1))
  expect_equal(params(autoplot(scale_aev()))$y.range, c(-0.5, 0.5))
  expect_identical(
    params(autoplot(scale_aev()))$y$get_labels(),
    c("61%", "67%", "74%", "82%", "90%", "100%",
      "111%", "122%", "135%", "149%", "165%")
  )
})

test_that("a wider range and a coarser step both take effect", {

  # Slides 6 and 7 ran 50% to 200%, which is a log range of plus and minus
  # log(2), in eighths of an octave.
  octave <- log(2)
  p <- autoplot(scale_aev(), log_range = c(-octave, octave), log_step = octave / 4)

  expect_equal(params(p)$y.range, c(-octave, octave))
  expect_identical(
    params(p)$y$get_labels(),
    c("50%", "59%", "71%", "84%", "100%", "119%", "141%", "168%", "200%")
  )
})

test_that("a coarser step gives fewer labels", {

  # 0.2 does not divide 0.5, so twenty per cent steps need a range that it does
  # divide. The error message says as much when they do not.
  fine <- autoplot(scale_aev())
  coarse <- autoplot(scale_aev(), log_range = c(-0.6, 0.6), log_step = 0.2)

  expect_length(params(fine)$y$get_breaks(), 11L)
  expect_length(params(coarse)$y$get_breaks(), 7L)
  expect_identical(
    params(coarse)$y$get_labels(),
    c("55%", "67%", "82%", "100%", "122%", "149%", "182%")
  )
})

test_that("an asymmetric range is allowed and still puts a break on 100%", {

  # Slide 3's older default was 60% to 160%, which is not symmetric in logs.
  p <- autoplot(scale_aev(), log_range = c(-0.5, 0.4), log_step = 0.1)

  expect_equal(params(p)$y.range, c(-0.5, 0.4))
  expect_true(0 %in% round(params(p)$y$get_breaks(), 10))
  expect_true("100%" %in% params(p)$y$get_labels())
})

test_that("the ends of the range are always labelled ticks", {

  # The axis has no expansion, so its ends are its outermost grid lines. A step
  # that fell short of them would leave the panel open at the top.
  for (step in c(0.1, 0.25, 0.5)) {
    breaks <- aev_log_breaks(c(-0.5, 0.5), step)
    expect_equal(range(breaks), c(-0.5, 0.5))
  }
})

#### The residual panel is untouched ####

test_that("changing the A/E axis leaves the residual axis alone", {

  # The whole point of the bound on `log_range`: the two panels are told apart
  # by how far they reach, and the residual panel must keep its own scale.
  p <- autoplot(scale_aev(), log_range = c(-0.9, 0.9), log_step = 0.3)

  expect_equal(params(p, 2L)$y.range, c(-3.5, 3.5))
  expect_equal(params(p, 2L)$y$get_breaks(), -3:3)
  expect_identical(
    vapply(params(p, 2L)$y$get_labels(), deparse, character(1)),
    c("-\"3\"", "-\"2\"", "-\"1\"", "\"0\"", "+\"1\"", "+\"2\"", "+\"3\"")
  )
})

test_that("the labels are chosen by how far the axis reaches, not by chance", {

  # The dispatch used to be `breaks %in% -3:3`, which an A/E axis could satisfy
  # by having whole-number breaks of its own.
  expect_identical(aev_panel_labels(c(-1, 0, 1)), c("37%", "100%", "272%"))
  # The residual side is plotmath, so compare the construction rather than a
  # string. That difference IS the dispatch: percentages are text, signed
  # residuals are expressions.
  expect_identical(deparse(aev_panel_labels(-3:3)[[4]]), "\"0\"")
  expect_true(is.expression(aev_panel_labels(-3:3)))
  expect_type(aev_panel_labels(c(-1, 0, 1)), "character")
})

#### What is refused ####

test_that("a range beyond the residual panel's is refused, with the reason", {

  expect_error(autoplot(scale_aev(), log_range = c(-1, 1)), "must lie within")
  expect_error(autoplot(scale_aev(), log_range = c(-2, 2)), "could not be told")
})

test_that("a step that does not divide the range is refused", {

  # 0.3 into 0.5 leaves the top of the axis without a tick, and would also miss
  # 100% on an asymmetric range.
  expect_error(autoplot(scale_aev(), log_step = 0.3), "must divide both ends")
  expect_no_error(autoplot(scale_aev(), log_step = 0.25))

  # And it says which range would have worked, because the commonest request
  # of all -- twenty per cent steps on the default range -- is one of these.
  expect_error(autoplot(scale_aev(), log_step = 0.2), "c(-0.6, 0.6) would work",
               fixed = TRUE)
})

test_that("a malformed range or step is refused", {

  expect_error(autoplot(scale_aev(), log_range = 0.5), "two increasing numbers")
  expect_error(autoplot(scale_aev(), log_range = c(0.5, -0.5)), "two increasing numbers")
  expect_error(autoplot(scale_aev(), log_range = c(NA, 0.5)), "two increasing numbers")
  expect_error(autoplot(scale_aev(), log_step = 0), "single positive number")
  expect_error(autoplot(scale_aev(), log_step = c(0.1, 0.2)), "single positive number")
})

#### Off scale follows the axis ####

test_that("what counts as off scale moves with the range", {

  # A/E of 150% is on a wide axis and off a narrow one, and the marker has to
  # agree with the axis it is drawn against.
  aev <- create_aev(A = 1500, E = 1000, V = 100)
  drawable <- function(p) {
    nrow(p$layers[[which(vapply(p$layers, function(l) class(l$geom)[1], "") == "GeomAevInterval")]]$data)
  }

  expect_identical(drawable(autoplot(aev, log_range = c(-0.6, 0.6), log_step = 0.2)), 1L)
  expect_identical(drawable(autoplot(aev, log_range = c(-0.2, 0.2), log_step = 0.1)), 0L)
})

#### The trap these arguments replace ####

test_that("replacing the y scale by hand still breaks the residual panel", {

  # Not a defect to fix -- a ggplot cannot stop a user replacing its scale. It
  # is why the arguments exist, and this records what happens without them so
  # the reason survives.
  hand_rolled <- suppressMessages(
    autoplot(scale_aev()) + ggplot2::scale_y_continuous(breaks = log(c(0.9, 1, 1.1)))
  )

  expect_false(isTRUE(all.equal(params(hand_rolled, 2L)$y.range, c(-3.5, 3.5))))
})


#### A drawable record that sits outside the axis ####

test_that("a record drawn from outside the axis does not drag the panel onto the residual scale", {

  # THE ONE THE REVIEW FOUND, 2026-08-26. A record is drawable when its inner
  # interval reaches the axis, which its centre need not do. Its centre used to
  # train the scale, taking the A/E panel's range past 1 -- and the two panels
  # are told apart by how far their ranges reach, so the A/E panel was read as
  # the residual one and drawn on the deviance axis, labels and all. No error,
  # no warning, and a confidently wrong chart.
  #
  # A/E of 300% with a sigma of 0.61: log A/E is 1.099, well outside the axis,
  # and the inner interval reaches back to 0.49 which is inside it.
  aev <- create_aev(A = c(1100, 12), E = c(1000, 4), V = c(2500, 6))

  expect_identical(as.character(aev_status(aev)), c("ok", "ok"))

  p <- autoplot(aev)

  expect_equal(params(p, 1L)$y.range, c(-0.5, 0.5))
  expect_true("100%" %in% params(p, 1L)$y$get_labels())
  expect_equal(params(p, 2L)$y.range, c(-3.5, 3.5))
})

test_that("the marker is still drawn where its data says, not where the scale was trained", {

  # The clamp is for the scale's benefit only. If it ever reached the drawing,
  # every off-axis marker would pile up on the panel edge instead of running
  # off it, which would look almost right and be wrong.
  aev <- create_aev(A = c(1100, 12), E = c(1000, 4), V = c(2500, 6))
  frame <- aev_plot_frame(aev, c(-0.5, 0.5))

  expect_equal(frame$log_ratio[2], log(3), tolerance = 1e-9)
  expect_equal(frame$clamped[2], 0.5)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  grob <- ggplot2::layer_grob(
    autoplot(aev),
    which(vapply(autoplot(aev)$layers, function(l) class(l$geom)[1], "") == "GeomAevInterval")
  )[[1]]

  # npc of log(3) on an axis running -0.5 to 0.5, which is well above the top.
  expect_gt(grob$y[2], 1)
})
