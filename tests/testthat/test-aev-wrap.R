# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Wrapping the record names under the A/E panel.
#
# THE WORK IS SPLIT IN TWO AND THE SPLIT IS THE POINT. ggplot2 sizes the axis
# row before the panel width is known and draws it afterwards, so the number of
# lines is settled at build time against an assumed width, and where the breaks
# fall is settled at draw time against the real one. Both halves are tested
# here, and so is the fact that they agree on the line count.

# A device, because measuring a string needs one even when the answer is in
# inches.
with_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

axis_text <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(p)
  cell <- which(table$layout$name == "axis-b-1-1")
  aev_text_within(table$grobs[[cell]])
}

# What is actually drawn, which is not what the scale was given.
#
# READ OUT OF THE FORCED WRAPPER, not by searching the whole tree for text:
# `forceGrob` keeps the unresolved original beside the resolved copy, so a
# search finds the build-time labels as readily as the drawn ones.
drawn_labels <- function(p, width = 8, height = 7) {

  grDevices::pdf(NULL, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(p))

  wrappers <- Filter(
    function(grob) inherits(grob, "forcedgrob") && inherits(grob, "aev_axis_labels"),
    grobs_in(forced)
  )

  if (length(wrappers) == 0L) return(NULL)

  as.character(aev_text_within(wrappers[[1]])$label)
}

#### Packing ####

test_that("words are taken until the next one would not fit", {

  expect_identical(aev_pack(c("aa", "bb", "cc"), 5, nchar), c("aa bb", "cc"))
  expect_identical(aev_pack(c("aa", "bb", "cc"), 8, nchar), c("aa bb cc"))
  # A word longer than the limit takes a line of its own rather than vanishing.
  expect_identical(aev_pack(c("aaaaaaaa", "bb"), 3, nchar), c("aaaaaaaa", "bb"))
})

test_that("blank space between words is not a word", {

  expect_identical(aev_words("  male   retirees "), c("male", "retirees"))
  expect_identical(aev_words("widows"), "widows")
})

#### Breaking into a fixed number of lines ####

test_that("one line, or one word, is left alone", {

  expect_identical(aev_wrap_to_lines("male retirees", 1L, nchar), "male retirees")
  expect_identical(aev_wrap_to_lines("widows", 2L, nchar), "widows")
})

test_that("the break is the most even one that fits", {

  # Not first-fit at some arbitrary limit: the narrowest limit that still fits
  # in two lines, which is what makes the two lines the same sort of length.
  expect_identical(
    aev_wrap_to_lines("female dependants under 65", 2L, nchar),
    "female dependants\nunder 65"
  )
  expect_identical(
    aev_wrap_to_lines("a bb ccc dddd", 2L, nchar),
    "a bb ccc\ndddd"
  )
})

test_that("a label is never broken into more lines than it was given", {

  for (lines in 1:3) {
    broken <- aev_wrap_to_lines("one two three four five six", lines, nchar)
    expect_lte(length(strsplit(broken, "\n", fixed = TRUE)[[1]]), lines)
  }
})

#### How many lines to reserve ####

test_that("the line count falls out of the assumed width and the record count", {

  # 119 characters across the axis, measured. Two records get sixty each and
  # need no break; twelve get ten each and cannot hold the same label.
  expect_identical(aev_label_lines(rep("deferred pensioners", 2)), 1L)
  expect_gt(aev_label_lines(rep("deferred pensioners", 12)), 1L)
})

test_that("the line count is capped, and is never below one", {

  expect_identical(aev_label_lines(character(0)), 1L)
  expect_identical(aev_label_lines(rep("a b c d e f g h i j k", 40)), aev_axis_max_lines)
})

#### What the scale is given, and why ####

test_that("the labels the scale gets carry real newlines", {

  # This is what makes ggplot2 measure the axis row deep enough. Without it the
  # row is one line tall and the second line is drawn outside it.
  aev <- create_aev(A = rep(1000, 8), E = rep(1000, 8), V = rep(2500, 8))
  names(aev) <- rep("female dependants under 65", 8)

  expect_true(all(grepl("\n", axis_text(autoplot(aev))$label)))
})

test_that("short names in a small chart are not wrapped at all", {

  aev <- create_aev(A = c(1000, 1000), E = c(1000, 1000), V = c(2500, 2500))
  names(aev) <- c("male", "female")

  expect_false(any(grepl("\n", axis_text(autoplot(aev))$label)))

  # And no wrapper is put in the way when there is nothing to wrap.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(autoplot(aev))
  cell <- which(table$layout$name == "axis-b-1-1")
  expect_false(inherits(table$grobs[[cell]], "aev_axis_labels"))
})

#### The draw-time pass ####

test_that("the break is chosen by drawn width, not by counting characters", {

  # `i` is narrow and `m` is wide, so the two disagree. Counting characters
  # puts the break after the first word; measuring puts it after the second.
  label <- "iiiiiiiiii mmmm mmmm"

  expect_identical(aev_wrap_to_lines(label, 2L, nchar), "iiiiiiiiii\nmmmm mmmm")

  by_width <- with_device({
    width_of <- function(string) {
      grid::convertWidth(grid::stringWidth(string), "inches", valueOnly = TRUE)
    }
    aev_wrap_to_lines(label, 2L, width_of)
  })

  expect_identical(by_width, "iiiiiiiiii mmmm\nmmmm")
})

test_that("what is drawn is the width-based break, not the one built in", {

  aev <- create_aev(A = rep(1000, 8), E = rep(1000, 8), V = rep(2500, 8))
  names(aev) <- rep("iiiiiiiiii mmmm mmmm", 8)

  built <- axis_text(autoplot(aev))$label
  drawn <- drawn_labels(autoplot(aev))

  expect_true(all(built == "iiiiiiiiii\nmmmm mmmm"))
  expect_true(all(drawn == "iiiiiiiiii mmmm\nmmmm"))
})

test_that("the line count is the same however wide the chart is drawn", {

  # The height was reserved for that many lines and there is no second chance
  # to change it, so the draw-time pass must never want more.
  aev <- create_aev(A = rep(1000, 8), E = rep(1000, 8), V = rep(2500, 8))
  names(aev) <- rep("female dependants under 65", 8)

  reserved <- aev_label_lines(names(aev))

  widths <- c(5, 8, 14)
  counts <- vapply(widths, function(width) {
    drawn <- drawn_labels(autoplot(aev), width = width)
    max(vapply(strsplit(drawn, "\n", fixed = TRUE), length, integer(1)))
  }, integer(1))

  # The cap holds. On its own this says little -- `aev_wrap_to_width()` cannot
  # return more lines than it is given, so this half of the test is mostly
  # checking the function against its own construction.
  expect_true(all(counts <= reserved))

  # THIS is the part that is about the width. Packing at the column's real
  # width means more room can never need MORE lines, which is a property of the
  # packing and not of the ceiling, and would catch a draw-time pass that had
  # stopped reading the column at all.
  expect_true(all(diff(counts) <= 0L))
})

test_that("the wrapper resolves when the chart is drawn", {

  aev <- create_aev(A = rep(1000, 8), E = rep(1000, 8), V = rep(2500, 8))
  names(aev) <- rep("female dependants under 65", 8)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  forced <- grid::forceGrob(ggplot2::ggplotGrob(autoplot(aev)))
  wrappers <- Filter(
    function(grob) inherits(grob, "forcedgrob") && inherits(grob, "aev_axis_labels"),
    grobs_in(forced)
  )

  expect_length(wrappers, 1L)

  # NOT `length(children) > 0`, which the wrapper satisfies from construction:
  # it holds the axis it wraps. What only exists once `makeContent` has run is a
  # label actually broken across lines.
  expect_true(any(grepl("\n", labels_in(wrappers[[1]]), fixed = TRUE)))
})


#### A name that already fits is left alone ####

test_that("reserving two lines does not force two lines on every name", {

  # FOUND BY REVIEW, 2026-08-26. The draw-time pass broke each name into the
  # reserved number of lines regardless of the column, because it bisected for
  # the NARROWEST width that fitted. One long name anywhere on the chart then
  # split every other multi-word name, however much room it had.
  chars <- function(string) nchar(string)

  expect_identical(aev_wrap_to_width("male lives", 20, 2L, chars), "male lives")
  expect_identical(aev_wrap_to_width("male lives", 6, 2L, chars), "male\nlives")
})

test_that("the reserved count is still a ceiling", {

  # A name that will not fit in the reserved lines is broken as evenly as it
  # can be into them, rather than into as many as it would like. The band was
  # made deep enough for that many and no more.
  chars <- function(string) nchar(string)

  expect_identical(
    length(strsplit(aev_wrap_to_width("a b c d e f", 3, 2L, chars), "\n")[[1]]),
    2L
  )
})
