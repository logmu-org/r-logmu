# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# NOTES ARE LINES OF TEXT AND NOTHING MORE. Tim's rule of 2026-08-28: the chart
# is not coupled to `aev()`, so its notes must not be either. Nothing in the
# chart knows what an overdispersion is, or a Z, or a test mortality -- it takes
# strings and prints them. That is the shape that ports to Python and to C#.

noted_aev <- function() {
  aev <- create_aev(A = c(1100, 970), E = c(1000, 1000), V = c(2500, 1600))
  names(aev) <- c("a", "b")
  aev
}

caption_of <- function(plot) plot$labels$caption

drawn_caption <- function(plot) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  table <- ggplot2::ggplotGrob(plot)
  cell <- which(table$layout$name == "caption")
  if (length(cell) == 0L) return(character(0))
  labels_in(table$grobs[[cell]])
}

#### What a note is ####

test_that("no notes means no caption", {

  expect_null(aev_notes_text(NULL))
  expect_null(aev_notes_text(character(0)))
  expect_null(caption_of(autoplot(noted_aev())))
})

test_that("one note per line, in the order given", {

  expect_identical(aev_notes_text("only one"), "only one")
  expect_identical(aev_notes_text(c("first", "second")), "first\nsecond")

  # The order is the caller's, not sorted and not reversed.
  expect_identical(aev_notes_text(c("z", "a", "m")), "z\na\nm")
})

test_that("an empty line is kept, so a caller can space a block", {

  expect_identical(aev_notes_text(c("a", "", "b")), "a\n\nb")
})

test_that("notes must be text, and must not be missing", {

  expect_error(aev_notes_text(2.0), "must be a character vector")
  expect_error(aev_notes_text(list("a")), "must be a character vector")
  expect_error(aev_notes_text(c("a", NA)), "cannot contain missing values")
})

#### On the chart ####

test_that("the notes become the chart's caption", {

  plot <- autoplot(noted_aev(), notes = c("Overdispersion 2.0", "Z = 1.31"))

  expect_identical(caption_of(plot), "Overdispersion 2.0\nZ = 1.31")
})

test_that("the notes are actually drawn", {

  # Not just set on the object: present in the table ggplot2 builds, which is
  # what the reader sees.
  drawn <- drawn_caption(autoplot(noted_aev(), notes = c("first note", "second note")))

  expect_true(any(grepl("first note", drawn, fixed = TRUE)))
  expect_true(any(grepl("second note", drawn, fixed = TRUE)))
})

test_that("a caption added afterwards wins, as any ggplot2 user would expect", {

  plot <- autoplot(noted_aev(), notes = "from notes") +
    ggplot2::labs(caption = "from labs")

  expect_identical(caption_of(plot), "from labs")
})

test_that("the caption is left-aligned against the plot, not ggplot2's corner", {

  # Slide 2 left-aligns everything of this kind, and the theme has said so since
  # the chart was built. This pins it now that something actually uses it.
  theme <- theme_aev()

  expect_identical(ggplot2::calc_element("plot.caption.position", theme), "plot")
  expect_equal(ggplot2::calc_element("plot.caption", theme)$hjust, 0)
})

test_that("notes know nothing about what they say", {

  # THE POINT OF THE DESIGN. A note is a string; the chart never parses it, so
  # anything a caller wants recorded goes through unchanged -- including text
  # that has nothing to do with mortality at all.
  # The pound sign is written as an escape so this file stays ASCII, which the
  # package requires. What is tested is that it survives the chart unchanged.
  odd <- c("Prepared by AN Other", "1 + 1 = 2", "\u00A3 amounts", "")

  plot <- autoplot(noted_aev(), notes = odd)

  expect_identical(caption_of(plot), paste(odd, collapse = "\n"))
})

test_that("plot() passes notes through as well", {

  # `plot()` forwards `...` to `autoplot()`, so anything that works on one
  # works on the other.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(plot(noted_aev(), notes = "through plot()"))
})

#### The three blocks of text, and how big they are ####

test_that("the title is twice the base size, as the prototype has it", {

  # `TitleFontSize` is 24 CSS px against the design's 12 px base text. ggplot2's
  # own default is rel(1.2), which is what this chart used until 2026-08-28 and
  # which made the title read as timid: the prototype's title is 2.7 times its
  # axis text, and rel(1.2) made it 1.5.
  title <- ggplot2::calc_element("plot.title", theme_aev())

  expect_equal(title$size, aev_reference_text_size * 2)
  expect_equal(title$hjust, 0)
})

test_that("the three blocks step down in size and in weight of colour", {

  theme <- theme_aev()
  size <- function(name) ggplot2::calc_element(name, theme)$size

  expect_gt(size("plot.title"), size("plot.subtitle"))
  expect_gt(size("plot.subtitle"), size("plot.caption"))

  # And each step down is a step lighter, so the order reads without rules or
  # indentation between the blocks.
  expect_identical(ggplot2::calc_element("plot.title", theme)$colour, aev_palette$title_text)
  expect_identical(ggplot2::calc_element("plot.subtitle", theme)$colour, aev_palette$axis_title_text)
  expect_identical(ggplot2::calc_element("plot.caption", theme)$colour, aev_palette$axis_scale_text)
})

test_that("the title and subtitle are set off from what follows them", {

  # `DefaultMargin` is 10 px against a 12 px base, so five sixths of the base
  # size. Without it the title sat straight on top of the legend.
  theme <- theme_aev()
  expected <- aev_reference_text_size * 10 / 12

  for (name in c("plot.title", "plot.subtitle")) {
    # Read in the units it was made in rather than converting: `convertHeight`
    # needs a device, and this suite refuses to open one implicitly.
    margin <- ggplot2::calc_element(name, theme)$margin
    expect_identical(as.character(grid::unitType(margin[3])), "points")
    expect_equal(as.numeric(margin[3]), expected)
  }
})

test_that("the headings follow the base size like everything else", {

  larger <- theme_aev(base_size = aev_reference_text_size * 2)

  expect_equal(
    ggplot2::calc_element("plot.title", larger)$size,
    aev_reference_text_size * 4
  )
  expect_equal(
    as.numeric(ggplot2::calc_element("plot.title", larger)$margin[3]),
    aev_reference_text_size * 2 * 10 / 12
  )
})

#### Arguments that are not arguments ####

test_that("an unrecognised argument is refused, not dropped", {

  # `...` exists so `plot()` can forward to `autoplot()`. Anything else landing
  # there used to be discarded in silence, so `plot(aev, main = "...")` drew a
  # chart with no title and said nothing about why.
  aev <- create_aev(A = 1100, E = 1000, V = 2500)

  expect_error(autoplot(aev, zzz = 1), "Unused argument")
  expect_error(autoplot(aev, zzz = 1), "It takes")
})

test_that("the base R names are translated by hand", {

  # The mistake most likely to be made, and the least likely to explain itself.
  aev <- create_aev(A = 1100, E = 1000, V = 2500)

  expect_error(autoplot(aev, main = "t"), "use title", fixed = TRUE)
  expect_error(autoplot(aev, cex = 2), "theme_aev(base_size = )", fixed = TRUE)
  expect_error(autoplot(aev, xlab = "x"), "names(aev)", fixed = TRUE)

  # `sub` needs no advice: it is a prefix of `subtitle`, so R matches it there
  # before the check runs. Base graphics habits get the right answer by luck.
  expect_identical(autoplot(aev, sub = "s")$labels$subtitle, "s")
})

test_that("plot() refuses them too, since it forwards", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(plot(create_aev(A = 1100, E = 1000, V = 2500), main = "t"), "use title")
})

test_that("everything the chart does take is still accepted", {

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  aev <- create_aev(A = 1100, E = 1000, V = 2500)

  expect_no_error(autoplot(
    aev, title = "t", subtitle = "s", notes = "n",
    residuals = FALSE, log_range = c(-0.6, 0.6), log_step = 0.2
  ))
})

test_that("R's own partial matching still reaches the real arguments", {

  # NOT a hole in the check. `note` and `residual` are prefixes of real
  # arguments, so R matches them before `...` ever sees them -- which is
  # ordinary R behaviour and worth knowing rather than defeating. Only names
  # that match nothing reach the check.
  aev <- create_aev(A = 1100, E = 1000, V = 2500)

  expect_identical(autoplot(aev, note = "n")$labels$caption, "n")
  expect_error(autoplot(aev, log = c(-1, 1)), "matches multiple formal arguments")
})
