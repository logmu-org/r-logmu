# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Shared by the chart tests, which all have to look inside `grid` objects.

# NO TEST MAY OPEN A DEVICE BY ACCIDENT.
#
# Anything that measures a grob needs a device, and with none open R starts the
# default one and writes an `Rplots.pdf` into `tests/testthat`. It is easy to
# miss -- it happened twice while the chart was being written, and once it
# reached a commit and would have shipped in the tarball.
#
# So the default device is made an error. Every test here opens its own with
# `grDevices::pdf(NULL)`, which is unaffected; a test that forgets fails by
# name instead of leaving a file behind.
options(device = function(...) {
  stop("a device was opened implicitly -- open one with grDevices::pdf(NULL)")
})

# CONVERTING A UNIT NEEDS A DEVICE EVEN WHEN THE UNIT IS ABSOLUTE. With none
# open, `grid` opens the default and leaves an `Rplots.pdf` in `tests/testthat`
# for `git status` to find. That one escaped a first sweep of the test helpers
# and only showed up when a whole file ran.
pixels <- function(length) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::convertWidth(length, "inches", valueOnly = TRUE) * 96
}

# Every grob anywhere in a tree, itself included. `labels_in` walks the same
# way; this hands back the grobs so a test can look at what they are drawn in.
#
# BOTH KINDS OF CHILD. A `gTree` keeps them in `children` and a `gtable` in
# `grobs`, and a plot is a nest of the two.
grobs_in <- function(grob) {
  children <- c(grob$children, grob$grobs)
  c(
    list(grob),
    unlist(lapply(children, grobs_in), recursive = FALSE, use.names = FALSE)
  )
}

# Every piece of text anywhere in a grob tree. An axis is an `absoluteGrob`
# wrapping a `gtable` wrapping a `titleGrob` wrapping the text, and the depth
# is not something a test should have to know.
labels_in <- function(grob) {
  if (!is.null(grob$label)) return(as.character(grob$label))
  children <- c(grob$children, grob$grobs)
  if (length(children) == 0L) return(character(0))
  unlist(lapply(children, labels_in), use.names = FALSE)
}
