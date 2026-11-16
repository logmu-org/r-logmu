# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

codes_of <- function(x) as.vector(unclass(x))

ab <- category_domain(labels = c("a", "b"))

test_that("as_categorical() turns a factor's levels into a domain", {
  x <- as_categorical(factor(c("def", "act", "def"), levels = c("def", "act", "pen")))
  expect_true(is_categorical(x))
  expect_identical(codes_of(x), c(1L, 2L, 1L))
  expect_identical(domain(x)$labels, c("def", "act", "pen"))
  expect_false(domain(x)$allows_undefined)
  expect_null(domain(x)$sorted_codes)
})

test_that("a factor names the level that stands for undefined", {
  x <- as_categorical(factor(c("a", "none", "b")), undefined_level = "none")
  expect_identical(codes_of(x), c(1L, 0L, 2L))
  expect_identical(domain(x)$labels, c("a", "b"))
  expect_true(domain(x)$allows_undefined)
  expect_identical(format(x), c("a", "none", "b"))
})

test_that("an ordered factor gives its level order, with undefined in its place", {
  f <- factor(c("high", "unknown", "low"), levels = c("low", "unknown", "high"), ordered = TRUE)
  x <- as_categorical(f, undefined_level = "unknown")
  expect_identical(domain(x)$labels, c("low", "high"))
  expect_identical(domain(x)$sorted_codes, c(1L, 0L, 2L))
  expect_identical(codes_of(x), c(2L, 0L, 1L))
  expect_identical(domain(as_categorical(factor("a", ordered = TRUE)))$sorted_codes, 1L)
  expect_error(as_categorical(factor("a", ordered = TRUE), undefined_level = "none"),
               "must be a level of an ordered factor")
})

test_that("as_categorical() refuses NA, so missing values need a level", {
  expect_error(as_categorical(factor(c("a", NA))), "a level of their own")
  expect_error(as_categorical(c("a", NA), ab), "cannot hold NA")
  expect_error(as_categorical(c(1, NA), category_domain(2)), "no NA")
  missing_as_level <- factor(c("a", "unknown"))
  x <- as_categorical(missing_as_level, undefined_level = "unknown")
  expect_identical(codes_of(x), c(1L, 0L))
})

test_that("text is matched to a domain, and unknown text is an error", {
  domain <- category_domain(labels = c("low", "high"), allows_undefined = TRUE)
  x <- as_categorical(c("high", "(undefined)", "low"), domain)
  expect_identical(codes_of(x), c(2L, 0L, 1L))
  expect_error(as_categorical("hgih", domain), "Not a category of this domain: \"hgih\"")
  expect_error(as_categorical(factor("hgih"), domain), "Not a category")
  expect_error(as_categorical(factor("low"), domain, undefined_level = "x"), "not both")
})

test_that("text needs a domain", {
  expect_error(as_categorical(c("b", "a")), "Text needs a `domain`")
})

test_that("numbers are codes in a domain", {
  domain <- category_domain(60, allows_undefined = TRUE)
  expect_identical(codes_of(as_categorical(c(12, 0, 60), domain)), c(12L, 0L, 60L))
  expect_error(as_categorical(61, domain), "from 0 to 60")
  expect_error(as_categorical(1.5, domain), "whole numbers")
  expect_error(as_categorical(0, category_domain(3)), "does not allow code 0")
  expect_error(as_categorical(1), "need a `domain`")
})

test_that("as_categorical() of a categorical checks its domain", {
  x <- as_categorical(c("a", "b"), ab)
  expect_identical(as_categorical(x), x)
  expect_identical(as_categorical(x, category_domain(labels = c("a", "b"))), x)
  expect_error(as_categorical(x, category_domain(2)), "different domain")
  expect_error(as_categorical(list(1)), "class `list`")
})

test_that("subsetting keeps the domain and never leaves NA", {
  allows <- as_categorical(c(1, 2, 3), category_domain(3, allows_undefined = TRUE))
  expect_identical(codes_of(allows[2:3]), c(2L, 3L))
  expect_identical(domain(allows[2]), domain(allows))
  expect_identical(codes_of(allows[c(1, 5)]), c(1L, 0L))
  expect_identical(codes_of(allows[[3]]), 3L)

  forbids <- as_categorical(c(1, 2), category_domain(2))
  expect_error(forbids[3], "does not allow code 0")
  expect_error(length(forbids) <- 3, "does not allow code 0")
  lengthened <- allows
  length(lengthened) <- 4
  expect_identical(codes_of(lengthened), c(1L, 2L, 3L, 0L))
})

test_that("assignment takes a categorical in the same domain or text", {
  domain <- category_domain(labels = c("a", "b", "c"))
  x <- as_categorical(c("a", "b"), domain)
  x[2] <- "c"
  expect_identical(format(x), c("a", "c"))
  x[[1]] <- as_categorical("b", domain)
  expect_identical(format(x), c("b", "c"))
  expect_error(x[1] <- "z", "Not a category")
  expect_error(x[1] <- 2L, "Assign a categorical")
  expect_error(x[1] <- as_categorical("a", ab), "different domains")
})

test_that("c() combines categorical vectors in the same domain only", {
  combined <- c(as_categorical("a", ab), as_categorical("b", ab))
  expect_identical(format(combined), c("a", "b"))
  expect_error(c(as_categorical("a", ab), as_categorical(1, category_domain(1))), "same domain")
})

test_that("comparison is by equality, and by order only in an ordered domain", {
  unordered <- as_categorical(c("a", "b", "a"), ab)
  expect_identical(unordered == "a", c(TRUE, FALSE, TRUE))
  expect_identical(unordered != unordered[1], c(FALSE, TRUE, FALSE))
  expect_identical(unordered == 2L, c(FALSE, TRUE, FALSE))
  expect_error(unordered < "b", "ordered domain")
  expect_error(unordered + 1, "not defined for categorical vectors")
  expect_error(max(unordered), "not defined for categorical vectors")
})

test_that("an ordered domain compares and sorts by its order, code 0 included", {
  # Codes 1, 2, 3 are low, mid, high; the order puts unknown between mid and high.
  domain <- category_domain(labels = c("low", "mid", "high"), undefined_label = "unknown",
                            sorted_codes = c(1, 2, 0, 3))
  x <- as_categorical(c("high", "unknown", "low", "mid"), domain)
  expect_identical(x < "unknown", c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(x >= "unknown", c(TRUE, TRUE, FALSE, FALSE))
  expect_identical(x > 2L, c(TRUE, TRUE, FALSE, FALSE))
  expect_error(x > 9L, "from 0 to 3")
  expect_identical(format(sort(x)), c("low", "mid", "unknown", "high"))
  expect_identical(levels(factor(x)), c("low", "mid", "unknown", "high"))
})

test_that("%in% and match() compare labels", {
  x <- as_categorical(c("a", "b", "a"), ab)
  expect_identical(x %in% "b", c(FALSE, TRUE, FALSE))
  expect_identical(match(x, x[2]), c(NA, 1L, NA))
  # Documented: numbers on the right are compared with the labels.
  expect_identical(x %in% 2L, c(FALSE, FALSE, FALSE))
})

test_that("base R tools see a categorical through its codes and labels", {
  x <- as_categorical(c("b", "a", "b"), ab)
  expect_identical(as.character(x), c("b", "a", "b"))
  expect_identical(as.integer(x), c(2L, 1L, 2L))
  expect_identical(format(sort(x)), c("a", "b", "b"))
  expect_identical(format(unique(x)), c("b", "a"))
  expect_identical(format(rep(x[2], 2)), c("a", "a"))
  expect_identical(as.vector(table(x)), c(1L, 2L))
  expect_identical(levels(factor(x)), c("a", "b"))
  frame <- data.frame(status = x)
  expect_true(is_categorical(frame$status))
  expect_identical(format(frame$status[2:3]), c("a", "b"))
})

test_that("forcats::as_factor() gives every category of the domain, in order", {
  skip_if_not_installed("forcats")
  unordered <- category_domain(labels = c("low", "high"), undefined_label = "none")
  f <- forcats::as_factor(as_categorical(c("high", "none"), unordered))
  expect_identical(levels(f), c("none", "low", "high"))
  expect_false(is.ordered(f))

  ordered <- category_domain(labels = c("low", "high"), undefined_label = "none",
                             sorted_codes = c(1, 0, 2))
  x <- as_categorical(c("high", "none"), ordered)
  g <- forcats::as_factor(x)
  expect_identical(levels(g), c("low", "none", "high"))
  expect_true(is.ordered(g))
  expect_identical(as.character(g), c("high", "none"))
  # And back again, exactly, with or without the domain.
  expect_identical(as_categorical(g, ordered), x)
  round_trip <- as_categorical(g, undefined_level = "none")
  expect_identical(codes_of(round_trip), codes_of(x))
  expect_true(logmu:::same_category_domain(domain(round_trip), ordered))
})

test_that("a categorical prints its labels", {
  x <- as_categorical(c(1, 0), category_domain(2, allows_undefined = TRUE, name = "types"))
  expect_output(print(x), "<categorical\\[2\\]> types\n\\[1\\] 1 0")
})

test_that("a naxp_encoded vector has its naxp as its domain", {
  skip_if_not_installed("naxp")
  postcode <- test_naxp("\\9\\9\\9\\9 \\A\\A")
  x <- test_naxp_encoded(c(5L, 0L), postcode)
  expect_identical(domain(x), postcode)
  expect_error(domain(1), "has a domain")
})

# `R/naxp_interop.R` reads naxp's objects directly instead of importing the
# package. If naxp ever changes how it stores them, this is the test that says so.
test_that("logmu reads naxp's own objects as naxp makes them", {
  skip_if_not_installed("naxp")
  postcode <- naxp::naxp_parse("\\9\\9\\9\\9 \\A\\A")
  encoded <- naxp::naxp_encoded(c(30L, 0L, 6760000L), postcode)
  expect_true(logmu:::is_naxp(postcode))
  expect_true(logmu:::is_naxp_encoded(encoded))
  expect_identical(logmu:::naxp_pattern(postcode), naxp::naxp_pattern(postcode))
  expect_identical(logmu:::naxp_of(encoded), postcode)
  expect_identical(logmu:::naxp_codes(encoded), c(30L, 0L, 6760000L))
  expect_identical(logmu:::new_naxp_encoded_codes(c(30L, 0L, 6760000L), postcode), encoded)
})
