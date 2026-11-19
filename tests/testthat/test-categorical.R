# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("categorical() takes a size, labels or both", {
  by_size <- categorical(3)
  expect_true(is_categorical(by_size))
  expect_identical(by_size$size, 3L)
  expect_null(by_size$labels)
  expect_null(by_size$sorted_codes)
  expect_identical(names(unclass(by_size)), c("size", "labels", "sorted_codes"))

  by_labels <- categorical(labels = c("a", "b"))
  expect_identical(by_labels$size, 2L)

  expect_identical(categorical(2, labels = c("a", "b"))$size, 2L)
  expect_error(categorical(3, labels = c("a", "b")), "one label per code")
  expect_error(categorical(), "needs a `size` or `labels`")
})

test_that("categorical() checks its arguments", {
  expect_error(categorical(0), "whole number from 1")
  expect_error(categorical(2.5), "whole number from 1")
  expect_error(categorical(NA_real_), "whole number from 1")
  expect_error(categorical(.Machine$integer.max), "whole number from 1")
  expect_error(categorical(labels = c("a", "a")), "must be unique; repeated: \"a\"")
  expect_error(categorical(labels = c("a", " b")), "valid names")
  expect_error(categorical(labels = character(0)), "at least one label")
  expect_error(categorical(2, name = ""), "valid name")
  # The arguments of the first design are gone.
  expect_error(categorical(2, allows_undefined = TRUE), "unused argument")
  expect_error(categorical(2, undefined_label = "none"), "unused argument")
})

test_that("sorted_codes holds every category once, and code 0 has no place", {
  ordered <- categorical(labels = c("low", "unknown", "mid", "high"), sorted_codes = c(1, 3, 4, 2))
  expect_identical(ordered$sorted_codes, c(1L, 3L, 4L, 2L))
  expect_true(logmu:::is_ordered_domain(ordered))
  # Ranks by code + 1: code 0 ranks lowest, below every category.
  expect_identical(logmu:::category_ranks(ordered), c(0L, 1L, 4L, 2L, 3L))

  expect_error(categorical(3, sorted_codes = c(1, 2)), "every code from 1 to 3 exactly once")
  expect_error(categorical(3, sorted_codes = c(1, 2, 2)), "exactly once")
  expect_error(categorical(3, sorted_codes = c(0, 1, 2, 3)), "exactly once")
  expect_error(categorical(3, sorted_codes = c(0, 1, 2)), "exactly once")
  expect_error(categorical(3, sorted_codes = c(1, 2.5, 3)), "exactly once")
  expect_error(categorical(3, sorted_codes = c("1", "2", "3")), "exactly once")
})

test_that("unlabelled codes are written as plain digits, and only those are read back", {
  domain <- categorical(100000)
  # Never "1e+05", whatever the scipen setting.
  old_options <- options(scipen = -9)
  on.exit(options(old_options), add = TRUE)
  expect_identical(logmu:::category_text(domain, c(0L, 7L, 100000L)), c(NA, "7", "100000"))
  expect_identical(logmu:::category_codes(domain, c(NA, "7", "100000")), c(0L, 7L, 100000L))
  expect_error(logmu:::category_codes(domain, "0"), "Not a category")
  expect_error(logmu:::category_codes(domain, "07"), "Not a category")
  expect_error(logmu:::category_codes(domain, "100001"), "Not a category")
  expect_error(logmu:::category_codes(domain, "7.0"), "Not a category")
})

test_that("text is read as codes, with NA as code 0 and unknown text an error", {
  domain <- categorical(labels = c("a", "b"))
  expect_identical(logmu:::category_codes(domain, c("b", NA, "a")), c(2L, 0L, 1L))
  expect_error(logmu:::category_codes(domain, c("a", "c", "d")), "Not a category of this domain: \"c\", \"d\"")
  expect_identical(logmu:::category_text(domain, c(2L, 0L)), c("b", NA))
})

test_that("domains are the same if their content is, whatever their names", {
  a <- categorical(labels = c("x", "y"), name = "first")
  b <- categorical(labels = c("x", "y"), name = "second")
  expect_true(logmu:::same_domain(a, b))
  expect_false(logmu:::same_domain(a, categorical(labels = c("y", "x"))))
  expect_false(logmu:::same_domain(a, categorical(2)))
  expect_false(logmu:::same_domain(a, categorical(labels = c("x", "y"), sorted_codes = c(1, 2))))
  expect_false(logmu:::same_domain(
    categorical(labels = c("x", "y"), sorted_codes = c(1, 2)),
    categorical(labels = c("x", "y"), sorted_codes = c(2, 1))
  ))
})

test_that("a naxp is a domain equal only to a naxp with the same pattern", {
  skip_if_not_installed("naxp")
  a <- categorical(labels = c("x", "y"))
  expect_false(logmu:::same_domain(a, test_naxp("\\9")))
  expect_false(logmu:::same_domain(test_naxp("\\9"), a))
  expect_true(logmu:::same_domain(test_naxp("\\9"), test_naxp("\\9")))
  expect_false(logmu:::same_domain(test_naxp("\\9"), test_naxp("\\A")))
})

test_that("a domain is immutable and has a name", {
  domain <- categorical(3, name = "types")
  expect_identical(name(domain), "types")
  name(domain) <- "renamed"
  expect_identical(name(domain), "renamed")
  expect_error(domain$size <- 4L, "immutable")
  expect_error(domain[["size"]] <- 4L, "immutable")
  expect_error(domain["size"] <- 4L, "immutable")
})

test_that("a domain prints its shape and order", {
  expect_output(print(categorical(60, name = "types")), "<categorical: types> 60 categories, unordered")
  expect_output(print(categorical(1)), "<categorical> 1 category, unordered")
  expect_output(
    print(categorical(labels = c("low", "unknown", "mid", "high"), sorted_codes = c(1, 3, 4, 2))),
    "4 categories, ordered\nLabels: low, unknown, mid, high\nOrder: low < mid < high < unknown"
  )
})
