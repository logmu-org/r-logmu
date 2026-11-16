# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("category_domain() takes a size, labels or both", {
  by_size <- category_domain(3)
  expect_true(is_category_domain(by_size))
  expect_identical(by_size$size, 3L)
  expect_null(by_size$labels)
  expect_false(by_size$allows_undefined)
  expect_null(by_size$undefined_label)
  expect_null(by_size$sorted_codes)

  by_labels <- category_domain(labels = c("a", "b"))
  expect_identical(by_labels$size, 2L)

  expect_identical(category_domain(2, labels = c("a", "b"))$size, 2L)
  expect_error(category_domain(3, labels = c("a", "b")), "one label per code")
  expect_error(category_domain(), "needs a `size` or `labels`")
})

test_that("category_domain() checks its arguments", {
  expect_error(category_domain(0), "whole number from 1")
  expect_error(category_domain(2.5), "whole number from 1")
  expect_error(category_domain(NA_real_), "whole number from 1")
  expect_error(category_domain(.Machine$integer.max), "whole number from 1")
  expect_error(category_domain(labels = c("a", "a")), "must be unique; repeated: \"a\"")
  expect_error(category_domain(labels = c("a", " b")), "valid names")
  expect_error(category_domain(labels = character(0)), "at least one label")
  expect_error(category_domain(2, allows_undefined = NA), "TRUE or FALSE")
  expect_error(category_domain(2, name = ""), "valid name")
})

test_that("an undefined label allows undefined, and contradicting it is an error", {
  labelled <- category_domain(labels = c("low", "high"), undefined_label = "unknown")
  expect_true(labelled$allows_undefined)
  expect_identical(labelled$undefined_label, "unknown")
  expect_error(category_domain(labels = c("low", "high"), allows_undefined = FALSE,
                               undefined_label = "unknown"),
               "`allows_undefined` is FALSE")
})

test_that("undefined without a label is written by default, and only where allowed", {
  expect_identical(category_domain(labels = c("a", "b"), allows_undefined = TRUE)$undefined_label,
                   "(undefined)")
  expect_identical(category_domain(3, allows_undefined = TRUE)$undefined_label, "0")
  expect_null(category_domain(3)$undefined_label)
})

test_that("the undefined label may not be the text of a code", {
  expect_error(category_domain(labels = c("a", "b"), undefined_label = "b"), "already the text")
  expect_error(category_domain(5, undefined_label = "5"), "already the text")
  # "05" is not how code 5 is written, and 6 is past the end.
  expect_no_error(category_domain(5, undefined_label = "05"))
  expect_no_error(category_domain(5, undefined_label = "6"))
})

test_that("sorted_codes holds every code once, 0 included where allowed", {
  ordered <- category_domain(labels = c("low", "mid", "high"), undefined_label = "unknown",
                             sorted_codes = c(1, 2, 0, 3))
  expect_identical(ordered$sorted_codes, c(1L, 2L, 0L, 3L))
  expect_true(logmu:::is_ordered_domain(ordered))
  expect_identical(logmu:::category_domain_ranks(ordered), c(3L, 1L, 2L, 4L))

  expect_error(category_domain(3, sorted_codes = c(1, 2)), "every code of the domain exactly once")
  expect_error(category_domain(3, sorted_codes = c(1, 2, 2)), "exactly once")
  expect_error(category_domain(3, sorted_codes = c(0, 1, 2, 3)), "exactly once")
  expect_error(category_domain(3, allows_undefined = TRUE, sorted_codes = c(1, 2, 3)), "and 0")
  expect_error(category_domain(3, sorted_codes = c(1, 2.5, 3)), "exactly once")
  expect_error(category_domain(3, sorted_codes = c("1", "2", "3")), "exactly once")
})

test_that("unlabelled codes are written as plain digits, and only those are read back", {
  domain <- category_domain(100000, allows_undefined = TRUE)
  # Never "1e+05", whatever the scipen setting.
  old_options <- options(scipen = -9)
  on.exit(options(old_options), add = TRUE)
  expect_identical(logmu:::category_domain_text(domain, c(0L, 7L, 100000L)), c("0", "7", "100000"))
  expect_identical(logmu:::category_domain_codes(domain, c("0", "7", "100000")), c(0L, 7L, 100000L))
  expect_error(logmu:::category_domain_codes(domain, "07"), "Not a category")
  expect_error(logmu:::category_domain_codes(domain, "100001"), "Not a category")
  expect_error(logmu:::category_domain_codes(domain, "7.0"), "Not a category")
})

test_that("undefined text is read only where the domain allows it", {
  allows <- category_domain(labels = c("a", "b"), allows_undefined = TRUE)
  forbids <- category_domain(labels = c("a", "b"))
  expect_identical(logmu:::category_domain_codes(allows, c("(undefined)", "b")), c(0L, 2L))
  expect_error(logmu:::category_domain_codes(forbids, "(undefined)"), "Not a category")
})

test_that("domains are the same if their content is, whatever their names", {
  a <- category_domain(labels = c("x", "y"), name = "first")
  b <- category_domain(labels = c("x", "y"), name = "second")
  expect_true(logmu:::same_category_domain(a, b))
  expect_false(logmu:::same_category_domain(a, category_domain(labels = c("y", "x"))))
  expect_false(logmu:::same_category_domain(a, category_domain(2)))
  expect_false(logmu:::same_category_domain(a, category_domain(labels = c("x", "y"), sorted_codes = c(1, 2))))
  expect_false(logmu:::same_category_domain(a, category_domain(labels = c("x", "y"), allows_undefined = TRUE)))
  # The undefined label plays no part.
  expect_true(logmu:::same_category_domain(
    category_domain(labels = c("x", "y"), undefined_label = "none"),
    category_domain(labels = c("x", "y"), undefined_label = "unknown")
  ))
  # Nor does the order of the arguments that give the same order.
  expect_false(logmu:::same_category_domain(
    category_domain(labels = c("x", "y"), sorted_codes = c(1, 2)),
    category_domain(labels = c("x", "y"), sorted_codes = c(2, 1))
  ))
})

test_that("a naxp is a domain equal only to a naxp with the same pattern", {
  skip_if_not_installed("naxp")
  a <- category_domain(labels = c("x", "y"))
  expect_false(logmu:::same_category_domain(a, test_naxp("\\9")))
  expect_false(logmu:::same_category_domain(test_naxp("\\9"), a))
  expect_true(logmu:::same_category_domain(test_naxp("\\9"), test_naxp("\\9")))
  expect_false(logmu:::same_category_domain(test_naxp("\\9"), test_naxp("\\A")))
})

test_that("a domain is immutable and has a name", {
  domain <- category_domain(3, name = "types")
  expect_identical(name(domain), "types")
  name(domain) <- "renamed"
  expect_identical(name(domain), "renamed")
  expect_error(domain$size <- 4L, "immutable")
  expect_error(domain[["size"]] <- 4L, "immutable")
})

test_that("a domain prints its shape and order", {
  expect_output(print(category_domain(60, allows_undefined = TRUE, name = "types")),
                "<category_domain: types> 60 codes, unordered, undefined allowed as \"0\"")
  expect_output(
    print(category_domain(labels = c("low", "mid", "high"), undefined_label = "unknown",
                          sorted_codes = c(1, 2, 0, 3))),
    "3 codes, ordered, undefined allowed as \"unknown\"\nLabels: low, mid, high\nOrder: low < mid < unknown < high"
  )
})
