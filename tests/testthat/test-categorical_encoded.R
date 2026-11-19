# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

codes_of <- function(x) as.vector(unclass(x))

ab <- categorical(labels = c("a", "b"))
sex <- categorical(labels = c("female", "male"), name = "sex")

test_that("encode() matches text to a domain, with NA as code 0", {
  x <- encode(c("male", NA, "female"), sex)
  expect_true(is_categorical_encoded(x))
  expect_false(is_categorical(x))
  expect_identical(codes_of(x), c(2L, 0L, 1L))
  expect_identical(domain(x), sex)
  expect_identical(codes_of(encode(NA, sex)), 0L)
  expect_identical(codes_of(encode(character(0), sex)), integer(0))
})

test_that("unknown text is an error, not code 0", {
  expect_error(encode(c("male", "Male"), sex), "Not a category of this domain: \"Male\"")
  expect_error(encode(factor("mlae"), sex), "Not a category")
})

test_that("a factor is read through its labels, and never makes its own domain", {
  status <- factor(c("def", NA, "act"), levels = c("act", "def", "pen"))
  x <- encode(status, categorical(labels = levels(status)))
  expect_identical(codes_of(x), c(2L, 0L, 1L))
  expect_error(encode(status), "needs a `domain`")
})

test_that("encode() refuses numbers and anything else that is not text", {
  types <- categorical(60)
  expect_error(encode(c(12, 60), types), "For codes, use categorical_encoded\\(\\)")
  expect_error(encode(12L, types), "For codes, use categorical_encoded\\(\\)")
  expect_error(encode(list("a"), ab), "not an object of class `list`")
  expect_error(encode("a", "not a domain"), "must be a `categorical` or a naxp")
})

test_that("encode() of an encoded vector checks its domain", {
  x <- encode(c("a", "b"), ab)
  expect_identical(encode(x, categorical(labels = c("a", "b"))), x)
  expect_error(encode(x, categorical(2)), "different domain")
})

test_that("categorical_encoded() takes codes, with 0 for invalid and no NA", {
  types <- categorical(60, name = "types")
  x <- categorical_encoded(c(12, 0, 60), types)
  expect_identical(codes_of(x), c(12L, 0L, 60L))
  expect_identical(domain(x), types)
  expect_error(categorical_encoded(61, types), "from 0 to 60")
  expect_error(categorical_encoded(1.5, types), "whole numbers")
  expect_error(categorical_encoded(-1, types), "whole numbers")
  expect_error(categorical_encoded(c(1, NA), types), "no NA")
  expect_error(categorical_encoded("12", types), "whole numbers")
  expect_error(categorical_encoded(1, ab$labels), "must be a `categorical`")
})

test_that("decode() and as.character() give the text, with NA for code 0", {
  x <- encode(c("male", NA, "female"), sex)
  expect_identical(decode(x), c("male", NA, "female"))
  expect_identical(as.character(x), c("male", NA, "female"))
  expect_identical(decode(encode(decode(x), sex)), decode(x))
  expect_error(decode("male"), "can be decoded")
})

test_that("is.na() and anyNA() are TRUE for code 0", {
  x <- encode(c("male", NA, "female"), sex)
  expect_identical(is.na(x), c(FALSE, TRUE, FALSE))
  expect_true(anyNA(x))
  expect_false(anyNA(x[c(1, 3)]))
})

test_that("subsetting keeps the domain and never leaves NA", {
  x <- categorical_encoded(c(1, 2, 3), categorical(3))
  expect_identical(codes_of(x[2:3]), c(2L, 3L))
  expect_identical(domain(x[2]), domain(x))
  expect_identical(codes_of(x[c(1, 5)]), c(1L, 0L))
  expect_identical(codes_of(x[c(1, NA)]), c(1L, 0L))
  expect_identical(codes_of(x[[3]]), 3L)
  lengthened <- x
  length(lengthened) <- 4
  expect_identical(codes_of(lengthened), c(1L, 2L, 3L, 0L))
})

test_that("assignment takes a vector in the same domain, text or NA", {
  domain <- categorical(labels = c("a", "b", "c"))
  x <- encode(c("a", "b"), domain)
  x[2] <- "c"
  expect_identical(decode(x), c("a", "c"))
  x[[1]] <- encode("b", domain)
  expect_identical(decode(x), c("b", "c"))
  x[1] <- NA
  expect_identical(codes_of(x), c(0L, 3L))
  x[3] <- "a"
  expect_identical(codes_of(x), c(0L, 3L, 1L))
  x[5] <- "a"
  expect_identical(codes_of(x), c(0L, 3L, 1L, 0L, 1L))
  expect_error(x[1] <- "z", "Not a category")
  expect_error(x[1] <- 2L, "Assign a categorical_encoded")
  expect_error(x[1] <- encode("a", ab), "different domains")
})

test_that("c() combines vectors in the same domain only", {
  combined <- c(encode("a", ab), encode(NA, ab), encode("b", ab))
  expect_identical(codes_of(combined), c(1L, 0L, 2L))
  expect_error(c(encode("a", ab), categorical_encoded(1, categorical(1))), "same domain")
})

test_that("comparison gives TRUE or FALSE only, with code 0 an ordinary code", {
  x <- encode(c("a", NA, "b"), ab)
  expect_identical(x == "a", c(TRUE, FALSE, FALSE))
  expect_identical(x != "a", c(FALSE, TRUE, TRUE))
  expect_identical(x == x, c(TRUE, TRUE, TRUE))
  expect_identical(x == x[2], c(FALSE, TRUE, FALSE))
  expect_identical("b" == x, c(FALSE, FALSE, TRUE))
})

test_that("comparing with NA, numbers or another domain is an error", {
  x <- encode(c("a", NA), ab)
  expect_error(x == NA, "`is.na\\(\\)` finds it")
  expect_error(x == NA_character_, "`is.na\\(\\)` finds it")
  expect_error(x == c("a", NA), "`is.na\\(\\)` finds it")
  expect_error(x == 1L, "only with one in the same domain, or with text")
  expect_error(x == encode("a", categorical(labels = c("a", "c"))), "different domains")
  expect_error(x == "c", "Not a category")
})

test_that("only comparison is defined, and ordering only in an ordered domain", {
  x <- encode(c("a", "b"), ab)
  expect_error(x < "b", "ordered domain")
  expect_error(x + 1, "not defined for categorical_encoded vectors")
  expect_error(!x, "not defined for categorical_encoded vectors")
  expect_error(max(x), "not defined for categorical_encoded vectors")
  expect_error(sqrt(x), "not defined for categorical_encoded vectors")
})

test_that("an ordered domain compares and sorts by its order, with code 0 lowest", {
  # Codes 1..4 are low, unknown, mid, high; unknown is ranked above high.
  domain <- categorical(labels = c("low", "unknown", "mid", "high"), sorted_codes = c(1, 3, 4, 2))
  x <- encode(c("high", "unknown", NA, "low", "mid"), domain)
  expect_identical(x < "mid", c(FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_identical(x >= "high", c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(x > x[3], c(TRUE, TRUE, FALSE, TRUE, TRUE))
  # Code 0 is not NA to sort(), so it is kept, and sorts first.
  expect_identical(decode(sort(x)), c(NA, "low", "mid", "high", "unknown"))
  expect_identical(decode(sort(x, decreasing = TRUE)), c("unknown", "high", "mid", "low", NA))
  expect_identical(levels(factor(x)), c("low", "mid", "high", "unknown"))
})

test_that("an unordered domain sorts by code, with code 0 first", {
  x <- encode(c("b", NA, "a"), ab)
  expect_identical(decode(x[order(x)]), c(NA, "a", "b"))
})

test_that("%in% and match() compare text", {
  x <- encode(c("a", "b", NA), ab)
  expect_identical(x %in% "b", c(FALSE, TRUE, FALSE))
  expect_identical(match(x, x[2]), c(NA, 1L, NA))
})

test_that("base R tools see a vector through its codes and labels", {
  x <- encode(c("b", "a", "b", NA), ab)
  expect_identical(as.integer(x), c(2L, 1L, 2L, 0L))
  expect_identical(decode(unique(x)), c("b", "a", NA))
  expect_identical(decode(rep(x[2], 2)), c("a", "a"))
  expect_identical(as.vector(table(x)), c(1L, 2L))
  expect_identical(levels(factor(x)), c("a", "b"))
  frame <- data.frame(status = x)
  expect_true(is_categorical_encoded(frame$status))
  expect_identical(decode(frame$status[2:3]), c("a", "b"))
})

test_that("forcats::as_factor() gives every category of the domain, in order, and NA for code 0", {
  skip_if_not_installed("forcats")
  unordered <- categorical(labels = c("low", "high"))
  f <- forcats::as_factor(encode(c("high", NA), unordered))
  expect_identical(levels(f), c("low", "high"))
  expect_identical(as.character(f), c("high", NA))
  expect_false(is.ordered(f))

  ordered <- categorical(labels = c("low", "high", "unknown"), sorted_codes = c(1, 3, 2))
  x <- encode(c("high", NA, "unknown"), ordered)
  g <- forcats::as_factor(x)
  expect_identical(levels(g), c("low", "unknown", "high"))
  expect_true(is.ordered(g))
  # And back again, exactly.
  expect_identical(encode(g, ordered), x)
})

test_that("a categorical_encoded vector prints its labels, with code 0 marked invalid", {
  x <- categorical_encoded(c(1, 0), categorical(2, name = "types"))
  expect_output(print(x), "<categorical_encoded\\[2\\]> types\n\\[1\\] 1 +<invalid \\(0\\)>")
  expect_identical(format(x), c("1", "<invalid (0)>"))
})

test_that("domain() is defined for encoded vectors only", {
  expect_error(domain(1), "has a domain")
  expect_error(domain(sex), "has a domain")
})

test_that("encode() and decode() hand a naxp to the naxp package", {
  skip_if_not_installed("naxp")
  postcode <- test_naxp("\\9\\9\\9\\9 \\A\\A")
  x <- encode(c("1234AB", NA, "not a postcode"), postcode)
  expect_true(logmu:::is_naxp_encoded(x))
  expect_identical(domain(x), postcode)
  expect_identical(codes_of(x)[2:3], c(0L, 0L))
  expect_true(codes_of(x)[1] > 0L)
  expect_identical(decode(x), c("1234AB", NA, NA))
  expect_identical(codes_of(encode(factor("1234AB"), postcode)), codes_of(x)[1])
  expect_identical(encode(x, postcode), x)
  expect_error(encode(x, test_naxp("\\A")), "different naxp")
  expect_error(encode(1234, postcode), "For codes, use naxp::naxp_encoded\\(\\)")
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
