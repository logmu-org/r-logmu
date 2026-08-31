# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests for the band constructor family and the `includes` collection.
#
# THE ONE EDGE RULE IS WHAT MOST OF THIS PINS. Band edges are `c(from, thresholds, to)` with a NULL
# bound meaning unbounded, and `by` generating the interior thresholds. Every difference between
# `ages(65, 95, by = 5)` excluding the outside and `bands(.i$pension, thresholds = ...)` being open
# at both ends falls out of which bounds the user supplied -- not out of which constructor was
# called. The tests are written to fail if that ever becomes two rules again.

test_that("singular constructors return an include, plural an includes", {
  expect_true(is_include(band(.i$pension, 0, 100)))
  expect_true(is_include(age(65, 95)))
  expect_true(is_include(period(2010, 2020)))

  expect_true(is_includes(bands(.i$pension, 0, 100)))
  expect_true(is_includes(ages(65, 95)))
  expect_true(is_includes(periods(2010, 2020)))
})

test_that("a plural constructor with no `by` still returns an includes of one", {
  # The return type does not depend on which arguments were supplied, which is the whole reason the
  # singular constructors were kept rather than letting `by = NULL` collapse to a scalar.
  one <- ages(65, 95)
  expect_true(is_includes(one))
  expect_equal(length(one), 1L)
})

test_that("`by` divides the range into bands, outside excluded", {
  a <- ages(65, 95, by = 5)
  expect_equal(length(a), 6L)
  expect_equal(names(a), c("65-70", "70-75", "75-80", "80-85", "85-90", "90-95"))
  expect_equal(group_names(a), rep("age", 6L))
})

test_that("thresholds with no bounds give bands open at both ends", {
  b <- bands(.i$pension, thresholds = c(5000, 10000, 20000))
  expect_equal(length(b), 4L)
  expect_equal(names(b), c("< 5000", "5000-10000", "10000-20000", ">= 20000"))
  expect_equal(group_names(b), rep("pension", 4L))
})

test_that("a bound closes the end it is given for", {
  b <- bands(.i$pension, from = 0, thresholds = c(5000, 10000))
  expect_equal(length(b), 3L)
  expect_equal(names(b), c("0-5000", "5000-10000", ">= 10000"))
})

test_that("`by` must divide `to - from` exactly", {
  # Refused rather than yielding a silent short final band, which would be a reporting trap once it
  # reached a chart axis.
  expect_error(ages(65, 96, by = 5), "divide")
  expect_silent(ages(65, 95, by = 5))
})

test_that("`by` and `thresholds` cannot both be given", {
  expect_error(ages(65, 95, by = 5, thresholds = 80), "not both")
})

test_that("`by` needs both bounds", {
  expect_error(ages(65, by = 5), "needs both")
  expect_error(ages(to = 95, by = 5), "needs both")
})

test_that("thresholds must increase", {
  expect_error(bands(.i$pension, thresholds = c(10000, 5000)), "strictly increasing")
  expect_error(bands(.i$pension, from = 6000, thresholds = 5000), "strictly increasing")
})

test_that("a band needs at least one edge", {
  expect_error(ages(), "at least one edge")
})

# ---- what the bands actually resolve to ------------------------------------

test_that("banding a time-invariant variable gates rather than clips", {
  b <- bands(.i$pension, thresholds = 5000)

  # Below the threshold the first band admits the whole of time; above it, none of it.
  expect_equal(period_included(b[[1]], .i = list(pension = 3000)), include_all_interval())
  expect_equal(period_included(b[[1]], .i = list(pension = 7000)), include_none_interval())
  expect_equal(period_included(b[[2]], .i = list(pension = 7000)), include_all_interval())
})

test_that("banding is clopen at every internal edge", {
  # An individual exactly on a threshold belongs to the UPPER band and not the lower one. This is
  # what stops adjacent bands double-counting, and it is why the gate is `>= lo & < hi`.
  b <- bands(.i$pension, thresholds = 5000)
  expect_equal(period_included(b[[1]], .i = list(pension = 5000)), include_none_interval())
  expect_equal(period_included(b[[2]], .i = list(pension = 5000)), include_all_interval())
})

test_that("an unbounded side does not narrow the interval", {
  after <- band(.t, 2010, NULL)
  expect_equal(
    period_included(after, .i = list()),
    datey::datey_interval(datey::datey(2010), datey::datey(datey::valid_years_end))
  )

  before <- band(.t, NULL, 2020)
  expect_equal(
    period_included(before, .i = list()),
    datey::datey_interval(datey::datey(datey::valid_years_start), datey::datey(2020))
  )
})

test_that("age bands resolve per individual", {
  a <- ages(65, 95, by = 15)
  .i <- list(birth = datey::datey(1950))
  expect_equal(period_included(a[[1]], .i = .i),
               datey::datey_interval(datey::datey(2015), datey::datey(2030)))
  expect_equal(period_included(a[[2]], .i = .i),
               datey::datey_interval(datey::datey(2030), datey::datey(2045)))
})

test_that("a NaN banded variable falls in no band", {
  # IEEE 754, deliberately: a NaN compares FALSE against every threshold, so the record is simply
  # absent rather than erroring. Note base R would give NA here, which needs a third logical state
  # logmu does not have.
  b <- bands(.i$pension, thresholds = c(5000, 10000))
  for (i in seq_len(length(b))) {
    expect_equal(period_included(b[[i]], .i = list(pension = NaN)), include_none_interval())
  }
})

# ---- duration, the general form of age -------------------------------------

# THE RULE THESE PIN is that `duration()` is sugar and not a second mechanism. Every one of them
# compares the shorthand against the general `band(.t - origin, ...)` form or against `age()`,
# which is the duration whose origin happens to be birth.

test_that("duration constructors follow the family's singular and plural rule", {
  expect_true(is_include(duration(.i$entry, 0, 10)))
  expect_true(is_includes(durations(.i$entry, 0, 10)))
  expect_equal(length(durations(.i$entry, 0, 10, by = 5)), 2L)
})

test_that("duration bands are named by their edges and grouped by their origin", {
  d <- durations(.i$entry, 0, 10, by = 5)
  expect_equal(names(d), c("0-5", "5-10"))
  # Two duration sets in one breakdown are told apart only by the origin, so it is part of the
  # group name rather than a bare "duration". A chart axis is where a collision would surface.
  expect_equal(group_names(d), rep("duration since entry", 2L))
  expect_equal(group_names(durations(.i$retirement, 0, 5)), "duration since retirement")
})

test_that("duration bands clip exposure from the origin", {
  d <- durations(.i$entry, 0, 10, by = 5)
  .i <- list(entry = datey::datey(2000))
  expect_equal(period_included(d[[1]], .i = .i),
               datey::datey_interval(datey::datey(2000), datey::datey(2005)))
  expect_equal(period_included(d[[2]], .i = .i),
               datey::datey_interval(datey::datey(2005), datey::datey(2010)))
})

test_that("the shorthand builds exactly what the general form builds", {
  expect_equal(duration(.i$entry, 0, 5)$terms, band(.t - .i$entry, 0, 5)$terms)
})

test_that("age is the duration whose origin is birth", {
  expect_equal(duration(.i$birth, 65, 95)$terms, age(65, 95)$terms)
})

test_that("the origin may be computed, not only a bare field", {
  d <- duration(min(.i$entry, .i$retirement), 0, 5)
  expect_equal(
    period_included(d, .i = list(entry = datey::datey(2000),
                                 retirement = datey::datey(2010))),
    datey::datey_interval(datey::datey(2000), datey::datey(2005))
  )
})

test_that("an origin that moves with time is refused at its own argument", {
  # Left to `band_variable`, this would reach the user as a complaint about banded variables when
  # what they got wrong is one argument of `durations()`.
  expect_error(duration(.t, 0, 5), "origin")
  expect_error(durations(.i$entry + .t, 0, 5), "origin")
})

test_that("duration bands take the family's edge rules unchanged", {
  expect_error(durations(.i$entry, 0, 11, by = 5), "divide")
  expect_equal(names(durations(.i$entry, thresholds = c(5, 10))),
               c("< 5", "5-10", ">= 10"))
})

# ---- the includes collection -----------------------------------------------

test_that("includes() flattens and keeps both levels of naming", {
  std <- includes(ages(65, 95, by = 15), periods(2000, 2020, by = 10))

  expect_equal(length(std), 4L)
  expect_equal(names(std), c("65-80", "80-95", "2000-2010", "2010-2020"))
  expect_equal(group_names(std), c("age", "age", "period", "period"))
})

test_that("naming an includes argument renames its group", {
  std <- includes(cohort = ages(65, 95, by = 15))
  expect_equal(group_names(std), c("cohort", "cohort"))
  expect_equal(names(std), c("65-80", "80-95"))
})

test_that("naming an include argument renames that element", {
  std <- includes(young = age(65, 80), old = age(80, 95))
  expect_equal(names(std), c("young", "old"))
  expect_equal(group_names(std), c("age", "age"))
})

test_that("an indicator can be an element and stays an indicator", {
  # Relabelling must not demote it to a plain include, which would lose the {0,1} guarantee behind
  # V = E.
  std <- includes(male = include(.i$sex == "male"))
  expect_true(is_indicator(std[[1]]))
  expect_equal(names(std), "male")
})

test_that("subsetting keeps the type", {
  std <- includes(ages(65, 95, by = 15), periods(2000, 2020, by = 10))
  expect_true(is_includes(std[1:2]))
  expect_equal(length(std[1:2]), 2L)
  expect_true(is_include(std[[1]]))
  expect_false(is_includes(std[[1]]))
})

test_that("renaming writes through to the elements", {
  # THE ELEMENTS ARE THE SOURCE OF TRUTH and the list names are a cache. A rename that touched only
  # the cache would revert the moment anything rebuilt from the elements, silently.
  std <- ages(65, 95, by = 15)
  names(std)[1] <- "young"

  expect_equal(names(std), c("young", "80-95"))
  expect_equal(std[[1]]$name, "young")
  expect_equal(names(includes(std)), c("young", "80-95"))
})

test_that("group_names<- writes through, and recycles a single value", {
  std <- bands(.i$pension, thresholds = 5000)
  group_names(std) <- "Pension amount"

  expect_equal(group_names(std), rep("Pension amount", 2L))
  expect_equal(std[[1]]$group_name, "Pension amount")
  expect_equal(group_names(includes(std)), rep("Pension amount", 2L))

  group_names(std) <- c("a", "b")
  expect_equal(group_names(std), c("a", "b"))
})

test_that("names<- does not recycle", {
  # `names(x)[1] <- "y"` on an includes with no names yet yields a length-1 vector, and recycling
  # that would rename every element instead of the first.
  std <- ages(65, 95, by = 15)
  expect_error({ names(std) <- "only one" }, "one name per include")
  expect_error({ group_names(std) <- c("a", "b", "c") }, "one per include")
})

test_that("renaming keeps the type and the other level", {
  std <- includes(ages(65, 95, by = 15), periods(2000, 2020, by = 10))
  names(std) <- c("w", "x", "y", "z")

  expect_true(is_includes(std))
  expect_equal(group_names(std), c("age", "age", "period", "period"))
})

test_that("includes() refuses anything that is not an include", {
  expect_error(includes(42), "must be an `include`")
  expect_error(includes(ages(65, 95), "x"), "must be an `include`")
})

test_that("include() builds an indicator and refuses time", {
  expect_true(is_indicator(include(.i$sex == "male")))
  expect_true(is_include(include(.i$sex == "male")))
  expect_error(include(.t > 2010), "must not depend on time")
})
