# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests for `category()` and `categories()`.
#
# THE RULE MOST OF THIS PINS is that each argument of `...` is ONE group. A character vector joins
# its values rather than splitting into several groups, which is what makes `c("def", "dep")` a
# single group and `"def", "dep"` two. The rest pins the three checks and, above all, the principle
# behind their defaults: every check runs whenever it has what it needs, so `.known` and
# `.exhaustive` come on with `.source` and `.disjoint` is always on.

status_levels <- function() {
  factor(c("act", "def", "pen"), levels = c("act", "def", "pen", "dep"))
}

test_that("the singular returns an indicator and the plural an includes", {
  one <- category(.i$status, "act")
  expect_true(is_indicator(one))
  expect_true(is_include(one))          # an indicator is an include

  many <- categories(.i$status, "act", "def")
  expect_true(is_includes(many))
  expect_equal(length(many), 2L)
})

test_that("every group is an indicator, so V = OmegaE holds within each", {
  # Not decoration: a plain include would lose the {0,1} guarantee that V = E rests on.
  groups <- categories(.i$status, "act", "def")
  expect_true(all(vapply(groups, is_indicator, logical(1L))))
})

test_that("one argument is one group, and a vector joins rather than splits", {
  expect_equal(length(categories(.i$status, "def", "dep")), 2L)
  expect_equal(length(categories(.i$status, c("def", "dep"))), 1L)
})

test_that("a group tests membership", {
  joined <- category(.i$status, c("def", "dep"))
  expect_true(logical_value(joined, .i = list(status = "dep")))
  expect_true(logical_value(joined, .i = list(status = "def")))
  expect_false(logical_value(joined, .i = list(status = "act")))

  single <- category(.i$status, "act")
  expect_true(logical_value(single, .i = list(status = "act")))
  expect_false(logical_value(single, .i = list(status = "def")))
})

test_that("an indicator group resolves to all of time or to nothing", {
  group <- category(.i$status, "act")
  expect_equal(period_included(group, .i = list(status = "act")), include_all_interval())
  expect_equal(period_included(group, .i = list(status = "def")), include_none_interval())
})

test_that("names default to the values and are overridden by the ones given", {
  expect_equal(names(categories(.i$status, "act", "def")), c("act", "def"))
  expect_equal(names(categories(.i$status, c("def", "dep"))), "def+dep")
  expect_equal(names(categories(.i$status, active = "act", other = c("def", "dep"))),
               c("active", "other"))
  # A name given to only one of them leaves the other on its default.
  expect_equal(names(categories(.i$status, "act", other = c("def", "dep"))),
               c("act", "other"))
})

test_that("the group name is the field, without the pronoun", {
  expect_equal(group_names(categories(.i$status, "act", "def")), c("status", "status"))
})

test_that("two groups may not share a name", {
  expect_error(categories(.i$status, x = "act", x = "def"), "own name")
})

# ---- `.source` -------------------------------------------------------------

test_that("a factor source gives its declared levels, unused ones included", {
  # The whole reason for declared rather than observed: a breakdown reused across datasets must
  # give the same groups in the same order, or a comparison between two `aev`s lines up the wrong
  # rows -- or fails the length check outright.
  groups <- categories(.i$status, .source = status_levels())
  expect_equal(names(groups), c("act", "def", "pen", "dep"))
})

test_that("a character source gives its distinct values, sorted", {
  groups <- categories(.i$status, .source = c("pen", "act", "pen", NA))
  expect_equal(names(groups), c("act", "pen"))
})

test_that("a dataset source finds the column named by the field", {
  data <- data.frame(status = status_levels(), scheme = c("A", "B", "A"))
  expect_equal(names(categories(.i$status, .source = data)), c("act", "def", "pen", "dep"))
  expect_equal(names(categories(.i$scheme, .source = data)), c("A", "B"))
})

test_that("a dataset source needs a plain field and a column of that name", {
  data <- data.frame(status = status_levels())
  expect_error(categories(.i$missing, .source = data), "no column")
  expect_error(categories(ifelse(.i$flag, .i$status, .i$other), .source = data), "plain field")
})

test_that("a source that is not categorical is refused", {
  expect_error(categories(.i$status, .source = c(1, 2, 3)), "factor")
  expect_error(categories(.i$status, .source = factor(character(0))), "no categories")
})

# ---- the three checks ------------------------------------------------------

test_that("`.known` catches a mistyped category, and is on by default with a source", {
  expect_error(categories(.i$status, "act", "depp", "def", "pen", .source = status_levels()),
               "no category")
  expect_silent(categories(.i$status, "act", "depp", .known = FALSE, .exhaustive = FALSE,
                           .source = status_levels()))
})

test_that("`.exhaustive` catches a value left out of every group", {
  # The failure it exists for: the scheme gains a status next year and its exposure silently
  # vanishes from the breakdown. Nothing downstream can catch that, since a breakdown is not
  # required to be a partition.
  expect_error(categories(.i$status, "act", "def", "pen", .source = status_levels()),
               "in no group")
  expect_silent(categories(.i$status, "act", "def", "pen",
                           .source = status_levels(), .exhaustive = FALSE))
})

test_that("`.known` is reported before `.exhaustive`", {
  # A typo trips both -- "depp" is unknown and "dep" is uncovered -- and only the first tells the
  # user what to do about it.
  expect_error(categories(.i$status, "act", "def", "pen", "depp", .source = status_levels()),
               "no category")
})

test_that("`.disjoint` is on without any source at all", {
  expect_error(categories(.i$status, c("act", "def"), "def"), "more than one group")
  expect_silent(categories(.i$status, c("act", "def"), "def", .disjoint = FALSE))
})

test_that("a value repeated within one group is refused", {
  # Across groups a repeat means overlapping groups, which a user can mean. Within one group it
  # means nothing.
  expect_error(categories(.i$status, c("def", "def")), "repeats")
})

test_that("the within-group check is not `.disjoint` in disguise", {
  # Kept as its own test deliberately. Folding the within-group check into `.disjoint` fails the
  # test above only on its message, so this line is the sole witness that turning `.disjoint` off
  # does not also turn this off.
  expect_error(categories(.i$status, c("def", "def"), .disjoint = FALSE), "repeats")
})

test_that("a check that has nothing to check against is an error, not a no-op", {
  expect_error(categories(.i$status, "act", .known = TRUE), "need one")
  expect_error(categories(.i$status, "act", .exhaustive = TRUE), "need one")
})

test_that("the checks are vacuous when the groups came from the source", {
  # `.known` and `.exhaustive` default TRUE here and hold trivially, so the source-only form needs
  # no special case.
  expect_silent(categories(.i$status, .source = status_levels()))
})

test_that("the flags must be TRUE or FALSE", {
  expect_error(categories(.i$status, "act", .source = status_levels(), .known = "yes"),
               "TRUE or FALSE")
  expect_error(categories(.i$status, "act", .source = status_levels(), .known = NA),
               "TRUE or FALSE")
  expect_error(categories(.i$status, "act", .disjoint = c(TRUE, TRUE)), "TRUE or FALSE")
})

# ---- shapes refused --------------------------------------------------------

test_that("a category variable may not depend on time", {
  expect_error(categories(.t, "act"), "must not depend on time")
  expect_error(categories(ifelse(.x > 65, .i$a, .i$b), "act"), "must not depend on time")
})

test_that("groups must be character", {
  expect_error(categories(.i$status, 1), "character vector")
  expect_error(categories(.i$status, c("act", NA)), "character vector")
  expect_error(categories(.i$status, character(0)), "character vector")
})

test_that("with no groups and no source there is nothing to build", {
  expect_error(categories(.i$status), "or a `.source`")
  expect_error(category(.i$status), "needs one category")
})

test_that("`category()` takes exactly one group", {
  expect_error(category(.i$status, "def", "dep"), "Use `categories\\(\\)`")
  expect_true(is_indicator(category(.i$status, c("def", "dep"))))
})

test_that("`category()` with a source checks the value without demanding coverage", {
  # `.exhaustive` is deliberately absent from the singular: under the dependent default it would
  # refuse this, the ordinary single-group call.
  expect_silent(category(.i$status, "act", .source = status_levels()))
  expect_error(category(.i$status, "depp", .source = status_levels()), "no category")
})

# ---- interop ---------------------------------------------------------------

test_that("categories collect into a wider breakdown", {
  breakdown <- includes(
    categories(.i$status, "act", "def", "pen", "dep"),
    ages(65, 85, by = 10)
  )
  expect_equal(length(breakdown), 6L)
  expect_equal(group_names(breakdown), c(rep("status", 4L), rep("age", 2L)))
})

test_that("a category intersects with another include", {
  both <- category(.i$status, "act") & age(65, 95)
  expect_true(is_include(both))
})

# ---- end to end over experience data ---------------------------------------

exp_status <- function(status) {
  exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955)),
      status    = status,
      E2R_start = datey::datey(rep(2015, 3)),
      E2R_end   = datey::datey(rep(2020, 3)),
      E2R_died  = rep(FALSE, 3)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
}

test_that("a categorical breakdown divides the exposure, and a declared level with no records is kept", {
  # The end-to-end witness, and the one that shows declared levels earning their keep: "pen" has no
  # records, so a breakdown built from observed values would be one row shorter and would not line
  # up with the same breakdown run over other data.
  data <- exp_status(factor(c("act", "def", "act"), levels = c("act", "def", "pen")))
  flat <- mortality_const(log_mu = -4)
  basis <- settings(overdispersion = 1)

  total <- unclass(aev(data, mortality = flat, settings = basis))$E
  groups <- unclass(aev(data, mortality = flat, settings = basis,
                        breakdown = categories(.i$status, .source = data)))$E

  expect_equal(length(groups), 3L)
  expect_equal(sum(groups), total)   # disjoint and exhaustive, so the parts are the whole
  expect_equal(groups[[3L]], 0)
})

test_that("a joined group is the sum of the values it joins", {
  data <- exp_status(factor(c("act", "def", "pen"), levels = c("act", "def", "pen")))
  flat <- mortality_const(log_mu = -4)
  basis <- settings(overdispersion = 1)

  separate <- unclass(aev(data, mortality = flat, settings = basis,
                          breakdown = categories(.i$status, "def", "pen")))$E
  joined <- unclass(aev(data, mortality = flat, settings = basis,
                        breakdown = category(.i$status, c("def", "pen"))))$E

  expect_equal(joined, sum(separate))
})
