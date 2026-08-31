# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests text handling end to end: a column of strings reaches veil as a factor, the binding merges
# every factor's levels into ONE numbering for the crossing, and `passEncodeText` resolves each text
# literal against that numbering before anything reasons about values.
#
# The two halves are tested separately on purpose. `cpp_veil_prepare` shows what the encoding pass did
# to the tree; `aev()` shows that the answer is right. A tree assertion alone would pass if the codes
# were consistent but meant the wrong strings.

cols <- list(
  birth  = datey::datey(1960),
  sex    = factor("male"),
  salary = 30000,
  smoker = TRUE
)

lit_types <- function(res) res$types[res$kinds == "lit"]
lit_values <- function(res) res$values[res$kinds == "lit"]
root_kind <- function(res) res$kinds[[res$root + 1L]]
root_value <- function(res) res$values[[res$root + 1L]]

# ---- the encoding pass -----------------------------------------------------

test_that("a text literal compared with a factor becomes a category code", {
  res <- cpp_veil_prepare(it_ast(~ .i$sex == "male"), cols)

  # No node is text any more: the literal is a category index, which is what the engine compares.
  expect_false("text" %in% res$types)
  expect_equal(lit_types(res), "category")
})

# The other direction through the mapping, and the only one the core uses: an index rendered back to
# the string it stands for. Without it the tree dump would account for `.i$sex == "male"` against a
# bare `0` whose meaning lives in a table the reader cannot see.
test_that("an encoded literal is reported as its string, not as its index", {
  expect_equal(lit_values(cpp_veil_prepare(it_ast(~ .i$sex == "male"), cols)), "male")
})

test_that("the reported string follows the string, not the column's own level order", {
  # "male" is the second level of `sex` and the first of `spouse`, so a literal rendered from either
  # column's own numbering would disagree with the other. Both report "male".
  two <- list(
    birth  = cols$birth,
    sex    = factor("male", levels = c("female", "male")),
    spouse = factor("male", levels = c("male", "female"))
  )

  expect_equal(lit_values(cpp_veil_prepare(it_ast(~ .i$sex == "male"), two)), "male")
  expect_equal(lit_values(cpp_veil_prepare(it_ast(~ .i$spouse == "male"), two)), "male")
})

test_that("a text literal naming no level folds the comparison to a constant", {
  # Not an error: a typo, or a level present in one dataset of a batch and absent from another, is
  # ordinary. No record can satisfy it, so the answer is settled before the engine runs.
  eq <- cpp_veil_prepare(it_ast(~ .i$sex == "mail"), cols)
  expect_equal(root_kind(eq), "lit")
  expect_equal(root_value(eq), "FALSE")

  ne <- cpp_veil_prepare(it_ast(~ .i$sex != "mail"), cols)
  expect_equal(root_kind(ne), "lit")
  expect_equal(root_value(ne), "TRUE")
})

test_that("`%in%` encodes every element of its set", {
  two <- list(birth = cols$birth, sex = factor(c("male", "female")))
  res <- cpp_veil_prepare(it_ast(~ .i$sex %in% c("male", "female")), two)

  expect_false("text" %in% res$types)
  expect_equal(lit_types(res), c("category", "category"))
})

test_that("a `%in%` set of levels that are all absent folds every element away", {
  res <- cpp_veil_prepare(it_ast(~ .i$sex %in% c("mail", "femail")), cols)

  # Each comparison settled on its own, so what is left is `FALSE | FALSE` and no comparison at all.
  expect_false("EQ" %in% res$labels)
  expect_equal(res$labels[[res$root + 1L]], "or")

  # The literals the folded comparisons used to read are still in the arena and still text. They are
  # unreachable -- lowering walks from the roots and the include gates -- so they are left alone
  # rather than rewritten, since a literal may have another use this pass is not looking at.
  expect_equal(sum(res$types == "text"), 2L)
})

test_that("a character column is refused, naming the fix", {
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex == "male"), list(sex = "male")),
               "must arrive as a factor")
})

# ---- ordering --------------------------------------------------------------

# Categories are unordered, and the refusal has to speak in the terms the user wrote. Someone writing
# `.i$sex < "male"` wrote strings and has never heard of a category -- it is the name of the encoded
# form they never see.
refusal_message <- function(expr) {
  tryCatch({expr; NA_character_}, error = conditionMessage)
}

test_that("ordering a text column is refused in the user's own terms", {
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex < "male"), cols), "cannot order text")
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex >= "male"), cols), "cannot order text")

  # Reversed, so the refusal follows the operand rather than which side it sits on.
  expect_error(cpp_veil_prepare(it_ast(~ "male" > .i$sex), cols), "cannot order text")

  # "category" is the name of the encoded form the user never sees, so it must not appear in a
  # message about an expression they wrote in strings.
  expect_false(grepl("category", refusal_message(cpp_veil_prepare(it_ast(~ .i$sex < "male"), cols))))
})

# Ordering two factor columns has no literal in it at all, so this is the case that the ban on
# Category rather than on Text is what catches.
test_that("ordering one text column against another is refused too", {
  two <- list(birth = cols$birth, sex = factor("male"), scheme = factor("A"))

  expect_error(cpp_veil_prepare(it_ast(~ .i$sex < .i$scheme), two), "cannot order text")

  equality <- cpp_veil_prepare(it_ast(~ .i$sex == .i$scheme), two)
  expect_equal(equality$types[[equality$root + 1L]], "bool")
})

# ---- select ----------------------------------------------------------------

# A select may produce a category, because both branches are then already indices into the shared
# mapping. It may not choose between written strings: the mapping is built from the columns before
# any expression is typed, so a string no column holds has no index to be chosen.
test_that("a select between two text columns is a category", {
  two <- list(birth = cols$birth, smoker = TRUE, sex = factor("male"), spouse = factor("female"))
  res <- cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, .i$sex, .i$spouse)), two)

  expect_equal(res$types[[res$root + 1L]], "category")
})

test_that("a select between written strings is refused, pointing at the rewrite", {
  expect_error(cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, "a", "b")), cols),
               "cannot choose between text values")

  # Mixed too, which is the form that reads most naturally and so is the one worth refusing clearly.
  expect_error(cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, .i$sex, "male")), cols),
               "cannot choose between text values")
})

# ---- the shared numbering --------------------------------------------------

exp_two_columns <- function(sex, spouse) {
  exp_data(
    list(
      birth     = datey::datey(rep(1950, length(sex))),
      sex       = sex,
      spouse    = spouse,
      E2R_start = datey::datey(rep(2015, length(sex))),
      E2R_end   = datey::datey(rep(2020, length(sex))),
      E2R_died  = rep(FALSE, length(sex))
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
}

basis <- settings(overdispersion = 1)
flat <- mortality_const(log_mu = -4)

exposure_of <- function(data, ...) unclass(aev(data, mortality = flat, settings = basis, ...))$E

# THE TEST THE GLOBAL NUMBERING EXISTS FOR. The two columns hold the same strings but number their
# own levels in opposite orders, so comparing their raw factor codes gets both rows wrong -- and gets
# them wrong the same way, which is why it would not look like a bug. Only a numbering shared across
# the crossing makes "male" mean one thing in both columns.
test_that("two factors with opposite level orders still compare by their strings", {
  data <- exp_two_columns(
    sex    = factor(c("male", "female"), levels = c("male", "female")),
    spouse = factor(c("male", "female"), levels = c("female", "male"))
  )

  both <- exposure_of(data)
  matched <- exposure_of(data, weight = .i$sex == .i$spouse)

  expect_equal(matched, both)
})

test_that("a text literal means the same thing in every column of the crossing", {
  data <- exp_two_columns(
    sex    = factor(c("male", "female"), levels = c("male", "female")),
    spouse = factor(c("female", "male"), levels = c("female", "male"))
  )

  # Row 1 is a male with a female spouse, row 2 the reverse, so each gate takes exactly one row and
  # the two together take the lot.
  by_sex <- exposure_of(data, include = include(.i$sex == "male"))
  by_spouse <- exposure_of(data, include = include(.i$spouse == "male"))

  expect_equal(by_sex + by_spouse, exposure_of(data))
})

# ---- end to end over experience data ---------------------------------------

exp_one_column <- function(sex) {
  exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955)),
      sex       = sex,
      E2R_start = datey::datey(rep(2015, 3)),
      E2R_end   = datey::datey(c(2020, 2020, 2018)),
      E2R_died  = c(FALSE, FALSE, TRUE)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
}

test_that("an include on a text column splits the exposure exactly", {
  data <- exp_one_column(c("male", "female", "male"))

  male <- exposure_of(data, include = include(.i$sex == "male"))
  female <- exposure_of(data, include = include(.i$sex == "female"))

  expect_gt(male, 0)
  expect_gt(female, 0)
  expect_equal(male + female, exposure_of(data))
})

test_that("`!=` is the complement of `==` over a text column", {
  data <- exp_one_column(c("male", "female", "male"))

  expect_equal(exposure_of(data, include = include(.i$sex != "male")),
               exposure_of(data, include = include(.i$sex == "female")))
})

test_that("a level that no record holds gives an empty result rather than an error", {
  data <- exp_one_column(c("male", "female", "male"))
  result <- aev(data, mortality = flat, include = include(.i$sex == "mail"), settings = basis)

  expect_equal(unclass(result)$A, 0)
  expect_equal(unclass(result)$E, 0)
})

test_that("`%in%` over every level takes the whole population", {
  data <- exp_one_column(c("male", "female", "male"))

  expect_equal(exposure_of(data, include = include(.i$sex %in% c("male", "female"))),
               exposure_of(data))
})

# The select half of the ruling above, run rather than merely typed. Both rows are male and both
# spouses female, so the select decides the answer on its own: `.i$smoker` picks the member's sex for
# row 1 and the spouse's for row 2. `.i$smoker` is therefore an exact oracle, and `.i$sex == "male"`
# -- which takes both rows -- is what the answer must NOT be.
test_that("a category-valued select runs and picks the right column per record", {
  data <- exp_data(
    list(
      birth     = datey::datey(c(1950, 1950)),
      smoker    = c(TRUE, FALSE),
      sex       = factor(c("male", "male")),
      spouse    = factor(c("female", "female")),
      E2R_start = datey::datey(rep(2015, 2)),
      E2R_end   = datey::datey(rep(2020, 2)),
      E2R_died  = c(FALSE, FALSE)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )

  chosen <- exposure_of(data, weight = ifelse(.i$smoker, .i$sex, .i$spouse) == "male")

  expect_equal(chosen, exposure_of(data, weight = .i$smoker))
  expect_lt(chosen, exposure_of(data, weight = .i$sex == "male"))
})

test_that("a text comparison works as a weight as well as a gate", {
  data <- exp_one_column(c("male", "female", "male"))

  expect_equal(exposure_of(data, weight = .i$sex == "male"),
               exposure_of(data, include = include(.i$sex == "male")))
})

# ---- missing values --------------------------------------------------------

# A CATEGORY HAS NO MISSING VALUE. Its storage is `int`, which carries no platform-independent NA --
# R's NA_INTEGER is INT_MIN, which no other platform repeats, and a sentinel invented inside veil
# would be no more portable for being ours. Only double, datey and durationy have a missing state,
# and they have it because their own representations provide one.
#
# So an NA factor is refused where it is read, beside the identical refusal for a logical, and the
# dataset says "unknown" with a level of its own instead.

test_that("a factor column carrying an NA is refused", {
  data <- exp_one_column(c("male", NA, "male"))

  expect_error(exposure_of(data, include = include(.i$sex == "male")),
               "must not contain NA")
})

test_that("the refusal names the fix", {
  data <- exp_one_column(c("male", NA, "male"))

  expect_error(exposure_of(data, include = include(.i$sex == "male")), "addNA")
})

# THE REFUSAL IS EAGER, NOT LAZY, and deliberately recorded here because it is the surprising half.
# `prepareColumns` reads every column of the dataset, so an NA in a factor stops a calculation that
# never mentions that column. Identical to the missing-logical refusal beside it, which behaves the
# same way for the same reason -- this is a property of when columns are READ, not of this rule.
#
# `exp_data()` itself still accepts the NA; the complaint comes from the engine.
test_that("an NA factor stops even a calculation that does not reference it", {
  data <- exp_one_column(c("male", NA, "male"))

  expect_no_error(exp_one_column(c("male", NA, "male")))
  expect_error(exposure_of(data), "must not contain NA")
})

test_that("an explicit level is the way to carry an unknown", {
  known <- exp_one_column(c("male", "unknown", "male"))

  male <- exposure_of(known, include = include(.i$sex == "male"))
  unknown <- exposure_of(known, include = include(.i$sex == "unknown"))

  # The unknown group is a group like any other: it has exposure, it can be gated on, and it adds
  # back into the total rather than quietly disappearing from it.
  expect_gt(unknown, 0)
  expect_equal(male + unknown, exposure_of(known))
})
