# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the R -> veil crossing: the it_node walk, the literal typing, and the object-leaf lowering.
# The risky code here is the reading of R objects -- attributes, class dispatch, the column-major
# matrix, term unpacking -- so these tests use real objects rather than hand-built lists.

cols <- c("birth", "sex", "pension")

test_that("a pronoun expression crosses as the right calls", {
  res <- cpp_veil_ingest_ast(it_ast(~ .i$pension > 0), cols)

  expect_equal(res$monikers, "GT")
  expect_equal(res$node_count, 3L) # field, lit, call
  expect_equal(res$root, 2L)       # a call builds its children first
})

test_that("`.x` desugars to time minus birth", {
  res <- cpp_veil_ingest_ast(it_ast(~ .x), cols)

  expect_equal(res$monikers, "sub")
  expect_true("time" %in% res$kinds)
  expect_true("field" %in% res$kinds)
})

test_that("unary and binary minus map to different ops", {
  expect_equal(cpp_veil_ingest_ast(it_ast(~ -.i$pension), cols)$monikers, "neg")
  expect_equal(cpp_veil_ingest_ast(it_ast(~ .i$pension - 1), cols)$monikers, "sub")
})

test_that("a bare integer literal is a count, not a datey", {
  # datey and durationy are both integer clicks; only the R class tells them apart, and an
  # unclassed integer is a count.
  res <- cpp_veil_ingest_ast(it_ast(~ .i$pension > 5L), cols)

  expect_equal(res$lit_types, "double")
})

test_that("datey and durationy literals keep their own types", {
  datey_res <- cpp_veil_ingest_ast(it_ast(~ .t > datey::datey(2010)), cols)
  expect_equal(datey_res$lit_types, "datey")

  durationy_res <- cpp_veil_ingest_ast(it_ast(~ .x > datey::durationy(65)), cols)
  expect_equal(durationy_res$lit_types, "durationy")
})

test_that("a mortality_const lowers to a double literal", {
  # The first argument is `q`, so log_mu must be named -- an unnamed -4.5 would be read as a
  # mortality rate and fold to log(-4.5).
  m <- mortality_const(log_mu = -4.5)
  res <- cpp_veil_ingest_ast(it_obj(m), cols)

  expect_equal(res$kinds, "lit")
  expect_equal(res$lit_types, "double")
  expect_length(res$objs, 0) # a const needs no object; it is just a scalar
})

test_that("a mortality_table lowers to vector_log_mu(table, birth)", {
  log_mu <- matrix(log(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)), nrow = 3, ncol = 2)
  tbl <- mortality_table(x0 = 60, t0 = 2010, log_mu = log_mu)

  res <- cpp_veil_ingest_ast(it_obj(tbl), cols)

  expect_equal(res$monikers, "vector_log_mu")
  expect_true("obj" %in% res$kinds)
  expect_true("field" %in% res$kinds) # birth

  expect_length(res$objs, 1)
  obj <- res$objs[[1]]
  expect_equal(obj$kind, "mortality_table")
  expect_equal(obj$age_count, 3L)
  expect_equal(obj$period_count, 2L)
})

test_that("a mortality_table's log_mu keeps its orientation", {
  # The core indexes log_mu[age + period * age_count], which is R's own column-major order. If the
  # two ever disagree the table is silently transposed, so compare the values directly.
  log_mu <- matrix(c(-4.6, -4.5, -4.4, -4.3, -4.2, -4.1), nrow = 3, ncol = 2)
  tbl <- mortality_table(x0 = 60, t0 = 2010, log_mu = log_mu)

  obj <- cpp_veil_ingest_ast(it_obj(tbl), cols)$objs[[1]]

  expect_equal(obj$log_mu, as.vector(log_mu))
  # Spot-check the indexing arithmetic itself: age 1, period 1 (0-based).
  expect_equal(obj$log_mu[[1 + 1 * obj$age_count + 1]], log_mu[2, 2])
})

test_that("age() lowers to a single offset term on birth", {
  res <- cpp_veil_ingest_ast(it_obj(age(65, 95)), cols)

  expect_length(res$objs, 1)
  inc <- res$objs[[1]]
  expect_equal(inc$kind, "include")
  expect_length(inc$terms, 1)

  term <- inc$terms[[1]]
  expect_equal(term$kind, "offset")
  expect_lt(term$from_clicks, term$to_clicks)

  # The offset is a node in the same tree, not a column index. `age()` is `.t - .i$birth`, so the
  # node it points at is a field -- and it is the only field the tree has, which is what stops this
  # passing for whatever node happened to land at that id.
  expect_equal(res$kinds[[term$offset + 1L]], "field")
  expect_equal(sum(res$kinds == "field"), 1L)
})

test_that("period() lowers to a single absolute term", {
  res <- cpp_veil_ingest_ast(it_obj(period(2010, 2020)), cols)

  inc <- res$objs[[1]]
  expect_length(inc$terms, 1)
  expect_equal(inc$terms[[1]]$kind, "absolute")
  expect_lt(inc$terms[[1]]$from_clicks, inc$terms[[1]]$to_clicks)
})

test_that("an intersected include carries both a gate and its bounds", {
  male <- indicator(.i$sex == "male")
  res <- cpp_veil_ingest_ast(it_obj(male & age(65, 95)), cols)

  inc <- res$objs[[1]]
  kinds <- vapply(inc$terms, function(t) t$kind, character(1))
  expect_setequal(kinds, c("gate", "offset"))

  # The gate's own expression is ingested into the same tree, so its condition is there too.
  expect_true("EQ" %in% res$monikers)
})

test_that("a mortality expression splices rather than nesting an object", {
  base <- mortality_const(log_mu = -4.5)
  res <- cpp_veil_ingest_ast(it_obj(mortality(base + 0.05)), cols)

  expect_equal(res$monikers, "add")
  expect_equal(sum(res$kinds == "lit"), 2L) # the const and the 0.05, both spliced in
})

test_that("is.na crosses as an op", {
  res <- cpp_veil_ingest_ast(it_ast(~ is.na(.i$pension)), cols)

  expect_equal(res$monikers, "is_na")
})

test_that("`%in%` desugars to a chain of equality tests", {
  # `x %in% c(a, b)` means `x == a | x == b`, so it needs no op and no set-valued operand.
  res <- cpp_veil_ingest_ast(it_ast(~ .i$sex %in% c("male", "female")), cols)

  expect_setequal(res$monikers, c("EQ", "or"))
  expect_equal(sum(res$monikers == "EQ"), 2L)
  expect_equal(sum(res$monikers == "or"), 1L)

  # The left-hand side is built once and shared by both comparisons -- a DAG, not a tree.
  expect_equal(sum(res$kinds == "field"), 1L)
  expect_equal(res$lit_types, c("text", "text"))
})

test_that("a single-element `%in%` is just an equality test", {
  res <- cpp_veil_ingest_ast(it_ast(~ .i$sex %in% "male"), cols)

  expect_equal(res$monikers, "EQ")
})

test_that("an empty `%in%` set is rejected by the parser", {
  # it_fold_set turns this away before it can cross, so the binding's own empty-set handling
  # (membership of nothing is false, never NA) is unreachable from R and stays only as a guard
  # for a front end that is less strict.
  expect_error(it_ast(~ .i$sex %in% character(0)), "character or numeric")
})

test_that("a numeric `%in%` set folds its elements to doubles", {
  res <- cpp_veil_ingest_ast(it_ast(~ .i$pension %in% c(1L, 2L)), cols)

  expect_equal(res$lit_types, c("double", "double"))
})

test_that("if and ifelse both cross as a selection", {
  both <- c("select", "GT")

  expect_setequal(cpp_veil_ingest_ast(it_ast(~ if (.i$pension > 0) 1 else 2), cols)$monikers, both)
  expect_setequal(cpp_veil_ingest_ast(it_ast(~ ifelse(.i$pension > 0, 1, 2)), cols)$monikers, both)
})

test_that("an if without an else is rejected", {
  # A selection needs a value on both paths.
  expect_error(cpp_veil_ingest_ast(it_ast(~ if (.i$pension > 0) 1), cols), "argument")
})

test_that("the crossing rejects what it cannot lower", {
  expect_error(cpp_veil_ingest_ast(it_ast(~ .i$missing > 0), cols), "not a column")
})
