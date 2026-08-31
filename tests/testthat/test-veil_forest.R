# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the forest and the sharing pass: several expressions compiled as one calculation, with
# sub-expressions they have in common computed once.
#
# THE ANSWERS ALONE CANNOT TELL YOU WHETHER SHARING HAPPENED -- they are the same either way. So the
# tests assert two things together: that every output is what R computes for that expression, and that
# the INSTRUCTION COUNT is what sharing implies. Checking only the second would pass for a pass that
# merged things it should not have, and checking only the first would pass for one that did nothing.

cols <- list(
  age    = c(30, 45, 60),
  weight = c(1.5, 2.5, 0.5),
  height = c(1.6, 1.8, 1.7),
  birth  = datey::datey(c(1960, 1965, 1975))
)

eval_multi <- function(...) cpp_veil_eval_multi(list(...), cols)

test_that("several expressions each get their own output, in order", {
  res <- eval_multi(it_ast(~ .i$age + 1), it_ast(~ .i$weight * 2), it_ast(~ .i$height))

  expect_equal(res$output_count, 3L)
  expect_equal(res$values[[1]], cols$age + 1)
  expect_equal(res$values[[2]], cols$weight * 2)
  expect_equal(res$values[[3]], cols$height)
  expect_equal(res$types, c("double", "double", "double"))
})

test_that("a sub-expression two outputs share is computed once", {
  # `.i$age * .i$weight` appears in both. Sharing it means one multiply, not two.
  shared <- eval_multi(
    it_ast(~ .i$age * .i$weight + 1),
    it_ast(~ .i$age * .i$weight + 2)
  )
  expect_equal(shared$values[[1]], cols$age * cols$weight + 1)
  expect_equal(shared$values[[2]], cols$age * cols$weight + 2)

  # The same two expressions with nothing in common, as the yardstick: four instructions there
  # against three here, the difference being the multiply that was shared.
  distinct <- eval_multi(
    it_ast(~ .i$age * .i$weight + 1),
    it_ast(~ .i$age * .i$height + 2)
  )

  expect_equal(shared$instruction_count, distinct$instruction_count - 1L)
  expect_gt(shared$shared_nodes, 0L)
})

test_that("two identical outputs collapse to one instruction each way", {
  res <- eval_multi(it_ast(~ .i$age * 2), it_ast(~ .i$age * 2))

  expect_equal(res$output_count, 2L)
  expect_equal(res$values[[1]], cols$age * 2)
  expect_equal(res$values[[2]], cols$age * 2)

  # Both roots became the same node, so the block holds one multiply and names it twice.
  expect_equal(res$instruction_count, 1L)
})

test_that("sharing does not merge expressions that only look alike", {
  # Same shape, different field. Nothing may be shared but the reads of `age`.
  res <- eval_multi(it_ast(~ .i$age - .i$weight), it_ast(~ .i$age - .i$height))

  expect_equal(res$values[[1]], cols$age - cols$weight)
  expect_equal(res$values[[2]], cols$age - cols$height)
  expect_equal(res$instruction_count, 2L)
})

test_that("subtraction is not treated as commutative", {
  res <- eval_multi(it_ast(~ .i$age - .i$weight), it_ast(~ .i$weight - .i$age))

  expect_equal(res$values[[1]], cols$age - cols$weight)
  expect_equal(res$values[[2]], cols$weight - cols$age)
  expect_equal(res$instruction_count, 2L)
})

test_that("association is not rewritten, so differently spelled products stay apart", {
  # `(a*b)*c` and `a*(b*c)` are the same number and different subtrees. The pass shares nodes; it
  # does not reassociate, and it is not meant to.
  res <- eval_multi(
    it_ast(~ (.i$age * .i$weight) * .i$height),
    it_ast(~ .i$age * (.i$weight * .i$height))
  )

  expect_equal(res$values[[1]], (cols$age * cols$weight) * cols$height)
  expect_equal(res$values[[2]], cols$age * (cols$weight * cols$height))
  expect_equal(res$instruction_count, 4L)
})

test_that("a literal of one type does not share with an equal literal of another", {
  # 1990 as a plain number and 1990 as a datey are not the same value, whatever they print as.
  res <- eval_multi(it_ast(~ .i$age + 1990), it_ast(~ .i$birth > 1990))

  expect_equal(res$values[[1]], cols$age + 1990)
  expect_equal(res$values[[2]], as.double(unclass(cols$birth) > 1990 * 534360))
})

test_that("outputs of different types travel together", {
  res <- eval_multi(it_ast(~ .i$age > 40), it_ast(~ .i$birth), it_ast(~ .i$weight))

  expect_equal(res$types, c("bool", "datey", "double"))
  expect_equal(res$values[[1]], as.double(cols$age > 40))
  expect_equal(res$values[[2]], as.double(unclass(cols$birth)))
  expect_equal(res$values[[3]], cols$weight)
})

test_that("a single root still works, and asking a forest for `the` root does not", {
  res <- eval_multi(it_ast(~ .i$age * 3))
  expect_equal(res$output_count, 1L)
  expect_equal(res$values[[1]], cols$age * 3)

  expect_error(cpp_veil_eval_multi(list(), cols), "at least one expression")
})

test_that("a NaN literal does not merge with another double literal", {
  # THE SHARING KEY IS A `std::map` KEY, SO IT NEEDS A TOTAL ORDER, and `<` on doubles does not
  # give one: NaN is unordered against everything, so a tuple comparison finds `!(NaN < x)` and
  # `!(x < NaN)` and calls them EQUIVALENT. Before the key was changed to compare BITS, a literal
  # NaN merged with whatever double literal it met on the way down the map -- undefined behaviour
  # by the container's own precondition, with a wrong answer attached (test-veil_aev.R holds it).
  #
  # BUILT AS RAW NODES rather than through `it_ast`, because the R front end now refuses a written
  # NaN at construction. The core has to be right on its own terms regardless: it is deliberately
  # R-free, and another front end need not make the same check.
  res <- cpp_veil_eval_multi(list(list(kind = "lit", value = NaN),
                                  list(kind = "lit", value = 1)), cols)

  # Asserted in BOTH directions. Whichever of the two the map happens to insert first is the one
  # that survives the bad merge, so a test looking at only one output passes half the time.
  expect_true(all(is.nan(res$values[[1]])))
  expect_equal(res$values[[2]], rep(1, length(cols$age)))

  # The direct witness: nothing was merged at all. Without it the assertions above would still
  # hold for a pass that merged the two into a third thing that happened to print correctly.
  expect_equal(res$shared_nodes, 0L)
})

test_that("two literals holding the same double still share", {
  # The other half of the same rule. Keying on bits must not stop ordinary sharing working, or the
  # test above would pass for a pass that had simply been switched off.
  res <- cpp_veil_eval_multi(list(list(kind = "lit", value = 2.5),
                                  list(kind = "lit", value = 2.5)), cols)

  expect_equal(res$values[[1]], rep(2.5, length(cols$age)))
  expect_gt(res$shared_nodes, 0L)
})
