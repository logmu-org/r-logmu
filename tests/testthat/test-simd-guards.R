# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# TWO LISTS THAT MUST MATCH, AND NOTHING ELSE CHECKS THAT THEY DO.
#
# The build says which CPU features each SIMD kernel is ALLOWED TO USE: the `-m`
# flags in `src/Makevars.win` and in the rules `configure` writes. The dispatch
# says which features the CPU MUST HAVE before that kernel is called: the
# `__builtin_cpu_supports` calls in `src/vec_ops/vec_op_dispatch.cpp`.
#
# A compiler may emit any instruction the flags permit. So if the build permits a
# feature the guard does not test for, the kernel can execute an instruction the
# processor does not implement, and the process dies on the spot -- no R error,
# no traceback, just a dead session.
#
# THIS IS NOT HYPOTHETICAL. The AVX-512 kernel was built with five features and
# guarded on one (`avx512f`) from the first version until 2026-08-31. A CPU with
# F but without CD, BW, DQ or VL would have taken that path and crashed.
#
# The source comments say the flags and the guard must be changed together. A
# comment is not a check; this is. It passes today and exists for the day someone
# adds `-mavx512vnni` without knowing the guard is there.
#
# ON READING BUILD FILES WITH REGULAR EXPRESSIONS. It is not elegant, and the
# failure mode is the right way round: a parsing change makes this shout, and
# somebody looks. The alternative -- no check -- fails silently on a user's
# machine, which is what happened.

# `-mfoo` flags on the compile line for one kernel.
simd_build_flags <- function(text, kernel) {
  line <- grep(kernel, text, value = TRUE, fixed = TRUE)
  line <- grep("-m", line, value = TRUE, fixed = TRUE)
  if (length(line) == 0L) return(character(0))
  sort(unique(sub("^-m", "", unlist(regmatches(line, gregexpr("-m[a-z0-9]+", line))))))
}

# Feature names inside `__builtin_cpu_supports("...")`, split into the two
# guarded blocks. The avx512 block returns `avx512_tbl`, the avx2 one `avx2_tbl`.
simd_guard_features <- function(text, table) {

  # ANCHOR ON THE `return`, NOT THE NAME. Each table is also DECLARED earlier in
  # the same function, and matching the declaration puts the start of the block
  # before any `if` at all.
  stop_at <- grep(paste0("return ", table), text, fixed = TRUE)
  expect_length(stop_at, 1L)

  opens <- grep("^\\s*if \\(", text[seq_len(stop_at)])
  expect_gt(length(opens), 0L)
  block <- text[max(opens):stop_at]

  sort(unique(unlist(regmatches(
    block, gregexpr('(?<=cpu_supports\\(")[a-z0-9]+(?=")', block, perl = TRUE)
  ))))
}

read_source <- function(path) readLines(test_path("..", "..", path), warn = FALSE)

test_that("every feature the AVX2 kernel is built with is also guarded", {

  flags <- simd_build_flags(read_source("src/Makevars.win"), "vec_kernel_avx2.cpp")
  guard <- simd_guard_features(read_source("src/vec_ops/vec_op_dispatch.cpp"), "avx2_tbl")

  expect_gt(length(flags), 0L)
  expect_identical(
    flags, guard,
    info = paste0(
      "built with: ", paste(flags, collapse = ", "),
      " | guarded on: ", paste(guard, collapse = ", ")
    )
  )
})

test_that("every feature the AVX-512 kernel is built with is also guarded", {

  # The one that was actually wrong: five flags, one guard, until 2026-08-31.
  flags <- simd_build_flags(read_source("src/Makevars.win"), "vec_kernel_avx512.cpp")
  guard <- simd_guard_features(read_source("src/vec_ops/vec_op_dispatch.cpp"), "avx512_tbl")

  expect_gt(length(flags), 0L)
  expect_identical(
    flags, guard,
    info = paste0(
      "built with: ", paste(flags, collapse = ", "),
      " | guarded on: ", paste(guard, collapse = ", ")
    )
  )
})

test_that("the Unix build uses the same flags as the Windows one", {

  # `configure` writes the Unix rules and `Makevars.win` holds the Windows ones.
  # They are separate files that must agree, or a tier is guarded correctly on
  # one platform and not the other.
  configure <- read_source("configure")
  windows <- read_source("src/Makevars.win")

  for (kernel in c("vec_kernel_avx2.cpp", "vec_kernel_avx512.cpp")) {
    expect_identical(
      simd_build_flags(configure, kernel),
      simd_build_flags(windows, kernel),
      info = kernel
    )
  }
})

test_that("the baseline kernel is built with no feature flags at all", {

  # It must run anywhere, so it carries no `-m` flags and needs no guard. If it
  # ever gains one, it needs a guard too -- and this test says so.
  expect_identical(
    simd_build_flags(read_source("src/Makevars.win"), "vec_kernel_baseline.cpp"),
    character(0)
  )
  expect_identical(
    simd_build_flags(read_source("configure"), "vec_kernel_baseline.cpp"),
    character(0)
  )
})

test_that("the tier the guards select is the tier that runs", {

  # The lists agreeing is the static half. This is the dynamic half: whatever
  # `select_tier()` chose is what `vec_active_tier()` reports, and its lane count
  # is the one that tier is built for.
  tier <- vec_active_tier()
  lanes <- vec_active_lanes()

  expect_true(tier %in% c("baseline", "avx2", "avx512"))
  expect_identical(
    lanes,
    switch(tier, baseline = 2L, avx2 = 4L, avx512 = 8L)
  )
})
