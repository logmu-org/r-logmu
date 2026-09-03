# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# WHAT THE BINARY CONTAINS, NOT WHAT THE BUILD FILES ASKED FOR.
#
# `test-simd-guards.R` checks our INTENT: that the features a kernel is built
# with are the features its runtime guard tests for. This file checks the
# OUTCOME, by disassembling the compiled library and looking for one class of
# instruction.
#
# THE BUG THIS EXISTS FOR. gcc on mingw-w64 emitted `vmovapd %ymm3,(%rcx)` in
# `tier::avx2::exp_V_V`, where `rcx` held a STACK address. An aligned 256-bit
# move faults unless the address is 32-byte aligned, and the Windows x64 ABI
# guarantees only 16 -- gcc will not even promise more, refusing
# `-mpreferred-stack-boundary=5` as "not between 3 and 4". So the store
# succeeded or killed the process depending on where the call chain happened to
# leave the stack: about half the time. It killed R with no error and no
# traceback, in a different place each run, and it took two days.
#
# `-Wa,-muse-unaligned-vector-move` in `src/Makevars.win` rewrites every
# aligned vector move as the unaligned encoding, which on AVX hardware costs
# nothing. That flag is a MITIGATION FOR A BUG WE DO NOT OWN -- GCC PR 49001
# (May 2011) and PR 54412, both still open. Rtools could change binutils,
# somebody could add a kernel without the flag, or gcc could start spilling
# somewhere new. In each case the build succeeds and the failure comes back as
# an intermittent crash on a user's machine.
#
# MATCH ANY MEMORY OPERAND, NOT JUST `rsp` AND `rbp`. The first version of this
# check looked only for stack-pointer displacements, on the reasoning that gcc
# rounds a register with `and $-32` before using it as an aligned address. IT
# DOES NOT ALWAYS, and the instruction that was actually crashing addressed the
# stack through `rcx`. That check ran clean against the very binary that was
# dying. With the flag on, the correct count of aligned wide memory accesses is
# ZERO, so that is what this asserts.
#
# WINDOWS ONLY, DELIBERATELY. An aligned stack access is not a defect on Linux
# or macOS: the SysV ABI lets gcc realign the stack, and it does.

# The library that is actually LOADED, which is the one whose instructions run.
# `system.file()` would find whatever happens to be installed, and that can be
# stale -- the copy on this machine was two weeks old when this was written.
loaded_library_path <- function() {
  dll <- getLoadedDLLs()[["logmu"]]
  if (is.null(dll)) return(NA_character_)
  dll[["path"]]
}

# Aligned 256/512-bit moves with any memory operand. `vmovdqa`, `vmovapd` and
# `vmovaps` all require natural alignment; the unaligned encodings (`vmovdqu`,
# `vmovupd`, `vmovups`) do not, and are what the build must produce.
aligned_wide_memory_moves <- function(disassembly) {
  wide <- grep("vmovdqa|vmovapd|vmovaps", disassembly, value = TRUE)
  wide <- grep("ymm|zmm", wide, value = TRUE)
  grep("(", wide, value = TRUE, fixed = TRUE)
}

test_that("the compiled library holds no aligned 256/512-bit memory access", {

  skip_on_os(c("mac", "linux", "solaris"))
  skip_if_not(
    grepl("x86|x64", R.version$arch, ignore.case = TRUE),
    "the fault is x86-64 specific"
  )

  objdump <- Sys.which("objdump")
  skip_if(objdump == "", "objdump is not on the path")

  dll <- loaded_library_path()
  skip_if(is.na(dll) || !file.exists(dll), "the compiled library is not loaded")

  disassembly <- suppressWarnings(
    system2(objdump, c("-d", shQuote(dll)), stdout = TRUE, stderr = FALSE)
  )
  skip_if(length(disassembly) < 100L, "objdump produced no disassembly")

  offenders <- aligned_wide_memory_moves(disassembly)

  expect_identical(
    length(offenders), 0L,
    info = paste0(
      "The build emitted ", length(offenders), " aligned wide memory access(es). ",
      "Any one of them faults on Windows if its address is not 32-byte aligned, ",
      "which the ABI does not guarantee. Check that `SIMD_SAFE_STACK` in ",
      "src/Makevars.win still reaches every SIMD kernel and that this binutils ",
      "still supports it. First offenders:\n",
      paste(utils::head(trimws(offenders), 5L), collapse = "\n")
    )
  )
})

test_that("the scan recognises an offending instruction when there is one", {

  # DISABLE-AND-RECHECK, BUILT IN. The test above passes by finding nothing,
  # which is also what a broken matcher does -- and a broken matcher is exactly
  # what happened here, so this is not hypothetical.
  offending <- c(
    # THE ONE THAT WAS ACTUALLY CRASHING, addressed through a register. The
    # first version of this file did not match it.
    "  0x7ffa1c161181 <exp_V_V+1105>:\tvmovapd %ymm3,(%rcx)",
    "     683:\tc5 fd 6f 44 24 40    \tvmovdqa 0x40(%rsp),%ymm0",
    "     ffd:\tc5 7d 6f 8c 24 20 01 \tvmovdqa 0x120(%rsp),%ymm9",
    "     6c4:\tc5 fd 28 0e          \tvmovapd (%rsi),%ymm1"
  )
  safe <- c(
    "     b5f:\tc5 fd 11 46 20       \tvmovupd %ymm0,0x20(%rsp)",
    "     4a1:\tc5 f9 6f 44 24 20    \tvmovdqa 0x20(%rsp),%xmm0",
    "     7c3:\tc5 fd 28 c1          \tvmovapd %ymm1,%ymm0"
  )

  expect_length(aligned_wide_memory_moves(offending), 4L)
  expect_length(aligned_wide_memory_moves(safe), 0L)
  expect_length(aligned_wide_memory_moves(c(safe, offending)), 4L)
})
