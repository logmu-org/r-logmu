# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests record chunking: how a dataset is divided for accumulation, and that the division is what
# fixes the order a total is added up in.
#
# WHY THE ORDER MATTERS. Floating point addition is not associative, so a sum is only as reproducible
# as its order. Accumulating per THREAD would make the total depend on which chunk landed on which
# worker, and so on timing -- the same binary, on the same machine, over the same data, disagreeing
# with itself between runs. Accumulating per CHUNK and folding the partials in chunk order takes the
# schedule out of the answer, and it does so whether or not anything is threaded.
#
# THE RULE IS TESTED DIRECTLY, not just through its consequences, because it decides what an answer
# means. The static assertions in RecordChunk.hpp pin a handful of points at compile time; these
# sweep the properties over many record counts, and pin WHICH equal-division rule it is -- which
# takes a record count that divides unevenly, since an even one looks the same under both.
#
# THE TWO-CHUNK CASE NEEDS OVER TEN THOUSAND INDIVIDUALS, so this file generates them rather than
# making the threshold a test parameter. A parameter would let a test pass against a threshold no
# release ever uses.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

records_per_chunk <- cpp_veil_record_chunks(1L)$records_per_chunk

test_that("the chunk rule divides a dataset into equal parts, in order, covering it exactly", {
  # Swept rather than spot-checked: every record must fall in exactly one chunk, the chunks must be
  # contiguous and ascending, and no two may differ in size by more than one.
  counts <- c(0L, 1L, 2L,
              records_per_chunk - 1L, records_per_chunk, records_per_chunk + 1L,
              2L * records_per_chunk, 2L * records_per_chunk + 1L,
              25000L, 99999L, 1000000L)

  for (n in counts) {
    chunks <- cpp_veil_record_chunks(n)
    sizes <- chunks$end_index - chunks$start_index

    if (n == 0L) {
      expect_length(chunks$start_index, 0L)
      next
    }

    expect_equal(chunks$start_index[1], 0L)
    expect_equal(chunks$end_index[length(chunks$end_index)], n)
    expect_equal(sum(sizes), n)

    # Contiguous: each chunk begins where the last one ended.
    if (length(sizes) > 1L) {
      expect_equal(chunks$start_index[-1], chunks$end_index[-length(chunks$end_index)])
    }

    expect_true(all(sizes > 0L))
    expect_lte(max(sizes) - min(sizes), 1L)
    expect_lte(max(sizes), records_per_chunk)
  }
})

test_that("the chunk count follows the record count, and the threshold is a maximum", {
  expect_length(cpp_veil_record_chunks(records_per_chunk)$start_index, 1L)

  # ONE MORE RECORD THAN THE THRESHOLD GIVES TWO CHUNKS, each about half the threshold. The number is
  # the largest a chunk may be, not the smallest.
  two <- cpp_veil_record_chunks(records_per_chunk + 1L)
  expect_length(two$start_index, 2L)
  expect_true(all((two$end_index - two$start_index) < records_per_chunk))
})

test_that("the division is the proportional form, not the remainder-to-the-front form", {
  # TWO RULES GIVE EQUAL-SIZED PARTS AND THEY ARE NOT THE SAME RULE. The proportional form used here
  # is `floor(i * n / k)` to `floor((i + 1) * n / k)`. The commoner alternative -- what numpy's
  # `array_split` and the MPI scatter helpers do -- divides evenly and hands the leftover records to
  # the FIRST few chunks. Both give parts differing by at most one, so the properties swept above
  # cannot tell them apart, and neither can an even division: fifteen thousand records split 7500 and
  # 7500 either way.
  #
  # IT MATTERS BECAUSE THE PARTITION DECIDES THE ORDER A TOTAL IS ACCUMULATED IN. Same sizes,
  # different membership, different last digits. So the rule is pinned by a case where the division
  # is UNEVEN, which is where the two forms part company.

  # An odd record count over two chunks. The extra record falls at the END here; the other form would
  # give 7501 and 7500.
  skip_if_not(records_per_chunk == 10000L, "the chunk threshold has been retuned")

  odd <- cpp_veil_record_chunks(15001L)
  expect_equal(odd$start_index, c(0L, 7500L))
  expect_equal(odd$end_index, c(7500L, 15001L))
  expect_equal(odd$end_index - odd$start_index, c(7500L, 7501L))

  # Three chunks, so the leftover is more obviously placed: 8333 / 8333 / 8334 here, against
  # 8334 / 8333 / 8333 the other way.
  three <- cpp_veil_record_chunks(25000L)
  expect_equal(three$end_index - three$start_index, c(8333L, 8333L, 8334L))
  expect_equal(three$start_index, c(0L, 8333L, 16666L))
})

test_that("an uneven division always gives the first chunk the smaller size", {
  # The same distinction stated so that it survives the threshold being retuned. Under the
  # proportional form the first chunk takes `floor(n / k)`; under the remainder-to-the-front form it
  # takes the ceiling. So wherever the division is uneven, the first chunk being the smaller one is
  # the discriminating fact, whatever `RecordsPerChunk` happens to be.
  #
  # It is NOT true that the sizes are non-decreasing in general -- ten records in four chunks divide
  # 2, 3, 2, 3 -- so the first chunk is what to assert, not the shape of the whole sequence.
  counts <- c(records_per_chunk + 1L, 15001L, 25000L, 99999L, 123457L)
  uneven <- 0L

  for (n in counts) {
    chunks <- cpp_veil_record_chunks(n)
    sizes <- chunks$end_index - chunks$start_index
    k <- length(sizes)

    if (n %% k == 0L) next
    uneven <- uneven + 1L

    expect_equal(sizes[1], n %/% k)
    expect_equal(min(sizes), n %/% k)
    expect_equal(max(sizes), n %/% k + 1L)
  }

  # The sweep is worthless if every count happened to divide evenly.
  expect_gt(uneven, 0L)
})

test_that("a total over two chunks is the sum of the contributions", {
  # THE ONLY TEST HERE THAT EXERCISES THE FOLDING, because everything smaller is a single chunk with
  # nothing to fold. Ten thousand and one individuals, so two partials, and the engine's total has to
  # agree with a plain sum of its own per-individual contributions. That is what catches a partial
  # never folded in, a chunk folded twice, or a boundary out by one -- none of which an analytic
  # answer would notice unless the error happened to be large.
  n <- records_per_chunk + 1L

  # Exposures are whole quarters of differing length, so individuals genuinely differ and a boundary
  # mistake cannot be masked by them all being identical.
  start <- rep(2010L * clicks_per_year, n)
  quarters <- rep(c(4L, 8L, 12L, 16L), length.out = n)

  big <- list(
    birth     = datey::datey(1940 + (seq_len(n) %% 30)),
    amount    = 100 + (seq_len(n) %% 7) * 250,
    E2R_start = datey_clicks(start),
    E2R_end   = datey_clicks(start + quarters * quarter),
    E2R_died  = (seq_len(n) %% 5L) == 0L
  )

  res <- cpp_veil_aev(it_obj(mortality_const(log_mu = -3.2)), it_ast(~ .i$amount), big,
                      quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(res$chunk_count, 2L)
  expect_equal(res$records_included, n)

  expect_equal(res$A, sum(res$contributions$A), tolerance = 1e-12)
  expect_equal(res$E, sum(res$contributions$E), tolerance = 1e-12)
  expect_equal(res$V, sum(res$contributions$V), tolerance = 1e-12)

  # And against the analytic answer, since a constant mortality makes each integrand constant.
  exposure_years <- quarters / 4
  expect_equal(res$E, sum(exp(-3.2) * big$amount * exposure_years), tolerance = 1e-9)
  expect_equal(res$A, sum(big$amount * as.double(big$E2R_died)), tolerance = 1e-9)
})

test_that("the chunk boundary does not decide which individual contributes what", {
  # The same individuals, with one extra appended so that the dataset falls either side of the
  # threshold and so divides into one chunk or two. Every original individual's own contribution must
  # be identical, because chunking decides the ORDER of the sum and nothing about the arithmetic.
  n <- records_per_chunk
  start <- rep(2010L * clicks_per_year, n + 1L)
  quarters <- rep(c(4L, 8L, 12L, 16L), length.out = n + 1L)

  make <- function(count) {
    list(
      birth     = datey::datey(1940 + (seq_len(count) %% 30)),
      amount    = 100 + (seq_len(count) %% 7) * 250,
      E2R_start = datey_clicks(start[seq_len(count)]),
      E2R_end   = datey_clicks(start[seq_len(count)] + quarters[seq_len(count)] * quarter),
      E2R_died  = (seq_len(count) %% 5L) == 0L
    )
  }

  one <- cpp_veil_aev(it_obj(mortality_const(log_mu = -3.2)), it_ast(~ .i$amount), make(n),
                      quarter_scale, NULL, no_overdispersion, 1L)
  two <- cpp_veil_aev(it_obj(mortality_const(log_mu = -3.2)), it_ast(~ .i$amount), make(n + 1L),
                      quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(one$chunk_count, 1L)
  expect_equal(two$chunk_count, 2L)

  expect_equal(two$contributions$A[seq_len(n)], one$contributions$A)
  expect_equal(two$contributions$E[seq_len(n)], one$contributions$E)
  expect_equal(two$contributions$V[seq_len(n)], one$contributions$V)
})
