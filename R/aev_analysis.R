# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Resolve the two calculation choices from whichever of the three places they
# were given: the call itself, a `settings` object, or the default.
#
# THE CALL WINS, so a single analysis can disagree with the settings it was
# otherwise handed. `overdispersion` has NO default at any level -- if neither
# the call nor a `settings` supplies it, that is an error rather than an
# assumption. See `settings()` for why.
resolve_settings <- function(settings, overdispersion, time_scale) {
  if (!is.null(settings)) ensure_is_settings(settings)

  overdispersion <- overdispersion %||% settings$overdispersion
  if (is.null(overdispersion)) {
    stop("`overdispersion` is required. Give it directly, or pass a `settings` ",
         "object that carries it.", call. = FALSE)
  }

  clicks <- if (!is.null(time_scale)) {
    time_scale_clicks(time_scale)
  } else if (!is.null(settings)) {
    settings$time_scale_clicks
  } else {
    time_scale_clicks(default_time_scale)
  }

  list(overdispersion = check_overdispersion(overdispersion),
       time_scale_clicks = clicks)
}

# The columns as the engine wants them: a plain named list, with the data.frame
# machinery (row names and the class) left behind.
#
# CHARACTER COLUMNS BECOME FACTORS HERE, which is the whole of R's side of text
# handling. The engine compares text as integer indices, and `factor()` is R's
# own C implementation of turning a column of strings into exactly that -- so
# the O(rows) work happens once, in the code that already exists for it, and
# what crosses is an integer column plus a handful of level strings. The engine
# merges those level tables into one numbering per crossing; nothing here has to
# know how the levels of two columns line up.
#
# A column that is already a factor is left alone, levels and all.
exp_data_columns <- function(exp_data) {
  columns <- unclass(exp_data)
  attributes(columns) <- list(names = names(columns))

  text <- vapply(columns, is.character, logical(1L))
  columns[text] <- lapply(columns[text], factor)
  columns
}

# The population include for one breakdown element.
#
# `subset_intersect` DROPS THE LABELS, deliberately -- an intersection of an age
# band and a population belongs to no single dimension. That is why the labels
# on the result are read from the breakdown itself rather than from these.
aev_element_include <- function(population, element) {
  if (is.null(population)) return(element)
  subset_intersect(population, element)
}

#' Actual, expected and variance over experience data
#'
#' @description
#' Calculates \eqn{A}, \eqn{E} and \eqn{V} over an `exp_data` under a given
#' mortality, optionally restricted to a sub-population and optionally broken
#' down into groups.
#'
#' `mortality` and `weight` are **pronoun expressions**, so they may be written
#' out in place (`.i$pension`), or reference a `mortality`, `variable` or
#' `indicator` object by name, or combine the two.
#'
#' @section Breakdowns:
#' With no `breakdown` the result is a length-1 `aev`. With one, the result is a
#' single `aev` with one record per element -- not a list -- carrying the
#' breakdown's `names()` and [group_names()] so a chart can tell which records
#' were ages and which were amounts.
#'
#' Every element is intersected with `include`, so `include` says who is in the
#' population and `breakdown` says how to divide them. Elements need not be
#' disjoint and need not cover everybody: examining A/E on intersecting subsets
#' is legitimate, and a record outside every element is simply absent.
#'
#' @section Overdispersion:
#' `overdispersion` is required, either directly or through `settings`. It
#' scales `V`, and so scales every confidence interval and residual read from
#' the result. There is no default: see [settings()].
#'
#' @section Similarity and distance:
#' A second weighting factor \eqn{s} says how much a record should count at all,
#' as opposed to how much of it counts. It differs from `weight` in where it
#' appears: the weight is squared in \eqn{V} and the similarity is not, so
#' \eqn{A} and \eqn{E} take \eqn{sw} while \eqn{V} takes \eqn{sw^2}. The effect
#' is that halving a record's similarity doubles the variance it implies, while
#' leaving the number of parameters a model may support unchanged.
#'
#' It may be written either way round, and they are one quantity:
#'
#' * `val_similarity` is \eqn{s} itself, intended to lie in \eqn{[0, 1]}, where
#'   1 counts a record fully and 0 not at all.
#' * `val_distance` is \eqn{d = -\log s}, so 0 counts a record fully, `log(2)`
#'   counts it at half, and larger counts it less.
#'
#' Give one or the other, never both. `val_distance` is usually the easier of
#' the two: a decay kernel is `val_distance = (2025 - .t) / 10` and a Gaussian
#' is `val_distance = (x / h)^2`, where the similarity form would need the
#' exponential written out each time. Distances also add, so several reasons to
#' discount a record combine by addition.
#'
#' A similarity of zero counts a record not at all, but [include()] is the way to
#' leave records out: it clips exposure rather than multiplying through it.
#'
#' @section Why the names begin with `val_`:
#' A similarity or a distance is ordinarily a relation between two things. These
#' arguments are functions of one record and one time, so the far end of the
#' relation is left implicit: it is the valuation the analysis is aimed at, taken
#' as a whole and at its as-at date. Leaving it implicit is what allows them to
#' be written as ordinary variables of the experience record.
#'
#' You write a `val_` argument yourself, and nothing derives it from valuation
#' data. The prefix keeps `similarity` and `distance` free for a later, more
#' general form that would take two records and two times.
#'
#' @section The bound on a similarity:
#' A similarity is a proportion and must lie in \eqn{[0, 1]}; a distance,
#' being \eqn{-\log s}, must not be negative.
#'
#' Neither is enforced value by value. Both are expressions of the record and of
#' time, so their values are not known until the walk, and testing each one
#' would put a comparison in the innermost loop to police something already
#' stated here.
#'
#' What **logmu** does instead is analytic, and costs nothing. It works out the
#' range each expression can take -- from the numbers written in it and the
#' range of each column in the data -- and refuses the calculation before
#' reading a record if that range lies **wholly** outside the bound. So
#' `val_similarity = 2` is refused, and so is a similarity given as a column
#' whose values run from 5 to 200.
#'
#' It refuses only what is certainly wrong. A range that merely permits a
#' violation is accepted, because the arithmetic that derives it is
#' conservative: `val_distance = (2025 - .t) / 10` would otherwise be rejected,
#' since the range of `.t` alone allows a negative distance that no exposure in
#' the data reaches. A check that fires on correct code would be worse than one that
#' occasionally stays quiet.
#'
#' @param exp_data The experience data.
#' @param mortality A pronoun expression for \eqn{\log\mu}, or a `mortality`.
#' @param include An `include` naming the population. `NULL` means everybody.
#' @param weight A pronoun expression for the weight \eqn{w}, or a `variable`.
#'   `NULL` means a weight of 1, i.e. a count of lives.
#' @param val_similarity,val_distance Two spellings of the second weighting
#'   factor \eqn{s}, related by \eqn{d = -\log s}. Give one or the other. A
#'   similarity belongs in \eqn{[0, 1]} and a distance must not be negative; see
#'   *Similarity and distance*, *Why the names begin with `val_`* and *The bound
#'   on a similarity*.
#' @param breakdown An `include` or [includes()] dividing the population into
#'   groups. `NULL` gives a single ungrouped result.
#' @param settings A [settings()] object supplying `overdispersion` and
#'   `time_scale`.
#' @param overdispersion,time_scale Given directly, these override `settings`.
#' @param threads Worker threads to use. `0` asks for as many as the machine
#'   reports. Cannot change any answer.
#' @returns An `aev` with one record per breakdown element, or one record if
#'   there is no breakdown.
#' @examples
#' data <- exp_data(
#'   list(
#'     birth     = datey::datey(c(1945, 1950, 1955)),
#'     pension   = c(5000, 12000, 30000),
#'     E2R_start = datey::datey(c(2015, 2015, 2015)),
#'     E2R_end   = datey::datey(c(2020, 2020, 2018)),
#'     E2R_died  = c(FALSE, FALSE, TRUE)
#'   ),
#'   exp_start = datey::datey(2015),
#'   exp_end   = datey::datey(2020)
#' )
#'
#' basis <- settings(overdispersion = 2)
#'
#' aev(data, mortality = mortality_const(log_mu = -4), settings = basis)
#'
#' # Weighted by pension, broken down by amount.
#' aev(data,
#'     mortality = mortality_const(log_mu = -4),
#'     weight    = .i$pension,
#'     breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
#'     settings  = basis)
#' @export
aev <- function(exp_data,
                mortality,
                include = NULL,
                weight = NULL,
                val_similarity = NULL,
                val_distance = NULL,
                breakdown = NULL,
                settings = NULL,
                overdispersion = NULL,
                time_scale = NULL,
                threads = cpp_veil_default_threads()) {

  ensure_is_exp_data(exp_data)

  if (missing(mortality)) {
    stop("`mortality` is required.", call. = FALSE)
  }

  resolved <- resolve_settings(settings, overdispersion, time_scale)

  # CAPTURED FROM THE CALLER'S FRAME, both of them, because these are pronoun
  # expressions rather than values. `missing()` rather than testing the value:
  # `substitute(weight)` on an unsupplied argument yields its default, so a
  # value test could not tell `weight = NULL` from no weight at all -- and both
  # mean the same thing here, which is why they are treated together.
  caller <- parent.frame()
  mortality_ast <- it_capture(substitute(mortality), caller)

  plan <- aev_plan(mortality_ast,
                   aev_optional_ast(substitute(weight), missing(weight), caller),
                   aev_optional_ast(substitute(val_similarity), missing(val_similarity), caller),
                   aev_optional_ast(substitute(val_distance), missing(val_distance), caller),
                   include, breakdown, resolved$overdispersion)

  run <- cpp_veil_run(
    plan$specs,
    exp_data_columns(exp_data),
    resolved$time_scale_clicks,
    FALSE,
    as.integer(threads)
  )

  plan$finalise(run$results)
}

# An optional pronoun argument, captured in the caller's frame. `missing()`
# rather than a value test, because `substitute()` on an unsupplied argument
# yields its default, so a value test could not tell `weight = NULL` from no
# weight at all -- and both mean the same thing here.
aev_optional_ast <- function(expr, absent, caller) {
  if (absent || is.null(expr)) NULL else it_capture(expr, caller)
}

# THE SEAM `batch()` RUNS THROUGH. An analysis is two halves: the specifications
# it wants run, and how to read the engine's answers back as its own type.
# Nothing between them touches the data, so several analyses can have their
# specifications compiled and scheduled together and still finish as themselves.
#
# `aev()` above is the degenerate batch of one: build the plan, run it, finalise
# it. It is deliberately written that way rather than kept as a separate path,
# so a batched aev and a lone one cannot drift apart.
#
# The plan takes ASTs rather than expressions because the capture must happen in
# the USER'S frame, and by here we are one or two frames away from it.
aev_plan <- function(mortality_ast, weight_ast, val_similarity_ast, val_distance_ast,
                     include, breakdown, overdispersion) {
  # `val_similarity` and `val_distance` are two spellings of ONE quantity,
  # related by d = -log s, and giving both is refused AT THE ENGINE BOUNDARY
  # rather than here. A check in this function was written first and turned out
  # to be invisible: breaking it failed no test, because the engine refuses the
  # same pair with the same words. The boundary is where the rule has to hold, since
  # it guards every front end; `batch()` adds its own only because it can say
  # WHICH analysis was wrong, which the engine cannot.
  #
  # Neither spelling is converted into the other on the way there. See
  # `SimilarityForm` in `AevRecipe.hpp` for why that matters.
  if (!is.null(include)) {
    if (is_includes(include)) {
      stop("`include` takes a single `include`. Use `breakdown` for several.",
           call. = FALSE)
    }
    ensure_is_include(include)
  }

  elements <- if (is.null(breakdown)) NULL else as_includes(breakdown)

  # ONE SPEC PER ELEMENT, which is what the batch entry point exists to run: a
  # breakdown of G groups is G blocks over one dataset, compiled together and
  # scheduled together, with the data read once.
  spec_includes <- if (is.null(elements)) {
    list(include)
  } else {
    lapply(unclass(elements), function(element) aev_element_include(include, element))
  }

  specs <- lapply(spec_includes, function(spec_include) {
    list(
      mortality      = mortality_ast,
      weight         = weight_ast,
      val_similarity = val_similarity_ast,
      val_distance   = val_distance_ast,
      include        = spec_include,
      overdispersion = overdispersion
    )
  })

  list(
    specs = specs,
    finalise = function(results) {
      totals <- function(what) {
        vapply(results, function(result) result[[what]], numeric(1L))
      }

      # UNCHECKED, DELIBERATELY. The invariants `create_aev()` enforces are for
      # a triple somebody typed; a computed one is whatever the arithmetic came
      # to, and it may legitimately be degenerate. A mortality underflowing to
      # exactly zero -- `exp()` gets there at a log mortality around -746 -- puts
      # E and V at zero while A still counts the deaths that happened, and that
      # is a true statement about an impossible model.
      #
      # THE COST OF CHECKING IS THE WHOLE RESULT, NOT THE BAD ELEMENT. This is
      # vectorised over the breakdown, so one underflowing band would take every
      # other band's answer with it, and `batch()` finalises in a plain loop, so
      # it would take the other elements of the batch too.
      result <- create_aev_unchecked(A = totals("A"), E = totals("E"), V = totals("V"))

      # The labels come from the BREAKDOWN, not from the includes that were run
      # -- intersecting with the population drops them, and rightly so.
      if (is.null(elements)) {
        return(result)
      }

      # A hand-built `includes` may carry neither level of naming, in which case
      # `group_names()` is all NA. Setting that would be worse than leaving it
      # off: `group_names()` returning NAs reads as "these have groups and I
      # lost them", where NULL says plainly that there were none.
      groups <- group_names(elements)
      if (all(is.na(groups))) groups <- NULL

      set_aev_labels(result, names(elements), groups)
    }
  )
}
