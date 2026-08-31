# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

###### Typing ######

include_class <- c("include", "logmu_function")

#' Includes: a single time interval per individual
#'
#' @description
#' An `include` maps each individual to a single clopen time interval
#' \eqn{[\nu, \tau)} of exposure to count. To guarantee the interval is always
#' a single, well-defined period, includes are built only through these
#' constructors -- never from raw `.t` expressions.
#'
#' The whole family is one construct with three shorthands:
#'
#' * [band()] and [bands()] band any permitted variable.
#' * [age()] and [ages()] are shorthand for banding `.x`, i.e. `.t - .b`.
#' * [duration()] and [durations()] are shorthand for banding `.t` minus an
#'   origin, e.g. `.t - .i$entry`. [age()] is the special case whose origin is
#'   birth.
#' * [period()] and [periods()] are shorthand for banding `.t`.
#'
#' Every singular constructor returns one `include`; every plural constructor
#' returns an [includes()], even when it holds a single band. So no
#' constructor's return type depends on which arguments were supplied.
#'
#' Banding cuts a variable at edges, so it suits a continuous one. A discrete
#' field is divided by [category()] and [categories()] instead, which follow
#' the same singular and plural rule.
#'
#' @section What may be banded:
#' The banded variable must be **increasing in `.t` at unit slope**, which
#' permits exactly three shapes:
#'
#' * a time-invariant variable, e.g. `.i$pension` -- this *gates* exposure,
#'   since the individual is in one band throughout;
#' * `.t` itself -- this *clips* exposure to calendar bounds;
#' * `.t` minus a time-invariant `datey` expression, e.g. `.t - .i$retirement`,
#'   `.t - min(.i$entry, .i$retirement)` (and `.x`, which is `.t - .b`) -- this
#'   clips exposure to bounds measured from that origin.
#'
#' A decreasing expression such as `2020 - .t` is refused: it would flip the
#' band to `(a, b]`, so adjacent bands would either double-count or drop an
#' instant.
#'
#' @section Band edges:
#' The edges are `c(from, thresholds, to)`, and a `NULL` bound means
#' unbounded. `by` generates the interior thresholds from `from` and `to`, and
#' must divide `to - from` exactly. `by` and `thresholds` cannot both be given.
#'
#' So `ages(65, 95, by = 5)` gives six bands and excludes the outside, because
#' both bounds were supplied, while
#' `bands(.i$pension, thresholds = c(5000, 10000, 20000))` gives four bands
#' open at both ends, because neither was.
#'
#' @section Missing values:
#' A `NaN` banded variable compares `FALSE` against every threshold, following
#' IEEE 754, so such a record falls in no band and is simply absent. Note that
#' base R differs here: `NaN < 5` is `NA` in R, which needs a third logical
#' state **logmu** does not have.
#'
#' @section Combining:
#' Includes (and indicators) combine with `&`, which intersects: the result
#' selects the time within *both* operands. Internally an include is a
#' conjunction of terms -- interval bounds plus indicator gates -- resolved
#' together by [period_included()].
#'
#' `period_included()` resolves the interval for a single individual's facts
#' `.i`. It is the plain-R reference path for testing and understanding, not
#' the performance path. An offset that is `NA` yields the empty interval.
#'
#' @param variable A pronoun expression to band. See *What may be banded*.
#' @param since The origin a duration is measured from: a time-invariant
#'   `datey` expression, e.g. `.i$entry` or `min(.i$entry, .i$retirement)`.
#' @param from,to Outer band edges, or `NULL` for unbounded. Dates (`datey` or
#'   year numbers) when banding `.t`; durations (`durationy` or year counts)
#'   when banding `.t` minus an origin; ordinary numbers otherwise.
#' @param by The width of each band. Requires both `from` and `to`, and must
#'   divide `to - from` exactly.
#' @param thresholds Interior band edges, given explicitly instead of `by`.
#' @param expr A time-invariant pronoun expression.
#' @param x An `include`.
#' @param .i A named list of scalar facts.
#' @returns
#' `band()`, `age()`, `duration()` and `period()` return an `include`;
#' `include()` returns an `indicator`, which is also an `include`.
#'
#' `bands()`, `ages()`, `durations()` and `periods()` return an `includes`.
#'
#' `is_include()` returns a scalar `logical`.
#'
#' `period_included()` returns a `datey_interval`.
#' @examples
#' period_included(period(2010, 2020), .i = list())
#' period_included(age(65, 95), .i = list(birth = datey::datey(1950)))
#' period_included(age(65, 95) & period(2010, 2040),
#'                 .i = list(birth = datey::datey(1950)))
#'
#' # Six age bands, outside excluded.
#' ages(65, 95, by = 5)
#'
#' # Four amount bands, open at both ends.
#' bands(.i$pension, thresholds = c(5000, 10000, 20000))
#'
#' # Duration since an origin other than birth.
#' durations(.i$retirement, 0, 10, by = 5)
#'
#' # The origin may be computed, not only a bare field.
#' band(.t - min(.i$entry, .i$retirement), 0, 5)
#' @name include
NULL

# An include is a conjunction of terms. Each term is one of:
#   absolute : calendar bounds        (from, to are datey)
#   offset   : bounds offset by an expression (from, to are durationy; the
#              offset is a time-invariant datey AST, not necessarily a field)
#   gate     : an indicator that must hold (ast is a pronoun AST)
# `period_included()` resolves them together: the interval is the intersection
# of the bounds, gated by the indicators. Either bound of an interval term may
# be NULL, meaning unbounded on that side.
#
# `group_name` and `name` are the two label levels: a 15-element breakdown of
# ages, periods and amounts has to tell a chart which six are the ages, so a
# single flat name is not enough. `name` becomes the list name once the include
# is collected into an `includes`. Both are NULL for a hand-built include.
new_include <- function(terms, group_name = NULL, name = NULL) {
  structure(list(terms = terms, group_name = group_name, name = name),
            class = include_class)
}

absolute_term  <- function(from, to)         list(kind = "absolute", from = from, to = to)
offset_term    <- function(offset, from, to) list(kind = "offset", offset = offset, from = from, to = to)
gate_term      <- function(ast)              list(kind = "gate", ast = ast)

# ---- what may be banded ----------------------------------------------------

# The AST for `.x`, which is what `age()` and `ages()` band. `it_build_ast`
# desugars a written `.x` to exactly this, so the two routes agree by
# construction rather than by coincidence.
band_age_ast <- function() it_call("-", list(it_time(), it_field("birth")))

# Classify a banded variable into the term kind it produces. The permitted
# shapes are the three that are increasing in `.t` at unit slope; anything else
# would not resolve to a single clopen interval.
band_variable <- function(ast) {
  if (!it_uses_t(ast)) return(list(kind = "static", ast = ast))

  if (identical(ast$kind, "time")) return(list(kind = "absolute"))

  # `.t` minus ANY time-invariant expression, not just a bare field. The origin
  # a duration is measured from is as much a computed thing as anything else --
  # `.t - min(.i$entry, .i$retirement)` is a perfectly ordinary duration -- and
  # a user who has written `.t - .i$entry` is entitled to be surprised that the
  # other is refused. Its `datey`-ness is checked where it is resolved, since
  # that needs the data.
  if (identical(ast$kind, "call") &&
      identical(ast$fn, "-") &&
      length(ast$args) == 2L &&
      identical(ast$args[[1L]]$kind, "time") &&
      !it_uses_t(ast$args[[2L]])) {
    return(list(kind = "offset", offset = ast$args[[2L]]))
  }

  stop("A banded variable must be time-invariant, `.t`, or `.t` minus a ",
       "time-invariant expression (such as `.x`, `.t - .i$retirement` or ",
       "`.t - min(.i$entry, .i$retirement)`).",
       call. = FALSE)
}

# ---- edges -----------------------------------------------------------------

# Band edges are computed in numeric years and converted back to the term's own
# type as each term is built. `as.numeric()` on a datey gives the year and on a
# durationy the count of years, so a user may supply either form.
band_bound_years <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (datey::is_datey(x) || datey::is_durationy(x)) x <- as.numeric(x)
  if (!is_single_pure_finite_numeric(x)) {
    stop("`", what, "` must be a single finite number, `datey` or `durationy`.",
         call. = FALSE)
  }
  as.numeric(x)
}

# The full edge list, NULL at either end meaning unbounded. One rule covers
# every constructor: `c(from, thresholds, to)`, with `by` generating the
# interior thresholds.
band_edges <- function(from, to, by, thresholds) {
  if (!is.null(by) && !is.null(thresholds)) {
    stop("Give either `by` or `thresholds`, not both.", call. = FALSE)
  }

  if (!is.null(by)) {
    if (is.null(from) || is.null(to)) {
      stop("`by` needs both `from` and `to`.", call. = FALSE)
    }
    by <- band_bound_years(by, "by")
    if (by <= 0) stop("`by` must be positive.", call. = FALSE)

    steps <- (to - from) / by
    # Refused rather than rounded: a silent short final band is a reporting trap
    # once it reaches a chart axis.
    if (!isTRUE(all.equal(steps, round(steps))) || round(steps) < 1L) {
      stop("`by` must divide `to - from` exactly.", call. = FALSE)
    }
    # length.out rather than a `by` step, so the last edge is exactly `to`.
    return(as.list(seq(from, to, length.out = round(steps) + 1L)))
  }

  if (is.null(thresholds)) return(list(from, to))

  thresholds <- vapply(thresholds, band_bound_years, numeric(1L), what = "thresholds")
  c(list(from), as.list(thresholds), list(to))
}

band_check_increasing <- function(edges) {
  known <- unlist(edges[!vapply(edges, is.null, logical(1L))])
  if (length(known) > 1L && any(diff(known) <= 0)) {
    stop("Band edges must be strictly increasing.", call. = FALSE)
  }
  if (length(known) == 0L) {
    stop("A band needs at least one edge.", call. = FALSE)
  }
}

# ---- terms and labels ------------------------------------------------------

band_term <- function(variable, lo, hi) {
  switch(variable$kind,
    absolute = absolute_term(
      if (is.null(lo)) NULL else get_single_valid_datey(lo),
      if (is.null(hi)) NULL else get_single_valid_datey(hi)),
    offset = offset_term(
      variable$offset,
      if (is.null(lo)) NULL else get_single_valid_durationy(lo),
      if (is.null(hi)) NULL else get_single_valid_durationy(hi)),
    static = gate_term(band_gate_ast(variable$ast, lo, hi))
  )
}

# A banded time-invariant variable gates exposure rather than clipping it, so
# it lowers to the clopen comparison written out. Only the bounds that exist
# appear, so an unbounded side costs nothing at all.
band_gate_ast <- function(ast, lo, hi) {
  lower <- if (is.null(lo)) NULL else it_call(">=", list(ast, it_lit(lo)))
  upper <- if (is.null(hi)) NULL else it_call("<", list(ast, it_lit(hi)))
  if (is.null(lower)) return(upper)
  if (is.null(upper)) return(lower)
  it_call("&", list(lower, upper))
}

band_number <- function(x) format(x, trim = TRUE, scientific = FALSE)

# The group name for a banded expression. `.i$pension` reads better as
# `pension` once it reaches a chart axis, and the pronoun carries no
# information a reader of the chart wants.
band_group_name <- function(q) sub("^\\.i\\$", "", deparse1(q))

band_name <- function(lo, hi) {
  if (is.null(lo) && is.null(hi)) return("all")
  if (is.null(lo)) return(paste0("< ", band_number(hi)))
  if (is.null(hi)) return(paste0(">= ", band_number(lo)))
  paste0(band_number(lo), "-", band_number(hi))
}

# ---- the general constructors ----------------------------------------------

build_band <- function(ast, group_name, from, to, by, thresholds) {
  variable <- band_variable(ast)
  from <- band_bound_years(from, "from")
  to <- band_bound_years(to, "to")

  edges <- band_edges(from, to, by, thresholds)
  band_check_increasing(edges)

  lapply(seq_len(length(edges) - 1L), function(i) {
    lo <- edges[[i]]
    hi <- edges[[i + 1L]]
    new_include(list(band_term(variable, lo, hi)), group_name, band_name(lo, hi))
  })
}

# One band, and exactly one: a singular constructor with thresholds would be a
# contradiction, so it takes only the outer edges.
build_one_band <- function(ast, group_name, from, to) {
  build_band(ast, group_name, from, to, NULL, NULL)[[1L]]
}

#' @rdname include
#' @export
band <- function(variable, from = NULL, to = NULL) {
  q <- substitute(variable)
  build_one_band(it_capture(q, parent.frame()), band_group_name(q), from, to)
}

#' @rdname include
#' @export
bands <- function(variable, from = NULL, to = NULL, by = NULL, thresholds = NULL) {
  q <- substitute(variable)
  new_includes(build_band(it_capture(q, parent.frame()), band_group_name(q),
                          from, to, by, thresholds))
}

#' @rdname include
#' @export
age <- function(from = NULL, to = NULL) build_one_band(band_age_ast(), "age", from, to)

#' @rdname include
#' @export
ages <- function(from = NULL, to = NULL, by = NULL, thresholds = NULL) {
  new_includes(build_band(band_age_ast(), "age", from, to, by, thresholds))
}

# `duration()` and `durations()` band the time since an origin, so they build
# the same `.t - origin` shape `band()` already takes. `age()` is the special
# case whose origin is birth.
#
# The origin is checked here rather than left to `band_variable`, which would
# answer a wrong origin with a message about banded variables when what the
# user got wrong is one argument.
duration_ast <- function(ast) {
  if (it_uses_t(ast)) {
    stop("The origin of a duration must be time-invariant, so it may not use `.t`.",
         call. = FALSE)
  }
  it_call("-", list(it_time(), ast))
}

# The origin is part of the group name because two duration sets measured from
# different origins would otherwise be indistinguishable in one breakdown, and
# a chart axis is where that would surface.
duration_group_name <- function(q) paste0("duration since ", band_group_name(q))

#' @rdname include
#' @export
duration <- function(since, from = NULL, to = NULL) {
  q <- substitute(since)
  build_one_band(duration_ast(it_capture(q, parent.frame())),
                 duration_group_name(q), from, to)
}

#' @rdname include
#' @export
durations <- function(since, from = NULL, to = NULL, by = NULL, thresholds = NULL) {
  q <- substitute(since)
  new_includes(build_band(duration_ast(it_capture(q, parent.frame())),
                          duration_group_name(q), from, to, by, thresholds))
}

#' @rdname include
#' @export
period <- function(from = NULL, to = NULL) build_one_band(it_time(), "period", from, to)

#' @rdname include
#' @export
periods <- function(from = NULL, to = NULL, by = NULL, thresholds = NULL) {
  new_includes(build_band(it_time(), "period", from, to, by, thresholds))
}

# `include(expr)` and `indicator(expr)` build the SAME object and differ only in
# what they signal at the call site. Returning a plain include would make
# `male & retired` degrade to gate terms and lose the {0,1} guarantee that V = E
# rests on.
#' @rdname include
#' @export
include <- function(expr) {
  ast <- it_capture(substitute(expr), parent.frame())
  if (it_uses_t(ast)) {
    stop("An `include` expression must not depend on time `.t`. Use `period()` ",
         "or `band()` for a time interval.", call. = FALSE)
  }
  new_indicator(ast)
}

#' @rdname include
#' @export
is_include <- function(x) inherits(x, "include")

ensure_is_include <- function(x) {
  if (!is_include(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be an `include`.", call. = FALSE)
  }
}

# ---- resolving -------------------------------------------------------------

# Resolve one interval term to its `datey` bounds, or NULL for no exposure. An
# absent bound stays NULL and simply does not narrow the interval.
include_term_bounds <- function(term, .i) {
  if (identical(term$kind, "absolute")) {
    return(list(start = term$from, end = term$to))
  }
  # The offset is an expression, so it is evaluated rather than looked up. A
  # missing field inside it still reports itself by name, from `it_eval`.
  offset <- it_eval(term$offset, .i = .i)
  if (!datey::is_datey(offset) || length(offset) != 1L) {
    stop(sprintf("The offset `%s` must be a single `datey`.", it_deparse(term$offset)),
         call. = FALSE)
  }
  if (is.na(offset)) return(NULL)   # undefined offset -> no exposure
  list(start = if (is.null(term$from)) NULL else offset + term$from,
       end   = if (is.null(term$to))   NULL else offset + term$to)
}

#' @rdname include
#' @export
period_included <- function(x, .i) UseMethod("period_included")

#' @rdname include
#' @export
period_included.default <- function(x, .i) {
  arg_name <- deparse(substitute(x))
  stop("S3 function `period_included` is not implemented for `", arg_name, "`.", call. = FALSE)
}

#' @rdname include
#' @export
period_included.include <- function(x, .i) {
  nu  <- datey::datey(datey::valid_years_start)
  tau <- datey::datey(datey::valid_years_end)
  for (term in x$terms) {
    if (identical(term$kind, "gate")) {
      if (!isTRUE(indicator_as_logical(it_eval(term$ast, .i)))) {
        return(include_none_interval())
      }
    } else {
      bounds <- include_term_bounds(term, .i)
      if (is.null(bounds)) return(include_none_interval())
      if (!is.null(bounds$start) && bounds$start > nu)  nu  <- bounds$start
      if (!is.null(bounds$end)   && bounds$end   < tau) tau <- bounds$end
    }
  }
  if (!(nu < tau)) return(include_none_interval())
  datey::datey_interval(nu, tau)
}

# The two special resolved intervals: all-of-time is the widest representable
# interval; the empty interval contributes no exposure.
include_all_interval <- function() {
  datey::datey_interval(datey::datey(datey::valid_years_start),
                        datey::datey(datey::valid_years_end))
}
include_none_interval <- function() datey::NA_datey_interval_

# ---- intersection (`&`) ----------------------------------------------------

# Coerce an operand of `&` to a subset (include or indicator). Already-built
# objects pass through; a `~` formula becomes an indicator (and so must be
# time-invariant).
coerce_to_subset <- function(x) {
  if (is_include(x)) return(x)   # includes, and indicators (indicator is-a include)
  if (inherits(x, "formula")) {
    ast <- it_build_ast(x[[length(x)]], environment(x))
    if (it_uses_t(ast)) {
      stop("A formula combined with `&` must be time-invariant (an indicator).", call. = FALSE)
    }
    return(new_indicator(ast))
  }
  stop("`&` operands must be includes, indicators or time-invariant pronoun formulas.",
       call. = FALSE)
}

# The conjunction terms an operand contributes.
subset_terms <- function(x) {
  if (is_indicator(x)) return(list(gate_term(x$ast)))
  x$terms
}

# Intersect two subsets. Two pure indicators stay an indicator (logical AND, so
# the {0,1}/V=E nature is kept); anything involving an interval is an include.
#
# The labels are DELIBERATELY DROPPED: an intersection of an age band and a
# period band belongs to no single dimension, and a wrong label is worse on a
# chart than none.
subset_intersect <- function(a, b) {
  a <- coerce_to_subset(a)
  b <- coerce_to_subset(b)
  if (is_indicator(a) && is_indicator(b)) {
    return(new_indicator(it_call("&", list(a$ast, b$ast))))
  }
  new_include(c(subset_terms(a), subset_terms(b)))
}

# ---- printing --------------------------------------------------------------

format_include_bounds <- function(from, to) {
  if (is.null(from) && is.null(to)) return("unbounded")
  if (is.null(from)) return(sprintf("before %s", format(to)))
  if (is.null(to)) return(sprintf("from %s", format(from)))
  sprintf("[%s, %s)", format(from), format(to))
}

format_include_term <- function(term) {
  switch(term$kind,
    absolute = sprintf("period %s", format_include_bounds(term$from, term$to)),
    offset = {
      # `age` is the name for the one origin that has one; everything else is
      # rendered as what the user wrote.
      label <- if (identical(it_deparse(term$offset), ".i$birth")) {
        "age"
      } else {
        sprintf("since %s", it_deparse(term$offset))
      }
      sprintf("%s %s", label, format_include_bounds(term$from, term$to))
    },
    gate = sprintf("where %s", it_deparse(term$ast))
  )
}

#' @export
#' @noRd
print.include <- function(x, ...) {
  if (is.null(x$name)) {
    cat("<include>\n")
  } else {
    cat("<include: ", x$group_name, " ", x$name, ">\n", sep = "")
  }
  for (term in x$terms) cat("  ", format_include_term(term), "\n", sep = "")
  invisible(x)
}
