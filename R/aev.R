# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

create_aev_unchecked <- function(A, E, V) {
  structure(list(A = A, E = E, V = V), class = "aev" )
}

#' Is an object an `aev`?
#'
#' Tests whether `x` is an `aev` -- the three-field record type [aev()] returns
#' and [create_aev()] builds. See [aev_properties].
#'
#' @param x An object to test.
#' @returns A single `TRUE` or `FALSE`. Never raises an error, whatever `x` is.
#' @examples
#' is_aev(create_aev(A = 1, E = 2, V = 3))
#' is_aev(1:3)
#' @export
# `inherits()` RATHER THAN `class(x) == "aev"`, WHICH THREW. A class vector of
# more than one element made that comparison length-2, and R refuses a condition
# of length other than one, so `is_aev(matrix(1))` and `is_aev(as.POSIXlt(t))`
# raised an error where a predicate has to answer FALSE. Found when this was
# exported, having been harmless while every caller held an aev already.
#
# `length()` is a METHOD on an aev and reports the number of records, so the
# field count has to be asked for through `unclass()`. Reading `length(x)`
# here would have silently started asking a different question.
is_aev <- function(x) inherits(x, "aev") && is.list(x) && length(unclass(x)) == 3L

validate_aev <- function(aev) {
  if (!is_aev(aev)) {
    stop("`", deparse(substitute(aev)), "` is not an `aev`.", call. = FALSE)
  }

  u <- unclass(aev)

  A <- u[["A"]]
  E <- u[["E"]]
  V <- u[["V"]]

  ensure_is_pure_double(A)
  ensure_is_pure_double(E)
  ensure_is_pure_double(V)

  cpp_validate_aev(A, E, V)
}

#' Properties of an aev
#'
#' @description
#'
#' An `aev` comprises the following properties:
#'
#' |Property|Formula|
#' |:-------|:------|
#' | `A`    | Actual deaths weighted by \eqn{w}
#' | `E`    | 'Expected' deaths weighted by \eqn{w}
#' | `V`    | 'Expected' deaths weighted by \eqn{w^2}
#'
#' where \eqn{w} is an arbitrary (non-negative) weight used to calculate them.
#'
#' [create_aev()] refuses a triple that breaks any of the following, so an `aev`
#' you build by hand is guaranteed to satisfy them:
#' (a) none of `A`, `E` or `V` are negative,
#' (b) `E` and `V` are either both zero or both non-zero, and
#' (c) `A` cannot be non-zero if `E` and `V` are zero.
#'
#' A **computed** `aev` is not checked against them, and (b) and (c) can fail
#' when a mortality underflows: `exp()` reaches exactly zero at a log mortality
#' of about -746, so `E` and `V` are then zero while `A` still counts the deaths
#' that happened. That is a legitimate statement about an impossible model, so it
#' is returned rather than raised -- an analysis running for hours, or a [batch()]
#' whose other elements are sound, must not be lost to one underflowing cell. The
#' calculated properties below take their ordinary IEEE values there, so
#' `A_over_E` reads `Inf` and `deviance_residual` reads `NaN`.
#'
#' These columns cannot be modified individually.
#'
#' In typical use, `aev`s are created by **logmu** analytic functions.
#' A `create_aev` function is provided for testing and illustration.
#'
#' `A`, `E` or `V` have the following statistical properties
#' *if the mortality used to calculate them is correct*:
#'
#' \deqn{\mathbb{E}\,(A - E) = 0}
#' \deqn{\mathrm{Var}(A-E) = \mathbb{E}\,V}
#'
#' An `aev` provides the following *calculated* properties:
#'
#' |Property|Formula|
#' |:-------|:------|
#' | `A_minus_E`           | \eqn{A - E}
#' | `A_minus_E_stddev`    | \eqn{\sqrt{V}}
#' | `A_over_E`            | \eqn{A / E}
#' | `log_A_over_E`        | \eqn{\log(A / E)}
#' | `log_A_over_E_stddev` | \eqn{\sqrt{V} / E}
#' | `log_A_over_E_95pc`   | \eqn{k\sqrt{V} / E} where \eqn{k \approx 1.96} -- used by **logmu** to display 95% confidence intervals
#' | `Pearson_residual`    | \eqn{(A - E) / \sqrt{V}}
#' | `deviance_residual`   | \eqn{\mathrm{sign}(A-E)\sqrt{2E/V\cdot\Big[A\cdot\log(A/E)-(A-E)\Big]}}
#'
#' `aev`s can be added, which simply means adding their A, E and V components.
#'
#' The addition of `aev`s is statistically legitimate
#' *provided they relate to independent experience data*,
#' i.e. data that does not intersect by time *and* individual.
#' This means that it *is* legitimate to add `aev`s
#' for the same individual provided they relate to non-overlapping time periods, or
#' for overlapping time periods provided they relate to different individuals.
#'
#' An `aev` is a **vector of records**: `length()` reports how many, and `[`
#' subsets records, keeping `A`, `E` and `V` together.
#'
#' When an `aev` comes from a breakdown it also carries two levels of labelling,
#' `names()` and [group_names()], exactly as the [includes] it came from did.
#' Both are absent from a hand-built `aev`.
#'
#' Two `aev`s may be added, which adds `A`, `E` and `V` record by record. The
#' labels travel with the sum. Adding record by record asserts that the records
#' correspond, so labels that disagree say they do not, and the addition is
#' refused. A missing set of labels is not a disagreement: adding a hand-built
#' `aev` to a labelled one keeps the labels.
#'
#' Both labels may be replaced. `names<-` takes exactly one label per record;
#' `group_names<-` takes one per record or a single value for all of them.
#'
#' @param x An `aev` object.
#' @param i,name The `aev` property being requested.
#' @param A,E,V Columns used to construct an `aev`.
#' @returns An `aev`.
#' @examples
#' aev <- create_aev(A = c(1100, 0), E = c(1000, 1), V = c(2500, 1))
#' aev
#'
#' # Data properties:
#' aev$A
#' aev$E
#' aev$V
#'
#' # Calculated properties:
#' aev$A_minus_E
#' aev$A_minus_E_stddev
#' aev$A_over_E
#' aev$log_A_over_E
#' aev$log_A_over_E_stddev
#' aev$log_A_over_E_95pc
#' aev$Pearson_residual
#' aev$deviance_residual
#'
#' # Addition:
#' aev + create_aev(A = c(100, 1), E = c(200, 10), V = c(1500, 10))
#'
#' # A vector of records: `length()` counts records and `[` subsets them.
#' length(aev)
#' aev[1]
#'
#' # Labels, as a breakdown would supply them:
#' names(aev) <- c("65-70", "70-75")
#' group_names(aev) <- "age"
#' aev
#' names(aev[2])
#'
#' # Validity:
#' # (a) none of `A`, `E` or `V` are negative,
#' try(create_aev(A = -1, E = 2, V = 3))
#' # (b) `E` and `V` are either both zero or both non-zero, and
#' try(create_aev(A = 0, E = 0, V = 3))
#' # (c) `A` cannot be non-zero if `E` and `V` are zero.
#' try(create_aev(A = 1, E = 0, V = 0))
#' @name aev_properties
NULL

#' @rdname aev_properties
#' @export
create_aev <- function(A, E, V) {
  aev <- create_aev_unchecked(A, E, V)
  validate_aev(aev)
  aev
}

#' @rdname aev_properties
#' @export
`[[.aev` <- function(x, i) {

  u <- unclass(x)
  if (i == "A_minus_E_stddev") return(sqrt(u[["V"]]))
  if (i == "log_A_over_E_stddev") return(sqrt(u[["V"]]) / u[["E"]])
  if (i == "log_A_over_E_95pc") return(sqrt(u[["V"]]) / u[["E"]] * 1.9599639845400540)

  A <- u[["A"]]
  E <- u[["E"]]

  if (i == "A_minus_E") return(A - E)
  if (i == "A_over_E") return(A / E)
  if (i == "log_A_over_E") return(log(A / E))

  V <- u[["V"]]

  if (i == "Pearson_residual") return((A - E) / sqrt(V))
  if (i == "deviance_residual") return(cpp_deviance_residual(A, E, V))

  u[[i]]
}

#' @rdname aev_properties
#' @export
`$.aev` <- function(x, name) x[[name]]

#### Element labels ####
#
# AN `aev` IS A LIST OF THREE PARALLEL FIELDS PRESENTED AS A VECTOR OF RECORDS,
# which is the shape `POSIXlt` has, and it inherits the same collision. The
# list's own `names` attribute holds "A", "E" and "V", and `$` and `[[` resolve
# through it, so the per-element labels a breakdown produces cannot live there.
#
# Base R's answer for `POSIXlt` is the one taken here: `names()` is an S3 method
# reading a separate attribute, so element labels and field names coexist. The
# price is that every method describing the vector -- `length`, `[`, `names`,
# `group_names` -- has to tell the SAME STORY about how long it is, so they are
# written together and tested together. A missing one is how this pattern fails:
# before `length.aev` existed, `length(aev)` answered 3 whatever the aev held.

# The number of records, read from a field rather than through `length()`, which
# would dispatch straight back to here.
aev_records <- function(x) length(unclass(x)[["A"]])

# Labels ride along as attributes and are absent unless something set them. A
# hand-built aev has none, which is why every reader below tolerates NULL.
aev_element_names_attr <- "element_names"
aev_group_names_attr <- "group_names"

set_aev_labels <- function(x, element_names = NULL, group_names = NULL) {
  if (!is.null(element_names)) attr(x, aev_element_names_attr) <- element_names
  if (!is.null(group_names)) attr(x, aev_group_names_attr) <- group_names
  x
}

# One label level, carried through an addition. Adding two aevs element-wise
# asserts that their rows correspond, so labels that disagree say they do not --
# the same misalignment the length check refuses, one level down. An absent set
# is not a disagreement: a hand-built aev added to a labelled one keeps the
# labels rather than erasing them.
combined_aev_labels <- function(a, b, what, operands) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (!identical(a, b)) {
    stop("The AEVs ", operands, " cannot be added because their `", what, "` differ.", call. = FALSE)
  }
  a
}

# Shared by both setters. Neither recycles: an aev's labels come from a
# breakdown that already has one per element, so a short vector is a mistake
# rather than a shorthand.
checked_aev_labels <- function(value, records, what) {
  value <- as.character(value)
  if (length(value) != records) {
    stop("`", what, "` must have one label per record (", records, ").", call. = FALSE)
  }
  if (anyNA(value)) stop("`", what, "` cannot be NA.", call. = FALSE)
  value
}

#' @rdname aev_properties
#' @usage NULL
#' @export
length.aev <- function(x) aev_records(x)

#' @rdname aev_properties
#' @usage NULL
#' @export
names.aev <- function(x) attr(x, aev_element_names_attr, exact = TRUE)

#' @rdname aev_properties
#' @usage NULL
#' @export
`names<-.aev` <- function(x, value) {
  if (is.null(value)) {
    attr(x, aev_element_names_attr) <- NULL
    return(x)
  }
  attr(x, aev_element_names_attr) <- checked_aev_labels(value, aev_records(x), "names")
  x
}

#' @rdname aev_properties
#' @usage NULL
#' @export
group_names.aev <- function(x) attr(x, aev_group_names_attr, exact = TRUE)

#' @rdname aev_properties
#' @usage NULL
#' @export
`group_names<-.aev` <- function(x, value) {
  if (is.null(value)) {
    attr(x, aev_group_names_attr) <- NULL
    return(x)
  }
  # A SINGLE VALUE IS RECYCLED here and not in `names<-`, matching `includes`:
  # naming a whole band set after one thing is the ordinary case.
  if (length(value) == 1L) value <- rep(value, aev_records(x))
  attr(x, aev_group_names_attr) <- checked_aev_labels(value, aev_records(x), "group_names")
  x
}

#' @rdname aev_properties
#' @export
# SUBSETS RECORDS, NOT FIELDS. It used to be `unclass(x)[i]`, which picked
# fields: `x[2]` returned an "aev" holding only E, still full length, failing
# `is_aev()`. Harmless while every aev was length 1 and wrong the moment a
# breakdown returns fifteen.
`[.aev` <- function(x, i) {
  u <- unclass(x)
  result <- create_aev_unchecked(u[["A"]][i], u[["E"]][i], u[["V"]][i])
  set_aev_labels(result, names(x)[i], group_names(x)[i])
}

#' @rdname aev_properties
#' @usage NULL
#' @export
Ops.aev <- function(e1, e2) {

  if (missing(e2)) {
    if (.Generic == "+") { return(e1) }
    stop("Unary operator `", .Generic, "` is undefined for `", deparse(substitute(e1)), "`.", call. = FALSE)
  }

  if (is_aev(e2) && .Generic == "+") {

    u1 <- unclass(e1)
    A1 <- u1[["A"]]
    E1 <- u1[["E"]]
    V1 <- u1[["V"]]

    u2 <- unclass(e2)
    A2 <- u2[["A"]]
    E2 <- u2[["E"]]
    V2 <- u2[["V"]]

    if (length(A1) != length(A2)) {
      stop("The AEVs `", deparse(substitute(e1)), "` and `", deparse(substitute(e2)), "` cannot be added because they have different lengths.", call. = FALSE)
    }

    # NOT VALIDATED, DELIBERATELY. Adding two sound aevs cannot break the
    # invariants, so a check here only ever fires on a DEGENERATE operand -- and
    # a computed aev is allowed to be degenerate, because `exp()` underflowing to
    # zero is a legitimate answer rather than a fault (see [aev_properties]).
    # Validating the sum would put the error back exactly where it hurts most:
    # summing results is what a long-running job does at the end, after the
    # expensive part is already paid for.
    aev_sum <- create_aev_unchecked(A1 + A2, E1 + E2, V1 + V2)

    # The labels are part of the result, not decoration: without this the sum of
    # two breakdowns comes back with its rows unnamed.
    operands <- paste0("`", deparse(substitute(e1)), "` and `", deparse(substitute(e2)), "`")

    return(set_aev_labels(
      aev_sum,
      combined_aev_labels(names(e1), names(e2), "names", operands),
      combined_aev_labels(group_names(e1), group_names(e2), "group_names", operands)
    ))
  }

  stop("Binary operator `", .Generic, "` is undefined for `", deparse(substitute(e1)), "` and `", deparse(substitute(e2)), "`.", call. = FALSE)
}

#### Invalidate modification ####
#' @rdname aev_properties
#' @usage NULL
#' @export
`[[<-.aev` <- function(x, i, j, value) invalid()
#' @rdname aev_properties
#' @usage NULL
#' @export
`$<-.aev` <- function(x, name, value) invalid()
#' @rdname aev_properties
#' @usage NULL
#' @export
`[<-.aev` <- function(x, i, value) invalid()

invalid <- function() stop("Modification of an aev is invalid", call. = FALSE)

#' @returns `x`, invisibly.
#' @export
print.aev <- function(x, ...) {
  u <- unclass(x)
  cat("<aev[", aev_records(x), "]>\n", sep = "")

  frame <- data.frame(A = u$A, E = u$E, V = u$V, `A/E` = x$A_over_E, `95% conf` = x$log_A_over_E_95pc, `dev resid` = x$deviance_residual, check.names = FALSE)

  # The labels lead, so a breakdown reads down the page as the groups it came
  # from. Both are absent on a hand-built aev, and neither is invented here.
  element_names <- names(x)
  groups <- group_names(x)
  if (!is.null(element_names)) frame <- cbind(name = element_names, frame)
  if (!is.null(groups)) frame <- cbind(group = groups, frame)

  print(frame, ...)
  invisible(x)
}
#' @returns A vector of `character`.
#' @export
format.aev <- function(x, ...) {
  A_over_E <- x[["A_over_E"]]
  log_A_over_E_95pc <- x[["log_A_over_E_95pc"]]
  deviance_residual <- x[["deviance_residual"]]
  result <- sprintf("A/E = %.1f%% \u00B1%.1f%% (dev resid = %f)", A_over_E * 100, log_A_over_E_95pc * 100, deviance_residual)
  gsub("-", "\u2212", result, fixed = TRUE)
}

#### Coercion to a data frame ####

#' Coerce an `aev` to a data frame
#'
#' @description
#' One row per record, carrying the labels, the raw triple and the three
#' calculated properties a chart reads. This is the frame **logmu**'s own
#' plotting is built on, and the way to take an `aev` into **ggplot2**,
#' **dplyr** or anything else that works on data frames.
#'
#' @section Columns:
#'
#' |Column|Contents|
#' |:-----|:-------|
#' | `name`                | the element label, `names()` on the `aev`
#' | `group`               | the group label, [group_names()] on the `aev`
#' | `A`, `E`, `V`         | the triple itself
#' | `A_over_E`            | \eqn{A / E}
#' | `log_A_over_E_stddev` | \eqn{\sqrt{V} / E}
#' | `deviance_residual`   | see [aev_properties]
#'
#' The calculated columns are named after the properties that produce them, so
#' `frame$A_over_E` and `aev$A_over_E` are the same word for the same quantity.
#' The five remaining properties are left out because each is a line of
#' arithmetic on `A`, `E` and `V`, and naming them here would fix five more
#' column names for no gain.
#'
#' `name` and `group` are always present, and are `NA` on an `aev` that carries
#' no labels. The columns therefore depend on the input's type and never on its
#' values, which is what lets frames from a broken-down `aev` and an ungrouped
#' one be stacked with `rbind()`.
#'
#' @param x An `aev` object.
#' @param row.names Row names for the result, or `NULL` for the default.
#' @param optional Ignored. Present because the generic has it; every column
#'   name here is already syntactic.
#' @param ... Ignored.
#' @returns A `data.frame` with one row per record of `x`.
#' @examples
#' aev <- create_aev(A = c(1100, 40), E = c(1000, 50), V = c(2500, 125))
#' names(aev) <- c("65-70", "70-75")
#' group_names(aev) <- "age"
#'
#' as.data.frame(aev)
#'
#' # An unlabelled `aev` gives the same columns, with the labels NA.
#' as.data.frame(create_aev(A = 1100, E = 1000, V = 2500))
#' @export
as.data.frame.aev <- function(x, row.names = NULL, optional = FALSE, ...) {

  records <- aev_records(x)
  u <- unclass(x)

  # NULL BECOMES A COLUMN OF NA rather than a missing column. See *Columns*:
  # the shape has to follow the type and not the labelling, or a frame from an
  # ungrouped `aev` could not be stacked on one from a breakdown.
  label_column <- function(value) {
    if (is.null(value)) rep(NA_character_, records) else as.character(value)
  }

  # No `stringsAsFactors`: its default has been FALSE since R 4.0 and the option
  # that could override it was removed in 4.1, so passing it would say nothing.
  data.frame(
    name                = label_column(names(x)),
    group               = label_column(group_names(x)),
    A                   = u[["A"]],
    E                   = u[["E"]],
    V                   = u[["V"]],
    A_over_E            = x[["A_over_E"]],
    log_A_over_E_stddev = x[["log_A_over_E_stddev"]],
    deviance_residual   = x[["deviance_residual"]],
    row.names           = row.names
  )
}
