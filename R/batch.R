# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The registry of batchable operations. PRIVATE, and deliberately: `batch()`
# parallelises **logmu**'s own analyses and is not an extension point, so there
# is no published plan/finalise contract and the seam stays free to change.
#
# An operation is recognised by the head of the call, which may be qualified --
# a user with a namespace conflict writes `logmu::aev(...)` and means the same
# thing.
batch_operations <- function() {
  list(aev = list(definition = aev, plan = batch_plan_aev))
}

# The arguments a batch may supply on an operation's behalf. Each is the dotted
# formal of `batch()` minus its dot, and each is a DEFAULT: the operation's own
# value wins, the batch's applies when it has none, and it is an error only when
# neither supplies one.
#
# `threads` is NOT among them. It belongs to the run rather than to an analysis,
# it cannot change an answer (reduction is per chunk, not per thread), and once
# several operations' specifications share one crossing there is a single pool
# and no per-operation meaning left.
batch_defaulted_arguments <- c("exp_data", "include", "weight", "val_similarity",
                               "val_distance", "breakdown", "settings",
                               "overdispersion", "time_scale")

# `val_similarity` and `val_distance` are two spellings of ONE quantity, so an
# analysis naming either of them has said its piece and NEITHER default applies.
# Injected argument by argument they would combine into the pair that is
# refused, which is a confusing way to be told that a default exists.
batch_paired_arguments <- c("val_similarity", "val_distance")

#' Run several analyses over one dataset in a single pass
#'
#' @description
#' `batch()` takes named calls to **logmu** analyses and runs them together:
#' their specifications are compiled and scheduled as one crossing, the data is
#' read once, and each analysis comes back as whatever it would have returned on
#' its own.
#'
#' @section Settings are defaults, not overrides:
#' Every dotted argument supplies a default for the analyses inside. An analysis
#' naming its own value wins; the batch's applies when it has none; it is an
#' error only when neither supplies one. So `overdispersion` remains required
#' without **logmu** ever assuming a value for it.
#'
#' The dot marks a batch setting, and it is what stops a setting colliding with
#' an analysis you have named. **An element of a batch may not be named with a
#' leading dot**, which reserves the whole dotted namespace so that settings
#' added in future cannot break existing code.
#'
#' `.threads` is the exception to the default rule: it belongs to the run rather
#' than to any one analysis, so the batch's value is used throughout and a
#' `threads` argument inside a batched call is ignored. It cannot change an
#' answer, only a duration.
#'
#' @section What a batch may not do:
#' **No element may use another element's result.** Every specification is
#' compiled before a record is read and the results exist only once the pass is
#' over, so a dependency between elements could not be honoured. `batch()`
#' refuses a call that mentions a sibling's name rather than letting it resolve
#' silently against something of the same name in your workspace.
#'
#' A batch runs over **one experience dataset**. Analyses may differ in
#' `time_scale`, which is worth doing to see whether the integration interval
#' moves the answer; a batch then makes one pass per distinct scale.
#'
#' @param ... Named calls to **logmu** analyses, currently [aev()]. Each name
#'   becomes a name in the result and may not begin with a dot.
#' @param .exp_data,.include,.weight,.val_similarity,.val_distance,.breakdown,.settings,.overdispersion,.time_scale
#'   Defaults for the analyses, each standing in for the argument of the same
#'   name wherever an analysis does not supply its own. `.val_similarity` and
#'   `.val_distance` are two spellings of one quantity, so only one may be given
#'   and an analysis naming either of them takes neither default.
#' @param .threads Worker threads for the whole batch. `0` asks for as many as
#'   the machine reports.
#' @returns A named list holding each analysis's own result, in the order
#'   written. Nothing else -- no class, no attributes.
#' @examples
#' data <- exp_data(
#'   list(
#'     birth     = datey::datey(c(1945, 1950, 1955)),
#'     pension   = c(5000, 12000, 30000),
#'     male      = c(TRUE, FALSE, TRUE),
#'     E2R_start = datey::datey(c(2015, 2015, 2015)),
#'     E2R_end   = datey::datey(c(2020, 2020, 2018)),
#'     E2R_died  = c(FALSE, FALSE, TRUE)
#'   ),
#'   exp_start = datey::datey(2015),
#'   exp_end   = datey::datey(2020)
#' )
#'
#' b <- batch(
#'   .exp_data       = data,
#'   .overdispersion = 2,
#'   .weight         = .i$pension,
#'   light = aev(mortality = mortality_const(log_mu = -4.5)),
#'   heavy = aev(mortality = mortality_const(log_mu = -4.0), overdispersion = 1)
#' )
#'
#' b$light
#' b$heavy
#' @export
batch <- function(...,
                  .exp_data = NULL,
                  .include = NULL,
                  .weight = NULL,
                  .val_similarity = NULL,
                  .val_distance = NULL,
                  .breakdown = NULL,
                  .settings = NULL,
                  .overdispersion = NULL,
                  .time_scale = NULL,
                  .threads = cpp_veil_default_threads()) {

  caller <- parent.frame()
  elements <- as.list(substitute(list(...)))[-1L]
  element_names <- batch_element_names(elements)

  if (!is.null(substitute(.val_similarity)) && !is.null(substitute(.val_distance))) {
    stop("Give either `.val_similarity` or `.val_distance`, not both.", call. = FALSE)
  }

  # THE DEFAULTS TRAVEL AS EXPRESSIONS, not values. A weight is a pronoun
  # expression, so it must reach `it_capture()` unevaluated and be captured in
  # the user's own frame. An absent default and one written as `NULL` are
  # indistinguishable here, which is right: both mean there is no default.
  defaults <- list(
    exp_data       = substitute(.exp_data),
    include        = substitute(.include),
    weight         = substitute(.weight),
    val_similarity = substitute(.val_similarity),
    val_distance   = substitute(.val_distance),
    breakdown      = substitute(.breakdown),
    settings       = substitute(.settings),
    overdispersion = substitute(.overdispersion),
    time_scale     = substitute(.time_scale)
  )

  planned <- lapply(seq_along(elements), function(i) {
    batch_plan_element(elements[[i]], element_names[[i]], element_names,
                       defaults, caller)
  })

  batch_run(planned, element_names, as.integer(.threads))
}

# ---- checking the elements -------------------------------------------------

batch_element_names <- function(elements) {
  if (length(elements) == 0L) {
    stop("`batch()` needs at least one analysis.", call. = FALSE)
  }

  given <- names(elements)
  if (is.null(given) || any(!nzchar(given))) {
    stop("Every analysis in a `batch()` must be named, since the names are how ",
         "the results are found.", call. = FALSE)
  }

  # The ban is what makes the dotted namespace safe to extend: a setting added
  # in a later release cannot collide with a name a user has already chosen.
  dotted <- startsWith(given, ".")
  if (any(dotted)) {
    stop("An analysis in a `batch()` may not be named with a leading dot (",
         paste0("`", given[dotted], "`", collapse = ", "),
         "). Names beginning with a dot are batch settings.", call. = FALSE)
  }

  repeated <- unique(given[duplicated(given)])
  if (length(repeated) > 0L) {
    stop("Two analyses in a `batch()` are named ",
         paste0("`", repeated, "`", collapse = ", "),
         ". Every one needs its own name.", call. = FALSE)
  }

  given
}

# The operation's name from the head of the call, `NULL` if the expression is
# not a call to a plain or namespace-qualified function.
batch_operation_name <- function(expr) {
  if (!is.call(expr)) return(NULL)

  head <- expr[[1L]]
  if (is.call(head) && identical(as.character(head[[1L]]), "::")) {
    head <- head[[3L]]
  }
  if (!is.symbol(head)) return(NULL)

  as.character(head)
}

batch_operation <- function(expr, name) {
  operations <- batch_operations()
  operation <- batch_operation_name(expr)

  if (!is.null(operation) && !is.null(operations[[operation]])) {
    return(operations[[operation]])
  }

  # `batch(overdispersion = 2, ...)` -- a setting written without its dot --
  # arrives here as an element whose expression is not an analysis at all, and
  # its name says exactly what was meant.
  hint <- if (name %in% batch_defaulted_arguments || identical(name, "threads")) {
    sprintf(" Did you mean `.%s`?", name)
  } else {
    ""
  }

  stop("`", name, "` is not a call to a logmu analysis (",
       paste0("`", names(operations), "()`", collapse = ", "), ").", hint,
       call. = FALSE)
}

# ---- planning one element --------------------------------------------------

batch_plan_element <- function(expr, name, all_names, defaults, caller) {
  operation <- batch_operation(expr, name)

  # NO ELEMENT MAY SEE ANOTHER'S RESULT. `batch()` binds nothing, so a mention
  # of a sibling would resolve in the caller's frame -- erroring when nothing is
  # there, and silently using a stale object of that name when something is.
  mentioned <- intersect(all.vars(expr), all_names)
  if (length(mentioned) > 0L) {
    stop("`", name, "` refers to ", paste0("`", mentioned, "`", collapse = ", "),
         ", which is another analysis in the same `batch()`. No analysis in a ",
         "batch can use another's result, because they are all run in one pass.",
         call. = FALSE)
  }

  given <- as.list(match.call(definition = operation$definition, call = expr))[-1L]

  paired <- vapply(batch_paired_arguments, function(argument) !is.null(given[[argument]]),
                   logical(1L))

  for (argument in batch_defaulted_arguments) {
    if (argument %in% batch_paired_arguments && any(paired)) next
    if (is.null(given[[argument]]) && !is.null(defaults[[argument]])) {
      given[[argument]] <- defaults[[argument]]
    }
  }

  operation$plan(given, caller, name)
}

# `aev()`'s half of the protocol: turn the resolved argument expressions into
# the plan `aev()` itself would have built. Values are evaluated in the user's
# frame and pronoun expressions are captured there, which is why that frame is
# carried this far rather than being taken from here.
batch_plan_aev <- function(given, caller, name) {
  value <- function(argument) {
    if (is.null(given[[argument]])) NULL else eval(given[[argument]], caller)
  }
  ast <- function(argument) {
    if (is.null(given[[argument]])) NULL else it_capture(given[[argument]], caller)
  }

  exp_data <- value("exp_data")
  if (is.null(exp_data)) {
    stop("`", name, "` has no `exp_data`. Give it in the call, or once as ",
         "`.exp_data` for the whole batch.", call. = FALSE)
  }
  ensure_is_exp_data(exp_data)

  if (is.null(given[["mortality"]])) {
    stop("`", name, "` has no `mortality`, which is required.", call. = FALSE)
  }

  # The engine refuses this pair too, and must, since it guards every front end.
  # This one exists to name the analysis, which the engine cannot: a batch of a
  # dozen would otherwise report the fault without saying whose it is.
  if (!is.null(given[["val_similarity"]]) && !is.null(given[["val_distance"]])) {
    stop("`", name, "` gives both `val_similarity` and `val_distance`. They are ",
         "two spellings of one quantity, so give either, not both.", call. = FALSE)
  }

  resolved <- resolve_settings(value("settings"), value("overdispersion"),
                               value("time_scale"))

  list(
    exp_data          = exp_data,
    time_scale_clicks = resolved$time_scale_clicks,
    plan              = aev_plan(ast("mortality"), ast("weight"),
                                 ast("val_similarity"), ast("val_distance"),
                                 value("include"), value("breakdown"),
                                 resolved$overdispersion)
  )
}

# ---- running ---------------------------------------------------------------

# ONE CROSSING PER DISTINCT `time_scale`, since the time grid is a property of
# the run rather than of a specification. Everything else about a batch is
# shared, which is the point of it: the columns are prepared once, and every
# specification in a group is compiled and scheduled together.
batch_run <- function(planned, element_names, threads) {
  exp_data <- planned[[1L]]$exp_data
  for (i in seq_along(planned)) {
    if (!identical(planned[[i]]$exp_data, exp_data)) {
      stop("`", element_names[[i]], "` uses a different `exp_data` from `",
           element_names[[1L]], "`. A batch runs over one experience dataset; ",
           "combine the data, or use separate batches.", call. = FALSE)
    }
  }

  columns <- exp_data_columns(exp_data)
  scales <- unlist(lapply(planned, function(one) one$time_scale_clicks))

  results <- vector("list", length(planned))
  for (scale in unique(scales)) {
    group <- which(scales == scale)
    specs <- unlist(lapply(planned[group], function(one) one$plan$specs),
                    recursive = FALSE)

    run <- cpp_veil_run(specs, columns, scale, FALSE, threads)

    # The specifications went in group order, so they come back in it. Each
    # analysis takes as many as it asked for and finalises them itself.
    taken <- 0L
    for (i in group) {
      wanted <- length(planned[[i]]$plan$specs)
      results[[i]] <- planned[[i]]$plan$finalise(run$results[taken + seq_len(wanted)])
      taken <- taken + wanted
    }
  }

  names(results) <- element_names
  results
}
