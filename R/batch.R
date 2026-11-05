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
  list(aev = list(definition = aev, plan = batch_plan_aev),
       fit = list(definition = fit, plan = batch_plan_fit),
       compare_models = list(definition = compare_models, plan = batch_plan_compare_models))
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
#
# A default an analysis has no use for is simply not read: `.breakdown` means
# nothing to a fit and `.test_mortality` nothing to an A/E, so one batch can
# hold both without either refusing the other's settings.
batch_defaulted_arguments <- c("exp_data", "include", "weight", "val_relevance",
                               "val_distance", "breakdown", "test_mortality", "Z",
                               "settings", "overdispersion", "time_scale",
                               "L_tolerance", "armijo", "max_iterations",
                               "max_halvings")

# PAIRS THAT ARE ONE QUANTITY SAID TWO WAYS, so an analysis naming either half
# has said its piece and NEITHER default applies. Injected argument by argument
# they would combine into the pair that is refused, which is a confusing way to
# be told that a default exists.
#
# `val_relevance` and `val_distance` are two spellings of the second weighting
# factor. `test_mortality` and `Z` are the yardstick and the number it would be
# used to work out.
batch_paired_arguments <- list(c("val_relevance", "val_distance"),
                               c("test_mortality", "Z"))

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
#' @section Fits in a batch:
#' [fit()] and [compare_models()] may be batched alongside [aev()]. The pass
#' each makes over the population, for \eqn{Z} and the figures its header
#' reports, joins the batch's shared pass. Once that pass is over, every model
#' in the batch -- each fit's, and each candidate of each comparison -- is
#' fitted together: each iteration's walk of the data serves every model still
#' iterating, and each converges on its own. The answers are exactly those of
#' the same calls made alone; only the time taken changes.
#'
#' **A fit that fails does not stop the batch.** Where a lone `fit()` would
#' raise an error because the model is not identifiable or did not converge,
#' its element holds that error instead -- a condition of class `fit_failure`,
#' carrying the message and the engine's `status` -- and one warning names
#' every element that failed. `stop()` on the element raises the error exactly
#' as the lone call would have. Anything else that goes wrong, such as a
#' malformed argument or an interrupt, still stops the whole batch.
#'
#' @param ... Named calls to **logmu** analyses: [aev()], [fit()] or
#'   [compare_models()]. Each name becomes a name in the result and may not
#'   begin with a dot.
#' @param .exp_data,.include,.weight,.val_relevance,.val_distance,.breakdown,.test_mortality,.Z,.settings,.overdispersion,.time_scale,.L_tolerance,.armijo,.max_iterations,.max_halvings
#'   Defaults for the analyses, each standing in for the argument of the same
#'   name wherever an analysis does not supply its own, and ignored by an
#'   analysis that has no such argument. `.val_relevance` and `.val_distance`
#'   are two spellings of one quantity, as `.test_mortality` and `.Z` are two
#'   ways of setting one scale, so only one of each pair may be given, and an
#'   analysis naming either half of a pair takes neither default.
#' @param .threads Worker threads for the whole batch. `0` asks for as many as
#'   the machine reports.
#' @returns A named list holding each analysis's own result, in the order
#'   written, or a `fit_failure` condition for a fit that failed. Nothing else
#'   -- no class, no attributes.
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
                  .val_relevance = NULL,
                  .val_distance = NULL,
                  .breakdown = NULL,
                  .test_mortality = NULL,
                  .Z = NULL,
                  .settings = NULL,
                  .overdispersion = NULL,
                  .time_scale = NULL,
                  .L_tolerance = NULL,
                  .armijo = NULL,
                  .max_iterations = NULL,
                  .max_halvings = NULL,
                  .threads = cpp_veil_default_threads()) {

  caller <- parent.frame()
  elements <- as.list(substitute(list(...)))[-1L]
  element_names <- batch_element_names(elements)

  # THE DEFAULTS TRAVEL AS EXPRESSIONS, not values. A weight is a pronoun
  # expression, so it must reach `it_capture()` unevaluated and be captured in
  # the user's own frame. An absent default and one written as `NULL` are
  # indistinguishable here, which is right: both mean there is no default.
  defaults <- list(
    exp_data       = substitute(.exp_data),
    include        = substitute(.include),
    weight         = substitute(.weight),
    val_relevance  = substitute(.val_relevance),
    val_distance   = substitute(.val_distance),
    breakdown      = substitute(.breakdown),
    test_mortality = substitute(.test_mortality),
    Z              = substitute(.Z),
    settings       = substitute(.settings),
    overdispersion = substitute(.overdispersion),
    time_scale     = substitute(.time_scale),
    L_tolerance    = substitute(.L_tolerance),
    armijo         = substitute(.armijo),
    max_iterations = substitute(.max_iterations),
    max_halvings   = substitute(.max_halvings)
  )

  for (pair in batch_paired_arguments) {
    if (!is.null(defaults[[pair[[1L]]]]) && !is.null(defaults[[pair[[2L]]]])) {
      stop(sprintf("Give either `.%s` or `.%s`, not both.", pair[[1L]], pair[[2L]]),
           call. = FALSE)
    }
  }

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
  ensure_are_valid_names(given, "names")

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
  mentioned <- intersect(batch_mentioned_names(expr), all_names)
  if (length(mentioned) > 0L) {
    stop("`", name, "` refers to ", paste0("`", mentioned, "`", collapse = ", "),
         ", which is another analysis in the same `batch()`. No analysis in a ",
         "batch can use another's result, because they are all run in one pass.",
         call. = FALSE)
  }

  given <- as.list(match.call(definition = operation$definition, call = expr))[-1L]

  # Asked before anything is injected, so that a default cannot make an
  # analysis look as though it had named half of a pair itself.
  pair_named <- function(argument) {
    for (pair in batch_paired_arguments) {
      if (argument %in% pair) {
        return(any(vapply(pair, function(half) !is.null(given[[half]]), logical(1L))))
      }
    }
    FALSE
  }
  skipped <- vapply(batch_defaulted_arguments, pair_named, logical(1L))

  for (argument in batch_defaulted_arguments[!skipped]) {
    if (is.null(given[[argument]]) && !is.null(defaults[[argument]])) {
      given[[argument]] <- defaults[[argument]]
    }
  }

  operation$plan(given, caller, name)
}

# The symbols an element's expression could resolve in the caller's frame.
#
# NOT `all.vars()`, which also returns the field after a `$`. So `.i$group`
# mentions `group`, and an element named after a column it reads -- a natural
# thing to call a fit by group -- was refused as though it used a sibling's
# result. What follows a `$` or an `@` is a name inside an object, never a
# variable, so only the object before it is walked. A function's own name is
# left out, as `all.vars()` leaves it out.
batch_mentioned_names <- function(expr) {
  if (is.symbol(expr)) return(as.character(expr))
  if (!is.call(expr)) return(character())

  head <- expr[[1L]]
  if (identical(head, quote(`$`)) || identical(head, quote(`@`))) {
    return(batch_mentioned_names(expr[[2L]]))
  }
  # `logmu::fit` names a function in a namespace, not two variables.
  if (identical(head, quote(`::`)) || identical(head, quote(`:::`))) {
    return(character())
  }

  parts <- as.list(expr)[-1L]
  if (is.call(head)) parts <- c(list(head), parts)
  unique(unlist(lapply(parts, batch_mentioned_names), use.names = FALSE))
}

# The engine and the front ends refuse both halves of a pair too, and must,
# since they guard every caller. This one exists to name the analysis, which
# they cannot: a batch of a dozen would otherwise report the fault without
# saying whose it is.
batch_refuse_both <- function(given, name, pair, why) {
  if (!is.null(given[[pair[[1L]]]]) && !is.null(given[[pair[[2L]]]])) {
    stop("`", name, "` gives both `", pair[[1L]], "` and `", pair[[2L]], "`. ", why,
         call. = FALSE)
  }
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

  batch_refuse_both(given, name, c("val_relevance", "val_distance"),
                    "They are two spellings of one quantity, so give either, not both.")

  resolved <- resolve_settings(value("settings"), value("overdispersion"),
                               value("time_scale"))

  list(
    exp_data          = exp_data,
    time_scale_clicks = resolved$time_scale_clicks,
    plan              = aev_plan(ast("mortality"), ast("weight"),
                                 ast("val_relevance"), ast("val_distance"),
                                 value("include"), value("breakdown"),
                                 resolved$overdispersion)
  )
}

# `fit()`'s half of the protocol, and `compare_models()`'s below it. Both take
# the same arguments apart from the model, so what they share is resolved here,
# in the order the lone calls resolve it, and each planner adds its own.
batch_fitting_inputs <- function(given, caller, name) {
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

  batch_refuse_both(given, name, c("val_relevance", "val_distance"),
                    "They are two spellings of one quantity, so give either, not both.")
  batch_refuse_both(given, name, c("test_mortality", "Z"),
                    "`Z` is the number a test mortality would be used to work out, so give either, not both.")

  list(
    value = value,
    exp_data = exp_data,
    resolved = resolve_settings(value("settings"), value("overdispersion"),
                                value("time_scale")),
    weight_ast = ast("weight"),
    relevance_ast = ast("val_relevance"),
    distance_ast = ast("val_distance"),
    test = fit_test_mortality(given[["test_mortality"]], is.null(given[["test_mortality"]]),
                              value("Z"), caller),
    test_expr = given[["test_mortality"]],
    include = value("include"),
    include_expr = given[["include"]],
    Z = value("Z"),

    # THE SIGNATURE'S OWN DEFAULTS where the analysis and the batch are both
    # silent, named once in `fit.R` so that the two cannot drift apart.
    method = function(start) {
      fit_method(start,
                 value("L_tolerance") %||% default_L_tolerance,
                 value("armijo") %||% default_armijo,
                 value("max_iterations") %||% default_max_iterations,
                 value("max_halvings") %||% default_max_halvings)
    }
  )
}

batch_plan_fit <- function(given, caller, name) {
  if (is.null(given[["model"]])) {
    stop("`", name, "` has no `model`, which is required.", call. = FALSE)
  }
  model <- eval(given[["model"]], caller)
  ensure_is_model(model)

  inputs <- batch_fitting_inputs(given, caller, name)

  list(
    exp_data          = inputs$exp_data,
    time_scale_clicks = inputs$resolved$time_scale_clicks,
    plan              = fit_plan(model, inputs$include, inputs$include_expr,
                                 inputs$weight_ast, inputs$relevance_ast,
                                 inputs$distance_ast, inputs$test, inputs$test_expr,
                                 inputs$Z, inputs$resolved,
                                 inputs$method(inputs$value("start")))
  )
}

batch_plan_compare_models <- function(given, caller, name) {
  if (is.null(given[["models"]])) {
    stop("`", name, "` has no `models`, which is required.", call. = FALSE)
  }
  models <- eval(given[["models"]], caller)
  labels <- comparison_labels(models)

  inputs <- batch_fitting_inputs(given, caller, name)

  list(
    exp_data          = inputs$exp_data,
    time_scale_clicks = inputs$resolved$time_scale_clicks,
    plan              = comparison_plan(models, labels, inputs$include,
                                        inputs$include_expr, inputs$weight_ast,
                                        inputs$relevance_ast, inputs$distance_ast,
                                        inputs$test, inputs$test_expr, inputs$Z,
                                        inputs$resolved, inputs$method(NULL))
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

  # Each group is run by `run_plans()`, the runner a lone analysis uses too:
  # one crossing for every specification in the group, then every fit in the
  # group together, sharing each walk of the data.
  #
  # A FAILED FIT IS HELD, NOT RAISED, and only a failed fit. One model that
  # will not converge must not take every other result in the batch with it,
  # so its `fit_failure` becomes its element. Any other error still stops the
  # batch: an interrupt has to, and a malformed argument is the user's to fix
  # before anything is worth keeping.
  results <- vector("list", length(planned))
  for (scale in unique(scales)) {
    group <- which(scales == scale)
    results[group] <- run_plans(lapply(planned[group], function(one) one$plan),
                                columns, scale, threads, hold_failures = TRUE)
  }

  names(results) <- element_names
  batch_warn_failures(results)
  results
}

# WARNED ABOUT ONCE, NAMING THE ELEMENTS, for the reason `compare_models()`
# gives about its candidates: a result silently replaced by an error is how
# somebody comes to believe they have ten fits when they have nine.
batch_warn_failures <- function(results) {
  failed <- vapply(results, inherits, logical(1L), what = "fit_failure")
  if (!any(failed)) return(invisible(NULL))

  statuses <- vapply(results[failed], function(failure) failure$status, character(1L))
  warning(sprintf(paste0("%d of %d analys%s in the batch could not be fitted: %s. ",
                         "Each holds the error it would have raised alone, and ",
                         "`stop()` on it raises that error."),
                  sum(failed), length(results),
                  if (length(results) == 1L) "is" else "es",
                  paste(sprintf("`%s` (%s)", names(results)[failed], statuses),
                        collapse = ", ")),
          call. = FALSE)
}
