# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Reserved `.i$` field names that may not appear in an `it_expr`.
it_reserved_fields <- c("E2R_start", "E2R_end", "E2R_died")

# Parse an `it()` pronoun expression.
#
# @param expr  An unevaluated R expression, typically from [substitute()].
# @param env   Environment in which to evaluate constant sub-expressions
#              (defaults to the caller's frame).
#
# @return A list:
#   \describe{
#     \item{`success`}{`TRUE` if parsing succeeded.}
#     \item{`expr`}{Transformed expression, or `NULL` on failure.}
#     \item{`fields`}{Character vector of `.i$field` names referenced.}
#     \item{`uses_t`}{`TRUE` if `.t` is referenced (directly or via `.x`).}
#     \item{`error`}{Error message string, or `NULL` on success.}
#   }
it_parse <- function(expr, env = parent.frame()) {
  expr <- substitute(expr)

  fields    <- character(0L)
  uses_t    <- FALSE
  error_msg <- NULL

  # -- helpers --------------------------------------------------------------

  # e is `.i$name`.
  is_dot_i_dollar <- function(e) {
    is.call(e) &&
      length(e) == 3L &&
      identical(e[[1L]], quote(`$`)) &&
      identical(e[[2L]], quote(.i))
  }

  # TRUE iff e (post-transform) contains any symbol beginning with `.`.
  depends_on_dot <- function(e) {
    if (is.atomic(e) || is.null(e))  return(FALSE)
    if (is.symbol(e))                return(startsWith(as.character(e), "."))
    if (is.call(e))                  return(any(vapply(as.list(e), depends_on_dot, logical(1L))))
    FALSE
  }

  # Record the first error only.
  set_error <- function(msg) {
    if (is.null(error_msg)) error_msg <<- msg
  }

  # -- recursive transform ---------------------------------------------------

  transform <- function(e) {
    if (!is.null(error_msg)) return(e)   # short-circuit after first error

    # Atomic literals are valid constant leaves.
    if (is.atomic(e) || is.null(e)) return(e)

    # Symbols ----------------------------------------------------------------
    if (is.symbol(e)) {
      nm <- as.character(e)
      switch(nm,
        ".b" = {
          fields <<- union(fields, "birth")
          return(quote(.i$birth))
        },
        ".x" = {
          # `.x` means `.t - .b` which means `.t - .i$birth`
          fields <<- union(fields, "birth")
          uses_t <<- TRUE
          return(quote(.t - .i$birth))
        },
        ".t" = {
          uses_t <<- TRUE
          return(e)
        },
        ".i" = {
          set_error("`.i` must be used as `.i$field_name`, not standalone")
          return(e)
        }
      )
      if (startsWith(nm, ".")) {
        set_error(sprintf(
          "Unknown pronoun `%s`; only `.i$field`, `.t`, `.x` and `.b` are allowed", nm
        ))
        return(e)
      }
      return(e)   # ordinary name; will be folded if inside a constant call
    }

    # Calls ------------------------------------------------------------------
    if (is.call(e)) {

      # `.i$field_name` -- handle before recursing so that `.i` is never
      # seen as a bare symbol.
      if (is_dot_i_dollar(e)) {
        field_nm <- as.character(e[[3L]])
        if (field_nm %in% it_reserved_fields) {
          set_error(sprintf(
            "`.i$%s` is a reserved field and cannot be used in an `it_expr`",
            field_nm
          ))
          return(e)
        }
        fields <<- union(fields, field_nm)
        return(e)
      }

      # Recurse into every part (function symbol + arguments).
      parts  <- lapply(as.list(e), transform)
      result <- as.call(parts)

      if (!is.null(error_msg)) return(result)

      # Fold sub-expressions that are free of pronoun dependencies.
      if (!depends_on_dot(result)) {
        val <- tryCatch(
          eval(result, envir = env),
          error = function(err)
            structure(list(msg = conditionMessage(err)), class = "it_eval_error")
        )
        if (inherits(val, "it_eval_error")) {
          set_error(paste0("Could not evaluate constant sub-expression: ", val$msg))
          return(result)
        }
        if (length(val) != 1L) {
          set_error(sprintf(
            "Constant sub-expressions must be scalar (length 1), not length %d",
            length(val)
          ))
          return(result)
        }
        if (!is.numeric(val) && !is.character(val) && !datey::is_datey(val)) {
          set_error(sprintf(
            "Constant sub-expressions must be numeric, character or datey, not `%s`",
            class(val)[[1L]]
          ))
          return(result)
        }
        return(val)
      }

      return(result)
    }

    e   # pairlists etc. -- leave untouched
  }

  # -- run ------------------------------------------------------------------

  result_expr <- tryCatch(
    transform(expr),
    error = function(e) {
      error_msg <<- conditionMessage(e)
      NULL
    }
  )

  result <- if (is.null(error_msg)) {
    list(
      success = TRUE,
      expr    = result_expr,
      fields  = fields,
      uses_t  = uses_t,
      error   = NULL
    )
  } else {
    list(
      success = FALSE,
      expr    = NULL,
      fields  = NULL,
      uses_t  = NULL,
      error   = error_msg
    )
  }


  return(result)
}
