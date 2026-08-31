# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# Pronoun-expression AST
# ============================================================================
#
# A pronoun expression is an R expression written with the data pronouns
# `.i$field`, `.t`, `.b` and `.x` -- for example `~ .i$pension > 0` or
# `.t - .b`.  `it_build_ast()` turns such an expression into a small, explicit
# abstract syntax tree (AST) that downstream code (typing, lowering to the
# veil engine) can walk without having to re-interpret raw R language
# objects.
#
# The AST has exactly four kinds of node, each a plain list carrying a `kind`
# tag and the S3 class "it_node":
#
#   lit(value)        a folded constant scalar, e.g. 0, "male", datey(1960)
#   field(name)       a fact access `.i$name`
#   time()            the time pronoun `.t`
#   call(fn, args)    an operation; `fn` is the function name as a string and
#                     `args` is a list of child it_nodes
#
# Two pieces of pronoun sugar are removed while building, so later passes never
# have to know about them:
#
#   .b   ->  field("birth")
#   .x   ->  call("-", list(time(), field("birth")))
#
# Parsing also folds: any maximal sub-expression that mentions no pronoun is
# evaluated once, in the supplied environment, and stored as a single `lit`.
# So `1 + 2 * 3` becomes `lit(7)`, while `.i$pension > 0` stays symbolic.

# `.i$` fields that user code may never reference -- the E2R columns are hidden.
it_reserved_fields <- c("E2R_start", "E2R_end", "E2R_died")

# ---- node constructors -----------------------------------------------------

it_lit <- function(value) {
  structure(list(kind = "lit", value = value), class = "it_node")
}

it_field <- function(name) {
  structure(list(kind = "field", name = name), class = "it_node")
}

it_time <- function() {
  structure(list(kind = "time"), class = "it_node")
}

it_call <- function(fn, args) {
  structure(list(kind = "call", fn = fn, args = args), class = "it_node")
}

# An opaque concept object spliced in by reference, e.g. a `mortality`. The
# engine (or the relevant S3 generic) knows how to evaluate it; the AST just
# carries it.
it_obj <- function(value) {
  structure(list(kind = "obj", value = value), class = "it_node")
}

is_it_node <- function(x) inherits(x, "it_node")

# ---- raw-expression helpers ------------------------------------------------

# TRUE if `e` must stay symbolic rather than being folded to a constant: it
# mentions a pronoun, or a free symbol that resolves to a concept object or a
# `~` formula (both of which we splice rather than fold). Needs `env` to tell a
# referenced object from an ordinary constant.
it_is_dynamic_expr <- function(e, env) {
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (nm %in% c(".t", ".i", ".b", ".x")) return(TRUE)
    val <- get0(nm, envir = env, inherits = TRUE)
    return(is_logmu_function(val) || inherits(val, "formula"))
  }
  if (is.call(e)) {
    return(any(vapply(as.list(e), it_is_dynamic_expr, logical(1L), env = env)))
  }
  FALSE
}

# TRUE if `e` is exactly `.i$something`.
it_is_field_access <- function(e) {
  is.call(e) &&
    length(e) == 3L &&
    identical(e[[1L]], quote(`$`)) &&
    identical(e[[2L]], quote(.i))
}

# Field name from a `.i$field` call, with the reserved-field guard.
it_field_name <- function(e) {
  rhs <- e[[3L]]
  name <- if (is.symbol(rhs)) {
    as.character(rhs)
  } else if (is.character(rhs) && length(rhs) == 1L) {
    rhs
  } else {
    stop("`.i$` must be followed by a field name, as in `.i$pension`.", call. = FALSE)
  }
  if (name %in% it_reserved_fields) {
    stop(sprintf("`.i$%s` is a reserved field and cannot be used in a pronoun expression.", name),
         call. = FALSE)
  }
  name
}

# Function name of a call as a string ("+", ">", "if", "pkg::fn", ...).
it_call_name <- function(e) paste(deparse(e[[1L]]), collapse = "")

# ---- allowed vocabulary ----------------------------------------------------
#
# Pronoun-dependent (dynamic) calls may only use this vocabulary; everything
# else is rejected at parse. Constant sub-expressions are exempt -- they are
# folded by ordinary R evaluation and may use any function.

# Canonical operators / functions (after alias normalisation).
it_allowed_ops <- c(
  "+", "-", "*", "/", "^", "%%", "%/%",                    # arithmetic
  "==", "!=", "<", "<=", ">", ">=",                        # comparison
  "&", "|", "!", "xor", "is.na",                           # logical
  "%in%",                                                  # membership
  "if", "ifelse",                                          # conditional
  "abs", "sqrt", "exp", "log", "log1p", "log10", "expm1",  # unary math
  "sin", "cos", "floor", "ceiling", "round", "trunc", "sign",
  "max", "min", "clamp"                                    # variadic / ternary
)

# Idiomatic R that we accept and normalise to a canonical head. Safe because
# expressions are scalar and pure: `&&`/`||` short-circuiting is irrelevant and
# `pmax`/`pmin` equal `max`/`min` on scalars.
it_op_aliases <- c("&&" = "&", "||" = "|", "pmax" = "max", "pmin" = "min")

# Namespaces a call head may carry (so `base::abs(.i$x)` works).
it_allowed_namespaces <- c("base", "stats", "logmu")

# Constructs that are never allowed anywhere in an expression -- assignment and
# control flow have side effects or no meaning in a pure scalar function.
it_forbidden_heads <- c("<-", "=", "<<-", "->", "->>",
                        "for", "while", "repeat", "break", "next", "function", "~")

# Split a (possibly namespaced) head into namespace + bare name.
it_split_op <- function(fn) {
  parts <- strsplit(fn, "::", fixed = TRUE)[[1L]]
  list(ns = if (length(parts) == 2L) parts[[1L]] else "", name = parts[[length(parts)]])
}

# Normalise an alias to its canonical head, keeping any namespace.
it_normalise_op <- function(fn) {
  p <- it_split_op(fn)
  name <- if (p$name %in% names(it_op_aliases)) it_op_aliases[[p$name]] else p$name
  if (nzchar(p$ns)) paste0(p$ns, "::", name) else name
}

# Error unless `fn` is a permitted (possibly namespaced) operator/function.
it_check_allowed <- function(fn) {
  p <- it_split_op(fn)
  if (nzchar(p$ns) && !p$ns %in% it_allowed_namespaces) {
    stop(sprintf("Namespace `%s::` is not allowed in a pronoun expression.", p$ns), call. = FALSE)
  }
  if (!p$name %in% it_allowed_ops) {
    stop(sprintf("`%s` is not an allowed operator or function in a pronoun expression.", fn),
         call. = FALSE)
  }
  invisible()
}

# A handful of arity rules.
it_check_arity <- function(fn, n) {
  name <- it_split_op(fn)$name
  if (name == "log" && n != 1L) {
    stop("`log(x, base)` is not supported; use `log()` (natural), `log10()` or `log1p()`.",
         call. = FALSE)
  }
  if (name == "clamp" && n != 3L) {
    stop("`clamp(x, lo, hi)` takes exactly three arguments.", call. = FALSE)
  }
  invisible()
}

# The condition of `if`/`ifelse` must be time-invariant: a time-varying shape
# must be written branchlessly (clamp / min / max) so it vectorises.
it_check_conditional <- function(node) {
  if (node$fn %in% c("if", "ifelse") && it_uses_t(node$args[[1L]])) {
    stop(sprintf(
      paste0("The condition of `%s` must not depend on time `.t`; use `clamp()`, ",
             "`min()` or `max()` for a time-varying shape."), node$fn),
      call. = FALSE)
  }
  invisible()
}

# Recursively reject any forbidden construct, including inside subtrees that
# would otherwise be folded (so `{ z <- 1; z }` fails rather than running).
it_check_forbidden <- function(e) {
  if (is.call(e)) {
    head <- e[[1L]]
    if (is.symbol(head)) {
      hn <- as.character(head)
      if (hn %in% it_forbidden_heads) {
        stop(sprintf("`%s` is not allowed in a pronoun expression.", hn), call. = FALSE)
      }
      if (hn == "{" && length(e) != 2L) {
        stop("A block `{ ... }` with multiple statements is not allowed in a pronoun expression.",
             call. = FALSE)
      }
    }
    for (part in as.list(e)) it_check_forbidden(part)
  }
  invisible()
}

# ---- folding ---------------------------------------------------------------

# Evaluate a pronoun-free expression and turn the result into a node. Usually
# this is a constant `lit`, but an inline constructor call may produce a concept
# object (-> obj leaf / spliced AST) or a `~` formula, so we classify it.
it_fold <- function(e, env) {
  value <- tryCatch(
    eval(e, envir = env),
    error = function(err) {
      stop(sprintf("Could not evaluate the constant sub-expression `%s`: %s",
                   it_short_deparse(e), conditionMessage(err)),
           call. = FALSE)
    }
  )
  it_classify_value(value, e)
}

# Fold the right-hand side of `%in%`, which (unlike an ordinary constant) may
# be a character or numeric vector -- the membership set.
it_fold_set <- function(e, env) {
  value <- tryCatch(
    eval(e, envir = env),
    error = function(err) {
      stop(sprintf("Could not evaluate the `%%in%%` set `%s`: %s",
                   it_short_deparse(e), conditionMessage(err)),
           call. = FALSE)
    }
  )
  if (!is.atomic(value) || length(value) < 1L || !(is.character(value) || is.numeric(value))) {
    stop("The right-hand side of `%in%` must be a character or numeric vector.", call. = FALSE)
  }
  it_lit(value)
}

# A folded constant must be a single scalar of a type the engine understands.
# Anything else (a vector, a list, a `mortality`/`include` object, a function)
# is rejected here.  NB: splicing named concept objects into the tree is a
# separate, later step; until then they correctly fail this check.
it_check_constant <- function(value, e) {
  ok_atomic <- is.atomic(value) &&
    length(value) == 1L &&
    (is.logical(value) || is.numeric(value) || is.character(value))
  ok_time <- (datey::is_datey(value) || datey::is_durationy(value)) &&
    length(value) == 1L
  if (!ok_atomic && !ok_time) {
    stop(sprintf(
      paste0("`%s` did not fold to a single constant. A pronoun expression may only ",
             "combine `.i$field`, `.t`, `.b`, `.x` and scalar constants (got %s of length %d)."),
      it_short_deparse(e), class(value)[[1L]], length(value)),
      call. = FALSE)
  }
  it_check_not_missing(value, e)
}

# A folded constant may not be missing. Two separate reasons arriving at one
# error: a NaN double is provably NaN before any data is read, which makes the
# expression malformed rather than merely useless, and an NA logical, integer or
# string has no representation in the engine at all -- only `double`, `datey` and
# `durationy` carry a portable missing value, so the others cannot express one.
#
# NaN ARRIVING FROM THE DATA IS A DIFFERENT THING AND STAYS LEGAL. A column may
# hold NA, a calculation over it yields NaN, and that NaN is a legitimate result
# carried all the way to the aev. What is refused here is a missing value WRITTEN
# INTO THE EXPRESSION, where no dataset exists that could make it mean anything.
#
# `-Inf` is deliberately NOT caught. It is clean IEEE, and in a mortality it is
# the neutral value -- log mu of -Inf is a mortality rate of zero, contributing
# no exposure. Anything added here must test for missingness specifically.
it_check_not_missing <- function(value, e) {
  if (!is.na(value)) return(value)

  reason <- if (is.nan(value)) {
    "is NaN whatever the data holds"
  } else {
    "is missing"
  }
  stop(sprintf(
    paste0("`%s` %s, so it cannot appear in a pronoun expression. A missing value in the ",
           "DATA is fine and flows through to a NaN result; one written into the expression ",
           "has no data that could give it a meaning."),
    it_short_deparse(e), reason),
    call. = FALSE)
}

# Resolve a free symbol to a node: splice a referenced concept object or `~`
# formula, otherwise fold it to a constant.
it_resolve_symbol <- function(e, env) {
  value <- tryCatch(
    eval(e, envir = env),
    error = function(err) {
      stop(sprintf("Could not resolve `%s`: %s", it_short_deparse(e), conditionMessage(err)),
           call. = FALSE)
    }
  )
  it_classify_value(value, e)
}

it_classify_value <- function(value, e) {
  # A referenced concept object: splice its AST if it has one (indicator,
  # variable), otherwise keep it as an opaque leaf (mortality, include).
  if (is_logmu_function(value)) {
    if (is.list(value) && !is.null(value$ast)) return(value$ast)
    return(it_obj(value))
  }
  # A referenced `~` formula: splice its body, parsed in its own environment.
  if (inherits(value, "formula")) {
    return(it_build_ast(value[[length(value)]], environment(value)))
  }
  # Otherwise an ordinary constant.
  it_lit(it_check_constant(value, e))
}

# ---- the builder ------------------------------------------------------------

# Public entry: validate the whole expression for forbidden constructs once,
# then build it. (Recursion uses `it_build`, which does not re-scan.)
it_build_ast <- function(e, env = parent.frame()) {
  it_check_forbidden(e)
  it_build(e, env)
}

# Turn an unevaluated expression `e` into an it_node, folding constants in `env`.
it_build <- function(e, env) {

  # Already-literal scalar (a number, string or logical written in the source).
  if (is.atomic(e)) {
    return(it_fold(e, env))
  }

  # A bare name: either a pronoun, or a free variable to splice/fold.
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (nm == ".t") return(it_time())
    if (nm == ".b") return(it_field("birth"))
    if (nm == ".x") return(it_call("-", list(it_time(), it_field("birth"))))
    if (nm == ".i") {
      stop("`.i` must be used as `.i$field`, not on its own.", call. = FALSE)
    }
    if (startsWith(nm, ".")) {
      stop(sprintf("Unknown pronoun `%s`; the pronouns are `.i$field`, `.t`, `.b` and `.x`.", nm),
           call. = FALSE)
    }
    # A free symbol: splice a referenced object/formula, or fold a constant.
    return(it_resolve_symbol(e, env))
  }

  # A call of some kind.
  if (is.call(e)) {
    head <- e[[1L]]

    # Redundant parentheses and single-expression braces are transparent.
    if (identical(head, quote(`(`)) && length(e) == 2L) return(it_build(e[[2L]], env))
    if (identical(head, quote(`{`)) && length(e) == 2L) return(it_build(e[[2L]], env))

    # `.i$field` -- handled before recursion so `.i` is never seen alone.
    if (it_is_field_access(e)) {
      return(it_field(it_field_name(e)))
    }

    # A maximal sub-expression with no pronouns and no spliced references:
    # evaluate it once, now. (Constant folding may use any function.)
    if (!it_is_dynamic_expr(e, env)) {
      return(it_fold(e, env))
    }

    # A genuine pronoun/object-bearing operation: it must use the allowed
    # vocabulary, and we recurse into its arguments.
    fn <- it_normalise_op(it_call_name(e))
    it_check_allowed(fn)

    # `%in%` is special: its right-hand side folds to a vector membership set.
    if (identical(fn, "%in%")) {
      if (length(e) != 3L) stop("`%in%` takes two arguments.", call. = FALSE)
      return(it_call("%in%", list(it_build(e[[2L]], env), it_fold_set(e[[3L]], env))))
    }

    args <- lapply(as.list(e)[-1L], it_build, env = env)
    it_check_arity(fn, length(args))
    node <- it_call(fn, args)
    it_check_conditional(node)
    return(node)
  }

  stop(sprintf("Cannot parse this expression: `%s`.", it_short_deparse(e)), call. = FALSE)
}

# ---- public entry -----------------------------------------------------------

# Build the AST from an already-captured (substituted) argument `q`. A `~ ...`
# formula is parsed on its right-hand side, in the formula's own environment;
# anything else is parsed as a bare expression in `env`. This is the single
# coercion path shared by `it_ast()`, `pronoun_expressions()` and the concept
# constructors, so every form folds the same way.
it_capture <- function(q, env) {
  if (is.call(q) && identical(q[[1L]], quote(`~`))) {
    f <- eval(q, env)
    return(it_build_ast(f[[length(f)]], environment(f)))
  }
  it_build_ast(q, env)
}

# Capture either a bare pronoun expression or a `~` formula and build its AST.
# Returns the root it_node.
it_ast <- function(expr, env = parent.frame()) {
  it_capture(substitute(expr), env)
}

# ============================================================================
# Analysis over the AST
# ============================================================================

# Does the expression depend on time `.t`?
it_uses_t <- function(node) {
  switch(node$kind,
    time = TRUE,
    call = any(vapply(node$args, it_uses_t, logical(1L))),
    FALSE
  )
}

# Which `.i$` fields does the expression reference?
it_fields <- function(node) {
  if (identical(node$kind, "field")) return(node$name)
  if (identical(node$kind, "call")) {
    found <- unlist(lapply(node$args, it_fields), use.names = FALSE)
    return(unique(if (is.null(found)) character(0L) else found))
  }
  character(0L)
}

# Operations whose result is logical.
it_logical_ops <- c(
  "==", "!=", "<", "<=", ">", ">=",
  "&", "|", "!", "xor", "is.na", "%in%"
)

# Does the top of the expression produce a logical value?
it_is_logical <- function(node) {
  switch(node$kind,
    lit = is.logical(node$value),
    call = node$fn %in% it_logical_ops,
    FALSE
  )
}

# An indicator is a logical-valued expression that does not depend on time.
# This is the structural check that lets `aev()` know weight is {0,1} so V = E.
it_is_indicator <- function(node) {
  it_is_logical(node) && !it_uses_t(node)
}

# The values of every obj leaf in the tree (the referenced concept objects).
it_obj_leaves <- function(node) {
  switch(node$kind,
    obj = list(node$value),
    call = unlist(lapply(node$args, it_obj_leaves), recursive = FALSE),
    list()
  )
}

# ============================================================================
# Reference evaluation
# ============================================================================
#
# `it_eval()` evaluates an AST against a single individual's facts `.i` (a
# named list) and a time vector `.t`. This is the plain-R reference evaluator
# used for testing and teaching -- NOT the performance path, which lowers the
# AST to the veil engine. `.t` may be omitted for a time-invariant tree.

# Resolve a call head (bare operator/function name, or `pkg::fn`) to a function.
# `get()` handles operators like `*`/`==`/`%in%` that `str2lang()` cannot parse.
it_resolve_fn <- function(fn) {
  if (grepl("::", fn, fixed = TRUE)) {
    parts <- strsplit(fn, "::", fixed = TRUE)[[1L]]
    return(getExportedValue(parts[[1L]], parts[[2L]]))
  }
  get(fn, mode = "function")
}

# `obj_fn(value, .i, .t)` controls how an obj leaf evaluates. The default
# returns the object itself (mortality selection yields the chosen object);
# a `mortality` evaluates obj leaves via `log_mu`, turning the tree into a
# scalar log-mu formula.
it_eval <- function(node, .i, .t = NULL, obj_fn = function(o, .i, .t) o) {
  switch(node$kind,
    lit = node$value,
    obj = obj_fn(node$value, .i, .t),
    time = {
      if (is.null(.t)) {
        stop("This expression uses `.t`, but no `.t` was supplied.", call. = FALSE)
      }
      .t
    },
    field = {
      if (!is.list(.i) || !node$name %in% names(.i)) {
        stop(sprintf("Field `.i$%s` was not supplied in `.i`.", node$name), call. = FALSE)
      }
      .i[[node$name]]
    },
    call = {
      # `if`/`ifelse` are a per-individual selection (their condition is
      # time-invariant, hence scalar), not a per-element vector op.
      if (node$fn %in% c("if", "ifelse")) {
        cond <- it_eval(node$args[[1L]], .i = .i, .t = .t, obj_fn = obj_fn)
        branch <- if (isTRUE(as.logical(cond))) node$args[[2L]] else node$args[[3L]]
        return(it_eval(branch, .i = .i, .t = .t, obj_fn = obj_fn))
      }
      args <- lapply(node$args, it_eval, .i = .i, .t = .t, obj_fn = obj_fn)
      # Scalar semantics: `max`/`min` are element-wise (== scalar-per-time),
      # and `clamp` has no base equivalent.
      if (identical(node$fn, "max")) return(do.call(pmax, args))
      if (identical(node$fn, "min")) return(do.call(pmin, args))
      if (identical(node$fn, "clamp")) return(pmin(pmax(args[[1L]], args[[2L]]), args[[3L]]))
      do.call(it_resolve_fn(node$fn), args)
    },
    stop(sprintf("Cannot evaluate an AST node of kind `%s`.", node$kind), call. = FALSE)
  )
}

# ============================================================================
# Inspection
# ============================================================================

it_format_value <- function(v) {
  if (length(v) != 1L) {
    return(paste0("c(", paste(vapply(v, it_format_value, character(1L)), collapse = ", "), ")"))
  }
  if (is.character(v)) return(encodeString(v, quote = "\""))
  paste(format(v), collapse = " ")
}

it_node_label <- function(node) {
  switch(node$kind,
    lit   = paste0("lit ", it_format_value(node$value)),
    obj   = paste0("obj <", class(node$value)[[1L]], ">"),
    field = paste0("field .i$", node$name),
    time  = "time .t",
    call  = paste0("call `", node$fn, "`"),
    paste0("?", node$kind)
  )
}

# Build an indented-tree rendering of the AST as a character vector of lines.
it_tree_lines <- function(node, indent = "") {
  lines <- paste0(indent, it_node_label(node))
  if (identical(node$kind, "call")) {
    for (child in node$args) {
      lines <- c(lines, it_tree_lines(child, paste0(indent, "  ")))
    }
  }
  lines
}

# Print the AST as an indented tree. (Also the S3 print method once the
# package NAMESPACE has been regenerated with `devtools::document()`.)
it_tree <- function(node) {
  cat(it_tree_lines(node), sep = "\n")
  invisible(node)
}

#' @export
print.it_node <- function(x, ...) it_tree(x)

# Operators rendered infix by `it_deparse()`.
it_infix_ops <- c(
  "+", "-", "*", "/", "^", "%%", "%/%",
  "==", "!=", "<", "<=", ">", ">=",
  "&", "|", "%in%"
)

# Reconstruct a readable, R-ish one-line string from an AST. Round-trips the
# pronouns, so `it_deparse(it_build_ast(quote(.i$pension > 0)))` is
# "(.i$pension > 0)".
it_deparse <- function(node) {
  switch(node$kind,
    lit   = it_format_value(node$value),
    obj   = paste0("<", class(node$value)[[1L]], ">"),
    field = paste0(".i$", node$name),
    time  = ".t",
    call  = {
      parts <- vapply(node$args, it_deparse, character(1L))
      if (identical(node$fn, "if") && length(parts) >= 2L) {
        els <- if (length(parts) >= 3L) paste0(" else ", parts[[3L]]) else ""
        paste0("if (", parts[[1L]], ") ", parts[[2L]], els)
      } else if (length(parts) == 2L && node$fn %in% it_infix_ops) {
        paste0("(", parts[1L], " ", node$fn, " ", parts[2L], ")")
      } else if (length(parts) == 1L && node$fn %in% c("-", "+", "!")) {
        paste0(node$fn, parts[1L])
      } else {
        paste0(node$fn, "(", paste(parts, collapse = ", "), ")")
      }
    }
  )
}

# Short, single-line deparse of a raw expression for error messages.
it_short_deparse <- function(e) {
  txt <- paste(deparse(e), collapse = " ")
  if (nchar(txt) > 60L) paste0(substr(txt, 1L, 57L), "...") else txt
}
