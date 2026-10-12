# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `name()` is the one reader of an object's own name and `name<-` the one
# writer, and every label in the package -- typed or generated -- passes the
# one rule in `is_valid_name()`. Each entry point gets its own witness here, so
# that dropping the check at any one of them fails a test that names it.

reference <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 1)

group_data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948, 1952, 1938, 1960, 1943)),
    group     = c("a", "a", "a", "b", "b", "b", "c", "c", "c"),
    pension   = c(5000, 12000, 30000, 8000, 15000, 22000, 4000, 9000, 40000),
    E2R_start = datey::datey(rep(2015, 9)),
    E2R_end   = datey::datey(c(2020, 2018, 2020, 2019, 2020, 2020, 2020, 2017, 2020)),
    E2R_died  = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

# ---- the accessor ------------------------------------------------------------

test_that("name() reads the name every constructor stores", {
  expect_identical(name(mortality_const(log_mu = -4, name = "flat")), "flat")
  expect_identical(name(gompertz_mortality(name = "gomp")), "gomp")
  expect_identical(name(mortality_table(x0 = datey::durationy(60), t0 = datey::datey(2000),
                                        log_mu = matrix(-4, 2, 2), name = "tbl")), "tbl")
  expect_identical(name(model(reference, covariates(level = 1), name = "m")), "m")
  expect_null(name(reference))
  expect_null(name(model(reference, covariates(level = 1))))
})

test_that("an include carries a name and a group name, and an indicator likewise", {
  std <- ages(65, 95, by = 15)
  expect_identical(name(std[[1]]), "65-80")
  expect_identical(group_name(std[[1]]), "age")
  expect_identical(name(std[[1]]), names(std)[[1]])
  expect_identical(group_name(std[[1]]), group_names(std)[[1]])

  cat_a <- categories(.i$group, "a", "b")[[1]]
  expect_true(is_indicator(cat_a))
  expect_identical(name(cat_a), "a")
  expect_identical(group_name(cat_a), "group")

  expect_null(name(include(.i$pension > 0)))
  expect_null(group_name(include(.i$pension > 0)))
})

test_that("a covariate term's name is its own name", {
  terms <- covariates(male = .i$group == "a", .i$group == "b")
  expect_identical(name(terms[[1]]), "male")
  expect_null(name(terms[[2]]))
  # And it survives the algebra, which is why it rides on the item.
  shaped <- terms * variable(.x)
  expect_identical(name(shaped[[1]]), "male")
  crossed <- covariates(a = .i$group == "a") * covariates(b = .i$group == "b")
  expect_identical(name(crossed[[1]]), "a.b")
})

test_that("name<- sets, validates and clears, on an immutable object too", {
  m <- model(reference, covariates(level = 1))
  name(m) <- "one level"
  expect_identical(name(m), "one level")
  name(m) <- NULL
  expect_null(name(m))

  # A `logmu_function` blocks `$<-`, which is what forced the attribute rule.
  inc <- include(.i$pension > 0)
  expect_error(inc$name <- "x", "immutable")
  name(inc) <- "pensioners"
  group_name(inc) <- "who"
  expect_identical(name(inc), "pensioners")
  expect_identical(group_name(inc), "who")

  expect_error(name(m) <- "", "must be a valid name")
  expect_error(name(m) <- " x", "must be a valid name")
  expect_error(group_name(inc) <- "a  b", "must be a valid name")
  expect_error(name(m) <- c("a", "b"), "must be a valid name")
})

test_that("name() does not partially match names()", {
  # `attr(x, "name")` would return the element names here; the accessor must not.
  x <- structure(list(a = 1, b = 2), class = "something")
  expect_null(name(x))
  expect_identical(names(x), c("a", "b"))
})

# ---- one rule at every entry point ------------------------------------------

test_that("includes() refuses a bad label", {
  std <- ages(65, 95, by = 15)
  expect_error(includes(" young" = std[[1]]), "must be valid names")
  expect_error(includes("a  b" = std), "must be valid names")
  expect_silent(includes(young = std[[1]], "pension amount" = std[[2]]))
})

test_that("names<- and group_names<- on an includes refuse a bad label", {
  std <- ages(65, 95, by = 15)
  expect_error(names(std) <- c("young", "old "), "`names` must be valid names")
  expect_error(names(std) <- c("", "old"), "`names` must be valid names")
  expect_error(group_names(std) <- "age\tband", "`group_names` must be valid names")
  expect_silent(group_names(std) <- "age band")
})

test_that("covariates() refuses a bad name", {
  expect_error(covariates(" male" = .i$group == "a"), "must be valid names")
  expect_error(covariates("male  " = .i$group == "a"), "must be valid names")
})

test_that("categories() refuses a bad group name, written or defaulted", {
  expect_error(categories(.i$group, "a  b" = c("a", "b")), "`group names` must be valid names")
  # A default label is the values joined, so a blank value gives a blank label.
  expect_error(categories(.i$group, "", "a"), "`group names` must be valid names")
  expect_silent(categories(.i$group, other = "", "a"))
})

test_that("batch() refuses a bad element name", {
  expect_error(
    batch("men " = aev(.exp_data, mortality = reference), .exp_data = group_data,
          .settings = basis),
    "`names` must be valid names"
  )
})

test_that("compare_models() refuses a bad list name", {
  level_model <- model(reference, covariates(level = 1))
  group_model <- model(reference, covariates(a = .i$group == "a"))
  expect_error(
    compare_models(group_data, models = list("level " = level_model, group_model),
                   settings = basis),
    "`names\\(models\\)` must be valid names"
  )
})

test_that("diagnose() refuses a bad breakdown name", {
  f <- fit(group_data, model(reference, covariates(level = 1)), settings = basis)
  expect_error(
    diagnose(group_data, f, breakdown = list(" amount" = bands(.i$pension, thresholds = 10000))),
    "`names\\(breakdown\\)` must be valid names"
  )
})

test_that("names<- and group_names<- on an aev refuse a bad label", {
  result <- aev(group_data, mortality = reference, settings = basis,
                breakdown = bands(.i$pension, thresholds = 10000))
  expect_error(names(result) <- c(" low", "high"), "`names` must be valid names")
  expect_error(group_names(result) <- "a  b", "`group_names` must be valid names")
  expect_silent(names(result) <- c("low", "high"))
})

# ---- the labels the package makes for itself --------------------------------

test_that("every generated label satisfies the rule", {
  pass <- function(labels) all(vapply(labels, is_valid_name, logical(1L)))

  expect_true(pass(names(bands(.i$pension, thresholds = c(10000, 20000)))))
  expect_true(pass(group_names(bands(.i$pension, thresholds = c(10000, 20000)))))
  expect_true(pass(names(ages(65, 95, by = 15))))
  expect_true(pass(names(bands(.i$pension, from = 5000, to = 20000, by = 5000))))
  expect_true(pass(names(durations(.i$birth, thresholds = 70))))
  expect_true(pass(group_names(durations(.i$birth, thresholds = 70))))
  expect_true(pass(names(categories(.i$group, "a", c("b", "c")))))
  expect_identical(names(categories(.i$group, "a", c("b", "c"))), c("a", "b+c"))
  expect_true(pass(names(covariates(a = .i$group == "a") * covariates(b = .i$group == "b"))))
})

test_that("a long expression deparses to a single-spaced label", {
  # `deparse()` wraps at sixty bytes and indents the continuation, so the plain
  # join carries runs of spaces past the first line. The label is cut at 57, so
  # it never reaches them; this guards the cut against being lengthened.
  long <- quote(.i$pension > 10000 & .i$pension < 20000 & .i$group == "a" &
                  .i$group == "b" & .i$group == "c" & .i$pension != 15000)
  expect_true(length(deparse(long)) > 1L)
  expect_true(grepl("  ", paste(deparse(long), collapse = " "), fixed = TRUE))
  label <- it_short_deparse(long)
  expect_false(grepl("  ", label, fixed = TRUE))
  expect_true(is_valid_name(label))
  expect_lte(nchar(label), 60L)
})

# ---- the name is what prints --------------------------------------------------

test_that("a named object deparses as its name inside a pronoun expression", {
  flat <- mortality_const(log_mu = -4, name = "flat")
  expect_identical(it_deparse(it_capture(quote(flat), environment())), "flat")
  expect_identical(it_deparse(it_capture(quote(reference), environment())),
                   "<mortality_const>")

  m <- model(flat, covariates(level = 1))
  expect_match(paste(capture.output(print(m)), collapse = "\n"),
               "reference: flat", fixed = TRUE)
})

test_that("the fit header names the test mortality where it has a name", {
  # The variable is NOT called `yardstick`, so that what was typed and the name
  # differ and a header reading the typed text would be caught.
  z_reference <- mortality_const(log_mu = -4, name = "yardstick")
  f <- fit(group_data, model(reference, covariates(level = 1)),
           test_mortality = z_reference, settings = basis)
  expect_identical(f$labels$test_mortality, "yardstick")
  expect_match(paste(capture.output(print(f)), collapse = "\n"),
               "test mortality  yardstick", fixed = TRUE)

  # Unnamed, the label is what was typed, as before.
  g <- fit(group_data, model(reference, covariates(level = 1)),
           test_mortality = reference, settings = basis)
  expect_identical(g$labels$test_mortality, "reference")
})

test_that("the header names a test mortality that was spliced, not kept as a leaf", {
  # A pronoun-expression mortality, which every Gompertz is, is SPLICED into the
  # AST and its name does not survive the splice, so the name has to be read off
  # the value before that happens. Again the variable is called something other
  # than the name it carries.
  steep <- gompertz_mortality(slope_x = 0.12, name = "steeper Gompertz")
  expect_false(identical(it_capture(quote(steep), environment())$kind, "obj"))

  model_1 <- model(reference, covariates(level = 1))
  by_symbol <- fit(group_data, model_1, test_mortality = steep, settings = basis)
  expect_identical(by_symbol$labels$test_mortality, "steeper Gompertz")

  # Redundant parentheses and a call evaluated in place read the same.
  in_parens <- fit(group_data, model_1, test_mortality = (steep), settings = basis)
  expect_identical(in_parens$labels$test_mortality, "steeper Gompertz")
  by_call <- fit(group_data, model_1, settings = basis,
                 test_mortality = gompertz_mortality(slope_x = 0.12, name = "steeper Gompertz"))
  expect_identical(by_call$labels$test_mortality, "steeper Gompertz")

  # Taking the name must not change the tree the engine is given.
  here <- environment()
  for (typed in list(quote(steep), quote((steep)),
                     quote(gompertz_mortality(slope_x = 0.12, name = "steeper Gompertz")))) {
    expect_identical(fit_test_mortality(typed, FALSE, NULL, here)$ast, it_capture(typed, here))
  }

  # Unnamed, the label is still what was typed.
  unnamed <- gompertz_mortality(slope_x = 0.12)
  plain <- fit(group_data, model_1, test_mortality = unnamed, settings = basis)
  expect_identical(plain$labels$test_mortality, "unnamed")

  # The default reads the same whether left out or passed, and so does a
  # comparison, which shares the helper.
  absent <- fit(group_data, model_1, settings = basis)
  passed <- fit(group_data, model_1, test_mortality = default_mortality(), settings = basis)
  expect_identical(absent$labels$test_mortality, "logmu default")
  expect_identical(passed$labels$test_mortality, "logmu default")
  compared <- compare_models(group_data, list(model_1), test_mortality = steep,
                             settings = basis)
  expect_identical(compared$labels$test_mortality, "steeper Gompertz")
})

test_that("mortalities print their name and a line about themselves", {
  flat <- paste(capture.output(print(mortality_const(log_mu = -4, name = "flat"))),
                collapse = "\n")
  expect_match(flat, "<mortality_const: flat> log mu -4, mu 0.018316", fixed = TRUE)
  expect_match(paste(capture.output(print(reference)), collapse = "\n"),
               "<mortality_const> log mu -4", fixed = TRUE)

  gomp <- paste(capture.output(print(gompertz_mortality(name = "gomp"))), collapse = "\n")
  expect_match(gomp, "<mortality_expr: gomp> ", fixed = TRUE)

  tbl <- mortality_table(x0 = datey::durationy(60), t0 = datey::datey(2000),
                         log_mu = matrix(-4, 3, 2), name = "tbl")
  printed <- paste(capture.output(print(tbl)), collapse = "\n")
  expect_match(printed, "<mortality_table: tbl> 3 ages from ", fixed = TRUE)
  expect_match(printed, ", 2 years from ", fixed = TRUE)
})
