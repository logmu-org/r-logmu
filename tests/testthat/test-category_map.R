# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

codes_of <- function(x) as.vector(unclass(x))

skip_if_not_installed("naxp")

postcode <- test_naxp("\\9\\9\\9\\9 \\A\\A")
types <- category_domain(3, allows_undefined = TRUE, name = "types")
clusters <- category_domain(labels = c("low", "high"), undefined_label = "none", name = "clusters")

test_that("a naxp-keyed map looks up codes, and code 0 and unlisted codes get undefined", {
  oac <- category_map(
    keys = test_naxp_encoded(c(30L, 10L, 20L), postcode),
    values = as_categorical(c(3, 1, 2), types),
    name = "oac"
  )
  expect_true(is_category_map(oac))
  expect_s3_class(oac, c("category_map", "function"))
  expect_identical(name(oac), "oac")
  # Keys are stored sorted, with their values alongside.
  expect_identical(attr(oac, "keys"), c(10L, 20L, 30L))
  expect_identical(attr(oac, "value_type"), "categorical")

  result <- oac(test_naxp_encoded(c(20L, 99L, 0L, 30L), postcode))
  expect_identical(codes_of(result), c(2L, 0L, 0L, 3L))
  expect_identical(domain(result), types)
})

test_that("undefined can be any value of the type, including a real category", {
  rural <- category_map(
    keys = test_naxp_encoded(10L, postcode),
    values = as_categorical(1, types),
    undefined = as_categorical(3, types)
  )
  expect_identical(codes_of(rural(test_naxp_encoded(c(10L, 0L, 11L), postcode))), c(1L, 3L, 3L))
})

test_that("code 0 is never a key", {
  expect_error(category_map(test_naxp_encoded(c(0L, 10L), postcode), c(1, 2)),
               "Code 0 cannot be a key")
  expect_error(category_map(as_categorical(c(0, 1), types), c(1, 2)), "Code 0 cannot be a key")
})

test_that("a zero undefined needs a value domain that allows it, where it is reached", {
  strict <- category_domain(labels = c("low", "high"))
  expect_error(category_map(test_naxp_encoded(1L, postcode), as_categorical("low", strict)),
               "give `undefined` explicitly")
  expect_no_error(category_map(test_naxp_encoded(1L, postcode), as_categorical("low", strict),
                               undefined = as_categorical("high", strict)))
  # A categorical source with every code listed, and no code 0, never reaches it.
  closed <- category_domain(2)
  expect_no_error(category_map(as_categorical(c(1, 2), closed), as_categorical(c("low", "high"), strict)))
  expect_error(category_map(as_categorical(1, closed), as_categorical("low", strict)),
               "give `undefined` explicitly")
})

test_that("undefined must match the type of the values", {
  keys <- test_naxp_encoded(5L, postcode)
  expect_error(category_map(keys, 2.5, undefined = 1L), "same type")
  expect_error(category_map(keys, 2.5, undefined = c(1, 2)), "must have 1 value")
  expect_error(category_map(keys, as_categorical(1, types), undefined = as_categorical("low", clusters)),
               "same type")
  expect_error(category_map(keys, data.frame(x = 1), undefined = data.frame(y = 0)), "same type")
  expect_no_error(category_map(keys, data.frame(x = 1), undefined = data.frame(x = -1)))
})

test_that("a naxp-keyed map refuses other naxps and other keys", {
  oac <- category_map(test_naxp_encoded(1L, postcode), 2.5)
  expect_error(oac(test_naxp_encoded(1L, test_naxp("\\A"))), "naxp_encoded` codes of the naxp")
  expect_error(oac(1L), "naxp_encoded` codes")
  expect_error(category_map(1L, 2.5), "must be a `naxp_encoded` vector or a `categorical`")
  expect_error(category_map(test_naxp_encoded(c(1L, 1L), postcode), c(1, 2)), "repeated: \"code 1\"")
  expect_error(category_map(test_naxp_encoded(1L, postcode), c(1, 2)), "must have 1 value")
})

test_that("every value type gives all-zero bits by default", {
  keys <- test_naxp_encoded(5L, postcode)
  look_up <- test_naxp_encoded(c(5L, 6L), postcode)
  expect_identical(category_map(keys, TRUE)(look_up), c(TRUE, FALSE))
  expect_identical(category_map(keys, 7L)(look_up), c(7L, 0L))
  expect_identical(category_map(keys, -2.5)(look_up), c(-2.5, 0))
  # A NaN value is kept for its own key: only the unlisted key gets undefined.
  expect_identical(category_map(keys, NaN)(look_up), c(NaN, 0))

  district <- test_naxp("\\9\\9")
  to_district <- category_map(keys, test_naxp_encoded(42L, district))(look_up)
  expect_identical(codes_of(to_district), c(42L, 0L))
  expect_identical(attr(to_district, "naxp"), district)
})

test_that("logical and integer values may not be NA", {
  keys <- test_naxp_encoded(5L, postcode)
  expect_error(category_map(keys, NA), "no portable missing value")
  expect_error(category_map(keys, NA_integer_), "no portable missing value")
  expect_no_error(category_map(keys, NA_real_))
  expect_error(category_map(keys, factor("a")), "as_categorical")
  expect_error(category_map(keys, list(1)), "must be logical")
  expect_error(category_map(keys, "a"), "must be logical")
})

test_that("a struct is a data frame of fields, with undefined field by field", {
  coords <- category_map(
    keys = test_naxp_encoded(c(2L, 1L), postcode),
    values = data.frame(easting = c(200, 100), northing = c(20, 10), type = as_categorical(c(2, 1), types)),
    undefined = data.frame(easting = -1, northing = 0, type = as_categorical(3, types))
  )
  expect_identical(attr(coords, "value_type"), "struct")
  result <- coords(test_naxp_encoded(c(1L, 3L, 2L), postcode))
  expect_s3_class(result, "data.frame")
  expect_identical(names(result), c("easting", "northing", "type"))
  expect_identical(result$easting, c(100, -1, 200))
  expect_identical(result$northing, c(10, 0, 20))
  expect_identical(codes_of(result$type), c(1L, 3L, 2L))
})

test_that("struct field names are checked, case-sensitively", {
  keys <- test_naxp_encoded(1L, postcode)
  expect_error(category_map(keys, data.frame(`grid ref` = 1, check.names = FALSE)), "start with a letter")
  expect_error(category_map(keys, data.frame(`1st` = 1, check.names = FALSE)), "start with a letter")
  expect_error(category_map(keys, data.frame(a = 1, a = 2, check.names = FALSE)), "repeated: \"a\"")
  expect_no_error(category_map(keys, data.frame(a = 1, A = 2)))
  nested <- data.frame(a = 1)
  nested$inner <- data.frame(b = 2)
  expect_error(category_map(keys, nested), "may not nest")
})

test_that("a categorical-keyed map stores a value for every code, undefined for the rest", {
  to_cluster <- category_map(
    keys = as_categorical(c(3, 1), types),
    values = as_categorical(c("high", "low"), clusters),
    undefined = as_categorical("low", clusters)
  )
  # Codes 0, 1, 2 and 3; 0 and 2 are not listed and so take undefined.
  expect_identical(codes_of(attr(to_cluster, "values")), c(1L, 1L, 1L, 2L))
  expect_null(attr(to_cluster, "undefined"))
  expect_identical(format(to_cluster(as_categorical(c(1, 2, 3, 0), types))),
                   c("low", "low", "high", "low"))
})

test_that("a categorical-keyed map accepts text and factors in its source's labels", {
  status <- category_domain(labels = c("act", "def", "pen"))
  weights <- category_map(as_categorical(c("act", "def", "pen"), status), c(1, 0.5, 2))
  expect_identical(weights(c("pen", "act")), c(2, 1))
  expect_identical(weights(factor("def")), 0.5)
  expect_error(weights("dep"), "Not a category")
  expect_error(weights(as_categorical(1, category_domain(1))), "looks up categories in the domain")
})

test_that("a map between ordered domains must keep the order", {
  bands <- category_domain(labels = c("60", "65", "70", "75"), sorted_codes = 1:4)
  groups <- category_domain(labels = c("young", "old"), undefined_label = "unknown",
                            sorted_codes = c(1, 2, 0))
  keys <- as_categorical(c("60", "65", "70", "75"), bands)
  expect_no_error(category_map(keys, as_categorical(c("young", "young", "old", "old"), groups)))
  expect_error(category_map(keys, as_categorical(c("young", "old", "young", "old"), groups)),
               "keep the order")
  # Undefined has a place in the order like any other value: last, here.
  expect_no_error(category_map(keys, as_categorical(c("young", "old", "old", "unknown"), groups)))
  expect_error(category_map(keys, as_categorical(c("young", "unknown", "old", "old"), groups)),
               "keep the order")
  # The source's order is followed, not its codes.
  reversed <- category_domain(labels = c("60", "65", "70", "75"), sorted_codes = 4:1)
  expect_no_error(category_map(as_categorical(c("60", "65", "70", "75"), reversed),
                               as_categorical(c("old", "old", "young", "young"), groups)))
  # An unordered side means no check.
  unordered_groups <- category_domain(labels = c("young", "old"))
  expect_no_error(category_map(keys, as_categorical(c("old", "young", "old", "young"), unordered_groups)))
})

test_that("applying a map to a map is an error", {
  oac <- category_map(test_naxp_encoded(10L, postcode), as_categorical(1, types))
  to_cluster <- category_map(as_categorical(1, types), as_categorical("high", clusters))
  expect_error(to_cluster(oac), "Use `compose_maps\\(\\)`")
})

test_that("compose_maps(f, g) is one map that does what f(g(x)) does", {
  oac <- category_map(test_naxp_encoded(c(10L, 20L, 30L), postcode), as_categorical(c(1, 2, 3), types))
  to_cluster <- category_map(as_categorical(c(1, 2, 3), types),
                             as_categorical(c("low", "low", "high"), clusters))
  composed <- compose_maps(to_cluster, oac, name = "oac clusters")
  expect_true(is_category_map(composed))
  expect_identical(name(composed), "oac clusters")
  look_up <- test_naxp_encoded(c(30L, 10L, 99L, 0L), postcode)
  expect_identical(composed(look_up), to_cluster(oac(look_up)))

  # Three deep, ending in a double.
  weights <- category_map(as_categorical(c("low", "high"), clusters), c(0.9, 1.1), undefined = 1)
  expect_identical(compose_maps(weights, composed)(look_up), weights(to_cluster(oac(look_up))))
  expect_error(compose_maps(oac, to_cluster), "values of `g` to be in the source of `f`")
  expect_error(compose_maps(oac, 1), "must be a `category_map`")
})

test_that("composition carries a non-zero undefined through to unlisted keys", {
  oac <- category_map(test_naxp_encoded(10L, postcode), as_categorical(1, types))
  # The outer map sends code 0 to a real cluster, so the composed map's
  # unlisted postcodes must get that cluster too.
  undefined_to_high <- category_map(as_categorical(1, types), as_categorical("low", clusters),
                                    undefined = as_categorical("high", clusters))
  composed <- compose_maps(undefined_to_high, oac)
  look_up <- test_naxp_encoded(c(10L, 11L, 0L), postcode)
  expect_identical(format(composed(look_up)), c("low", "high", "high"))
  expect_identical(composed(look_up), undefined_to_high(oac(look_up)))

  # A categorical-keyed inner map composes the same way.
  dense_inner <- category_map(as_categorical(1, types), as_categorical(2, types))
  dense_composed <- compose_maps(undefined_to_high, dense_inner)
  every_code <- as_categorical(c(0, 1, 2, 3), types)
  expect_identical(dense_composed(every_code), undefined_to_high(dense_inner(every_code)))
})

test_that("a map prints its source and value type", {
  oac <- category_map(test_naxp_encoded(c(10L, 20L), postcode), as_categorical(c(1, 2), types), name = "oac")
  expect_output(print(oac), "<category_map: oac> 2 keys from naxp \\\\9\\\\9\\\\9\\\\9 \\\\A\\\\A to categorical \\(types\\)")
  coords <- category_map(as_categorical(1, types), data.frame(x = 1, y = 2))
  expect_output(print(coords), "<category_map> 4 codes from types to struct \\(x, y\\)")
})
