# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

codes_of <- function(x) as.vector(unclass(x))

skip_if_not_installed("naxp")

postcode <- test_naxp("\\9\\9\\9\\9 \\A\\A")
types <- categorical(3, name = "types")
clusters <- categorical(labels = c("low", "high", "unknown"), name = "clusters")

test_that("a naxp-keyed map looks up codes, and code 0 and unlisted codes get undefined", {
  oac <- category_map(
    keys = test_naxp_encoded(c(30L, 10L, 20L), postcode),
    values = categorical_encoded(c(3, 1, 2), types),
    name = "oac"
  )
  expect_true(is_category_map(oac))
  expect_s3_class(oac, c("category_map", "function"))
  expect_identical(name(oac), "oac")
  # Keys are stored sorted, with their values alongside.
  expect_identical(attr(oac, "keys"), c(10L, 20L, 30L))
  expect_identical(attr(oac, "value_type"), "categorical")

  # A categorical's default undefined is code 0, which is invalid.
  result <- oac(test_naxp_encoded(c(20L, 99L, 0L, 30L), postcode))
  expect_identical(codes_of(result), c(2L, 0L, 0L, 3L))
  expect_identical(is.na(result), c(FALSE, TRUE, TRUE, FALSE))
  expect_identical(domain(result), types)
})

test_that("undefined can be a declared category, to model unknowns", {
  to_cluster <- category_map(
    keys = test_naxp_encoded(10L, postcode),
    values = encode("low", clusters),
    undefined = encode("unknown", clusters)
  )
  expect_identical(decode(to_cluster(test_naxp_encoded(c(10L, 0L, 11L), postcode))),
                   c("low", "unknown", "unknown"))
})

test_that("code 0 is never a key", {
  expect_error(category_map(test_naxp_encoded(c(0L, 10L), postcode), c(1, 2)),
               "Code 0 is invalid and cannot be a key")
  expect_error(category_map(categorical_encoded(c(0, 1), types), c(1, 2)), "cannot be a key")
})

test_that("each value type has its own default undefined, and logical and integer have none", {
  keys <- test_naxp_encoded(5L, postcode)
  look_up <- test_naxp_encoded(c(5L, 6L, 0L), postcode)
  expect_identical(category_map(keys, -2.5)(look_up), c(-2.5, NaN, NaN))
  expect_error(category_map(keys, TRUE), "A logical map needs `undefined`")
  expect_error(category_map(keys, 7L), "An integer map needs `undefined`")
  expect_identical(category_map(keys, TRUE, undefined = FALSE)(look_up), c(TRUE, FALSE, FALSE))
  expect_identical(category_map(keys, 7L, undefined = -1L)(look_up), c(7L, -1L, -1L))
  # A NaN value is kept for its own key.
  expect_identical(category_map(keys, NaN, undefined = 0)(look_up), c(NaN, 0, 0))

  district <- test_naxp("\\9\\9")
  to_district <- category_map(keys, test_naxp_encoded(42L, district))(look_up)
  expect_identical(codes_of(to_district), c(42L, 0L, 0L))
  expect_identical(attr(to_district, "naxp"), district)
})

test_that("a categorical source reaches undefined too, through code 0", {
  # Every category is listed, but code 0 can still arrive.
  every <- categorical_encoded(c(1, 2, 3), types)
  expect_error(category_map(every, c(TRUE, FALSE, TRUE)), "needs `undefined`")
  weights <- category_map(every, c(1.1, 1.2, 1.3))
  expect_identical(weights(categorical_encoded(c(3, 0), types)), c(1.3, NaN))
})

test_that("undefined must match the type of the values", {
  keys <- test_naxp_encoded(5L, postcode)
  expect_error(category_map(keys, 2.5, undefined = 1L), "same type")
  expect_error(category_map(keys, 2.5, undefined = c(1, 2)), "must have 1 value")
  expect_error(category_map(keys, categorical_encoded(1, types), undefined = encode("low", clusters)),
               "same type")
  expect_error(category_map(keys, data.frame(x = 1), undefined = data.frame(y = 0)), "same type")
  expect_error(category_map(keys, TRUE, undefined = NA), "may not be NA")
  expect_no_error(category_map(keys, data.frame(x = 1), undefined = data.frame(x = -1)))
})

test_that("a naxp-keyed map refuses other naxps and other keys", {
  oac <- category_map(test_naxp_encoded(1L, postcode), 2.5)
  expect_error(oac(test_naxp_encoded(1L, test_naxp("\\A"))), "naxp_encoded` codes of the naxp")
  expect_error(oac(1L), "naxp_encoded` codes")
  expect_error(category_map(1L, 2.5), "must be a `naxp_encoded` or a `categorical_encoded`")
  expect_error(category_map(test_naxp_encoded(c(1L, 1L), postcode), c(1, 2)), "repeated: \"code 1\"")
  expect_error(category_map(test_naxp_encoded(1L, postcode), c(1, 2)), "must have 1 value")
})

test_that("logical and integer values may not be NA", {
  keys <- test_naxp_encoded(5L, postcode)
  expect_error(category_map(keys, NA, undefined = FALSE), "no portable missing value")
  expect_error(category_map(keys, NA_integer_, undefined = 0L), "no portable missing value")
  expect_no_error(category_map(keys, NA_real_))
  expect_error(category_map(keys, factor("a")), "encode\\(\\)")
  expect_error(category_map(keys, list(1)), "must be logical")
  expect_error(category_map(keys, "a"), "must be logical")
})

test_that("a struct is a data frame of fields, with undefined field by field", {
  coords <- category_map(
    keys = test_naxp_encoded(c(2L, 1L), postcode),
    values = data.frame(easting = c(200, 100), northing = c(20, 10),
                        type = categorical_encoded(c(2, 1), types))
  )
  expect_identical(attr(coords, "value_type"), "struct")
  result <- coords(test_naxp_encoded(c(1L, 3L, 2L), postcode))
  expect_s3_class(result, "data.frame")
  expect_identical(names(result), c("easting", "northing", "type"))
  expect_identical(result$easting, c(100, NaN, 200))
  expect_identical(result$northing, c(10, NaN, 20))
  expect_identical(codes_of(result$type), c(1L, 0L, 2L))

  given <- category_map(
    keys = test_naxp_encoded(1L, postcode),
    values = data.frame(easting = 100, type = categorical_encoded(1, types)),
    undefined = data.frame(easting = -1, type = categorical_encoded(3, types))
  )
  result <- given(test_naxp_encoded(c(1L, 0L), postcode))
  expect_identical(result$easting, c(100, -1))
  expect_identical(codes_of(result$type), c(1L, 3L))
})

test_that("a struct with a logical or integer field needs undefined", {
  keys <- test_naxp_encoded(1L, postcode)
  expect_error(category_map(keys, data.frame(weight = 1, flagged = TRUE, count = 2L)),
               "logical or integer field \\(`flagged`, `count`\\) needs `undefined`")
  expect_no_error(category_map(keys, data.frame(weight = 1, flagged = TRUE),
                               undefined = data.frame(weight = NaN, flagged = FALSE)))
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

test_that("a categorical-keyed map stores a value for code 0 and every category", {
  to_cluster <- category_map(
    keys = categorical_encoded(c(3, 1), types),
    values = encode(c("high", "low"), clusters),
    undefined = encode("unknown", clusters)
  )
  # Codes 0, 1, 2 and 3; 0 and 2 are not listed and so take undefined.
  expect_identical(codes_of(attr(to_cluster, "values")), c(3L, 1L, 3L, 2L))
  expect_null(attr(to_cluster, "undefined"))
  expect_identical(decode(to_cluster(categorical_encoded(c(1, 2, 3, 0), types))),
                   c("low", "unknown", "high", "unknown"))
})

test_that("a categorical-keyed map accepts text and factors in its source's labels", {
  status <- categorical(labels = c("act", "def", "pen"))
  weights <- category_map(encode(c("act", "def", "pen"), status), c(1, 0.5, 2))
  expect_identical(weights(c("pen", "act", NA)), c(2, 1, NaN))
  expect_identical(weights(factor("def")), 0.5)
  expect_error(weights("dep"), "Not a category")
  expect_error(weights(categorical_encoded(1, categorical(1))), "looks up categories in the domain")
})

test_that("a map between ordered domains must keep the order", {
  bands <- categorical(labels = c("60", "65", "70", "75"), sorted_codes = 1:4)
  groups <- categorical(labels = c("young", "old"), sorted_codes = c(1, 2))
  keys <- encode(c("60", "65", "70", "75"), bands)
  expect_no_error(category_map(keys, encode(c("young", "young", "old", "old"), groups)))
  expect_error(category_map(keys, encode(c("young", "old", "young", "old"), groups)),
               "keep the order")
  # Invalid values take no part in the check, wherever they fall.
  expect_no_error(category_map(keys, encode(c("young", NA, "old", "old"), groups)))
  expect_no_error(category_map(encode(c("60", "75"), bands), encode(c("young", "old"), groups)))
  expect_error(category_map(keys, encode(c("old", NA, "young", NA), groups)), "keep the order")
  # The source's order is followed, not its codes.
  reversed <- categorical(labels = c("60", "65", "70", "75"), sorted_codes = 4:1)
  expect_no_error(category_map(encode(c("60", "65", "70", "75"), reversed),
                               encode(c("old", "old", "young", "young"), groups)))
  # An unordered side means no check.
  unordered_groups <- categorical(labels = c("young", "old"))
  expect_no_error(category_map(keys, encode(c("old", "young", "old", "young"), unordered_groups)))
})

test_that("applying a map to a map is an error", {
  oac <- category_map(test_naxp_encoded(10L, postcode), categorical_encoded(1, types))
  to_cluster <- category_map(categorical_encoded(1, types), encode("high", clusters))
  expect_error(to_cluster(oac), "Use `compose_maps\\(\\)`")
})

test_that("compose_maps(f, g) is one map that does what f(g(x)) does", {
  oac <- category_map(test_naxp_encoded(c(10L, 20L, 30L), postcode), categorical_encoded(c(1, 2, 3), types))
  to_cluster <- category_map(categorical_encoded(c(1, 2, 3), types),
                             encode(c("low", "low", "high"), clusters))
  composed <- compose_maps(to_cluster, oac, name = "oac_clusters")
  expect_true(is_category_map(composed))
  expect_identical(name(composed), "oac_clusters")
  look_up <- test_naxp_encoded(c(30L, 10L, 99L, 0L), postcode)
  expect_identical(composed(look_up), to_cluster(oac(look_up)))

  # Three deep, ending in a double.
  weights <- category_map(encode(c("low", "high"), clusters), c(0.9, 1.1), undefined = 1)
  expect_identical(compose_maps(weights, composed)(look_up), weights(to_cluster(oac(look_up))))
  expect_error(compose_maps(oac, to_cluster), "values of `g` to be in the source of `f`")
  expect_error(compose_maps(oac, 1), "must be a `category_map`")
})

test_that("composition carries the outer undefined through to invalid and unlisted keys", {
  oac <- category_map(test_naxp_encoded(10L, postcode), categorical_encoded(1, types))
  # The inner map sends unlisted postcodes to code 0; the outer map sends code
  # 0 to a declared category, so the composed map must do the same.
  to_cluster <- category_map(categorical_encoded(1, types), encode("low", clusters),
                             undefined = encode("unknown", clusters))
  composed <- compose_maps(to_cluster, oac)
  look_up <- test_naxp_encoded(c(10L, 11L, 0L), postcode)
  expect_identical(decode(composed(look_up)), c("low", "unknown", "unknown"))
  expect_identical(composed(look_up), to_cluster(oac(look_up)))

  # A categorical-keyed inner map composes the same way.
  dense_inner <- category_map(categorical_encoded(1, types), categorical_encoded(2, types))
  dense_composed <- compose_maps(to_cluster, dense_inner)
  every_code <- categorical_encoded(c(0, 1, 2, 3), types)
  expect_identical(dense_composed(every_code), to_cluster(dense_inner(every_code)))
})

test_that("a map prints its source and value type", {
  oac <- category_map(test_naxp_encoded(c(10L, 20L), postcode), categorical_encoded(c(1, 2), types), name = "oac")
  expect_output(print(oac), "<category_map: oac> 2 keys from naxp \\\\9\\\\9\\\\9\\\\9 \\\\A\\\\A to categorical \\(types\\)")
  coords <- category_map(categorical_encoded(1, types), data.frame(x = 1, y = 2))
  expect_output(print(coords), "<category_map> 3 categories from types to struct \\(x, y\\)")
})
