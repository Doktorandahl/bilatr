test_that("get_root extracts the two-digit root code", {
  expect_equal(get_root("19"), "19")
  expect_equal(get_root("190"), "19")
  expect_equal(get_root("0211"), "02")
  expect_equal(get_root(19), "19")
})

test_that("assign_quad maps root codes to the four quad classes", {
  expect_equal(assign_quad(c("01", "05")), c(1L, 1L))
  expect_equal(assign_quad(c("06", "09")), c(2L, 2L))
  expect_equal(assign_quad(c("10", "14")), c(3L, 3L))
  expect_equal(assign_quad(c("15", "20")), c(4L, 4L))
  expect_true(is.na(assign_quad("99")))
})

test_that("assign_penta carves out verbal cooperation, protest, and reduce-relations", {
  expect_equal(assign_penta(c("01", "02"), assign_quad(c("01", "02"))), c(0L, 0L))
  expect_equal(assign_penta("14", assign_quad("14")), 4L)
  expect_equal(assign_penta("16", assign_quad("16")), 3L)
  # everything else inherits quad
  expect_equal(assign_penta("19", assign_quad("19")), assign_quad("19"))
})

test_that("assign_eventrootcode2 splits root 04 and folds 09/14/15/18/20", {
  # most roots map to themselves
  expect_equal(assign_eventrootcode2(c("01", "0211", "19")), c("01", "02", "19"))
  # root 04 splits three ways
  expect_equal(
    assign_eventrootcode2(c("04", "040", "041", "044", "045", "046")),
    c("040", "040", "044", "044", "046", "046")
  )
  # folded roots
  expect_equal(assign_eventrootcode2(c("09", "093")), c("10", "10"))
  expect_equal(assign_eventrootcode2(c("14", "150")), c("13", "13"))
  expect_equal(assign_eventrootcode2(c("180", "204")), c("19", "19"))
})

test_that("assign_bilatr_class follows EventRootCode2 with per-code refinements", {
  # regrouped-root defaults
  expect_equal(assign_bilatr_class(c("010", "190", "071", "0862")), c(0L, 10L, 5L, 6L))
  # per-code refinements that override the root default
  expect_equal(assign_bilatr_class("016"), 8L)
  expect_equal(assign_bilatr_class("018"), 4L)
  expect_equal(assign_bilatr_class("019"), 1L)
  expect_equal(assign_bilatr_class(c("0241", "0252")), c(0L, 0L))
  expect_equal(assign_bilatr_class(c("026", "028")), c(3L, 3L))
  expect_equal(assign_bilatr_class("041"), 0L)
  expect_true(is.na(assign_bilatr_class("not-a-code")))
})

test_that("bilatr_class_name labels the 0-10 classes and NAs anything else", {
  expect_equal(bilatr_class_name(0L), "Neutral / low-intensity statement")
  expect_equal(bilatr_class_name(10L), "Assault, fight, or mass violence")
  expect_true(is.na(bilatr_class_name(11L)))
})

test_that("cameo_lookup has no missing recodes or duplicate codes", {
  expect_false(any(is.na(cameo_lookup$QuadClass)))
  expect_false(any(is.na(cameo_lookup$PentaClass)))
  expect_false(any(is.na(cameo_lookup$EventRootCode2)))
  expect_false(any(is.na(cameo_lookup$BilatrClass)))
  expect_false(any(is.na(cameo_lookup$BilatrClassName)))
  expect_true(all(cameo_lookup$BilatrClass %in% 0:10))
  expect_false(any(duplicated(cameo_lookup$CAMEOEVENTCODE)))
})

test_that("cameo_lookup anchors match the project's established reference categories", {
  coop <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "044", ]
  hostile <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "19", ]
  expect_equal(coop$QuadClass, 1L)
  expect_equal(hostile$QuadClass, 4L)
  expect_equal(hostile$PentaClass, 4L)
})

test_that("recode_cameo attaches QuadClass/PentaClass by joining on a named code column", {
  events <- tibble::tibble(EventCode = c("044", "19", "010"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_equal(out$QuadClass, c(1L, 4L, 1L))
  expect_true(all(c("PentaClass", "PentaClass_modified", "GoldsteinScore") %in% names(out)))
})

test_that("recode_cameo also attaches EventRootCode2 and BilatrClass/Name", {
  events <- tibble::tibble(EventCode = c("044", "190", "180", "093"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_true(all(c("EventRootCode2", "BilatrClass", "BilatrClassName") %in% names(out)))
  expect_equal(out$EventRootCode2, c("044", "19", "19", "10"))
  expect_equal(out$BilatrClass, c(2L, 10L, 10L, 7L))
  expect_equal(out$BilatrClassName[2], "Assault, fight, or mass violence")
})

test_that("recode_cameo leaves unmatched codes as NA rather than erroring", {
  events <- tibble::tibble(EventCode = c("044", "not-a-code"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_equal(out$QuadClass, c(1L, NA_integer_))
})
