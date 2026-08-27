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

test_that("cameo_lookup has no missing recodes or duplicate codes", {
  expect_false(any(is.na(cameo_lookup$QuadClass)))
  expect_false(any(is.na(cameo_lookup$PentaClass)))
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

test_that("recode_cameo leaves unmatched codes as NA rather than erroring", {
  events <- tibble::tibble(EventCode = c("044", "not-a-code"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_equal(out$QuadClass, c(1L, NA_integer_))
})
