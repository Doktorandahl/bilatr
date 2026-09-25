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

test_that("quadclass_name labels the 1-4 classes and NAs anything else", {
  expect_equal(quadclass_name(1L), "Verbal cooperation")
  expect_equal(quadclass_name(2L), "Material cooperation")
  expect_equal(quadclass_name(3L), "Verbal conflict")
  expect_equal(quadclass_name(4L), "Material conflict")
  expect_true(is.na(quadclass_name(5L)))
})

test_that("quadclass_eventcodes lists the root codes feeding each class and NAs anything else", {
  expect_equal(quadclass_eventcodes(1L), "01, 02, 03, 04, 05")
  expect_equal(quadclass_eventcodes(2L), "06, 07, 08, 09")
  expect_equal(quadclass_eventcodes(3L), "10, 11, 12, 13, 14")
  expect_equal(quadclass_eventcodes(4L), "15, 16, 17, 18, 19, 20")
  expect_true(is.na(quadclass_eventcodes(5L)))
})

test_that("assign_penta carves out verbal cooperation, protest, and reduce-relations", {
  expect_equal(assign_penta(c("01", "02"), assign_quad(c("01", "02"))), c(0L, 0L))
  expect_equal(assign_penta("14", assign_quad("14")), 4L)
  expect_equal(assign_penta("16", assign_quad("16")), 3L)
  # everything else inherits quad
  expect_equal(assign_penta("19", assign_quad("19")), assign_quad("19"))
})

test_that("pentaclass_name labels the 0-4 classes and NAs anything else", {
  expect_equal(pentaclass_name(0L), "Make statement")
  expect_equal(pentaclass_name(1L), "Verbal cooperation")
  expect_equal(pentaclass_name(2L), "Material cooperation")
  expect_equal(pentaclass_name(3L), "Verbal conflict")
  expect_equal(pentaclass_name(4L), "Material conflict")
  expect_true(is.na(pentaclass_name(5L)))
})

test_that("pentaclass_eventcodes lists the root codes feeding each class and NAs anything else", {
  expect_equal(pentaclass_eventcodes(0L), "01, 02")
  expect_equal(pentaclass_eventcodes(1L), "03, 04, 05")
  expect_equal(pentaclass_eventcodes(2L), "06, 07, 08, 09")
  expect_equal(pentaclass_eventcodes(3L), "10, 11, 12, 13, 16")
  expect_equal(pentaclass_eventcodes(4L), "14, 15, 17, 18, 19, 20")
  expect_true(is.na(pentaclass_eventcodes(5L)))
})

test_that("cameo_lookup has no missing recodes or duplicate codes", {
  expect_false(any(is.na(cameo_lookup$QuadClass)))
  expect_false(any(is.na(cameo_lookup$QuadClassName)))
  expect_false(any(is.na(cameo_lookup$QuadClassEventCodes)))
  expect_false(any(is.na(cameo_lookup$PentaClass)))
  expect_false(any(is.na(cameo_lookup$PentaClassName)))
  expect_false(any(is.na(cameo_lookup$PentaClassEventCodes)))
  expect_false(any(is.na(cameo_lookup$PentaClass_modified)))
  expect_true(all(cameo_lookup$ModifiedRootCode %in% c(1:18, NA)))
  expect_false(any(duplicated(cameo_lookup$CAMEOEVENTCODE)))
})

test_that("cameo_lookup anchors match the project's established reference categories", {
  coop <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "044", ]
  hostile <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "19", ]
  expect_equal(coop$QuadClass, 1L)
  expect_equal(hostile$QuadClass, 4L)
  expect_equal(hostile$PentaClass, 4L)
})

test_that("cameo_lookup blanks ModifiedRootCode only on non-homogeneous top-level roots", {
  # roots "01" and "04" split their own sub-codes across several
  # ModifiedRootCode classes, so the top-level (two-digit) row is blank
  root01 <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "01", ]
  root04 <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "04", ]
  expect_true(is.na(root01$ModifiedRootCode))
  expect_true(is.na(root01$ModifiedRootCodeName))
  expect_true(is.na(root01$ModifiedRootCodeEventCodes))
  expect_true(is.na(root04$ModifiedRootCode))
  expect_true(is.na(root04$ModifiedRootCodeName))
  expect_true(is.na(root04$ModifiedRootCodeEventCodes))

  # but their sub-codes still get their individual classes
  sub <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE %in% c("010", "016", "018", "040", "041", "046"), ]
  expect_false(any(is.na(sub$ModifiedRootCode)))
  expect_false(any(is.na(sub$ModifiedRootCodeName)))
  expect_false(any(is.na(sub$ModifiedRootCodeEventCodes)))

  # a homogeneous root (e.g. "02") keeps its top-level value
  root02 <- cameo_lookup[cameo_lookup$CAMEOEVENTCODE == "02", ]
  expect_equal(root02$ModifiedRootCode, 2L)
  expect_equal(root02$ModifiedRootCodeName, "Appeal for action")

  # root 01's eventcodes note the 016/018/019 exception
  expect_equal(
    cameo_lookup$ModifiedRootCodeEventCodes[cameo_lookup$CAMEOEVENTCODE == "010"],
    "01 except 016, 018, 019"
  )

  # root 04's own bare code is dropped from "Consult, unspecified"'s
  # eventcodes, since root 04 itself is blanked at the top level
  expect_equal(
    cameo_lookup$ModifiedRootCodeEventCodes[cameo_lookup$CAMEOEVENTCODE == "040"],
    "040"
  )
})

test_that("recode_cameo attaches QuadClass/PentaClass by joining on a named code column", {
  events <- tibble::tibble(EventCode = c("044", "19", "010"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_equal(out$QuadClass, c(1L, 4L, 1L))
  expect_true(all(c("PentaClass", "PentaClass_modified", "GoldsteinScore") %in% names(out)))
})

test_that("recode_cameo also attaches QuadClassName/PentaClassName/EventCodes columns", {
  events <- tibble::tibble(EventCode = c("044", "190"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_true(all(c(
    "QuadClassName", "QuadClassEventCodes", "PentaClassName", "PentaClassEventCodes"
  ) %in% names(out)))
  expect_equal(out$QuadClassName, c("Verbal cooperation", "Material conflict"))
  expect_equal(out$PentaClassName, c("Verbal cooperation", "Material conflict"))
})

test_that("recode_cameo does not attach the retired classification schemes", {
  events <- tibble::tibble(EventCode = c("044", "190"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_false(any(c(
    "EventRootCode2", "EventRootCode2Name",
    "EventRootCode3", "EventRootCode3Name", "EventRootCode3RootCodes",
    "EventRootCode4", "EventRootCode4Name", "EventRootCode4RootCodes",
    "BilatrClass", "BilatrClassName", "BilatrClass2", "BilatrClass2Name"
  ) %in% names(out)))
})

test_that("assign_modified_root_code merges investigate and demand and otherwise follows EventRootCode3", {
  expect_equal(assign_modified_root_code(c("09", "091", "10", "100")), c(11L, 11L, 11L, 11L))
  expect_equal(assign_modified_root_code(c("01", "08")), c(1L, 10L))
  expect_equal(assign_modified_root_code(c("11", "12", "016")), c(12L, 13L, 13L))
  expect_equal(assign_modified_root_code(c("13", "15", "14")), c(14L, 14L, 15L))
  expect_equal(assign_modified_root_code(c("16", "17", "18", "19", "20")), c(16L, 17L, 18L, 18L, 18L))
  expect_true(is.na(assign_modified_root_code("99")))

  # Covers exactly classes 1-18 over every real CAMEO code.
  codes <- cameo_lookup$CAMEOEVENTCODE
  modroot <- assign_modified_root_code(codes)
  expect_equal(sort(unique(modroot)), 1:18)
})

test_that("modified_root_code_name labels every assign_modified_root_code() output value and NAs anything else", {
  expect_equal(modified_root_code_name(11L), "Investigate or demand")
  expect_equal(modified_root_code_name(18L), "Assault, fight, or mass violence")
  expect_true(is.na(modified_root_code_name(19L)))
  expect_true(is.na(modified_root_code_name(0L)))

  produced <- unique(assign_modified_root_code(cameo_lookup$CAMEOEVENTCODE))
  expect_false(any(is.na(modified_root_code_name(produced))))
})

test_that("modified_root_code_eventcodes lists the codes feeding each class, notes the 01 exception, and NAs anything else", {
  expect_equal(modified_root_code_eventcodes(1L), "01 except 016, 018, 019")
  expect_equal(modified_root_code_eventcodes(4L), "040")
  expect_equal(modified_root_code_eventcodes(5L), "041, 042, 043, 044")
  expect_equal(modified_root_code_eventcodes(11L), "09, 10")
  expect_equal(modified_root_code_eventcodes(13L), "12, 016")
  expect_equal(modified_root_code_eventcodes(18L), "18, 19, 20")
  expect_true(is.na(modified_root_code_eventcodes(19L)))

  # every class 1-18 has a non-NA eventcodes string
  expect_false(any(is.na(modified_root_code_eventcodes(1:18))))
})

test_that("recode_cameo also attaches ModifiedRootCode and its name/eventcodes columns", {
  events <- tibble::tibble(EventCode = c("09", "10", "11", "20"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_true(all(c(
    "ModifiedRootCode", "ModifiedRootCodeName", "ModifiedRootCodeEventCodes"
  ) %in% names(out)))
  expect_equal(out$ModifiedRootCode, c(11L, 11L, 12L, 18L))
  expect_equal(out$ModifiedRootCodeName, c(
    "Investigate or demand", "Investigate or demand", "Disapprove",
    "Assault, fight, or mass violence"
  ))
  expect_equal(out$ModifiedRootCodeEventCodes, c("09, 10", "09, 10", "11", "18, 19, 20"))
})

test_that("recode_cameo leaves unmatched codes as NA rather than erroring", {
  events <- tibble::tibble(EventCode = c("044", "not-a-code"))
  out <- recode_cameo(events, code_col = "EventCode")
  expect_equal(out$QuadClass, c(1L, NA_integer_))
})

test_that("recode_cameo skips recode columns already present in the data, with a warning", {
  events <- tibble::tibble(
    EventCode = c("044", "19"),
    QuadClass = c("keep-me", "and-me"),
    ModifiedRootCode = c(99L, 99L)
  )
  expect_warning(
    out <- recode_cameo(events, code_col = "EventCode"),
    "already present"
  )
  # pre-existing columns are untouched (no .x/.y, values unchanged)
  expect_equal(out$QuadClass, c("keep-me", "and-me"))
  expect_equal(out$ModifiedRootCode, c(99L, 99L))
  expect_false(any(grepl("\\.(x|y)$", names(out))))
  # non-clashing columns still attached
  expect_true(all(c("PentaClass", "PentaClassName", "ModifiedRootCodeName") %in% names(out)))
  expect_equal(out$PentaClass, c(1L, 4L))
})

test_that("recode_cameo does not warn when no recode columns pre-exist", {
  events <- tibble::tibble(EventCode = c("044", "19"))
  expect_no_warning(recode_cameo(events, code_col = "EventCode"))
})

test_that("assign_modified_root_code() reproduces its pre-0.10.1 snapshot exactly (0.10.1, audit §2.2)", {
  # 0.10.1 inlined the EventRootCode3 logic that assign_modified_root_code()
  # used to depend on (assign_eventrootcode3(), now deleted along with the
  # other 11 retired scheme functions) directly into a single case_when(),
  # renumbering straight to the 18 classes rather than via `erc3 - 1`. The
  # fixture below was captured from the pre-0.10.1 implementation, over
  # every cameo_lookup$CAMEOEVENTCODE plus every distinct EventCode in
  # original_code/gdelt_bilatr.rds (310 distinct codes total; see
  # dev/summary_0.10.1_api_cleanup.md), so this is a regression test against
  # actual production-relevant input, not just a handful of hand-picked
  # codes.
  snapshot <- readRDS(testthat::test_path("fixtures", "modified_root_code_snapshot.rds"))
  expect_equal(assign_modified_root_code(snapshot$code), snapshot$modified_root_code)
})
