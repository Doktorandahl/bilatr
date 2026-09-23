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

test_that("eventrootcode2_name labels every assign_eventrootcode2() output value and NAs anything else", {
  expect_equal(eventrootcode2_name("01"), "Make a public statement")
  expect_equal(eventrootcode2_name("044"), "Consult: meet, discuss, or visit")
  expect_equal(eventrootcode2_name("046"), "Consult: negotiate or mediate")
  expect_equal(eventrootcode2_name("19"), "Assault, fight, or mass violence")
  expect_true(is.na(eventrootcode2_name("99")))

  # every value assign_eventrootcode2() produces for a real CAMEO code
  # has a label (unlike arbitrary 2-digit strings, which may not
  # correspond to any real CAMEO root and are correctly unlabeled)
  produced <- unique(assign_eventrootcode2(cameo_lookup$CAMEOEVENTCODE))
  expect_false(any(is.na(eventrootcode2_name(produced))))

  # "044"/"046" deliberately match bilatr_class_name()'s text for the
  # same underlying category (see eventrootcode2_name()'s docs)
  expect_equal(eventrootcode2_name("044"), bilatr_class_name(2L))
  expect_equal(eventrootcode2_name("046"), bilatr_class_name(3L))
})

test_that("assign_eventrootcode3 relocates 016/018/019 and splits the 10/13 groupings", {
  # most roots still map to the same single class as EventRootCode2
  expect_equal(assign_eventrootcode3(c("01", "0211", "17")), c(1L, 2L, 18L))
  # 016 moves to Reject (with root 12), not root 01
  expect_equal(assign_eventrootcode3("016"), 14L)
  expect_equal(assign_eventrootcode3("12"), 14L)
  # 018/019 move to diplomatic cooperation (with root 05), not root 01
  expect_equal(assign_eventrootcode3("018"), 7L)
  expect_equal(assign_eventrootcode3("019"), 7L)
  expect_equal(assign_eventrootcode3("05"), 7L)
  # investigate (09) and demand (10) are split apart again
  expect_equal(assign_eventrootcode3(c("09", "093", "10", "100")), c(11L, 11L, 12L, 12L))
  # threaten (13) and exhibit force posture (15) stay merged...
  expect_equal(assign_eventrootcode3(c("13", "150")), c(15L, 15L))
  # ...but protest (14) becomes its own class
  expect_equal(assign_eventrootcode3("141"), 16L)
  expect_true(is.na(assign_eventrootcode3("99")))
})

test_that("eventrootcode3_name labels every assign_eventrootcode3() output value and NAs anything else", {
  expect_equal(eventrootcode3_name(14L), "Reject")
  expect_equal(eventrootcode3_name(15L), "Threaten or exhibit force posture")
  expect_equal(eventrootcode3_name(16L), "Protest")
  expect_true(is.na(eventrootcode3_name(20L)))

  produced <- unique(assign_eventrootcode3(cameo_lookup$CAMEOEVENTCODE))
  expect_false(any(is.na(eventrootcode3_name(produced))))
  expect_equal(sort(unique(produced)), 1:19)
})

test_that("eventrootcode3_rootcodes lists the codes feeding each class and NAs anything else", {
  expect_equal(eventrootcode3_rootcodes(7L), "05, 018, 019")
  expect_equal(eventrootcode3_rootcodes(14L), "12, 016")
  expect_equal(eventrootcode3_rootcodes(15L), "13, 15")
  expect_equal(eventrootcode3_rootcodes(5L), "041, 042, 043, 044")
  expect_true(is.na(eventrootcode3_rootcodes(20L)))
})

test_that("assign_eventrootcode4 splits 041, threaten/force-posture, and assault/fight/mass-violence", {
  # most roots still map to the same class as EventRootCode3
  expect_equal(assign_eventrootcode4(c("01", "0211", "17")), c(1L, 2L, 20L))
  # 016 still moves to Reject (with root 12)
  expect_equal(assign_eventrootcode4("016"), 15L)
  expect_equal(assign_eventrootcode4("12"), 15L)
  # 018/019 still move to diplomatic cooperation (with root 05)
  expect_equal(assign_eventrootcode4("018"), 8L)
  expect_equal(assign_eventrootcode4("019"), 8L)
  expect_equal(assign_eventrootcode4("05"), 8L)
  # 041 (discuss by telephone) is split out of the meet/discuss/visit class
  expect_equal(assign_eventrootcode4("041"), 5L)
  expect_equal(assign_eventrootcode4(c("042", "043", "044")), c(6L, 6L, 6L))
  # threaten (13), protest (14), and exhibit force posture (15) are all separate
  expect_equal(assign_eventrootcode4(c("13", "141", "150")), c(16L, 17L, 18L))
  # assault (18), fight (19), and mass violence (20) are all separate
  expect_equal(assign_eventrootcode4(c("180", "190", "200")), c(21L, 22L, 23L))
  expect_true(is.na(assign_eventrootcode4("99")))
})

test_that("eventrootcode4_name labels every assign_eventrootcode4() output value and NAs anything else", {
  expect_equal(eventrootcode4_name(5L), "Discuss by telephone")
  expect_equal(eventrootcode4_name(16L), "Threaten")
  expect_equal(eventrootcode4_name(18L), "Exhibit force posture")
  expect_equal(eventrootcode4_name(21L), "Assault")
  expect_equal(eventrootcode4_name(23L), "Use unconventional mass violence")
  expect_true(is.na(eventrootcode4_name(24L)))

  produced <- unique(assign_eventrootcode4(cameo_lookup$CAMEOEVENTCODE))
  expect_false(any(is.na(eventrootcode4_name(produced))))
  expect_equal(sort(unique(produced)), 1:23)
})

test_that("eventrootcode4_rootcodes lists the codes feeding each class and NAs anything else", {
  expect_equal(eventrootcode4_rootcodes(5L), "041")
  expect_equal(eventrootcode4_rootcodes(6L), "042, 043, 044")
  expect_equal(eventrootcode4_rootcodes(8L), "05, 018, 019")
  expect_equal(eventrootcode4_rootcodes(15L), "12, 016")
  expect_equal(eventrootcode4_rootcodes(16L), "13")
  expect_equal(eventrootcode4_rootcodes(18L), "15")
  expect_true(is.na(eventrootcode4_rootcodes(24L)))
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

test_that("assign_bilatr_class2 merges BilatrClass 7-8 and 9-10, passes 0-6 through", {
  expect_equal(assign_bilatr_class2(0:10), c(0:6, 7L, 7L, 8L, 8L))
  expect_true(is.na(assign_bilatr_class2(11L)))
  expect_true(is.na(assign_bilatr_class2(NA_integer_)))
})

test_that("bilatr_class2_name labels the 0-8 classes and NAs anything else", {
  expect_equal(bilatr_class2_name(0L), "Neutral / low-intensity statement")
  expect_equal(bilatr_class2_name(7L), "Disapprove, demand, reject, or reduce relations")
  expect_equal(bilatr_class2_name(8L), "Threaten, coerce, or use force")
  expect_true(is.na(bilatr_class2_name(9L)))
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

  # Only 09 and 10 collapse; everything else is a 1:1 relabelling of EventRootCode3.
  codes <- cameo_lookup$CAMEOEVENTCODE
  erc3 <- assign_eventrootcode3(codes)
  modroot <- assign_modified_root_code(codes)
  expect_equal(length(unique(modroot)), length(unique(erc3)) - 1L)
  expect_equal(nrow(unique(data.frame(erc3, modroot))), length(unique(erc3)))
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

  # Consistent with EventRootCode3's root codes, apart from the 09/10
  # merge, the 01 exception note, and the 04 top-level omission.
  erc3_codes <- eventrootcode3_rootcodes(c(2:3, 5:10, 13:19))
  expect_equal(modified_root_code_eventcodes(c(2:3, 5:10, 12:18)), erc3_codes)
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
