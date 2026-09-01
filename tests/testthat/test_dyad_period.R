test_that("order_event_classes places reference_category first, rest sorted after it", {
  expect_equal(
    order_event_classes(c("3", "1", "0", "2", "4"), reference_category = "2"),
    c("2", "0", "1", "3", "4")
  )
  # falls back to plain sort when no reference is supplied
  expect_equal(order_event_classes(c("b", "a", "c")), c("a", "b", "c"))
})

test_that("validate_reference_class warns and drops values absent from the data", {
  expect_warning(out <- validate_reference_class("z", c("a", "b"), "reference_category"))
  expect_null(out)
  expect_equal(validate_reference_class("a", c("a", "b"), "reference_category"), "a")
  expect_null(validate_reference_class(NULL, c("a", "b"), "reference_category"))
})

test_that("grouped_events_to_dyad_period produces one row per dyad-year with correctly ordered class columns", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")

  agg <- grouped_events_to_dyad_period(
    events,
    resolution = "yearly",
    grouping_var = "PentaClass",
    directed = TRUE,
    reference_category = 0
  )

  expect_true(all(c("dyad", "year", "total_events") %in% names(agg)))
  class_cols <- grep("^EventClass_", names(agg), value = TRUE)
  expect_equal(class_cols, c("EventClass_0", "EventClass_1", "EventClass_2", "EventClass_3", "EventClass_4"))
  expect_equal(nrow(agg), dplyr::n_distinct(paste(agg$dyad, agg$year)))
  expect_equal(rowSums(agg[class_cols]), agg$total_events)
})

test_that("directed = FALSE collapses ordered pairs to a single undirected dyad key", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "CHN"),
    Actor2CountryCode = c("CHN", "USA"),
    SQLDATE = c(20200101L, 20200201L),
    EventCode = c("044", "19")
  )
  events <- recode_cameo(events, code_col = "EventCode")

  directed <- grouped_events_to_dyad_period(events, resolution = "yearly", grouping_var = "PentaClass", directed = TRUE)
  undirected <- grouped_events_to_dyad_period(events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE)

  expect_equal(dplyr::n_distinct(directed$dyad), 2L)
  expect_equal(dplyr::n_distinct(undirected$dyad), 1L)
})

test_that("fill_dyad_period_skeleton fills gaps with zero counts and flags is_obs correctly", {
  agg <- tibble::tibble(
    dyad = "USA_CHN",
    year = 2016,
    EventClass_0 = 3L,
    total_events = 3L
  )
  filled <- fill_dyad_period_skeleton(agg, years = 2015:2017, resolution = "yearly")

  expect_equal(nrow(filled), 3L)
  expect_equal(filled$is_obs, c(0L, 1L, 0L))
  expect_equal(filled$total_events, c(0L, 3L, 0L))
})
