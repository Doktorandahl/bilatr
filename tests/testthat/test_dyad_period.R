test_that("order_event_classes places reference_category first, rest sorted after it", {
  expect_equal(
    order_event_classes(c("3", "1", "0", "2", "4"), reference_category = "2"),
    c("2", "0", "1", "3", "4")
  )
  # falls back to plain sort when no reference is supplied
  expect_equal(order_event_classes(c("b", "a", "c")), c("a", "b", "c"))
})

test_that("validate_reference_class errors on a value absent from the data (0.9.0; was a warning + NULL)", {
  expect_error(
    validate_reference_class("z", c("a", "b"), "reference_category"),
    "not present"
  )
  expect_equal(validate_reference_class("a", c("a", "b"), "reference_category"), "a")
  expect_null(validate_reference_class(NULL, c("a", "b"), "reference_category"))
})

# The blocks below were ported in 0.12.0 from grouped_events_to_dyad_period()
# and fill_dyad_period_skeleton() to aggregate_events() (which now also
# fills the skeleton) and prepare_stan_data() (which builds is_obs).

test_that("aggregate_events() gives one row per dyad-year with C-locale numerically ordered class columns", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")

  agg <- suppressMessages(aggregate_events(events, "PentaClass", "yearly", directed = TRUE))

  class_cols <- attr(agg, "action_cols")
  expect_equal(class_cols, c("0", "1", "2", "3", "4"))
  expect_named(agg, c("dyad", "actor1", "actor2", "year", "w_send", class_cols))
  expect_equal(nrow(agg), dplyr::n_distinct(paste(agg$dyad, agg$year)))
  # every event is counted exactly once
  expect_equal(sum(as.matrix(agg[class_cols])), nrow(events))
})

test_that("directed = FALSE collapses ordered pairs to a single undirected dyad key", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "CHN"),
    Actor2CountryCode = c("CHN", "USA"),
    SQLDATE = c(20200101L, 20200201L),
    EventCode = c("044", "19")
  )
  events <- recode_cameo(events, code_col = "EventCode")

  directed <- aggregate_events(events, "PentaClass", "yearly", directed = TRUE)
  undirected <- aggregate_events(events, "PentaClass", "yearly", directed = FALSE)

  expect_equal(dplyr::n_distinct(directed$dyad), 2L)
  expect_equal(dplyr::n_distinct(undirected$dyad), 1L)
})

test_that("0d: aggregate_events(start =, end =) reports how many rows it dropped", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "USA"),
    Actor2CountryCode = c("RUS", "RUS", "RUS"),
    SQLDATE = c(20100101L, 20200101L, 20200201L),
    PentaClass = c(0, 1, 2)
  )
  expect_message(
    aggregate_events(events, "PentaClass", "yearly", start = "2020-01-01", end = "2020-12-31"),
    "dropping 1 row"
  )
  # no window: no filtering, no message
  expect_no_message(aggregate_events(events, "PentaClass", "yearly"))
})

test_that("gaps in the window are filled with zero counts, and is_obs flags the observed periods", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "USA"),
    Actor2CountryCode = c("CHN", "CHN", "CHN"),
    SQLDATE = c(20160101L, 20160201L, 20160301L),
    PentaClass = c(0, 0, 1)
  )
  agg <- aggregate_events(events, "PentaClass", "yearly", start = "2015-01-01", end = "2017-12-31")

  expect_equal(nrow(agg), 3L)
  expect_equal(agg$year, 2015:2017)
  expect_equal(agg$`0`, c(0L, 2L, 0L))
  expect_equal(agg$`1`, c(0L, 1L, 0L))
  expect_equal(agg$actor1, rep("USA", 3L))
  expect_equal(agg$actor2, rep("CHN", 3L))

  sd <- suppressMessages(prepare_stan_data(agg))
  expect_equal(sd$is_obs, matrix(c(0L, 1L, 0L), nrow = 1))
  expect_equal(sd$T, 3)
})
