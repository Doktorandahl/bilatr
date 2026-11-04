# Tests for aggregate_events(). No CmdStan needed.

# Four events with hand-computable IDs. Days since 1946-01-01: 1946-2019
# is 74 years with 18 leap days (1948, ..., 2016), so 2020-01-01 is day
# 74 * 365 + 18 = 27028. Hence
#   2019-12-31 = day 27027 = week 3861 (27027 = 3861 * 7; a Tuesday, like 1946-01-01)
#   2020-01-02 = day 27029 = week 3861
#   2020-01-15 = day 27042 = week 3863
#   2020-03-01 = day 27028 + 31 + 29 = 27088 = week 3869
# and month_id = year * 12 + month - 1: 2019-12 = 24239, 2020-01 = 24240,
# 2020-03 = 24242.
hand_events <- function() {
  tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "RUS", "CHN"),
    Actor2CountryCode = c("RUS", "RUS", "USA", "USA"),
    SQLDATE = c(20191231L, 20200102L, 20200115L, 20200301L),
    action = c("a", "b", "a", "b")
  )
}

agg <- function(...) aggregate_events(hand_events(), "action", ...)

test_that("yearly: year IDs, the full skeleton, and counts", {
  out <- agg("yearly")
  expect_named(out, c("dyad", "actor1", "actor2", "year", "w_send", "a", "b"))
  expect_identical(out$dyad, rep(c("CHN_USA", "RUS_USA", "USA_RUS"), each = 2))
  expect_identical(out$year, rep(2019:2020, 3))
  expect_identical(out$a, c(0L, 0L, 0L, 1L, 1L, 0L))
  expect_identical(out$b, c(0L, 1L, 0L, 0L, 0L, 1L))
})

test_that("monthly: month_id arithmetic across a year boundary, with month labels", {
  out <- agg("monthly")
  expect_identical(unique(out$month_id), 24239:24242)
  usa_rus <- out[out$dyad == "USA_RUS", ]
  expect_identical(usa_rus$a, c(1L, 0L, 0L, 0L))
  expect_identical(usa_rus$b, c(0L, 1L, 0L, 0L))
  labels <- attr(out, "time_labels")
  expect_identical(labels, tibble::tibble(
    month_id = 24239:24242,
    year = c(2019L, 2020L, 2020L, 2020L),
    month = c(12L, 1L, 2L, 3L)
  ))
})

test_that("weekly and daily IDs count from 1946-01-01", {
  weekly <- agg("weekly")
  expect_identical(range(weekly$week_id), c(3861L, 3869L))
  expect_identical(attr(weekly, "time_origin"), as.Date("1946-01-01"))
  rus_usa <- weekly[weekly$dyad == "RUS_USA" & weekly$a > 0, ]
  expect_identical(rus_usa$week_id, 3863L)
  week_labels <- attr(weekly, "time_labels")
  expect_identical(week_labels$week_start[week_labels$week_id == 3861L], as.Date("2019-12-31"))
  expect_identical(weekdays(week_labels$week_start[1]), weekdays(as.Date("1946-01-01")))

  daily <- agg("daily")
  expect_identical(range(daily$day_id), c(27027L, 27088L))
  expect_identical(nrow(daily), 3L * (27088L - 27027L + 1L))
  chn <- daily[daily$dyad == "CHN_USA" & daily$b > 0, ]
  expect_identical(chn$day_id, 27088L)
  day_labels <- attr(daily, "time_labels")
  expect_identical(day_labels$date[day_labels$day_id == 27088L], as.Date("2020-03-01"))
})

test_that("pre-1946 data moves the origin to 1 January of the earliest year, with a warning", {
  events <- hand_events()
  events$SQLDATE[1] <- 19400601L
  expect_warning(
    out <- aggregate_events(events, "action", "daily"),
    "before 1946-01-01.*1940-01-01.*not comparable"
  )
  expect_identical(attr(out, "time_origin"), as.Date("1940-01-01"))
  expect_identical(min(out$day_id), as.integer(as.Date("1940-06-01") - as.Date("1940-01-01")))

  expect_warning(
    out <- agg("weekly", start = "1930-05-05"),
    "1930-01-01"
  )
  expect_identical(attr(out, "time_origin"), as.Date("1930-01-01"))
  expect_identical(min(out$week_id), as.integer(as.Date("1930-05-05") - as.Date("1930-01-01")) %/% 7L)

  # yearly/monthly need no origin: no warning, origin NA
  expect_no_warning(out <- agg("yearly", start = "1930-05-05"))
  expect_true(is.na(attr(out, "time_origin")))
  expect_s3_class(attr(out, "time_origin"), "Date")
})

test_that("start/end: NULL is the full range, one-sided windows, and start > end errors", {
  full <- agg("yearly")
  expect_identical(agg("yearly", start = NULL, end = NULL), full)

  expect_message(left <- agg("monthly", start = "2020-01-01"), "dropping 1 row")
  expect_identical(range(left$month_id), c(24240L, 24242L))
  expect_false("USA_RUS" %in% left$dyad[left$a > 0])

  expect_message(right <- agg("monthly", end = as.Date("2020-01-31")), "dropping 1 row")
  expect_identical(range(right$month_id), c(24239L, 24240L))
  # CHN_USA's only event is outside the window: the dyad is gone
  expect_false("CHN_USA" %in% right$dyad)

  # the skeleton covers the window, not just the data
  wide <- agg("yearly", start = "2015-01-01", end = "2022-12-31")
  expect_identical(unique(wide$year), 2015:2022)

  expect_error(agg("yearly", start = "2020-06-01", end = "2020-01-01"), "`start` \\(2020-06-01\\) is after `end`")
  expect_error(agg("yearly", start = "not a date"), "`start` must be NULL, a single Date")
})

test_that("empty data and empty windows error, reporting the data's date range", {
  expect_error(aggregate_events(hand_events()[0, ], "action"), "`data` is empty")
  expect_error(
    suppressMessages(agg("yearly", start = "2021-01-01")),
    "no events fall within.*the data's dates span 2019-12-31 to 2020-03-01"
  )
  expect_error(
    suppressMessages(agg("yearly", end = "2010-01-01")),
    "no events fall within.*the data's dates span 2019-12-31 to 2020-03-01"
  )
})

test_that("a class with no in-window events gets no column, with a message", {
  msgs <- capture_messages(out <- agg("yearly", end = "2020-01-01"))
  expect_match(msgs, "zero in-window events and get no column: b", all = FALSE)
  expect_false("b" %in% names(out))
  expect_true("a" %in% names(out))
})

test_that("class columns default to C-locale numeric order; class_order must be a permutation", {
  events <- hand_events()
  events$action <- c("10", "2", "B", "a")
  expect_named(aggregate_events(events, "action")[-(1:5)], c("2", "10", "B", "a"))

  out <- agg("yearly", class_order = c("b", "a"))
  expect_named(out, c("dyad", "actor1", "actor2", "year", "w_send", "b", "a"))

  expect_error(agg("yearly", class_order = "a"), "permutation.*missing: b")
  expect_error(agg("yearly", class_order = c("a", "b", "c")), "not in the \\(in-window\\) data: c")
  expect_error(agg("yearly", class_order = c("a", "b", "b")), "duplicated: b")
})

test_that("class values that would clash with the fixed columns error", {
  events <- hand_events()
  events$action <- c("a", "year", "a", "b")
  expect_error(aggregate_events(events, "action"), "Offending class\\(es\\): \"year\"")
})

test_that("the output columns are exactly dyad, actor1, actor2, <time>, w_send, <classes>", {
  events <- hand_events()
  events$extra <- "ignored"
  for (res in c("yearly", "monthly", "weekly", "daily")) {
    out <- aggregate_events(events, "action", res)
    time_col <- c(yearly = "year", monthly = "month_id", weekly = "week_id", daily = "day_id")[[res]]
    expect_named(out, c("dyad", "actor1", "actor2", time_col, "w_send", "a", "b"))
    expect_type(out[[time_col]], "integer")
    expect_type(out$a, "integer")
    expect_identical(out, dplyr::arrange(out, .data$dyad, .data[[time_col]]))
  }
})

test_that("w_send is 1 when directed and the hand-computed sender share when undirected", {
  directed <- agg("yearly", directed = TRUE)
  expect_true(all(directed$w_send == 1))

  undirected <- agg("yearly", directed = FALSE)
  expect_identical(unique(undirected$dyad), c("CHN_USA", "RUS_USA"))
  # RUS_USA (side A = RUS): USA sends events 1 and 2, RUS sends event 3 -> 1/3
  expect_identical(unique(undirected$w_send[undirected$dyad == "RUS_USA"]), 1 / 3)
  expect_identical(unique(undirected$w_send[undirected$dyad == "CHN_USA"]), 1)
  expect_identical(unique(undirected$actor1), c("CHN", "RUS"))
  # the pair's events pooled: RUS_USA 2020 has events 2 (b) and 3 (a)
  expect_identical(undirected$a[undirected$dyad == "RUS_USA"], c(1L, 1L))
  expect_identical(undirected$b[undirected$dyad == "RUS_USA"], c(0L, 1L))
})

test_that("attributes are present", {
  out <- agg("weekly", directed = FALSE)
  expect_s3_class(attr(out, "time_labels"), "tbl_df")
  expect_named(attr(out, "time_labels"), c("week_id", "week_start"))
  expect_identical(attr(out, "time_labels")$week_id, 3861:3869)
  expect_identical(attr(out, "time_col"), "week_id")
  expect_identical(attr(out, "time_origin"), as.Date("1946-01-01"))
  expect_identical(attr(out, "resolution"), "weekly")
  expect_identical(attr(out, "grouping_var"), "action")
  expect_false(attr(out, "directed"))

  yearly <- agg("yearly")
  expect_identical(attr(yearly, "time_labels"), tibble::tibble(year = 2019:2020))
  expect_true(is.na(attr(yearly, "time_origin")))
})

test_that("events are validated against every row, including out-of-window rows", {
  events <- hand_events()
  events$SQLDATE[1] <- 20191399L
  expect_error(aggregate_events(events, "action", start = "2020-01-01"), "could not be parsed as a date")
})
