# Tests for the 0.9.0 event-data contract: ?bilatr_event_data,
# validate_bilatr_events(), and the D1-D4 fixes in
# dev/audit_2026-09-23.md / dev/claude_code_prompt_0.9.0_data_contract.md.
# No CmdStan needed.

test_that("D1: non-3-letter actor codes work for directed and undirected data (COW numeric codes)", {
  # Before 0.9.0: directed data with non-3-letter codes errored
  # ("Internal invariant violated: computed w_send != 1 for directed
  # data ..."), because side A/B were parsed back out of the `dyad`
  # string with str_sub(dyad, 1, 3)/(5, 7). Undirected data silently got
  # garbage country codes/dyad keys from the same parsing (audit D1).
  events <- tibble::tibble(
    Actor1CountryCode = c("2", "2", "365", "710", "710"),
    Actor2CountryCode = c("365", "710", "2", "2", "365"),
    SQLDATE = 20200101L,
    PentaClass = c(0, 1, 2, 3, 4)
  )

  directed <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = TRUE
  )
  expect_equal(attr(directed, "w_send")$w_send, rep(1, 5))

  undirected <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE
  )
  w_send <- attr(undirected, "w_send")
  expect_setequal(w_send$dyad, c("2_365", "2_710", "365_710"))
  # hand count: "2_365" pools event 1 (actor1 = "2" = side A) and event 3
  # (actor1 = "365" = side B) -> 1 of 2 as sender.
  expect_equal(w_send$w_send[w_send$dyad == "2_365"], 0.5)
  # "2_710" pools event 2 (actor1 = "2" = side A) and event 4
  # (actor1 = "710" = side B) -> 1 of 2 as sender.
  expect_equal(w_send$w_send[w_send$dyad == "2_710"], 0.5)
  # "365_710" is only event 5, actor1 = "710" = side B -> 0 of 1.
  expect_equal(w_send$w_send[w_send$dyad == "365_710"], 0)

  sd <- assemble_stan_data(
    events,
    years = 2020, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )
  expect_setequal(attr(sd, "country_codes"), c("2", "365", "710"))
  expect_true(all(sd$w_send == 1))

  dyad_ids <- attr(sd, "dyad_ids")
  expect_true(all(c("actor_a", "actor_b") %in% names(dyad_ids)))
  by_dyad <- dplyr::distinct(dyad_ids, dyad, dyad2, actor_a, actor_b)
  expect_equal(by_dyad$actor_a[by_dyad$dyad == "2_365"], "2")
  expect_equal(by_dyad$actor_b[by_dyad$dyad == "2_365"], "365")
  # dyad "365_2" (row 3's directed key) carries the same undirected dyad2
  # as "2_365" -- the reverse-direction event of the same pair.
  expect_equal(by_dyad$dyad2[by_dyad$dyad == "365_2"], "2_365")
})

test_that("D1: ISO2 and mixed-length actor labels work, directed and undirected", {
  events <- tibble::tibble(
    Actor1CountryCode = c("US", "US", "Freedonia"),
    Actor2CountryCode = c("CN", "Freedonia", "CN"),
    SQLDATE = 20200101L,
    PentaClass = c(0, 1, 2)
  )

  directed <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = TRUE
  )
  expect_setequal(directed$dyad, c("US_CN", "US_Freedonia", "Freedonia_CN"))
  expect_equal(attr(directed, "w_send")$w_send, c(1, 1, 1))

  undirected <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE
  )
  expect_setequal(undirected$dyad, c("CN_US", "Freedonia_US", "CN_Freedonia"))
})

test_that("custom actor1/actor2/date column names and a Date column give the same result as the default-named equivalent", {
  base <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "CHN", "RUS"),
    Actor2CountryCode = c("RUS", "CHN", "RUS", "USA"),
    SQLDATE = c(20200101L, 20200201L, 20200301L, 20200401L),
    PentaClass = c(0, 1, 2, 3)
  )
  renamed <- base %>%
    dplyr::rename(sender = Actor1CountryCode, target = Actor2CountryCode, action = PentaClass) %>%
    dplyr::mutate(when = as.Date(as.character(SQLDATE), format = "%Y%m%d")) %>%
    dplyr::select(-SQLDATE)

  sd_default <- assemble_stan_data(
    base,
    years = 2020, resolution = "monthly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )
  sd_custom <- assemble_stan_data(
    renamed,
    years = 2020, resolution = "monthly", grouping_var = "action",
    reference_category = 0, min_n_events = 1,
    actor1 = "sender", actor2 = "target", date = "when"
  )

  compare_fields <- c("D", "T", "A", "Y", "is_obs", "ctry_a", "ctry_b", "w_send")
  expect_equal(sd_default[compare_fields], sd_custom[compare_fields])
  expect_equal(attr(sd_default, "event_classes"), attr(sd_custom, "event_classes"))
  expect_equal(attr(sd_default, "country_codes"), attr(sd_custom, "country_codes"))
})

test_that("a user's own year/date/dyad/event_type columns are never read or overwritten (no clobbering)", {
  base <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "CHN"),
    Actor2CountryCode = c("RUS", "CHN", "RUS"),
    SQLDATE = c(20200101L, 20200201L, 20200301L),
    PentaClass = c(0, 1, 2)
  )
  clobber_prone <- base %>%
    dplyr::mutate(
      year = "not-a-year",
      date = "not-a-date",
      dyad = "not-a-dyad",
      event_type = "not-a-class"
    )

  compare_fields <- c("D", "T", "A", "Y", "is_obs", "ctry_a", "ctry_b", "w_send")
  sd_plain <- assemble_stan_data(
    base,
    years = 2020, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )
  sd_clobber_prone <- assemble_stan_data(
    clobber_prone,
    years = 2020, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )

  expect_equal(sd_plain[compare_fields], sd_clobber_prone[compare_fields])
})

test_that("validate_bilatr_events() collects every broken rule into one error", {
  # Breaks four rules at once, across different rows: an NA actor1 (row
  # 1), an actor1 containing "_" (row 3), a self-dyad (row 4), and an
  # unparseable date (row 4).
  events <- tibble::tibble(
    Actor1CountryCode = c(NA, "USA", "US_A", "CHN"),
    Actor2CountryCode = c("RUS", "RUS", "RUS", "CHN"),
    SQLDATE = c(20200101L, 20200102L, 20200103L, NA_integer_),
    PentaClass = c(0, 1, 2, 3)
  )

  err <- tryCatch(
    validate_bilatr_events(events, grouping_var = "PentaClass"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "missing \\(NA\\)")
  expect_match(err, "\"_\"")
  expect_match(err, "self-dyad")
  expect_match(err, "date")
})

test_that("D2: an NA grouping_var value errors, naming the count (bare ModifiedRootCode root, NEWS 0.8.0 case)", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA"),
    Actor2CountryCode = c("RUS", "RUS"),
    SQLDATE = 20200101L,
    EventCode = c("04", "010")
  )
  events <- suppressMessages(recode_cameo(events, code_col = "EventCode"))
  expect_true(is.na(events$ModifiedRootCode[1]))

  expect_error(
    grouped_events_to_dyad_period(events, resolution = "yearly", grouping_var = "ModifiedRootCode"),
    "1 row"
  )
  expect_error(
    grouped_events_to_dyad_period(events, resolution = "yearly", grouping_var = "ModifiedRootCode"),
    "ModifiedRootCode"
  )
})

test_that("D3: a class present only outside `years` is absent from event_classes/Y, and named in the message", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "USA"),
    Actor2CountryCode = c("RUS", "RUS", "RUS"),
    SQLDATE = c(19900101L, 20200101L, 20200201L),
    PentaClass = c(9, 0, 1)
  )

  msgs <- character(0)
  sd <- withCallingHandlers(
    assemble_stan_data(
      events,
      years = 2020, resolution = "yearly", grouping_var = "PentaClass", min_n_events = 1
    ),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_false("9" %in% attr(sd, "event_classes"))
  expect_true(any(grepl("9", msgs, fixed = TRUE)))
})

test_that("D4: a reference_category absent from the in-window data errors (present only outside the window too, D3+D4)", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "USA"),
    Actor2CountryCode = c("RUS", "RUS", "RUS"),
    SQLDATE = c(19900101L, 20200101L, 20200201L),
    PentaClass = c(9, 0, 1)
  )

  # never present at all
  expect_error(
    suppressMessages(assemble_stan_data(
      events,
      years = 2020, resolution = "yearly", grouping_var = "PentaClass",
      reference_category = 42, min_n_events = 1
    )),
    "not present"
  )
  # present, but only outside the analysis window
  expect_error(
    suppressMessages(assemble_stan_data(
      events,
      years = 2020, resolution = "yearly", grouping_var = "PentaClass",
      reference_category = 9, min_n_events = 1
    )),
    "not present"
  )
})

test_that("recode_cameo() errors on numeric EventCode, and reports unmatched/bare-root counts", {
  events_numeric <- tibble::tibble(EventCode = c(10L, 190L))
  expect_error(recode_cameo(events_numeric), "leading zero")

  events <- tibble::tibble(EventCode = c("044", "not-a-code", "04"))
  expect_message(recode_cameo(events, code_col = "EventCode"), "did not match")
  expect_message(recode_cameo(events, code_col = "EventCode"), "bare top-level")
})
