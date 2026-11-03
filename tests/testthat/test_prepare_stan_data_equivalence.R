# Equivalence of the 0.12.0 data path with 0.11.0: for event data,
# prepare_stan_data(aggregate_events(ev, ...), ...) must give exactly the
# Stan input assemble_stan_data(ev, ...) gives (every list field, and the
# event_classes/country_codes/grouping_var/dyad_ids attributes). This is
# what establishes that the replication results are unchanged (plan D15),
# so every comparison is expect_identical(). No CmdStan needed.

# Builds the Stan data both ways and compares. `years = a:b` maps to
# `start = "a-01-01", end = "b-12-31"`. Messages are not compared (the
# function-name prefixes differ by design).
expect_paths_identical <- function(events, grouping_var, years, resolution, directed = TRUE, ...,
                                   actor1 = "Actor1CountryCode", actor2 = "Actor2CountryCode",
                                   date = "SQLDATE") {
  old <- suppressMessages(assemble_stan_data(
    events,
    years = years, resolution = resolution, grouping_var = grouping_var,
    directed = directed, actor1 = actor1, actor2 = actor2, date = date, ...
  ))
  agg <- suppressMessages(aggregate_events(
    events,
    grouping_var = grouping_var, resolution = resolution, directed = directed,
    start = paste0(min(years), "-01-01"), end = paste0(max(years), "-12-31"),
    actor1 = actor1, actor2 = actor2, date = date
  ))
  # 0.12.0b (D16): the action_cols attribute from aggregate_events() is
  # used silently.
  expect_no_warning(
    new <- suppressMessages(prepare_stan_data(agg, ...))
  )

  expect_identical(names(new), names(old))
  for (field in names(old)) {
    expect_identical(new[[field]], old[[field]], label = paste0("new$", field))
  }
  for (a in c("event_classes", "country_codes", "grouping_var")) {
    expect_identical(attr(new, a), attr(old, a), label = paste0("attr(new, \"", a, "\")"))
  }
  old_ids <- dplyr::rename(attr(old, "dyad_ids"), actor1 = "actor_a", actor2 = "actor_b")
  new_ids <- attr(new, "dyad_ids")
  if (resolution == "monthly") {
    new_ids$month_id <- NULL
  }
  expect_identical(new_ids, old_ids)
  invisible(list(old = old, new = new))
}

fake_events_recoded <- function(...) {
  events <- make_fake_events(...)
  # A sparse pair, so a min_n_events threshold of 5 drops dyads.
  events <- dplyr::bind_rows(events, tibble::tibble(
    Actor1CountryCode = c("FRA", "FRA", "USA"),
    Actor2CountryCode = c("USA", "USA", "FRA"),
    SQLDATE = c(20160315L, 20170704L, 20180101L),
    EventCode = c("040", "190", "040")
  ))
  suppressMessages(recode_cameo(events, code_col = "EventCode"))
}

test_that("the full option grid on make_fake_events() is identical to assemble_stan_data()", {
  events <- fake_events_recoded()
  grid <- tidyr::expand_grid(
    directed = c(TRUE, FALSE),
    resolution = c("yearly", "monthly"),
    reference_category = list(NULL, "3"),
    min_n_events = c(1, 5),
    theta_filtered = c("off", "all", "subset")
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    filter_dyads <- if (g$theta_filtered == "subset") c("RUS_USA", "CHN_USA")
    expect_paths_identical(
      events, "PentaClass", 2015:2019, g$resolution,
      directed = g$directed,
      reference_category = g$reference_category[[1]],
      min_n_events = g$min_n_events,
      compute_theta_filtered = as.numeric(g$theta_filtered != "off"),
      filter_dyads = filter_dyads
    )
  }
})

test_that("min_n_events = 5 really drops dyads (the grid's dropping case is not vacuous)", {
  events <- fake_events_recoded()
  res1 <- expect_paths_identical(events, "PentaClass", 2015:2019, "yearly", min_n_events = 1)
  res5 <- expect_paths_identical(events, "PentaClass", 2015:2019, "yearly", min_n_events = 5)
  expect_lt(res5$new$D, res1$new$D)
  expect_false(any(grepl("FRA", attr(res5$new, "dyad_ids")$dyad)))
})

test_that("ModifiedRootCode and QuadClass, other chunk_size/prior settings, are identical", {
  events <- fake_events_recoded()
  # bare root codes have no ModifiedRootCode (see ?recode_cameo)
  mrc_events <- dplyr::filter(events, !is.na(.data$ModifiedRootCode))
  for (directed in c(TRUE, FALSE)) {
    expect_paths_identical(
      mrc_events, "ModifiedRootCode", 2015:2019, "yearly",
      directed = directed, reference_category = "5", chunk_size = 1200L,
      compute_theta_filtered = 1, rho_prior_a = 2, rho_prior_b = 3,
      compute_log_lik = 1, prior_only = 1
    )
    expect_paths_identical(events, "QuadClass", 2016:2018, "monthly", directed = directed)
  }
})

test_that("a window that drops rows and drops a class is identical", {
  events <- fake_events_recoded(n = 800, years = 2010:2019)
  events$action <- as.character(events$PentaClass)
  # class "9" occurs only outside the 2012-2014 window
  outside <- which(substr(events$SQLDATE, 1, 4) %in% c("2010", "2019"))[1:3]
  events$action[outside] <- "9"

  for (directed in c(TRUE, FALSE)) {
    for (resolution in c("yearly", "monthly")) {
      msgs <- capture_messages(aggregate_events(
        events, "action", resolution, directed,
        start = "2012-01-01", end = "2014-12-31"
      ))
      expect_match(msgs, "dropping [0-9]+ row\\(s\\) outside", all = FALSE)
      expect_match(msgs, "zero in-window events and get no column: 9", all = FALSE)
      res <- expect_paths_identical(events, "action", 2012:2014, resolution, directed = directed)
      expect_false("9" %in% attr(res$new, "event_classes"))
    }
  }
})

test_that("the COW-code fixture (test_event_data.R) is identical, with default and custom column names", {
  events <- tibble::tibble(
    Actor1CountryCode = c("2", "2", "365", "710", "710"),
    Actor2CountryCode = c("365", "710", "2", "2", "365"),
    SQLDATE = 20200101L,
    PentaClass = c(0, 1, 2, 3, 4)
  )
  renamed <- events %>%
    dplyr::rename(sender = "Actor1CountryCode", target = "Actor2CountryCode", action = "PentaClass") %>%
    dplyr::mutate(when = as.Date(as.character(.data$SQLDATE), format = "%Y%m%d"), .keep = "unused")

  for (directed in c(TRUE, FALSE)) {
    for (resolution in c("yearly", "monthly")) {
      expect_paths_identical(events, "PentaClass", 2020:2020, resolution, directed = directed, reference_category = 0)
      expect_paths_identical(
        renamed, "action", 2020:2020, resolution,
        directed = directed, actor1 = "sender", actor2 = "target", date = "when"
      )
    }
  }
})

test_that("the gamma fixtures (test_gamma_offset.R) are identical", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN"),
    Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA"),
    SQLDATE = 20180101L,
    PentaClass = c(0, 1, 0, 1, 0)
  )
  for (directed in c(TRUE, FALSE)) {
    res <- expect_paths_identical(events, "PentaClass", 2018:2018, "yearly", directed = directed)
  }
  # undirected: a non-trivial sender share made it through unchanged
  expect_true(any(res$new$w_send != 1))

  set.seed(1)
  iso <- tibble::tibble(
    Actor1CountryCode = c(rep("USA", 20), rep("RUS", 20), "ISO"),
    Actor2CountryCode = c(rep("RUS", 20), rep("USA", 20), "USA"),
    SQLDATE = 20180101L,
    EventCode = sample(bilatr::cameo_lookup$CAMEOEVENTCODE, 41, replace = TRUE)
  )
  iso <- suppressMessages(recode_cameo(iso, code_col = "EventCode"))
  for (directed in c(TRUE, FALSE)) {
    res <- expect_paths_identical(
      iso, "PentaClass", 2018:2018, "yearly",
      directed = directed, reference_category = 0, min_n_events = 5
    )
    expect_false("ISO" %in% attr(res$new, "country_codes"))
  }
})

test_that("mixed-case and numeric actor codes give the same dyad_id order on both paths", {
  events <- tibble::tibble(
    Actor1CountryCode = c("usa", "USA", "Abc", "2", "10", "abc", "B", "usa", "10"),
    Actor2CountryCode = c("USA", "abc", "2", "10", "usa", "B", "Abc", "2", "B"),
    SQLDATE = c(20200101L, 20200301L, 20200501L, 20210101L, 20210601L, 20200101L, 20210201L, 20211201L, 20200801L),
    PentaClass = c(0, 1, 2, 0, 1, 2, 0, 1, 2)
  )
  for (directed in c(TRUE, FALSE)) {
    for (resolution in c("yearly", "monthly")) {
      res <- expect_paths_identical(events, "PentaClass", 2020:2021, resolution, directed = directed)
      dyads <- unique(attr(res$new, "dyad_ids")$dyad)
      expect_identical(dyads, sort(dyads, method = "radix"))
    }
  }
})
