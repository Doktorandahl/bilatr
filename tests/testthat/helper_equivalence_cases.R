# The case grid of test_prepare_stan_data_equivalence.R. The stored
# reference, fixtures/assemble_stan_data_0.11.0_reference.rds, is a frozen
# output of an earlier, independent implementation of the data
# preparation, run on these same cases. Both sides build the events from
# these functions, so the inputs are the same by construction.
#
# Each case is a list:
#   events      the event table
#   grouping_var, years, resolution, directed, actor1, actor2, date
#   prep        named list of the Stan-data arguments shared by the
#               reference implementation and prepare_stan_data()
#               (reference_category, min_n_events, chunk_size, ...);
#               chunk_size is 100 unless a case sets it
# `years = a:b` maps to aggregate_events(start = "a-01-01", end = "b-12-31").
#
# Do not edit a case: the stored outputs only mean something for these
# inputs, and the reference cannot be regenerated with this version.

# This is a regression test of aggregate_events()/prepare_stan_data(),
# not of the CAMEO schemes. So the QuadClass/PentaClass columns of the
# recoded cases come from the frozen copies below (the mappings the
# reference was built with, keyed by root code), not from cameo_lookup:
# a later change to the schemes must not change these inputs.
.frozen_quad_class <- c(
  "01" = 1L, "02" = 1L, "03" = 1L, "04" = 1L, "05" = 1L,
  "06" = 2L, "07" = 2L, "08" = 2L, "09" = 2L,
  "10" = 3L, "11" = 3L, "12" = 3L, "13" = 3L, "14" = 3L,
  "15" = 4L, "16" = 4L, "17" = 4L, "18" = 4L, "19" = 4L, "20" = 4L
)
.frozen_penta_class <- c(
  "01" = 0L, "02" = 0L, "03" = 1L, "04" = 1L, "05" = 1L,
  "06" = 2L, "07" = 2L, "08" = 2L, "09" = 2L,
  "10" = 3L, "11" = 3L, "12" = 3L, "13" = 3L, "14" = 4L,
  "15" = 4L, "16" = 3L, "17" = 4L, "18" = 4L, "19" = 4L, "20" = 4L
)

# recode_cameo(), then QuadClass/PentaClass from the frozen mappings.
recode_with_frozen_classes <- function(events) {
  events <- suppressMessages(recode_cameo(events, code_col = "EventCode"))
  root <- substr(events$EventCode, 1, 2)
  events$QuadClass <- unname(.frozen_quad_class[root])
  events$PentaClass <- unname(.frozen_penta_class[root])
  events
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
  recode_with_frozen_classes(events)
}

# Events with an "action" column whose class "9" occurs only in 2010 and
# 2019, outside the 2012-2014 window.
window_events <- function() {
  events <- fake_events_recoded(n = 800, years = 2010:2019)
  events$action <- as.character(events$PentaClass)
  outside <- which(substr(events$SQLDATE, 1, 4) %in% c("2010", "2019"))[1:3]
  events$action[outside] <- "9"
  events
}

cow_events <- function() {
  tibble::tibble(
    Actor1CountryCode = c("2", "2", "365", "710", "710"),
    Actor2CountryCode = c("365", "710", "2", "2", "365"),
    SQLDATE = 20200101L,
    PentaClass = c(0, 1, 2, 3, 4)
  )
}

cow_events_renamed <- function() {
  cow_events() %>%
    dplyr::rename(sender = "Actor1CountryCode", target = "Actor2CountryCode", action = "PentaClass") %>%
    dplyr::mutate(when = as.Date(as.character(.data$SQLDATE), format = "%Y%m%d"), .keep = "unused")
}

gamma_events <- function() {
  tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN"),
    Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA"),
    SQLDATE = 20180101L,
    PentaClass = c(0, 1, 0, 1, 0)
  )
}

# ISO appears in one event only, so min_n_events = 5 drops its dyad.
iso_events <- function() {
  set.seed(1)
  iso <- tibble::tibble(
    Actor1CountryCode = c(rep("USA", 20), rep("RUS", 20), "ISO"),
    Actor2CountryCode = c(rep("RUS", 20), rep("USA", 20), "USA"),
    SQLDATE = 20180101L,
    EventCode = sample(bilatr::cameo_lookup$CAMEOEVENTCODE, 41, replace = TRUE)
  )
  recode_with_frozen_classes(iso)
}

mixed_case_events <- function() {
  tibble::tibble(
    Actor1CountryCode = c("usa", "USA", "Abc", "2", "10", "abc", "B", "usa", "10"),
    Actor2CountryCode = c("USA", "abc", "2", "10", "usa", "B", "Abc", "2", "B"),
    SQLDATE = c(20200101L, 20200301L, 20200501L, 20210101L, 20210601L, 20200101L, 20210201L, 20211201L, 20200801L),
    PentaClass = c(0, 1, 2, 0, 1, 2, 0, 1, 2)
  )
}

equivalence_case <- function(events, grouping_var, years, resolution, directed = TRUE, ...,
                             actor1 = "Actor1CountryCode", actor2 = "Actor2CountryCode",
                             date = "SQLDATE") {
  prep <- list(...)
  # The reference was built when chunk_size defaulted to 100; pass it
  # explicitly now that the default is 1.
  if (!"chunk_size" %in% names(prep)) {
    prep <- c(list(chunk_size = 100), prep)
  }
  list(
    events = events, grouping_var = grouping_var, years = years, resolution = resolution,
    directed = directed, actor1 = actor1, actor2 = actor2, date = date,
    prep = prep
  )
}

dir_label <- function(directed) if (directed) "directed" else "undirected"

# Named list of every case, grouped by the test_that() block that uses it
# (the name prefix). 74 cases.
equivalence_cases <- function() {
  cases <- list()
  events <- fake_events_recoded()

  # grid: 48 cases
  grid <- tidyr::expand_grid(
    directed = c(TRUE, FALSE),
    resolution = c("yearly", "monthly"),
    reference_category = list(NULL, "3"),
    min_n_events = c(1, 5),
    theta_filtered = c("off", "all", "subset")
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    name <- sprintf(
      "grid/%s/%s/ref_%s/min_%d/tf_%s",
      dir_label(g$directed), g$resolution, g$reference_category[[1]] %||% "NULL",
      g$min_n_events, g$theta_filtered
    )
    cases[[name]] <- equivalence_case(
      events, "PentaClass", 2015:2019, g$resolution,
      directed = g$directed,
      reference_category = g$reference_category[[1]],
      min_n_events = g$min_n_events,
      compute_theta_filtered = as.numeric(g$theta_filtered != "off"),
      filter_dyads = if (g$theta_filtered == "subset") c("RUS_USA", "CHN_USA")
    )
  }

  # min_n_events: 2 cases
  cases[["min_n/1"]] <- equivalence_case(events, "PentaClass", 2015:2019, "yearly", min_n_events = 1)
  cases[["min_n/5"]] <- equivalence_case(events, "PentaClass", 2015:2019, "yearly", min_n_events = 5)

  # other schemes and settings: 4 cases
  # bare root codes have no ModifiedRootCode (see ?recode_cameo)
  mrc_events <- dplyr::filter(events, !is.na(.data$ModifiedRootCode))
  for (directed in c(TRUE, FALSE)) {
    cases[[paste0("schemes/mrc/", dir_label(directed))]] <- equivalence_case(
      mrc_events, "ModifiedRootCode", 2015:2019, "yearly",
      directed = directed, reference_category = "5", chunk_size = 1200L,
      compute_theta_filtered = 1, rho_prior_a = 2, rho_prior_b = 3,
      compute_log_lik = 1, prior_only = 1
    )
    cases[[paste0("schemes/quad/", dir_label(directed))]] <- equivalence_case(
      events, "QuadClass", 2016:2018, "monthly", directed = directed
    )
  }

  # windowing: 4 cases
  w_events <- window_events()
  for (directed in c(TRUE, FALSE)) {
    for (resolution in c("yearly", "monthly")) {
      cases[[sprintf("window/%s/%s", dir_label(directed), resolution)]] <- equivalence_case(
        w_events, "action", 2012:2014, resolution, directed = directed
      )
    }
  }

  # COW codes, default and custom column names: 8 cases
  for (directed in c(TRUE, FALSE)) {
    for (resolution in c("yearly", "monthly")) {
      cases[[sprintf("cow/default/%s/%s", dir_label(directed), resolution)]] <- equivalence_case(
        cow_events(), "PentaClass", 2020:2020, resolution, directed = directed, reference_category = 0
      )
      cases[[sprintf("cow/renamed/%s/%s", dir_label(directed), resolution)]] <- equivalence_case(
        cow_events_renamed(), "action", 2020:2020, resolution,
        directed = directed, actor1 = "sender", actor2 = "target", date = "when"
      )
    }
  }

  # gamma fixtures: 4 cases
  for (directed in c(TRUE, FALSE)) {
    cases[[paste0("gamma/plain/", dir_label(directed))]] <- equivalence_case(
      gamma_events(), "PentaClass", 2018:2018, "yearly", directed = directed
    )
  }
  iso <- iso_events()
  for (directed in c(TRUE, FALSE)) {
    cases[[paste0("gamma/iso/", dir_label(directed))]] <- equivalence_case(
      iso, "PentaClass", 2018:2018, "yearly",
      directed = directed, reference_category = 0, min_n_events = 5
    )
  }

  # mixed-case and numeric actor codes: 4 cases
  for (directed in c(TRUE, FALSE)) {
    for (resolution in c("yearly", "monthly")) {
      cases[[sprintf("mixed/%s/%s", dir_label(directed), resolution)]] <- equivalence_case(
        mixed_case_events(), "PentaClass", 2020:2021, resolution, directed = directed
      )
    }
  }

  cases
}

# One case through aggregate_events() then prepare_stan_data().
new_path_stan_data <- function(case) {
  agg <- suppressMessages(aggregate_events(
    case$events,
    grouping_var = case$grouping_var, resolution = case$resolution, directed = case$directed,
    start = paste0(min(case$years), "-01-01"), end = paste0(max(case$years), "-12-31"),
    actor1 = case$actor1, actor2 = case$actor2, date = case$date
  ))
  suppressMessages(do.call(prepare_stan_data, c(list(agg), case$prep)))
}
