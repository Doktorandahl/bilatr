# Dyad-periods outside the actors' periods (1.0.0f, Part 3f):
# aggregate_events(actor_periods = ), the `active` column of the
# dyad-period contract, prepare_stan_data()'s dyad_ids, and the
# downstream active_only/Tier 3/icc_curves handling. No CmdStan needed.

spells <- function(...) dplyr::rename(gdelt_state_spells(...), actor = iso3c)

ev_table <- function(a1, a2, date, cls = "1") {
  tibble::tibble(Actor1CountryCode = a1, Actor2CountryCode = a2, SQLDATE = date, Cls = cls)
}

active_of <- function(agg, dyad) agg$active[agg$dyad == dyad]

test_that("without actor_periods there is no active column and the table is unchanged", {
  ev <- ev_table(c("USA", "SWE"), c("CHN", "NOR"), c(20000101L, 20010101L), c("1", "2"))
  agg <- aggregate_events(ev, "Cls")
  expect_false("active" %in% names(agg))
  expect_null(attr(agg, "actor_periods"))
})

test_that("single-spell actors are active throughout; the active column sits before the counts", {
  ev <- ev_table(c("USA", "SWE"), c("CHN", "NOR"), c(20000101L, 20010101L), c("1", "2"))
  agg <- aggregate_events(ev, "Cls", actor_periods = spells())
  expect_named(agg, c("dyad", "actor1", "actor2", "year", "w_send", "active", "1", "2"))
  expect_true(all(agg$active))
  expect_identical(
    dplyr::select(agg, -"active"),
    aggregate_events(ev, "Cls"),
    ignore_attr = TRUE
  )
  expect_named(attr(agg, "actor_periods"), c("actor", "start", "end"))
})

test_that("an actor that leaves (CSK) and one that enters (SVK) are inactive outside their spells", {
  ev <- ev_table(c("USA", "SVK"), c("CSK", "CZE"), c(19910601L, 19940601L))
  agg <- aggregate_events(ev, "Cls", start = "1991-01-01", end = "1995-12-31", actor_periods = spells())
  expect_identical(active_of(agg, "USA_CSK"), c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(active_of(agg, "SVK_CZE"), c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_true(all(agg[[ "1" ]][!agg$active] == 0L))
})

test_that("monthly resolution applies the at-any-point-in-the-period rule (HRV enters 1992-01-15)", {
  ev <- ev_table("HRV", "USA", 19920210L)
  agg <- aggregate_events(ev, "Cls", resolution = "monthly", start = "1991-11-01", end = "1992-02-29",
                          actor_periods = spells())
  # Nov 1991, Dec 1991 inactive; Jan 1992 (from the 15th) and Feb 1992 active
  expect_identical(agg$active, c(FALSE, FALSE, TRUE, TRUE))
  # CSK's last month is December 1992; SVK's first is January 1993
  ev2 <- ev_table(c("USA", "SVK"), c("CSK", "USA"), c(19921201L, 19930101L))
  agg2 <- aggregate_events(ev2, "Cls", resolution = "monthly", actor_periods = spells())
  expect_identical(active_of(agg2, "USA_CSK"), c(TRUE, FALSE))
  expect_identical(active_of(agg2, "SVK_USA"), c(FALSE, TRUE))
})

test_that("a spell with a gap, and an open end, are handled", {
  # Kuwait has one post-1946 spell in COW v2024
  expect_identical(nrow(gdelt_state_spells()[gdelt_state_spells()$iso3c == "KWT", ]), 1L)
  periods <- tibble::tibble(
    actor = c("AAA", "AAA", "BBB"),
    start = c("2000-01-01", "2003-02-01", "1990-01-01"),
    end = c("2001-03-15", NA, NA)
  )
  ev <- ev_table("AAA", "BBB", c(20000601L, 20010210L, 20030310L, 20040101L))
  agg <- aggregate_events(ev, "Cls", start = "2000-01-01", end = "2005-12-31", actor_periods = periods)
  expect_identical(agg$active, c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
  expect_message(
    m <- aggregate_events(ev, "Cls", resolution = "monthly", start = "2001-02-01", end = "2003-03-31",
                          actor_periods = periods),
    "dropping 2 row"
  )
  expect_identical(m$active[1:3], c(TRUE, TRUE, FALSE)) # Feb, Mar (to the 15th), Apr 2001
  expect_identical(utils::tail(m$active, 3), c(FALSE, TRUE, TRUE)) # Jan, Feb, Mar 2003
  # Date columns work too
  periods_d <- dplyr::mutate(periods, start = as.Date(start), end = as.Date(end))
  expect_identical(
    aggregate_events(ev, "Cls", start = "2000-01-01", end = "2005-12-31", actor_periods = periods_d)$active,
    agg$active
  )
})

test_that("weekly and daily periods use their own bounds", {
  periods <- tibble::tibble(actor = c("AAA", "BBB"), start = c("2020-01-08", "1990-01-01"), end = NA)
  ev <- ev_table("AAA", "BBB", 20200110L)
  d <- aggregate_events(ev, "Cls", resolution = "daily", start = "2020-01-06", end = "2020-01-10",
                        actor_periods = periods)
  expect_identical(d$active, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  w <- aggregate_events(ev, "Cls", resolution = "weekly", start = "2019-12-31", end = "2020-01-13",
                        actor_periods = periods)
  # weeks from Tuesday 1946-01-01: 2019-12-31..2020-01-06 (inactive), 2020-01-07..13 (active)
  expect_identical(w$active, c(FALSE, TRUE))
})

test_that("an event in an inactive period is an error naming the dyad-periods", {
  ev <- ev_table(c("SVK", "SVK"), c("USA", "USA"), c(19920601L, 19930601L))
  expect_error(
    aggregate_events(ev, "Cls", actor_periods = spells()),
    "1 event\\(s\\) in 1 dyad-period\\(s\\).*SVK_USA year 1992 \\(1\\).*filter_gdelt_states"
  )
  # monthly after a year-level filter: South Sudan in February 2011
  ev2 <- ev_table("SSD", "SDN", c(20110215L, 20110801L))
  expect_error(aggregate_events(ev2, "Cls", resolution = "monthly", actor_periods = spells()),
               "SSD_SDN month_id")
  m <- aggregate_events(ev2, "Cls", resolution = "monthly", actor_periods = spells(whole_years = TRUE))
  expect_true(all(m$active))
})

test_that("actors missing from actor_periods, or none covered, are errors", {
  ev <- ev_table(c("USA", "ROM"), c("CHN", "HUN"), c(20000101L, 20000101L))
  expect_error(aggregate_events(ev, "Cls", actor_periods = spells()),
               "1 actor\\(s\\) in the data are not in `actor_periods`: ROM")
  ev2 <- ev_table("ROM", "XXX", 20000101L)
  expect_error(aggregate_events(ev2, "Cls", actor_periods = spells()), "covers none of the 2 actors")
  expect_error(aggregate_events(ev, "Cls", actor_periods = list()), "must be NULL or a data frame")
  expect_error(aggregate_events(ev, "Cls", actor_periods = data.frame(actor = "USA", from = "2000-01-01")),
               "missing column\\(s\\) `start`, `end`")
  bad <- tibble::tibble(actor = "USA", start = "not a date", end = NA)
  expect_error(aggregate_events(ev, "Cls", actor_periods = bad), "not dates")
  bad2 <- tibble::tibble(actor = "USA", start = "2001-01-01", end = "2000-01-01")
  expect_error(aggregate_events(ev, "Cls", actor_periods = bad2), "end before they start")
})

test_that("undirected aggregation marks activity on the ordered pair", {
  ev <- ev_table(c("USA", "SVK"), c("SVK", "USA"), c(19930601L, 19940601L))
  agg <- aggregate_events(ev, "Cls", directed = FALSE, start = "1992-01-01", end = "1994-12-31",
                          actor_periods = spells())
  expect_identical(agg$dyad, rep("SVK_USA", 3))
  expect_identical(agg$active, c(FALSE, TRUE, TRUE))
})

# --- the contract ---------------------------------------------------------------

dp_active <- function() {
  tibble::tibble(
    dyad = c("a_b", "a_b", "a_c"), period = 1:3 %% 2L + c(0L, 0L, 5L),
    x = c(1, 0, 2), y = c(0, 0, 1), active = c(TRUE, FALSE, TRUE)
  )
}

test_that("the validator accepts `active` and enforces its rules", {
  res <- validate_bilatr_dyad_period(dp_active(), time = "period", action_cols = c(x, y))
  expect_identical(res$active, "active")
  # `active` is a role column: not picked up as an action column
  expect_warning(res2 <- validate_bilatr_dyad_period(dp_active(), time = "period"), "x, y")
  expect_identical(res2$action_cols, c("x", "y"))
  bad <- dp_active()
  bad$x[2] <- 3
  expect_error(validate_bilatr_dyad_period(bad, time = "period", action_cols = c(x, y)),
               "rows with `active` FALSE have non-zero action counts.*\\(1 row\\(s\\); e.g. row\\(s\\) 2\\)")
  bad <- dp_active()
  bad$active[1] <- NA
  expect_error(validate_bilatr_dyad_period(bad, time = "period", action_cols = c(x, y)), "`active` has missing")
  bad <- dp_active()
  bad$active <- c(1, 0, 1)
  expect_error(validate_bilatr_dyad_period(bad, time = "period", action_cols = c(x, y)), "must be logical")
})

test_that("prepare_stan_data() carries active into dyad_ids; the Stan list is unchanged", {
  ev <- ev_table(c("USA", "SVK", "USA", "CZE"), c("CSK", "CZE", "CSK", "USA"),
                 c(19910601L, 19940601L, 19920601L, 19950601L), c("1", "2", "2", "1"))
  agg <- aggregate_events(ev, "Cls", start = "1991-01-01", end = "1995-12-31", actor_periods = spells())
  sd <- suppressMessages(prepare_stan_data(agg, reference_category = "1"))
  plain <- suppressMessages(prepare_stan_data(dplyr::select(agg, -"active"), reference_category = "1"))
  strip <- function(x) {
    attr(x, "dyad_ids") <- NULL
    x
  }
  expect_identical(strip(sd), strip(plain))
  expect_identical(dplyr::select(attr(sd, "dyad_ids"), -"active"), dplyr::select(attr(plain, "dyad_ids"), -"active"))
  ids <- attr(sd, "dyad_ids")
  expect_identical(names(ids)[8], "active")
  expect_identical(ids$active[ids$dyad == "USA_CSK"], c(TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_identical(ids$active[ids$dyad == "SVK_CZE"], c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_true(all(attr(plain, "dyad_ids")$active))
  # inactive periods are unobserved
  expect_true(all(t(sd$is_obs)[!ids$active] == 0L))
})

test_that("periods the grid adds take activity from the actor periods, else TRUE", {
  ev <- ev_table(c("USA", "SVK"), c("CSK", "CZE"), c(19910601L, 19940601L), c("1", "2"))
  agg <- aggregate_events(ev, "Cls", start = "1991-01-01", end = "1994-12-31", actor_periods = spells())
  sd <- suppressMessages(prepare_stan_data(agg, time_range = c(1990, 1996), action_cols = c("1", "2"),
                                           reference_category = "1"))
  ids <- attr(sd, "dyad_ids")
  expect_identical(ids$active[ids$dyad == "USA_CSK"], c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_identical(ids$active[ids$dyad == "SVK_CZE"], c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  # without the attribute, added periods are TRUE
  bare <- agg
  attr(bare, "actor_periods") <- NULL
  sd2 <- suppressMessages(prepare_stan_data(bare, time_range = c(1990, 1996), action_cols = c("1", "2"),
                                            reference_category = "1"))
  ids2 <- attr(sd2, "dyad_ids")
  expect_identical(ids2$active[ids2$dyad == "USA_CSK"], c(TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE))
})

test_that("min_n_events is unaffected by inactive periods", {
  ev <- ev_table(c("USA", "SVK", "SVK"), c("CSK", "CZE", "CZE"), c(19910601L, 19940601L, 19950601L), c("1", "2", "1"))
  agg <- aggregate_events(ev, "Cls", start = "1991-01-01", end = "1995-12-31", actor_periods = spells())
  sd <- suppressMessages(prepare_stan_data(agg, reference_category = "1", min_n_events = 2))
  expect_identical(unique(attr(sd, "dyad_ids")$dyad), "SVK_CZE")
})

# --- downstream, on the shipped example fit -------------------------------------

example_with_inactive <- function() {
  csv <- bilatr_example_fit()
  sd <- readRDS(list.files(dirname(csv[1]), "_stan_data\\.rds$", full.names = TRUE))
  ids <- attr(sd, "dyad_ids")
  # pretend dyad 1 is inactive in its first two periods, dyad 2 in its last
  ids$active <- !((ids$dyad_id == 1 & ids$time_index <= 2) | (ids$dyad_id == 2 & ids$time_index == max(ids$time_index)))
  attr(sd, "dyad_ids") <- ids
  list(csv = csv, stan_data = sd, n_inactive = sum(!ids$active))
}

test_that("extract_theta(active_only = TRUE) drops inactive dyad-periods (CSV path)", {
  ex <- example_with_inactive()
  all_rows <- suppressMessages(extract_theta(ex$csv, ex$stan_data, active_only = FALSE))
  kept <- suppressMessages(extract_theta(ex$csv, ex$stan_data))
  expect_true("active" %in% names(kept))
  expect_identical(nrow(all_rows) - nrow(kept), ex$n_inactive)
  expect_true(all(kept$active))
  expect_false(any(kept$dyad_id == 1 & kept$time_index <= 2))
  expect_identical(
    dplyr::filter(all_rows, active),
    kept
  )
  # stan data that predates `active`: every period counts as active
  old <- suppressMessages(extract_theta(ex$csv, max_memory_mb = 4096))
  expect_true(all(old$active))
  expect_identical(nrow(old), nrow(all_rows))
  expect_error(extract_theta(ex$csv, ex$stan_data, active_only = NA), "TRUE or FALSE")
})

test_that("diagnose_and_extract_bilatr() drops inactive periods and Tier 3 leaves them out", {
  ex <- example_with_inactive()
  res <- suppressMessages(diagnose_and_extract_bilatr(ex$csv, stan_data = ex$stan_data, max_memory_mb = 4096))
  res_all <- suppressMessages(diagnose_and_extract_bilatr(ex$csv, stan_data = ex$stan_data,
                                                          active_only = FALSE, max_memory_mb = 4096))
  expect_identical(nrow(res_all$theta) - nrow(res$theta), ex$n_inactive)
  expect_true(all(res$theta$active))
  expect_true(res$diagnostics$summary$tier3_active_only)
  expect_identical(res$diagnostics$summary$n_tier3_inactive_excluded, ex$n_inactive)
  n_theta <- res$diagnostics$tier3$n_theta
  expect_identical(n_theta[1], ex$stan_data$T - 2L)
  expect_identical(n_theta[2], ex$stan_data$T - 1L)
  expect_output(print(res$diagnostics), "active periods only, 3 inactive left out")
  # the example fit's own stan data has no `active`: no restriction
  plain <- suppressMessages(diagnose_and_extract_bilatr(ex$csv, max_memory_mb = 4096))
  expect_false(plain$diagnostics$summary$tier3_active_only)
  expect_true(all(plain$diagnostics$tier3$n_theta == ex$stan_data$T))
})

test_that("diagnose_convergence() restricts Tier 3 when the Stan data next to the CSVs marks inactive periods", {
  ex <- example_with_inactive()
  dir <- withr::local_tempdir()
  file.copy(ex$csv, dir)
  saveRDS(ex$stan_data, file.path(dir, "example_stan_data.rds"))
  csv <- file.path(dir, basename(ex$csv))
  diag <- suppressMessages(diagnose_convergence(csv, max_memory_mb = 4096))
  expect_true(diag$summary$tier3_active_only)
  expect_identical(sum(diag$tier3$n_theta), ex$stan_data$D * ex$stan_data$T - ex$n_inactive)
})

test_that("icc_curves()' fallback theta range reads active periods only", {
  ex <- example_with_inactive()
  ids <- attr(ex$stan_data, "dyad_ids")
  ids$active <- ids$dyad_id == 1 & ids$time_index == 3
  sd <- ex$stan_data
  attr(sd, "dyad_ids") <- ids
  out <- suppressMessages(icc_curves(ex$csv, stan_data = sd, n_grid = 5))
  theta1 <- suppressMessages(extract_theta(ex$csv, sd))
  expect_identical(nrow(theta1), 1L)
  # with one active dyad-period, both quantiles are its posterior mean
  expect_equal(range(out$theta), rep(theta1$mean, 2), tolerance = 1e-6)
})
