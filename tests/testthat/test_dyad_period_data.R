# Tests for the 0.12.0 dyad-period data contract: ?bilatr_dyad_period_data
# and validate_bilatr_dyad_period(). No CmdStan needed.

valid_dp <- function() {
  tibble::tibble(
    dyad = c("A_B", "A_B", "A_C", "B_C"),
    actor1 = c("A", "A", "A", "B"),
    actor2 = c("B", "B", "C", "C"),
    period = c(1L, 2L, 1L, 3L),
    w_send = c(0.5, 0.5, 1, 0.25),
    coop = c(1, 0, 2, 3),
    conf = c(0L, 4L, 1L, 0L)
  )
}

validate <- function(data, ...) {
  validate_bilatr_dyad_period(data, time = "period", ...)
}

test_that("a valid table passes and returns the resolved column names invisibly", {
  expect_invisible(res <- validate(valid_dp(), action_cols = c(coop, conf)))
  expect_identical(res, list(
    dyad = "dyad", time = "period", action_cols = c("coop", "conf"),
    actor1 = "actor1", actor2 = "actor2", w_send = "w_send"
  ))
})

test_that("absent optional columns resolve to NULL; dyad may be absent if both actors are present", {
  dp <- dplyr::select(valid_dp(), -"actor1", -"actor2", -"w_send")
  res <- validate(dp, action_cols = c(coop, conf))
  expect_null(res$actor1)
  expect_null(res$actor2)
  expect_null(res$w_send)

  no_dyad <- dplyr::select(valid_dp(), -"dyad")
  expect_identical(validate(no_dyad, action_cols = c(coop, conf))$dyad, "dyad")
})

test_that("missing structural columns are reported up front", {
  expect_error(validate(dplyr::select(valid_dp(), -"period"), action_cols = c(coop, conf)), "missing time column `period`")
  expect_error(
    validate(dplyr::select(valid_dp(), -"dyad", -"actor2"), action_cols = c(coop, conf)),
    "only one actor column.*missing dyad column"
  )
  expect_error(
    validate(valid_dp(), action_cols = c(coop, nope)),
    "doesn't exist"
  )
})

test_that("each rule fails on its own", {
  check_rule <- function(modify, pattern, ...) {
    dp <- modify(valid_dp())
    expect_error(validate(dp, action_cols = c(coop, conf), ...), pattern)
  }
  check_rule(function(d) { d$dyad[2] <- NA; d }, "`dyad` has missing")
  check_rule(function(d) { d$dyad[2] <- ""; d }, "`dyad` has empty-string")
  check_rule(function(d) { d$actor1[4] <- NA; d }, "`actor1` has missing")
  check_rule(function(d) { d$actor2[4] <- ""; d }, "`actor2` has empty-string")
  check_rule(function(d) { d$actor1[3] <- "A_X"; d }, "`actor1` contains \"_\"")
  check_rule(function(d) { d$actor2[4] <- "B"; d }, "`actor1` equals `actor2`")
  check_rule(function(d) { d$actor2[2] <- "D"; d }, "not constant within a dyad.*A_B")
  check_rule(function(d) { d$period <- as.character(d$period); d }, "`period` must be integer-valued")
  check_rule(function(d) { d$period <- c(1, 2.5, 1, 3); d }, "`period` is not integer-valued")
  check_rule(function(d) { d$period[1] <- NA; d }, "`period` has missing")
  check_rule(function(d) { d$w_send[3] <- 1.5; d }, "`w_send` is outside \\[0, 1\\]")
  check_rule(function(d) { d$w_send[3] <- NA; d }, "`w_send` has missing")
  check_rule(function(d) { d$w_send[2] <- 0.75; d }, "`w_send` is not constant within a dyad")
  check_rule(function(d) { d$w_send <- "half"; d }, "`w_send` must be numeric")
  check_rule(function(d) { d$coop[1] <- -1; d }, "`coop` has values that are not non-negative whole")
  check_rule(function(d) { d$coop[1] <- 0.5; d }, "`coop` has values that are not non-negative whole")
  check_rule(function(d) { d$conf[1] <- NA; d }, "`conf` has missing")
  check_rule(function(d) { d$conf <- as.character(d$conf); d }, "`conf` must be numeric")
  check_rule(function(d) { d$period[2] <- 1L; d }, "duplicate dyad x `period` rows.*A_B")
})

test_that("K >= 2 action columns are required", {
  expect_error(validate(valid_dp(), action_cols = coop), "at least two action-count columns.*selected 1 \\(coop\\)")
})

test_that("action_cols may not include a role column", {
  expect_error(validate(valid_dp(), action_cols = c(coop, conf, period)), "includes the dyad/actor/time/w_send column")
})

test_that("several broken rules are reported in a single error", {
  dp <- valid_dp()
  dp$dyad[1] <- NA
  dp$w_send[3] <- 2
  dp$coop[4] <- -3
  dp$period[4] <- 1.5
  err <- expect_error(validate(dp, action_cols = c(coop, conf)))
  expect_match(conditionMessage(err), "`dyad` has missing")
  expect_match(conditionMessage(err), "`w_send` is outside")
  expect_match(conditionMessage(err), "`coop` has values")
  expect_match(conditionMessage(err), "`period` is not integer-valued")
  expect_length(strsplit(conditionMessage(err), "\n- ")[[1]], 5)
})

test_that("default action_cols warns and lists the columns it chose", {
  expect_warning(
    res <- validate(valid_dp()),
    "`action_cols` not set.*coop, conf.*Set `action_cols` explicitly"
  )
  expect_identical(res$action_cols, c("coop", "conf"))
})

test_that("a default-selected non-count column errors and says it was picked up by the default", {
  dp <- valid_dp()
  dp$region <- c("eu", "eu", "as", "af")
  dp$share <- c(0.5, 0.5, 0, 1)
  err <- expect_error(suppressWarnings(validate(dp)))
  expect_match(conditionMessage(err), "`region` was picked up by the default `action_cols`")
  expect_match(conditionMessage(err), "`share` was picked up by the default `action_cols`")
  expect_no_match(conditionMessage(err), "`coop` was picked up")
  # explicitly selected bad columns get the ordinary rule errors instead
  err <- expect_error(validate(dp, action_cols = c(coop, share)))
  expect_match(conditionMessage(err), "action column `share` has values that are not non-negative whole")
  expect_no_match(conditionMessage(err), "picked up by the default")
})

test_that("tidyselect and character action_cols give identical results", {
  dp <- valid_dp()
  cols <- c("coop", "conf")
  by_tidy <- validate(dp, action_cols = c(coop, conf))
  expect_identical(validate(dp, action_cols = c("coop", "conf")), by_tidy)
  expect_identical(validate(dp, action_cols = cols), by_tidy)
  expect_identical(validate(dp, action_cols = dplyr::all_of(cols)), by_tidy)
  expect_identical(validate(dp, action_cols = coop:conf), by_tidy)
  expect_identical(validate(dp, action_cols = -c(dyad, actor1, actor2, period, w_send)), by_tidy)
})

test_that("a character variable named like no column does not trigger tidyselect's deprecation", {
  dp <- valid_dp()
  my_cols <- c("conf", "coop")
  expect_no_warning(res <- validate(dp, action_cols = my_cols))
  expect_identical(res$action_cols, c("conf", "coop"))
})

test_that("custom role column names work", {
  dp <- dplyr::rename(valid_dp(), pair = "dyad", from = "actor1", to = "actor2", ws = "w_send")
  res <- validate(dp, action_cols = c(coop, conf), dyad = "pair", actor1 = "from", actor2 = "to", w_send = "ws")
  expect_identical(res$actor1, "from")
  expect_identical(res$w_send, "ws")
  # absent under the default names: actors/w_send are not validated, and
  # the dyad column is required
  expect_error(validate(dp, action_cols = c(coop, conf)), "missing dyad column `dyad`")
})

test_that("leading/trailing whitespace in dyad/actor1/actor2 is an error (0.12.0b, D18)", {
  dp <- valid_dp()
  dp$actor1[1] <- " A"
  dp$actor2[3] <- "C "
  dp$dyad[4] <- "B_C "
  err <- expect_error(validate(dp, action_cols = c(coop, conf)))
  msg <- conditionMessage(err)
  expect_match(msg, "the following rule(s) were broken", fixed = TRUE)
  expect_match(msg, "`actor1` has leading/trailing whitespace (e.g. \" USA\" and \"USA\" would otherwise become two actors) (1 row(s); e.g. row(s) 1)", fixed = TRUE)
  expect_match(msg, "`actor2` has leading/trailing whitespace (e.g. \" USA\" and \"USA\" would otherwise become two actors) (1 row(s); e.g. row(s) 3)", fixed = TRUE)
  expect_match(msg, "`dyad` has leading/trailing whitespace", fixed = TRUE)
  # inner whitespace is fine
  ok <- valid_dp()
  ok$actor1 <- c("A x", "A x", "A x", "B")
  ok$dyad <- paste(ok$actor1, ok$actor2, sep = "_")
  expect_no_error(validate(ok, action_cols = c(coop, conf)))
})
