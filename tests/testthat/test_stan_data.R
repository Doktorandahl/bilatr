test_that("assemble_stan_data produces correctly shaped D x T x A arrays", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")

  sd <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 0,
    min_n_events = 1
  )

  expect_equal(dim(sd$Y), c(sd$D, sd$T, sd$A))
  expect_equal(dim(sd$is_obs), c(sd$D, sd$T))
  expect_equal(sum(sd$Y), nrow(events))
  expect_equal(attr(sd, "event_classes"), as.character(0:4))
})

test_that("assemble_stan_data() defaults rho_prior_a/b, compute_log_lik and anchor_scale, reproducing current behaviour when unset", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")

  sd_default <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 0,
    min_n_events = 1
  )
  expect_equal(sd_default$rho_prior_a, 8)
  expect_equal(sd_default$rho_prior_b, 2)
  expect_equal(sd_default$compute_log_lik, 0)
  expect_equal(sd_default$anchor_scale, 0.1)

  sd_custom <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 0,
    min_n_events = 1,
    rho_prior_a = 3,
    rho_prior_b = 3,
    compute_log_lik = 1,
    anchor_scale = 0.25
  )
  expect_equal(sd_custom$rho_prior_a, 3)
  expect_equal(sd_custom$rho_prior_b, 3)
  expect_equal(sd_custom$compute_log_lik, 1)
  expect_equal(sd_custom$anchor_scale, 0.25)
})

test_that("reference_category is reordered to action_index 1 (alpha[1] is the anchor position, not a raw event-class code)", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")

  sd <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 2,
    min_n_events = 1
  )
  expect_equal(attr(sd, "event_classes")[1], "2")
})

test_that("is_obs matches whether any events were observed in that dyad-period", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    min_n_events = 1
  )
  totals <- apply(sd$Y, c(1, 2), sum)
  expect_equal(sd$is_obs == 1L, totals > 0)
})

test_that("directed = FALSE yields fewer or equal dyads than directed = TRUE", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd_directed <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass")
  sd_undirected <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE)
  expect_lte(sd_undirected$D, sd_directed$D)
})

test_that("weighted = FALSE (default) and weighted = 'none' are both accepted, silent no-ops", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd_false <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", weighted = FALSE)
  sd_none <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", weighted = "none")
  expect_equal(sd_false, sd_none)
  expect_false(any(c("dyad_weight", "period_weight", "action_weight") %in% names(sd_false)))
})

test_that("weighted rejects anything other than FALSE/'none' (likelihood weighting removed in 0.4.6)", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  for (bad in list(TRUE, "all", "dyad", "dyad-period")) {
    expect_error(
      assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", weighted = bad),
      "removed in 0.4.6"
    )
  }
})

test_that("min_n_events drops low-activity dyads, and errors clearly if it drops all of them", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd_all <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", min_n_events = 1)
  expect_gt(sd_all$D, 0L)
  expect_error(
    assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", min_n_events = 10000),
    "No dyads"
  )
})

test_that("0e: an empty window gives a clear error naming the data's year range, not a min_n_events error", {
  events <- make_fake_events(years = 2015:2019)
  events <- recode_cameo(events, code_col = "EventCode")

  err <- tryCatch(
    assemble_stan_data(events, years = 2050, resolution = "yearly", grouping_var = "PentaClass"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "no events fall within")
  expect_match(err, "2015-2019")
  expect_false(grepl("No dyads", err))
})

test_that("0e: an empty `data` errors clearly instead of via the min_n_events path", {
  empty <- tibble::tibble(
    Actor1CountryCode = character(0),
    Actor2CountryCode = character(0),
    SQLDATE = integer(0),
    PentaClass = character(0)
  )
  err <- tryCatch(
    assemble_stan_data(empty, years = 2020, resolution = "yearly", grouping_var = "PentaClass"),
    error = function(e) conditionMessage(e)
  )
  expect_match(err, "`data` is empty")
})

test_that("dyad_ids attribute reattaches dyad_id to the dyad string for every observed dyad", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass")
  ids <- attr(sd, "dyad_ids")
  expect_equal(dplyr::n_distinct(ids$dyad_id), sd$D)
  expect_equal(max(ids$time_index), sd$T)
})
