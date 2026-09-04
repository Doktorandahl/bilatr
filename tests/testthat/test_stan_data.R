test_that("assemble_stan_data produces correctly shaped D x T x A arrays with 1s weights by default", {
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
  expect_equal(sd$dyad_weight, rep(1, sd$D))
  expect_equal(sd$period_weight, rep(1, sd$T))
  expect_equal(sd$action_weight, rep(1, sd$A))
  expect_equal(attr(sd, "event_classes"), as.character(0:4))
})

test_that("assemble_stan_data() defaults rho_prior_a/b and compute_log_lik, reproducing current behaviour when unset", {
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

  sd_custom <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 0,
    min_n_events = 1,
    rho_prior_a = 3,
    rho_prior_b = 3,
    compute_log_lik = 1
  )
  expect_equal(sd_custom$rho_prior_a, 3)
  expect_equal(sd_custom$rho_prior_b, 3)
  expect_equal(sd_custom$compute_log_lik, 1)
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

test_that("weighted = 'all' computes non-trivial, mean-1-normalized weights for all three components", {
  events <- make_fake_events(n = 800)
  events <- recode_cameo(events, code_col = "EventCode")
  sd <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    weighted = "all"
  )
  expect_false(all(sd$dyad_weight == 1))
  expect_false(all(sd$period_weight == 1))
  expect_false(all(sd$action_weight == 1))
  expect_equal(mean(sd$dyad_weight), 1, tolerance = 1e-8)
  expect_equal(mean(sd$period_weight), 1, tolerance = 1e-8)
  expect_equal(mean(sd$action_weight), 1, tolerance = 1e-8)
})

test_that("weighted accepts single components and hyphenated combinations, leaving the rest at 1s", {
  events <- make_fake_events(n = 800)
  events <- recode_cameo(events, code_col = "EventCode")

  sd_dyad <- assemble_stan_data(
    events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
    weighted = "dyad"
  )
  expect_false(all(sd_dyad$dyad_weight == 1))
  expect_equal(sd_dyad$period_weight, rep(1, sd_dyad$T))
  expect_equal(sd_dyad$action_weight, rep(1, sd_dyad$A))

  sd_dyad_period <- assemble_stan_data(
    events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
    weighted = "dyad-period"
  )
  expect_false(all(sd_dyad_period$dyad_weight == 1))
  expect_false(all(sd_dyad_period$period_weight == 1))
  expect_equal(sd_dyad_period$action_weight, rep(1, sd_dyad_period$A))

  # order within the hyphenated string shouldn't matter
  sd_period_action <- assemble_stan_data(
    events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
    weighted = "action-period"
  )
  expect_equal(sd_period_action$dyad_weight, rep(1, sd_period_action$D))
  expect_false(all(sd_period_action$period_weight == 1))
  expect_false(all(sd_period_action$action_weight == 1))
})

test_that("weighted = TRUE errors and points users at weighted = 'all'", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  expect_error(
    assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", weighted = TRUE),
    "no longer supported"
  )
})

test_that("weighted rejects unrecognized component names", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  expect_error(
    assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass", weighted = "dyad-bogus"),
    "Unrecognized weighting component"
  )
})

test_that("explicit weight vectors override the 1s default and weighted='all' computation", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    action_weight = c(1, 2, 3, 4, 5)
  )
  expect_equal(sd$action_weight, c(1, 2, 3, 4, 5))
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

test_that("dyad_ids attribute reattaches dyad_id to the dyad string for every observed dyad", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd <- assemble_stan_data(events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass")
  ids <- attr(sd, "dyad_ids")
  expect_equal(dplyr::n_distinct(ids$dyad_id), sd$D)
  expect_equal(max(ids$time_index), sd$T)
})

test_that("parse_weighted_arg parses FALSE, 'all', single components, and hyphenated combinations", {
  expect_equal(parse_weighted_arg(FALSE), character(0))
  expect_equal(parse_weighted_arg("all"), c("dyad", "period", "action"))
  expect_equal(parse_weighted_arg("dyad"), "dyad")
  expect_equal(parse_weighted_arg("dyad-period"), c("dyad", "period"))
  expect_equal(parse_weighted_arg("period-action"), c("period", "action"))
  expect_equal(parse_weighted_arg("dyad-period-action"), c("dyad", "period", "action"))
  # de-duplicates repeated components
  expect_equal(parse_weighted_arg("dyad-dyad"), "dyad")
})

test_that("parse_weighted_arg errors on TRUE, invalid components, and non-string input", {
  expect_error(parse_weighted_arg(TRUE), "no longer supported")
  expect_error(parse_weighted_arg("dyad-bogus"), "Unrecognized weighting component")
  expect_error(parse_weighted_arg(123), "must be `FALSE` or a single string")
  expect_error(parse_weighted_arg(c("dyad", "period")), "must be `FALSE` or a single string")
})

test_that("compute_default_weights only returns the requested components", {
  events_array <- array(sample(0:5, 2 * 3 * 4, replace = TRUE), dim = c(2, 3, 4))

  dyad_only <- compute_default_weights(events_array, components = "dyad")
  expect_named(dyad_only, "dyad_weight")

  dyad_period <- compute_default_weights(events_array, components = c("dyad", "period"))
  expect_named(dyad_period, c("dyad_weight", "period_weight"))

  all_three <- compute_default_weights(events_array, components = c("dyad", "period", "action"))
  expect_named(all_three, c("dyad_weight", "period_weight", "action_weight"))
})
