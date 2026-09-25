# --- diagnose_and_extract_bilatr(): fused output matches the four ---------
# --- separate calls exactly, sequential and parallel -----------------------

test_that("diagnose_and_extract_bilatr() matches calling the four functions separately", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # stable model (0.4.2+: no reflection symmetry, no basin to land in)

  diag_ref <- suppressWarnings(diagnose_convergence(fx$csv_files, n_dt = fx$n_dt, chunk_size = 3))
  theta_ref <- suppressWarnings(suppressMessages(extract_theta(fx$csv_files, fx$stan_data, chunk_size = 3)))
  alpha_ref <- suppressWarnings(extract_alpha(fx$csv_files, probs = c(0.05, 0.5, 0.95)))
  mu_ref <- suppressWarnings(extract_mu_intercept(fx$csv_files, probs = c(0.05, 0.5, 0.95)))

  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, chunk_size = 3
  )))

  expect_equal(fused$diagnostics$tier1, diag_ref$tier1)
  expect_equal(fused$diagnostics$tier2, diag_ref$tier2)
  expect_equal(fused$diagnostics$tier3, diag_ref$tier3)
  expect_equal(fused$diagnostics$summary, diag_ref$summary)

  expect_equal(
    dplyr::arrange(fused$theta, dyad_id, time_index),
    dplyr::arrange(theta_ref, dyad_id, time_index)
  )
  expect_equal(
    dplyr::arrange(fused$alpha, action_index)[c("action_index", "mean", "5%", "95%")],
    dplyr::arrange(alpha_ref, action_index)[c("action_index", "mean", "5%", "95%")]
  )
  expect_equal(
    dplyr::arrange(fused$mu_intercept, action_index)[c("action_index", "mean", "5%", "95%")],
    dplyr::arrange(mu_ref, action_index)[c("action_index", "mean", "5%", "95%")]
  )

  # parallel path gives the same result
  fused_par <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, chunk_size = 3, parallel = TRUE, n_workers = 2
  )))
  expect_equal(
    dplyr::arrange(fused_par$theta, dyad_id, time_index),
    dplyr::arrange(theta_ref, dyad_id, time_index)
  )
  expect_equal(fused_par$diagnostics$tier3, diag_ref$tier3)
})

test_that("diagnose_and_extract_bilatr() returns theta_filtered/theta_filtered_sd when compute_theta_filtered = 1, matching a direct read of the fit's own draws", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture(extra_data = list(
    compute_theta_filtered = 1, n_filter_dyads = 6L, filter_dyads = 1:6
  ))

  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt
  )))

  expect_true(nrow(fused$theta_filtered) > 0)
  expect_true(nrow(fused$theta_filtered_sd) > 0)
  expect_setequal(fused$theta_filtered$dyad_id, 1:6)

  direct <- posterior::summarise_draws(
    posterior::subset_draws(cmdstanr::read_cmdstan_csv(fx$csv_files, variables = "theta_filtered")$post_warmup_draws, variable = "theta_filtered"),
    mean = mean
  )
  direct <- dplyr::mutate(
    direct,
    dyad_id = as.integer(stringr::str_match(variable, "\\[(\\d+),")[, 2]),
    time_index = as.integer(stringr::str_match(variable, ",(\\d+)\\]")[, 2])
  )
  cmp <- merge(
    dplyr::select(fused$theta_filtered, dyad_id, time_index, mean),
    dplyr::select(direct, dyad_id, time_index, mean_direct = mean),
    by = c("dyad_id", "time_index")
  )
  expect_gt(nrow(cmp), 0)
  expect_equal(cmp$mean, cmp$mean_direct, tolerance = 1e-8)
})

test_that("diagnose_and_extract_bilatr()'s theta_filtered translates the filter_dyads subset position back to the true dyad_id, not the position itself", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # A deliberately non-contiguous, non-identity subset (dyads 3 and 5 out
  # of 6): if the position (1, 2) were joined against dyad_ids directly
  # instead of being translated via stan_data$filter_dyads, this would
  # silently mislabel the filtered rows as dyads 1/2.
  fx <- make_csv_diagnostics_fixture(extra_data = list(
    compute_theta_filtered = 1, n_filter_dyads = 2L, filter_dyads = c(3L, 5L)
  ))

  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt
  )))

  expect_setequal(unique(fused$theta_filtered$dyad_id), c(3L, 5L))
  expect_setequal(unique(fused$theta_filtered_sd$dyad_id), c(3L, 5L))
})

test_that("diagnose_and_extract_bilatr()'s theta_filtered/theta_filtered_sd are empty tibbles when compute_theta_filtered = 0 (default)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # compute_theta_filtered = 0 by default

  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt
  )))

  expect_equal(nrow(fused$theta_filtered), 0)
  expect_equal(nrow(fused$theta_filtered_sd), 0)
})

test_that("diagnose_and_extract_bilatr() folds Tier 1/2 and Tier 3 into one read when Tier 3 fits a single chunk (default chunk_size)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # stable model (0.4.2+: no reflection symmetry, no basin to land in)

  diag_ref <- suppressWarnings(diagnose_convergence(fx$csv_files, n_dt = fx$n_dt, chunk_size = 3))
  theta_ref <- suppressWarnings(suppressMessages(extract_theta(fx$csv_files, fx$stan_data, chunk_size = 3)))
  alpha_ref <- suppressWarnings(extract_alpha(fx$csv_files, probs = c(0.05, 0.5, 0.95)))

  # no chunk_size/max_memory_mb override: the fixture's whole Tier 3 set
  # is tiny relative to the 8192 MB default, so chunk_size_used ends up
  # >= length(tier3_vars) and .read_diagnostics_summary_from_csv() folds
  # Tier 1/2 and Tier 3 into a single combined read rather than a
  # separate Tier 1/2 pass plus a one-chunk Tier 3 "sweep"
  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt
  )))

  expect_equal(fused$diagnostics$tier1, diag_ref$tier1)
  expect_equal(fused$diagnostics$tier2, diag_ref$tier2)
  expect_equal(fused$diagnostics$tier3, diag_ref$tier3)
  expect_equal(
    dplyr::arrange(fused$theta, dyad_id, time_index),
    dplyr::arrange(theta_ref, dyad_id, time_index)
  )
  expect_equal(
    dplyr::arrange(fused$alpha, action_index)[c("action_index", "mean", "5%", "95%")],
    dplyr::arrange(alpha_ref, action_index)[c("action_index", "mean", "5%", "95%")]
  )
})

test_that("diagnose_and_extract_bilatr() attaches event_class labels the same as extract_alpha()", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  event_classes <- c("neutral", "verbal_coop", "material_coop", "conflict")

  alpha_ref <- suppressWarnings(extract_alpha(fx$csv_files, event_classes = event_classes))
  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, chunk_size = 3, event_classes = event_classes
  )))
  expect_equal(
    dplyr::arrange(fused$alpha, action_index)$event_class,
    dplyr::arrange(alpha_ref, action_index)$event_class
  )
})

test_that("diagnose_and_extract_bilatr() rejects an in-memory fit", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  expect_error(
    diagnose_and_extract_bilatr(fx$fit, fx$stan_data, n_dt = fx$n_dt),
    "only supports raw CmdStan CSV file paths"
  )
})

test_that("diagnose_and_extract_bilatr() requires n_dt when tiers includes 2 or 3", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  expect_error(
    diagnose_and_extract_bilatr(fx$csv_files, fx$stan_data, tiers = 1:3),
    "`n_dt` is required"
  )
  expect_no_error(
    suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(fx$csv_files, fx$stan_data, tiers = 1)))
  )
})
