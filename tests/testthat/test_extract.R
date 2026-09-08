# --- extract_theta(): CSV-path matches in-memory, tiers, sequential/parallel ---

test_that("extract_theta() from CSV files matches the in-memory path exactly, sequential and parallel", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # stable model, ordinary (right-basin) init
  theta_mem <- suppressWarnings(extract_theta(fx$fit, fx$stan_data))

  theta_csv <- suppressWarnings(suppressMessages(extract_theta(fx$csv_files, fx$stan_data, chunk_size = 3)))
  expect_equal(
    dplyr::arrange(theta_mem, dyad_id, time_index),
    dplyr::arrange(theta_csv, dyad_id, time_index)
  )

  theta_csv_par <- suppressWarnings(suppressMessages(extract_theta(
    fx$csv_files, fx$stan_data, chunk_size = 3, parallel = TRUE, n_workers = 2
  )))
  expect_equal(
    dplyr::arrange(theta_mem, dyad_id, time_index),
    dplyr::arrange(theta_csv_par, dyad_id, time_index)
  )

  # default max_memory_mb (large enough that this tiny fixture fits in one chunk)
  theta_csv_default <- suppressWarnings(suppressMessages(extract_theta(fx$csv_files, fx$stan_data)))
  expect_equal(
    dplyr::arrange(theta_mem, dyad_id, time_index),
    dplyr::arrange(theta_csv_default, dyad_id, time_index)
  )

  # shape: dyad identifiers correctly reattached, not just numeric equality
  expect_true(all(c("dyad_id", "time_index", "dyad", "dyad2", "year", "mean", "5%", "50%", "95%") %in% names(theta_csv)))
})

test_that("extract_theta() from CSV files rejects non-default probs", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  expect_error(
    suppressWarnings(extract_theta(fx$csv_files, fx$stan_data, probs = c(0.1, 0.5, 0.9))),
    "only supports the default"
  )
})

test_that("extract_theta() from CSV files messages about the max_memory_mb default only when unset", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  expect_message(
    suppressWarnings(extract_theta(fx$csv_files, fx$stan_data, chunk_size = 3)),
    "Using the default max_memory_mb"
  )
  expect_no_message(
    suppressWarnings(extract_theta(fx$csv_files, fx$stan_data, chunk_size = 3, max_memory_mb = 8192)),
    message = "Using the default max_memory_mb"
  )
})

# --- extract_alpha()/extract_mu_intercept(): CSV-path matches in-memory ------

test_that("extract_alpha()/extract_mu_intercept() from CSV files match the in-memory path, across multiple chain files", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # 2 chains -> 2 csv files, the roll-up scenario
  expect_gte(length(fx$csv_files), 2)

  alpha_mem <- suppressWarnings(extract_alpha(fx$fit))
  alpha_csv <- suppressWarnings(extract_alpha(fx$csv_files))
  expect_equal(
    dplyr::arrange(alpha_mem, action_index),
    dplyr::arrange(alpha_csv, action_index)
  )

  mu_mem <- suppressWarnings(extract_mu_intercept(fx$fit))
  mu_csv <- suppressWarnings(extract_mu_intercept(fx$csv_files))
  expect_equal(
    dplyr::arrange(mu_mem, action_index),
    dplyr::arrange(mu_csv, action_index)
  )
})

test_that("extract_alpha()/extract_mu_intercept() from CSV files carry event_classes labels the same as in-memory", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  event_classes <- c("neutral", "verbal_coop", "material_coop", "conflict")

  alpha_mem <- suppressWarnings(extract_alpha(fx$fit, event_classes = event_classes))
  alpha_csv <- suppressWarnings(extract_alpha(fx$csv_files, event_classes = event_classes))
  expect_equal(alpha_csv$event_class, alpha_mem$event_class)
})

# --- reflection-symmetry flip: CSV path must apply the same raw-draws flip --

test_that("extract_theta()/extract_alpha() from CSV files apply bilatr_orient()'s flip identically to the in-memory path", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 6
  Tn <- 4
  A <- 4

  # deliberately wrong-basin init, pinned near it (see test_orient.R for
  # the full rationale): both chains get the identical, deterministic
  # negative-alpha_raw[1] init, so this also exercises the multi-chain-
  # file case for the flip specifically.
  bad_init <- function() {
    a <- stats::rnorm(A, 0, 0.5)
    a <- a - mean(a)
    if (a[1] > 0) a <- -a
    list(
      theta_raw = matrix(0, D, Tn), mu_intercept = rep(0, A), alpha_raw = a,
      sigma_theta0 = 0.5, z_theta0 = stats::rnorm(D, 0, 0.5),
      log_process_noise_raw = rep(0, D), mu_log_noise = log(0.2), sigma_log_noise = 0.3,
      phi = rep(1, D), mu_log_phi = 0, sigma_log_phi = 0.5
    )
  }

  fx <- make_csv_diagnostics_fixture(
    stan_model = "stable",
    init = bad_init,
    extra_data = list(compute_log_lik = 0, anchor_scale = 0.1),
    # pin the chain near its (deliberately wrong-basin) init, per
    # test_orient.R's rationale: this fixture's small D/T dataset has a
    # likelihood barrier too small to reliably survive ordinary
    # (adapting) warmup otherwise
    adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
    iter_warmup = 2, iter_sampling = 5
  )

  # confirm the fixture actually landed in the wrong (negative) basin,
  # otherwise this test wouldn't be exercising the flip at all
  alpha1_raw <- posterior::extract_variable(fx$fit$draws("alpha[1]"), "alpha[1]")
  expect_lt(stats::median(alpha1_raw), 0)

  theta_mem <- suppressWarnings(extract_theta(fx$fit, fx$stan_data, stan_model = "stable"))
  theta_csv <- suppressWarnings(suppressMessages(extract_theta(
    fx$csv_files, fx$stan_data, stan_model = "stable", chunk_size = 3
  )))
  expect_equal(
    dplyr::arrange(theta_mem, dyad_id, time_index),
    dplyr::arrange(theta_csv, dyad_id, time_index)
  )
  # both must have been reoriented to the canonical positive alpha[1]
  alpha_csv <- suppressWarnings(extract_alpha(fx$csv_files, stan_model = "stable"))
  expect_gt(alpha_csv$mean[1], 0)

  alpha_mem <- suppressWarnings(extract_alpha(fx$fit, stan_model = "stable"))
  expect_equal(alpha_mem$mean[1], alpha_csv$mean[1])
})
