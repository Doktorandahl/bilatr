# --- extract_theta(): CSV-path matches in-memory, tiers, sequential/parallel ---

test_that("extract_theta() from CSV files matches the in-memory path exactly, sequential and parallel", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # stable model (0.4.2+: no reflection symmetry, no basin to land in)
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

