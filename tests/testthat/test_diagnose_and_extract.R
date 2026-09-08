# --- diagnose_and_extract_bilatr(): fused output matches the four ---------
# --- separate calls exactly, sequential and parallel -----------------------

test_that("diagnose_and_extract_bilatr() matches calling the four functions separately", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # stable model, ordinary (right-basin) init

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

# --- reflection-symmetry flip: fused path must orient Tier 1/2 AND Tier 3 -

test_that("diagnose_and_extract_bilatr() orients alpha/theta/mu_intercept identically to the separate calls", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 6
  Tn <- 4
  A <- 4

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
    adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
    iter_warmup = 2, iter_sampling = 5
  )

  alpha1_raw <- posterior::extract_variable(fx$fit$draws("alpha[1]"), "alpha[1]")
  expect_lt(stats::median(alpha1_raw), 0)

  theta_ref <- suppressWarnings(suppressMessages(extract_theta(
    fx$csv_files, fx$stan_data, stan_model = "stable", chunk_size = 3
  )))
  alpha_ref <- suppressWarnings(extract_alpha(
    fx$csv_files, stan_model = "stable", probs = c(0.05, 0.5, 0.95)
  ))
  mu_ref <- suppressWarnings(extract_mu_intercept(
    fx$csv_files, stan_model = "stable", probs = c(0.05, 0.5, 0.95)
  ))

  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, stan_model = "stable", chunk_size = 3
  )))

  expect_equal(
    dplyr::arrange(fused$theta, dyad_id, time_index),
    dplyr::arrange(theta_ref, dyad_id, time_index)
  )
  expect_gt(fused$alpha$mean[fused$alpha$action_index == 1], 0)
  expect_equal(
    dplyr::arrange(fused$alpha, action_index)[c("action_index", "mean", "5%", "95%")],
    dplyr::arrange(alpha_ref, action_index)[c("action_index", "mean", "5%", "95%")]
  )
  # mu_intercept is never flipped (alpha .* theta invariant under the
  # joint negation) -- confirm the fused path leaves it identical too
  expect_equal(
    dplyr::arrange(fused$mu_intercept, action_index)[c("action_index", "mean", "5%", "95%")],
    dplyr::arrange(mu_ref, action_index)[c("action_index", "mean", "5%", "95%")]
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
