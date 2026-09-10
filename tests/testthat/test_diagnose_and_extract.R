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

test_that("diagnose_and_extract_bilatr() folds Tier 1/2 and Tier 3 into one read when Tier 3 fits a single chunk (default chunk_size)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture() # stable model, ordinary (right-basin) init

  diag_ref <- suppressWarnings(diagnose_convergence(fx$csv_files, n_dt = fx$n_dt, chunk_size = 3))
  theta_ref <- suppressWarnings(suppressMessages(extract_theta(fx$csv_files, fx$stan_data, chunk_size = 3)))
  alpha_ref <- suppressWarnings(extract_alpha(fx$csv_files, probs = c(0.05, 0.5, 0.95)))

  # no chunk_size/max_memory_mb override: the fixture's whole Tier 3 set
  # is tiny relative to the 8192 MB default, so chunk_size_used ends up
  # >= length(tier3_vars) and .read_and_orient_draws_summary() is used
  # for a single combined Tier-1/2-plus-Tier-3 read rather than a
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

test_that("diagnose_and_extract_bilatr() with the pre-0.4.0 'alphanorm' alias orients identically to stan_model = 'stable' (B1)", {
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

  stable_result <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, stan_model = "stable", chunk_size = 3
  )))

  .reset_bilatr_alias_messaged()
  expect_message(
    alias_result <- suppressWarnings(diagnose_and_extract_bilatr(
      fx$csv_files, fx$stan_data, n_dt = fx$n_dt, stan_model = "alphanorm", chunk_size = 3
    )),
    "pre-0.4.0 name of 'stable'"
  )

  expect_equal(alias_result$theta, stable_result$theta)
  expect_equal(alias_result$alpha, stable_result$alpha)
  expect_equal(alias_result$mu_intercept, stable_result$mu_intercept)
  expect_gt(alias_result$alpha$mean[alias_result$alpha$action_index == 1], 0)
})

test_that("diagnose_and_extract_bilatr() orients correctly in the folded Tier-1/2-plus-Tier-3 read (default chunk_size)", {
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

  # default chunk_size: Tier 3 fits in one chunk, so this exercises the
  # combined Tier-1/2-plus-Tier-3 read in .read_and_orient_draws_summary(),
  # not the separate-reads branch the other orientation test above uses
  # (chunk_size = 3)
  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, stan_model = "stable"
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
})

test_that(".chunked_summarise_csv() flips only theta/theta_raw within a chunk that also holds log_lik, never the whole chunk (B5)", {
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
    extra_data = list(compute_log_lik = 1, anchor_scale = 0.1),
    adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
    iter_warmup = 2, iter_sampling = 5
  )
  alpha1_raw <- posterior::extract_variable(fx$fit$draws("alpha[1]"), "alpha[1]")
  expect_lt(stats::median(alpha1_raw), 0)

  prepared <- .prepare_fast_csv_read(fx$csv_files)
  var_tiers <- .classify_bilatr_tier(prepared$variables)
  tier3_vars <- var_tiers$variable[var_tiers$tier == 3L]
  theta_vars <- grep("^theta\\[", tier3_vars, value = TRUE)
  log_lik_vars <- grep("^log_lik\\[", tier3_vars, value = TRUE)
  # confirm the fixture actually exercises the mixed case B5 fixed
  expect_gt(length(theta_vars), 0)
  expect_gt(length(log_lik_vars), 0)

  flip_vars <- .bilatr_flip_variables("stable")

  # a single chunk holding every Tier 3 variable, so theta/theta_raw
  # (flip_vars) and log_lik (not in flip_vars) get negated within the
  # SAME chunk -- exactly the case B5 fixed (previously the whole chunk
  # flipped uniformly whenever it contained any flip_vars column)
  summ <- .chunked_summarise_csv(
    prepared, tier3_vars,
    chunk_size = length(tier3_vars), n_cores = 1L, flip_vars = flip_vars
  )

  raw <- cmdstanr::read_cmdstan_csv(fx$csv_files, variables = c("theta", "theta_raw", "log_lik"))
  raw_draws <- posterior::as_draws_df(raw$post_warmup_draws)
  raw_mean <- function(v) mean(posterior::extract_variable(raw_draws, v))

  summ_mean <- function(v) summ$mean[match(v, summ$variable)]

  expect_equal(
    summ_mean(theta_vars),
    -vapply(theta_vars, raw_mean, numeric(1)),
    tolerance = 1e-8, ignore_attr = TRUE
  )
  expect_equal(
    summ_mean(log_lik_vars),
    vapply(log_lik_vars, raw_mean, numeric(1)),
    tolerance = 1e-8, ignore_attr = TRUE
  )
})

test_that("diagnose_and_extract_bilatr(tiers = 3) still orients theta via the first-chunk path, with no alpha[1] leaking into alpha", {
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

  # tiers = 3 alone: no Tier 1/2 read at all, so orientation must come
  # from the first-Tier3-chunk path (.chunked_summarise_csv_with_orientation())
  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, stan_model = "stable",
    tiers = 3, chunk_size = 3
  )))

  expect_equal(
    dplyr::arrange(fused$theta, dyad_id, time_index),
    dplyr::arrange(theta_ref, dyad_id, time_index)
  )
  # tier 1 was never requested: alpha/mu_intercept must be empty, not a
  # silently-incomplete single-row table from the alpha[1] injected only
  # for orientation
  expect_equal(nrow(fused$alpha), 0)
  expect_equal(nrow(fused$mu_intercept), 0)
})

test_that("diagnose_and_extract_bilatr(tiers = c(2, 3)) orients theta via the injected-into-Tier-1/2-read path", {
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

  # tiers = c(2, 3): tier12_vars is Tier-2-only (non-empty), so alpha[1]
  # is injected into THAT read for orientation, then dropped afterward
  fused <- suppressWarnings(suppressMessages(diagnose_and_extract_bilatr(
    fx$csv_files, fx$stan_data, n_dt = fx$n_dt, stan_model = "stable",
    tiers = c(2, 3), chunk_size = 3
  )))

  expect_equal(
    dplyr::arrange(fused$theta, dyad_id, time_index),
    dplyr::arrange(theta_ref, dyad_id, time_index)
  )
  expect_equal(nrow(fused$alpha), 0)
  expect_equal(nrow(fused$mu_intercept), 0)
  expect_true(is.null(fused$diagnostics$tier1))
  expect_false(is.null(fused$diagnostics$tier2))
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
