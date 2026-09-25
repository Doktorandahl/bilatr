test_that("the bilatr Stan model compiles", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  mod <- compile_bilatr_model(opt_level = 1)
  expect_s3_class(mod, "CmdStanModel")
})

test_that("every registered Stan model compiles", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # `stable`/`ou`/`stable_gamma`, all current (alpha[1] > 0 by
  # construction, no legacy soft-anchor entries -- those were fully
  # retired in 0.10.0, see NEWS.md).
  expect_setequal(
    names(.bilatr_stan_models),
    c("stable", "ou", "stable_gamma")
  )

  for (name in names(.bilatr_stan_models)) {
    mod <- .compile_stan_model(name, opt_level = 1)
    expect_s3_class(mod, "CmdStanModel")
  }
})

test_that("every experimental model runs a short fixed-seed sample on real assembled data without erroring", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  stan_data <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 0,
    min_n_events = 1
  )

  experimental_models <- names(Filter(
    function(m) identical(m$status, "experimental"), .bilatr_stan_models
  ))
  expect_gt(length(experimental_models), 0)

  for (name in experimental_models) {
    fit <- suppressWarnings(fit_panel(
      stan_data,
      chains = 1,
      parallel_chains = 1,
      threads_per_chain = 1,
      iter_warmup = 25,
      iter_sampling = 5,
      seed = 1,
      opt_level = 1,
      output_dir = tempdir(),
      stan_model = name,
      refresh = 0,
      show_messages = FALSE
    ))
    expect_s3_class(fit, "CmdStanMCMC")
    expect_equal(posterior::ndraws(fit$draws()), 5)
  }
})

# The following three blocks moved from tests/testthat/test_fit_dev.R
# (0.10.1: R/fit_dev.R and its test file were deleted once fit_panel()/
# fit_dyad_ts() gained stan_model directly; these test bilatr_init_fn()/
# .alpha_raw_sum0_init() live behavior, with no dependency on the removed
# _dev wrappers).

test_that("bilatr_init_fn() initializes phi for the stable model, not phi_logn params", {
  stan_data <- list(D = 2, T = 3, A = 4)

  stable_init <- bilatr_init_fn(stan_data, stan_model = "stable")()
  expect_true("phi" %in% names(stable_init))
  expect_length(stable_init$phi, stan_data$D)
  # phi_logn-only params are gone
  expect_false(any(c("log_phi0_raw", "beta_logn") %in% names(stable_init)))
})

test_that(".alpha_raw_sum0_init() sums to exactly 0, is bounded away from the dot_self()=0 degeneracy, and starts in the anchored (alpha[1] > 0) basin", {
  for (A in c(2, 4, 5, 9)) {
    for (draw in 1:20) { # repeat: it's a random draw, the guarantees must hold every time
      v <- .alpha_raw_sum0_init(A)
      expect_length(v, A)
      expect_equal(sum(v), 0)
      expect_gt(sum(v^2), 0)
      expect_gt(v[1], 0)
    }
  }
})

test_that("bilatr_init_fn() builds correctly-shaped inits for every registered model", {
  stan_data <- list(D = 3, T = 5, A = 4)

  stable_init <- bilatr_init_fn(stan_data, stan_model = "stable")()
  expect_setequal(
    names(stable_init),
    c(
      "theta_raw", "mu_intercept", "alpha_raw",
      "sigma_theta0", "z_theta0",
      "log_process_noise_raw", "mu_log_noise", "sigma_log_noise", "phi",
      "mu_log_phi", "sigma_log_phi"
    )
  )
  # process_noise is non-centered: log_process_noise_raw, not process_noise
  expect_false("process_noise" %in% names(stable_init))
  expect_false("mu_theta0" %in% names(stable_init)) # hard-pinned, removed
  expect_length(stable_init$mu_intercept, stan_data$A) # sum_to_zero_vector[A], not A - 1
  expect_length(stable_init$alpha_raw, stan_data$A) # sum_to_zero_vector[A], not A - 1
  expect_equal(sum(stable_init$alpha_raw), 0)
  expect_gt(stable_init$alpha_raw[1], 0) # starts in the anchored basin
  expect_length(stable_init$z_theta0, stan_data$D)
  expect_gt(stats::sd(stable_init$z_theta0), 0) # real initial spread, not rep(0, D)

  ou_init <- bilatr_init_fn(stan_data, stan_model = "ou")()
  expect_setequal(
    names(ou_init),
    c(
      "theta_raw", "mu_intercept", "alpha_raw",
      "sigma_mu", "mu_dyad_raw",
      "rho", "mu_log_sd_stat", "sigma_log_sd_stat", "log_sd_stat_raw", "phi",
      "mu_log_phi", "sigma_log_phi"
    )
  )
  expect_false(any(c("mu_theta0", "mu_theta_bar") %in% names(ou_init))) # location pinned hard
  expect_length(ou_init$mu_intercept, stan_data$A)
  expect_length(ou_init$alpha_raw, stan_data$A)
  expect_equal(sum(ou_init$alpha_raw), 0)
  expect_gt(ou_init$alpha_raw[1], 0)
  expect_length(ou_init$mu_dyad_raw, stan_data$D)
  expect_gt(stats::sd(ou_init$mu_dyad_raw), 0) # real initial spread, not rep(0, D)
})

test_that("bilatr_init_fn()'s inits pass cmdstanr's init validation for every registered model", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # Specifically verifies what R/fit.R's bilatr_init_fn() assumes but
  # cannot check for itself at the R level: that cmdstanr's `init`
  # argument accepts alpha_raw's/mu_intercept's sum_to_zero_vector[A] as
  # their length-A CONSTRAINED representation, not the length-(A - 1)
  # unconstrained one. If cmdstanr instead required the unconstrained
  # form, mod$sample(init = ...) below would error.
  set.seed(1)
  D <- 2
  Tn <- 3
  A <- 4
  Y <- array(sample(0:5, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0, n_filter_dyads = 0, filter_dyads = integer(0), rho_prior_a = 8, rho_prior_b = 2
  )

  for (name in c("stable", "ou")) {
    mod <- .compile_stan_model(name, opt_level = 1)
    init_fn <- bilatr_init_fn(list(D = D, T = Tn, A = A), stan_model = name)
    fit <- suppressWarnings(suppressMessages(mod$sample(
      data = data_list, chains = 1, iter_warmup = 20, iter_sampling = 5,
      seed = 1, refresh = 0, threads_per_chain = 1,
      init = init_fn, output_dir = tempdir(), show_messages = FALSE
    )))
    expect_s3_class(fit, "CmdStanMCMC")
    expect_equal(posterior::ndraws(fit$draws()), 5)
  }
})

test_that(".compile_stan_model() caches compiled models by file + opt_level", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  mod1 <- .compile_stan_model("stable", opt_level = 1)
  mod2 <- .compile_stan_model("stable", opt_level = 1)
  expect_identical(mod1, mod2)

  mod3 <- .compile_stan_model("stable", opt_level = 1, force_recompile = TRUE)
  expect_s3_class(mod3, "CmdStanModel")
})

# A rigorous, mechanism-level test of the orientation fold (deliberately
# opposite alpha_raw inits, pinned so neither chain can adapt out of its
# basin, asserting large Rhat on the raw parameters as proof the fold --
# not coincidence -- is doing the work) lives in test_sign_ambiguity.R.
# See "the orientation fold reports alpha[1] > 0 and agreeing alpha/theta
# from two chains pinned in opposite alpha_raw basins" there.
