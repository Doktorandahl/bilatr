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
    fit <- suppressWarnings(fit_panel_dev(
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
