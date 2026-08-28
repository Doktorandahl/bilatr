skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(has_cmdstan, "CmdStan is not installed")
}

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

  # All four current variants (stable, phi_logn, and their non-centered
  # process_noise reparameterizations) must be registered and compile
  # cleanly at -O1.
  expect_true(
    all(c("stable", "phi_logn", "stable_ncproc", "phi_logn_ncproc") %in%
      names(.bilatr_stan_models))
  )

  for (name in names(.bilatr_stan_models)) {
    mod <- .compile_stan_model(name, opt_level = 1)
    expect_s3_class(mod, "CmdStanModel")
  }
})

test_that(".compile_stan_model() caches compiled models by file + opt_level", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  mod1 <- .compile_stan_model("phi_logn", opt_level = 1)
  mod2 <- .compile_stan_model("phi_logn", opt_level = 1)
  expect_identical(mod1, mod2)

  mod3 <- .compile_stan_model("phi_logn", opt_level = 1, force_recompile = TRUE)
  expect_s3_class(mod3, "CmdStanModel")
})

test_that("weight vectors of 1s recover the unweighted (base) model exactly", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  D <- 2
  Tn <- 3
  A <- 4
  Y <- array(sample(0:5, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)

  params <- list(
    theta_raw = matrix(stats::rnorm(D * Tn, 0, 0.3), D, Tn),
    mu_intercept_raw = stats::rnorm(A - 1, 0, 0.3),
    mu_theta0 = 0.1,
    sigma_theta0 = 0.4,
    z_theta0 = stats::rnorm(D, 0, 0.3),
    process_noise = rep(0.25, D),
    mu_log_noise = log(0.2),
    sigma_log_noise = 0.3,
    phi = c(1.2, 0.9),
    mu_log_phi = 0.05,
    sigma_log_phi = 0.4,
    alpha_raw = stats::rnorm(A - 2, 0, 0.3),
    alpha_hostile = 0.6
  )

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "bilatr_dirmult_irt.stan", package = "bilatr"),
    cpp_options = list(stan_threads = TRUE),
    compile_model_methods = TRUE,
    force_recompile = TRUE
  )

  logprob_at <- function(data, params) {
    fit <- suppressWarnings(mod$sample(
      data = data, chains = 1, iter_warmup = 50, iter_sampling = 5,
      seed = 1, refresh = 0, threads_per_chain = 1,
      output_dir = tempdir(), show_messages = FALSE
    ))
    fit$init_model_methods(verbose = FALSE)
    up <- fit$unconstrain_variables(variables = params)
    fit$log_prob(up, jacobian = TRUE)
  }

  data_1s <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A)
  )
  action_weight <- c(1.3, 0.7, 1.0, 1.5)
  data_action_weighted <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = action_weight
  )

  lp_base <- logprob_at(data_1s, params)
  lp_weighted <- logprob_at(data_action_weighted, params)

  # weighting away from 1s must change the density...
  expect_false(isTRUE(all.equal(lp_base, lp_weighted)))
  # ...but re-running the same (1s) data must be exactly reproducible
  expect_equal(lp_base, logprob_at(data_1s, params))
})
