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

  # `stable` plus the three experimental variants added alongside it
  # (phi_logn was retired to inst/stan/legacy/ in 0.3.2 and is not
  # registered).
  expect_setequal(
    names(.bilatr_stan_models),
    c("stable", "alphanorm", "ou", "alphanorm_ou")
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

  experimental_models <- setdiff(names(.bilatr_stan_models), .BILATR_DEFAULT_MODEL)
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
    log_process_noise_raw = stats::rnorm(D, 0, 0.3),
    mu_log_noise = log(0.2),
    sigma_log_noise = 0.3,
    phi = c(1.2, 0.9),
    mu_log_phi = 0.05,
    sigma_log_phi = 0.4,
    alpha_raw = stats::rnorm(A - 1, 0, 0.3)
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

test_that("alphanorm's soft sign anchor exactly accounts for the log-prob gap between mirror-image parameter states", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  D <- 2
  Tn <- 3
  A <- 5
  anchor_scale <- 0.1
  Y <- array(sample(0:5, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)

  sum0 <- function(x) x - mean(x)

  alpha_raw_pos <- sum0(c(2, stats::rnorm(A - 1, 0, 0.5)))
  if (alpha_raw_pos[1] < 0) alpha_raw_pos <- -alpha_raw_pos

  params_pos <- list(
    theta_raw = matrix(stats::rnorm(D * Tn, 0, 0.3), D, Tn),
    mu_intercept = sum0(stats::rnorm(A, 0, 0.3)),
    alpha_raw = alpha_raw_pos,
    sigma_theta0 = 0.6,
    z_theta0 = stats::rnorm(D, 0, 0.3),
    log_process_noise_raw = stats::rnorm(D, 0, 0.3),
    mu_log_noise = log(0.2),
    sigma_log_noise = 0.3,
    phi = c(1.2, 0.9),
    mu_log_phi = 0.05,
    sigma_log_phi = 0.4
  )

  # exact mirror image: flip alpha_raw and every theta-side quantity;
  # mu_intercept and every process/dispersion parameter are untouched,
  # per the model header's FLIP/UNCHANGED lists
  params_neg <- params_pos
  params_neg$alpha_raw <- -params_pos$alpha_raw
  params_neg$z_theta0 <- -params_pos$z_theta0
  params_neg$theta_raw <- -params_pos$theta_raw

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "bilatr_alphanorm.stan", package = "bilatr"),
    cpp_options = list(stan_threads = TRUE),
    compile_model_methods = TRUE,
    force_recompile = TRUE
  )

  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
    compute_log_lik = 0, anchor_scale = anchor_scale
  )

  logprob_at <- function(params) {
    fit <- suppressWarnings(mod$sample(
      data = data_list, chains = 1, iter_warmup = 50, iter_sampling = 5,
      seed = 1, refresh = 0, threads_per_chain = 1,
      output_dir = tempdir(), show_messages = FALSE
    ))
    fit$init_model_methods(verbose = FALSE)
    up <- fit$unconstrain_variables(variables = params)
    fit$log_prob(up, jacobian = TRUE)
  }

  lp_pos <- logprob_at(params_pos)
  lp_neg <- logprob_at(params_neg)

  alpha1_pos <- alpha_raw_pos[1] * sqrt(A / sum(alpha_raw_pos^2))
  expected_gap <- stats::plogis(alpha1_pos / anchor_scale, log.p = TRUE) -
    stats::plogis(-alpha1_pos / anchor_scale, log.p = TRUE)

  # the base likelihood, every prior, and the sum_to_zero_vector Jacobian
  # are all exactly invariant under this joint negation (see the model's
  # header, "REFLECTION SYMMETRY"), so the entire log-prob gap must come
  # from the anchor term alone
  expect_equal(lp_pos - lp_neg, expected_gap, tolerance = 1e-6)
  # and the anchor must actually favor the positive-alpha[1] mode
  expect_gt(lp_pos, lp_neg)
})

test_that("alphanorm_ou's soft sign anchor exactly accounts for the log-prob gap between mirror-image parameter states", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  set.seed(2)
  D <- 2
  Tn <- 4
  A <- 5
  anchor_scale <- 0.1
  Y <- array(sample(0:5, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)

  sum0 <- function(x) x - mean(x)

  alpha_raw_pos <- sum0(c(2, stats::rnorm(A - 1, 0, 0.5)))
  if (alpha_raw_pos[1] < 0) alpha_raw_pos <- -alpha_raw_pos

  params_pos <- list(
    theta_raw = matrix(stats::rnorm(D * Tn, 0, 0.3), D, Tn),
    mu_intercept = sum0(stats::rnorm(A, 0, 0.3)),
    alpha_raw = alpha_raw_pos,
    sigma_mu = 0.6,
    mu_dyad_raw = stats::rnorm(D, 0, 0.3),
    rho = 0.7,
    mu_log_sd_stat = log(0.8),
    sigma_log_sd_stat = 0.3,
    log_sd_stat_raw = stats::rnorm(D, 0, 0.3),
    phi = c(1.2, 0.9),
    mu_log_phi = 0.05,
    sigma_log_phi = 0.4
  )

  # exact mirror image: flip alpha_raw and every theta-side quantity
  # (mu_dyad_raw, theta_raw); mu_intercept and every process/dispersion/
  # ratio parameter are untouched, per the model header's FLIP/UNCHANGED
  # lists
  params_neg <- params_pos
  params_neg$alpha_raw <- -params_pos$alpha_raw
  params_neg$mu_dyad_raw <- -params_pos$mu_dyad_raw
  params_neg$theta_raw <- -params_pos$theta_raw

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "bilatr_alphanorm_ou.stan", package = "bilatr"),
    cpp_options = list(stan_threads = TRUE),
    compile_model_methods = TRUE,
    force_recompile = TRUE
  )

  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
    compute_log_lik = 0, rho_prior_a = 8, rho_prior_b = 2, anchor_scale = anchor_scale
  )

  logprob_at <- function(params) {
    fit <- suppressWarnings(mod$sample(
      data = data_list, chains = 1, iter_warmup = 50, iter_sampling = 5,
      seed = 1, refresh = 0, threads_per_chain = 1,
      output_dir = tempdir(), show_messages = FALSE
    ))
    fit$init_model_methods(verbose = FALSE)
    up <- fit$unconstrain_variables(variables = params)
    fit$log_prob(up, jacobian = TRUE)
  }

  lp_pos <- logprob_at(params_pos)
  lp_neg <- logprob_at(params_neg)

  alpha1_pos <- alpha_raw_pos[1] * sqrt(A / sum(alpha_raw_pos^2))
  expected_gap <- stats::plogis(alpha1_pos / anchor_scale, log.p = TRUE) -
    stats::plogis(-alpha1_pos / anchor_scale, log.p = TRUE)

  expect_equal(lp_pos - lp_neg, expected_gap, tolerance = 1e-6)
  expect_gt(lp_pos, lp_neg)
})
