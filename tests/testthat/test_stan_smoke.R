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

  # `stable`/`ou` (current, alpha[1] > 0 by construction) plus the
  # retired `stable_soft_anchor`/`ou_soft_anchor` (pre-0.4.2, kept
  # registered so their CmdStan output stays readable -- see NEWS.md).
  # phi_logn was retired to inst/stan/legacy/ in 0.3.2, and the pre-0.4.0
  # stable/ou were retired there in 0.4.0 when alphanorm/alphanorm_ou
  # were promoted; neither of those two is registered.
  expect_setequal(
    names(.bilatr_stan_models),
    c("stable", "ou", "stable_soft_anchor", "ou_soft_anchor")
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

test_that("every legacy model still runs a short fixed-seed sample (readable, not fit going forward)", {
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

  legacy_models <- names(Filter(
    function(m) identical(m$status, "legacy"), .bilatr_stan_models
  ))
  expect_setequal(legacy_models, c("stable_soft_anchor", "ou_soft_anchor"))

  # assemble_stan_data() no longer supplies dyad_weight/period_weight/
  # action_weight (retired in 0.4.6, see NEWS.md), but the legacy programs
  # still declare them -- add unit weights by hand here, since this test's
  # only job is confirming the retired .stan files still compile/run, not
  # exercising a real production data-assembly path (none of the current
  # runscripts fit these models).
  stan_data <- utils::modifyList(stan_data, list(
    dyad_weight = rep(1, stan_data$D),
    period_weight = rep(1, stan_data$T),
    action_weight = rep(1, stan_data$A)
  ))

  for (name in legacy_models) {
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

test_that("fit_panel_dev() reaches sampling with the pre-0.4.0 'alphanorm' alias (Step 6, verification section 3)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # before fit_bilatr() canonicalised stan_model itself, this reached
  # .compile_stan_model() (via .resolve_stan_model(), which
  # canonicalises) and .warn_if_wrong_basin() (via
  # .bilatr_flip_variables(), which also canonicalises) fine, but died
  # in bilatr_init_fn(), which switches on the raw name and has no
  # "alphanorm" entry of its own
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
  # "alphanorm" resolves to the legacy stable_soft_anchor program, which
  # still declares dyad_weight/period_weight/action_weight (retired from
  # assemble_stan_data() in 0.4.6, see NEWS.md) -- add unit weights by
  # hand, same reasoning as the legacy-model test above.
  stan_data <- utils::modifyList(stan_data, list(
    dyad_weight = rep(1, stan_data$D),
    period_weight = rep(1, stan_data$T),
    action_weight = rep(1, stan_data$A)
  ))

  .reset_bilatr_alias_messaged()
  fit <- suppressWarnings(suppressMessages(fit_panel_dev(
    stan_data,
    chains = 1,
    parallel_chains = 1,
    threads_per_chain = 1,
    iter_warmup = 25,
    iter_sampling = 5,
    seed = 1,
    opt_level = 1,
    output_dir = tempdir(),
    stan_model = "alphanorm",
    refresh = 0,
    show_messages = FALSE
  )))
  expect_s3_class(fit, "CmdStanMCMC")
  expect_equal(posterior::ndraws(fit$draws()), 5)
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
  # form, mod$sample(init = ...) below would error. `stable`/`ou` and
  # their retired `_soft_anchor` counterparts share this exact parameter
  # shape (they differ only in what the .stan program does with
  # alpha_raw's sign, not in shape -- see bilatr_init_fn()'s docs), so one
  # loop covers all four.
  set.seed(1)
  D <- 2
  Tn <- 3
  A <- 4
  Y <- array(sample(0:5, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
    compute_log_lik = 0, anchor_scale = 0.1, rho_prior_a = 8, rho_prior_b = 2
  )

  for (name in c("stable", "ou", "stable_soft_anchor", "ou_soft_anchor")) {
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

test_that("legacy stable_soft_anchor's soft sign anchor exactly accounts for the log-prob gap between mirror-image parameter states", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # Coverage for the retired program's still-relevant anchor math: the
  # current `stable` shares this exact reflection symmetry in its raw
  # parameters (same free sum_to_zero_vector alpha_raw), but folds the
  # sign into the REPORTED alpha/theta instead of anchoring it, so this
  # log-prob-gap mechanism is specific to the legacy soft-anchor program
  # -- see inst/stan/bilatr_alphanorm.stan's header, "IDENTIFICATION:
  # ORIENTATION FOLD".
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
    system.file("stan", "legacy/bilatr_stable_soft_anchor.stan", package = "bilatr"),
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

test_that("legacy ou_soft_anchor's soft sign anchor exactly accounts for the log-prob gap between mirror-image parameter states", {
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
    system.file("stan", "legacy/bilatr_ou_soft_anchor.stan", package = "bilatr"),
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

# A rigorous, mechanism-level test of the orientation fold (deliberately
# opposite alpha_raw inits, pinned so neither chain can adapt out of its
# basin, asserting large Rhat on the raw parameters as proof the fold --
# not coincidence -- is doing the work) lives in test_orient.R alongside
# the analogous legacy wrong-basin test, since both rely on the same
# pinning trick. See "the orientation fold reports alpha[1] > 0 and
# agreeing alpha/theta from two chains pinned in opposite alpha_raw
# basins" there.
