test_that("assemble_stan_data()'s output at weighted = FALSE differs from before 0.4.6 by exactly the three dropped weight fields", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  sd <- assemble_stan_data(
    events,
    years = 2015:2019,
    resolution = "yearly",
    grouping_var = "PentaClass",
    reference_category = 0,
    min_n_events = 1
  )
  # pre-0.4.6, this same call additionally returned dyad_weight/
  # period_weight/action_weight (all vectors of 1s); every other field
  # (and every value of every other field) is unchanged.
  expect_false(any(c("dyad_weight", "period_weight", "action_weight") %in% names(sd)))
  expect_setequal(
    names(sd),
    c(
      "D", "T", "A", "C", "is_obs", "Y", "rho_prior_a", "rho_prior_b",
      "compute_log_lik", "prior_only", "compute_theta_filtered",
      "n_filter_dyads", "filter_dyads", "anchor_scale"
    )
  )
})

test_that("the 0.4.6 stable program is bit-identical to the pre-0.4.6 one at unit weights", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # Proof this is a pure deletion, not a behavior change: dyad_weight/
  # period_weight/action_weight all being 1 must leave the target exactly
  # (not approximately) unchanged, since removing the elementwise
  # multiplies is IEEE-754-exact at unit weights. Compiles a pre-0.4.6
  # copy of the retired-weighting stable program (git blob ff3e7e1, the
  # commit immediately before this task) standalone -- NOT via
  # .compile_stan_model(), which always reads the package's live,
  # already-changed inst/stan/ -- and fits it alongside the current one
  # with identical data, a fixed (non-random) init, seed, and adaptation.
  old_stan_src <- system2("git", c("show", "ff3e7e1:inst/stan/bilatr_alphanorm.stan"), stdout = TRUE)
  old_stan_file <- tempfile(fileext = ".stan")
  writeLines(old_stan_src, old_stan_file)

  mod_old <- cmdstanr::cmdstan_model(
    old_stan_file,
    cpp_options = list(stan_threads = TRUE),
    force_recompile = TRUE
  )
  mod_new <- .compile_stan_model("stable", opt_level = 1, force_recompile = TRUE)

  set.seed(1)
  D <- 4
  Tn <- 5
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)

  data_new <- list(T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y, compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0, n_filter_dyads = 0, filter_dyads = integer(0))
  data_old <- utils::modifyList(
    data_new,
    list(dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A))
  )

  sum0 <- function(x) x - mean(x)
  set.seed(2)
  fixed_init <- list(
    theta_raw = matrix(stats::rnorm(D * Tn, 0, 0.3), D, Tn),
    mu_intercept = sum0(stats::rnorm(A, 0, 0.3)),
    alpha_raw = sum0(c(2, stats::rnorm(A - 1, 0, 0.3))),
    sigma_theta0 = 0.5,
    z_theta0 = stats::rnorm(D, 0, 0.3),
    log_process_noise_raw = stats::rnorm(D, 0, 0.3),
    mu_log_noise = log(0.2),
    sigma_log_noise = 0.3,
    phi = stats::rlnorm(D, 0, 0.2),
    mu_log_phi = 0.05,
    sigma_log_phi = 0.4
  )

  fit_at <- function(mod, data) {
    suppressWarnings(mod$sample(
      data = data, chains = 2, parallel_chains = 1,
      iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
      threads_per_chain = 1, init = function() fixed_init,
      output_dir = tempdir(), show_messages = FALSE
    ))
  }

  fit_old <- fit_at(mod_old, data_old)
  fit_new <- fit_at(mod_new, data_new)

  lp_old <- posterior::extract_variable_matrix(fit_old$draws("lp__"), "lp__")
  lp_new <- posterior::extract_variable_matrix(fit_new$draws("lp__"), "lp__")

  shared_vars <- intersect(
    posterior::variables(fit_old$draws()),
    posterior::variables(fit_new$draws())
  )
  draws_old <- posterior::as_draws_matrix(fit_old$draws(variables = shared_vars))
  draws_new <- posterior::as_draws_matrix(fit_new$draws(variables = shared_vars))

  bit_identical_lp <- identical(lp_old, lp_new)
  bit_identical_draws <- identical(draws_old, draws_new)

  diag_old <- fit_old$sampler_diagnostics()
  diag_new <- fit_new$sampler_diagnostics()
  bit_identical_diagnostics <- identical(diag_old, diag_new)

  cat("\nbit-identical lp__:        ", bit_identical_lp, "\n")
  cat("bit-identical draws matrix:", bit_identical_draws, "\n")
  cat("bit-identical diagnostics: ", bit_identical_diagnostics, "\n")
  if (!bit_identical_lp) {
    cat("max abs lp__ diff:", max(abs(lp_old - lp_new)), "\n")
  }
  if (!bit_identical_draws) {
    cat("max abs draws diff:", max(abs(draws_old - draws_new)), "\n")
  }

  expect_true(bit_identical_lp)
  expect_true(bit_identical_draws)
  expect_true(bit_identical_diagnostics)
})

test_that("weighted = 'dyad-period'/'all'/TRUE all error with the 0.4.6 removal message", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")
  for (bad in list("dyad-period", "all", TRUE)) {
    expect_error(
      assemble_stan_data(
        events,
        years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
        weighted = bad
      ),
      "removed in 0.4.6"
    )
  }
})
