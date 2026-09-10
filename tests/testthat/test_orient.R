.make_synthetic_draws <- function(extra) {
  n <- length(extra[[1]])
  df <- data.frame(
    .chain = 1L, .iteration = seq_len(n), .draw = seq_len(n),
    as.data.frame(extra, check.names = FALSE),
    check.names = FALSE
  )
  posterior::as_draws_df(df)
}

test_that("bilatr_orient() flips exactly the stable flip-list and leaves the rest unchanged", {
  set.seed(1)
  n <- 20
  extra <- list(
    `alpha[1]` = stats::rnorm(n, -2, 0.1), # deliberately negative median
    `alpha[2]` = stats::rnorm(n, 1, 0.1),
    `alpha_raw[1]` = stats::rnorm(n, -1, 0.1),
    `theta[1,1]` = stats::rnorm(n, -0.5, 0.1),
    `theta[1,2]` = stats::rnorm(n, -0.6, 0.1),
    `theta0[1]` = stats::rnorm(n, -0.3, 0.1),
    `z_theta0[1]` = stats::rnorm(n, -0.4, 0.1),
    `theta_raw[1,1]` = stats::rnorm(n, -0.1, 0.1),
    `mu_intercept[1]` = stats::rnorm(n, 0.2, 0.1),
    `phi[1]` = stats::rnorm(n, 1, 0.1),
    `sigma_theta0` = stats::rnorm(n, 0.5, 0.1)
  )
  draws <- .make_synthetic_draws(extra)

  flip_vars <- c("alpha", "alpha_raw", "theta", "theta0", "z_theta0", "theta_raw")
  unchanged_vars <- c("mu_intercept", "phi", "sigma_theta0")
  oriented <- bilatr_orient(draws, stan_model = "stable", variables = c(flip_vars, unchanged_vars))

  for (v in c("alpha[1]", "alpha[2]", "alpha_raw[1]", "theta[1,1]", "theta[1,2]", "theta0[1]", "z_theta0[1]", "theta_raw[1,1]")) {
    expect_equal(
      posterior::extract_variable(oriented, v), -extra[[v]],
      info = paste("expected", v, "to be negated")
    )
  }
  for (v in c("mu_intercept[1]", "phi[1]", "sigma_theta0")) {
    expect_equal(
      posterior::extract_variable(oriented, v), extra[[v]],
      info = paste("expected", v, "to be unchanged")
    )
  }

  # the invariant the whole symmetry rests on: alpha .* theta unaffected
  alpha1_theta11_before <- extra[["alpha[1]"]] * extra[["theta[1,1]"]]
  alpha1_theta11_after <- posterior::extract_variable(oriented, "alpha[1]") *
    posterior::extract_variable(oriented, "theta[1,1]")
  expect_equal(alpha1_theta11_after, alpha1_theta11_before)
})

test_that("bilatr_orient() flips exactly the ou flip-list and leaves the rest unchanged", {
  set.seed(2)
  n <- 20
  extra <- list(
    `alpha[1]` = stats::rnorm(n, -1.5, 0.1),
    `alpha_raw[1]` = stats::rnorm(n, -0.8, 0.1),
    `theta[1,1]` = stats::rnorm(n, -0.4, 0.1),
    `mu_dyad[1]` = stats::rnorm(n, -0.2, 0.1),
    `mu_dyad_raw[1]` = stats::rnorm(n, -0.3, 0.1),
    `theta_raw[1,1]` = stats::rnorm(n, -0.1, 0.1),
    `mu_intercept[1]` = stats::rnorm(n, 0.2, 0.1),
    `sigma_mu` = stats::rnorm(n, 0.5, 0.1),
    `sd_stat[1]` = stats::rnorm(n, 0.4, 0.1),
    `rho` = stats::runif(n, 0.5, 0.9),
    `within_between_ratio` = stats::rnorm(n, 1, 0.1)
  )
  draws <- .make_synthetic_draws(extra)

  flip_vars <- c("alpha", "alpha_raw", "theta", "mu_dyad", "mu_dyad_raw", "theta_raw")
  unchanged_vars <- c("mu_intercept", "sigma_mu", "sd_stat", "rho", "within_between_ratio")
  oriented <- bilatr_orient(draws, stan_model = "ou", variables = c(flip_vars, unchanged_vars))

  for (v in c("alpha[1]", "alpha_raw[1]", "theta[1,1]", "mu_dyad[1]", "mu_dyad_raw[1]", "theta_raw[1,1]")) {
    expect_equal(posterior::extract_variable(oriented, v), -extra[[v]], info = paste("expected", v, "to be negated"))
  }
  for (v in c("mu_intercept[1]", "sigma_mu", "sd_stat[1]", "rho", "within_between_ratio")) {
    expect_equal(posterior::extract_variable(oriented, v), extra[[v]], info = paste("expected", v, "to be unchanged"))
  }

  alpha1_theta11_before <- extra[["alpha[1]"]] * extra[["theta[1,1]"]]
  alpha1_theta11_after <- posterior::extract_variable(oriented, "alpha[1]") *
    posterior::extract_variable(oriented, "theta[1,1]")
  expect_equal(alpha1_theta11_after, alpha1_theta11_before)
})

test_that("bilatr_orient() does not flip when alpha[1]'s median is already positive", {
  set.seed(3)
  n <- 20
  extra <- list(
    `alpha[1]` = stats::rnorm(n, 2, 0.1),
    `theta[1,1]` = stats::rnorm(n, 0.5, 0.1)
  )
  draws <- .make_synthetic_draws(extra)
  oriented <- bilatr_orient(draws, stan_model = "stable", variables = c("alpha", "theta"))
  expect_equal(posterior::extract_variable(oriented, "alpha[1]"), extra[["alpha[1]"]])
  expect_equal(posterior::extract_variable(oriented, "theta[1,1]"), extra[["theta[1,1]"]])
})

test_that(".bilatr_flip_variables()/bilatr_orient() error on an unrecognized stan_model rather than silently no-op (B1)", {
  # B1: an unrecognized stan_model used to make .bilatr_flip_variables()
  # return character(0) -- indistinguishable from "this model has no
  # reflection symmetry" -- silently disabling orientation for a
  # wrong-basin fit with no error or warning. Both now go through
  # .canonical_stan_model() first, which errors instead.
  expect_error(.bilatr_flip_variables("some_unregistered_model"), "Unknown stan_model")

  set.seed(4)
  n <- 10
  extra <- list(`alpha[1]` = stats::rnorm(n, -3, 0.1), `alpha[2]` = stats::rnorm(n, 1, 0.1))
  draws <- .make_synthetic_draws(extra)
  expect_error(
    bilatr_orient(draws, stan_model = "some_unregistered_model", variables = "alpha"),
    "Unknown stan_model"
  )
})

test_that(".bilatr_flip_variables() accepts the pre-0.4.0 alphanorm/alphanorm_ou aliases and matches stable/ou exactly", {
  .reset_bilatr_alias_messaged()
  expect_message(
    alphanorm_flip <- .bilatr_flip_variables("alphanorm"),
    "pre-0.4.0 name of 'stable'"
  )
  expect_identical(alphanorm_flip, .bilatr_flip_variables("stable"))

  expect_message(
    alphanorm_ou_flip <- .bilatr_flip_variables("alphanorm_ou"),
    "pre-0.4.0 name of 'ou'"
  )
  expect_identical(alphanorm_ou_flip, .bilatr_flip_variables("ou"))
})

test_that("bilatr_orient() errors informatively if alpha[1] is not present in draws", {
  n <- 5
  draws <- .make_synthetic_draws(list(`theta[1,1]` = stats::rnorm(n)))
  expect_error(
    bilatr_orient(draws, stan_model = "stable", variables = "theta"),
    "alpha\\[1\\]"
  )
})

test_that(".warn_if_wrong_basin() fires when alpha[1] is negative, and bilatr_orient() recovers a positive orientation from that fit", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  D <- 3
  Tn <- 5
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
    compute_log_lik = 0, anchor_scale = 0.1
  )

  mod <- .compile_stan_model("stable", opt_level = 1)

  # deliberately seed into the WRONG basin (opposite sign convention from
  # .alpha_raw_sum0_init()), and pin the sampler near its init
  # (adapt_engaged = FALSE, tiny step_size, capped max_treedepth) so a
  # short run can't cross the (data-scale-dependent) likelihood barrier
  # between the two modes regardless of how large it happens to be for
  # this particular synthetic dataset -- this is what makes the test
  # deterministic rather than dependent on the barrier actually holding
  # for a tiny dataset.
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

  fit_wrong <- suppressWarnings(mod$sample(
    data = data_list, chains = 1, iter_warmup = 2, iter_sampling = 5,
    seed = 1, refresh = 0, threads_per_chain = 1,
    adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
    init = bad_init, output_dir = tempdir(), show_messages = FALSE
  ))
  expect_lt(stats::median(posterior::extract_variable(fit_wrong$draws("alpha[1]"), "alpha[1]")), 0)

  expect_warning(
    expect_message(.warn_if_wrong_basin(fit_wrong, "stable"), "Posterior median of alpha\\[1\\]"),
    "wrong-sign basin"
  )

  oriented <- bilatr_orient(fit_wrong$draws(variables = "alpha"), stan_model = "stable", variables = "alpha")
  expect_gt(stats::median(posterior::extract_variable(oriented, "alpha[1]")), 0)

  # and the right-basin case (bilatr_init_fn()'s actual, anchored init)
  # should report but not warn
  good_init <- bilatr_init_fn(list(D = D, T = Tn, A = A), stan_model = "stable")
  fit_right <- suppressWarnings(mod$sample(
    data = data_list, chains = 1, iter_warmup = 20, iter_sampling = 5,
    seed = 1, refresh = 0, threads_per_chain = 1,
    init = good_init, output_dir = tempdir(), show_messages = FALSE
  ))
  expect_no_warning(expect_message(.warn_if_wrong_basin(fit_right, "stable"), "Posterior median of alpha\\[1\\]"))
})

test_that(".warn_if_wrong_basin() errors on an unrecognized stan_model rather than silently no-op (B1)", {
  # both currently-registered models (stable, ou) DO have a reflection
  # symmetry since 0.4.0's promotion (see NEWS.md); .warn_if_wrong_basin()
  # now decides via .bilatr_flip_variables(), which itself validates
  # stan_model through .canonical_stan_model() -- an unrecognized name
  # errors rather than being silently treated as "no reflection
  # symmetry, nothing to check" (B1). In practice this stan_model has
  # already been validated by .compile_stan_model() earlier in
  # fit_bilatr(), so this case shouldn't arise from a real call, but the
  # guard must not paper over it if it somehow did.
  fake_fit <- list(draws = function(...) stop("should not be called"))
  expect_error(.warn_if_wrong_basin(fake_fit, "some_unregistered_model"), "Unknown stan_model")
})

test_that(".warn_if_wrong_basin() accepts the pre-0.4.0 alphanorm/alphanorm_ou aliases", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # a fake fit whose $draws() would error if ever called -- confirms
  # .warn_if_wrong_basin() treats the alias exactly like its canonical
  # name (both are flip-variable models, so this does NOT return early;
  # it must reach $draws())
  fake_fit <- list(draws = function(...) stop("should not be called"))
  expect_error(
    suppressMessages(.warn_if_wrong_basin(fake_fit, "alphanorm")),
    "should not be called"
  )
})
