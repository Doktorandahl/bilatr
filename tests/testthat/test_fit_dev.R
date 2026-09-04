test_that("fit_dyad_ts_dev() rejects an unknown stan_model before touching Stan", {
  fake_data <- list(D = 1, T = 3, A = 4)
  expect_error(
    fit_dyad_ts_dev(fake_data, stan_model = "not_a_real_model"),
    "Unknown stan_model"
  )
})

test_that("fit_panel_dev() rejects an unknown stan_model before touching Stan", {
  fake_data <- list(D = 3, T = 3, A = 4)
  expect_error(
    fit_panel_dev(fake_data, stan_model = "not_a_real_model"),
    "Unknown stan_model"
  )
})

test_that("fit_dyad_ts_dev() still validates D == 1 for a registered stan_model", {
  fake_data <- list(D = 3, T = 3, A = 4)
  expect_error(
    fit_dyad_ts_dev(fake_data, stan_model = "stable"),
    "expects a single-dyad"
  )
})

test_that("fit_panel_dev() still validates D >= 2 for a registered stan_model", {
  fake_data <- list(D = 1, T = 3, A = 4)
  expect_error(
    fit_panel_dev(fake_data, stan_model = "stable"),
    "expects multiple dyads"
  )
})

test_that("fit_dyad_ts_dev()/fit_panel_dev() default stan_model matches the exported functions' default", {
  expect_equal(formals(fit_dyad_ts_dev)$stan_model, quote(.BILATR_DEFAULT_MODEL))
  expect_equal(formals(fit_panel_dev)$stan_model, quote(.BILATR_DEFAULT_MODEL))
})

test_that("bilatr_init_fn() initializes phi for the stable model, not phi_logn params", {
  stan_data <- list(D = 2, T = 3, A = 4)

  stable_init <- bilatr_init_fn(stan_data, stan_model = "stable")()
  expect_true("phi" %in% names(stable_init))
  expect_length(stable_init$phi, stan_data$D)
  # phi_logn-only params are gone
  expect_false(any(c("log_phi0_raw", "beta_logn") %in% names(stable_init)))
})

test_that("bilatr_init_fn() matches the current parameterization (non-centered process noise, free alpha)", {
  stan_data <- list(D = 2, T = 3, A = 4)
  init <- bilatr_init_fn(stan_data, stan_model = "stable")()

  # process_noise is non-centered: log_process_noise_raw, not process_noise
  expect_true("log_process_noise_raw" %in% names(init))
  expect_false("process_noise" %in% names(init))
  expect_length(init$log_process_noise_raw, stan_data$D)

  # alpha[A] is freely estimated: no alpha_hostile, alpha_raw has length A - 1
  expect_false("alpha_hostile" %in% names(init))
  expect_length(init$alpha_raw, stan_data$A - 1)
})

test_that(".alpha_raw_sum0_init() sums to exactly 0 and is bounded away from the dot_self()=0 degeneracy", {
  for (A in c(2, 4, 5, 9)) {
    v <- .alpha_raw_sum0_init(A)
    expect_length(v, A)
    expect_equal(sum(v), 0)
    expect_gt(sum(v^2), 0)
  }
})

test_that("bilatr_init_fn() builds correctly-shaped inits for each experimental variant", {
  stan_data <- list(D = 3, T = 5, A = 4)

  alphanorm_init <- bilatr_init_fn(stan_data, stan_model = "alphanorm")()
  expect_setequal(
    names(alphanorm_init),
    c(
      "theta_raw", "mu_intercept", "alpha_raw", "sigma_theta0", "z_theta0",
      "log_process_noise_raw", "mu_log_noise", "sigma_log_noise", "phi",
      "mu_log_phi", "sigma_log_phi"
    )
  )
  expect_false("mu_theta0" %in% names(alphanorm_init)) # hard-pinned, removed
  expect_length(alphanorm_init$mu_intercept, stan_data$A) # sum_to_zero_vector[A], not A - 1
  expect_length(alphanorm_init$alpha_raw, stan_data$A)
  expect_equal(sum(alphanorm_init$alpha_raw), 0)

  ou_init <- bilatr_init_fn(stan_data, stan_model = "ou")()
  expect_setequal(
    names(ou_init),
    c(
      "theta_raw", "mu_intercept_raw", "alpha_raw", "mu_theta_bar", "sigma_mu",
      "mu_dyad_raw", "rho", "mu_log_sd_stat", "sigma_log_sd_stat",
      "log_sd_stat_raw", "phi", "mu_log_phi", "sigma_log_phi"
    )
  )
  expect_false("mu_theta0" %in% names(ou_init)) # role taken by mu_theta_bar
  expect_length(ou_init$mu_intercept_raw, stan_data$A - 1) # identification unchanged from stable
  expect_length(ou_init$mu_dyad_raw, stan_data$D)
  expect_true(ou_init$rho > 0 && ou_init$rho < 1)

  alphanorm_ou_init <- bilatr_init_fn(stan_data, stan_model = "alphanorm_ou")()
  expect_setequal(
    names(alphanorm_ou_init),
    c(
      "theta_raw", "mu_intercept", "alpha_raw", "sigma_mu", "mu_dyad_raw",
      "rho", "mu_log_sd_stat", "sigma_log_sd_stat", "log_sd_stat_raw", "phi",
      "mu_log_phi", "sigma_log_phi"
    )
  )
  expect_false(any(c("mu_theta0", "mu_theta_bar") %in% names(alphanorm_ou_init))) # location pinned hard
  expect_length(alphanorm_ou_init$mu_intercept, stan_data$A)
  expect_equal(sum(alphanorm_ou_init$alpha_raw), 0)
})
