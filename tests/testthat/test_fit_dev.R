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
