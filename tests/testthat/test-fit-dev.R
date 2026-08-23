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
    fit_dyad_ts_dev(fake_data, stan_model = "phi_logn"),
    "expects a single-dyad"
  )
})

test_that("fit_panel_dev() still validates D >= 2 for a registered stan_model", {
  fake_data <- list(D = 1, T = 3, A = 4)
  expect_error(
    fit_panel_dev(fake_data, stan_model = "phi_logn"),
    "expects multiple dyads"
  )
})

test_that("fit_dyad_ts_dev()/fit_panel_dev() default stan_model matches the exported functions' default", {
  expect_equal(formals(fit_dyad_ts_dev)$stan_model, quote(.BILATR_DEFAULT_MODEL))
  expect_equal(formals(fit_panel_dev)$stan_model, quote(.BILATR_DEFAULT_MODEL))
})

test_that("bilatr_init_fn() initializes phi for the stable model and log_phi0_raw/beta_logn for phi_logn", {
  stan_data <- list(D = 2, T = 3, A = 4)

  stable_init <- bilatr_init_fn(stan_data, stan_model = "stable")()
  expect_true("phi" %in% names(stable_init))
  expect_false("beta_logn" %in% names(stable_init))

  phi_logn_init <- bilatr_init_fn(stan_data, stan_model = "phi_logn")()
  expect_false("phi" %in% names(phi_logn_init))
  expect_true(all(c("log_phi0_raw", "beta_logn") %in% names(phi_logn_init)))
  expect_length(phi_logn_init$log_phi0_raw, stan_data$D)
})
