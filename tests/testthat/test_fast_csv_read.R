# --- .fast_read_post_warmup_draws() matches cmdstanr::read_cmdstan_csv() --
# --- exactly, and .prepare_fast_csv_read() cleans up after itself ---------

test_that(".fast_read_post_warmup_draws() matches cmdstanr::read_cmdstan_csv() for exact and base-name variables", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()

  prepared <- .prepare_fast_csv_read(fx$csv_files)
  on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)

  # exact indexed name
  ref <- cmdstanr::read_cmdstan_csv(fx$csv_files, variables = "theta[3,2]")$post_warmup_draws
  got <- .fast_read_post_warmup_draws(prepared, "theta[3,2]")
  expect_equal(as.numeric(got), as.numeric(ref))
  expect_equal(posterior::variables(got), posterior::variables(ref))

  # bare base name, expands to every indexed element (cmdstanr's own
  # matching_variables() prefix-fallback behavior)
  ref_alpha <- cmdstanr::read_cmdstan_csv(fx$csv_files, variables = "alpha")$post_warmup_draws
  got_alpha <- .fast_read_post_warmup_draws(prepared, "alpha")
  expect_equal(posterior::variables(got_alpha), posterior::variables(ref_alpha))
  expect_equal(as.numeric(got_alpha), as.numeric(ref_alpha))

  # multi-variable, mixed exact/base-name request
  ref_mix <- cmdstanr::read_cmdstan_csv(fx$csv_files, variables = c("alpha[1]", "mu_intercept"))$post_warmup_draws
  got_mix <- .fast_read_post_warmup_draws(prepared, c("alpha[1]", "mu_intercept"))
  expect_equal(posterior::variables(got_mix), posterior::variables(ref_mix))
  expect_equal(as.numeric(got_mix), as.numeric(ref_mix))
})

test_that(".fast_read_post_warmup_draws() flip = TRUE negates draws before assembly", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  prepared <- .prepare_fast_csv_read(fx$csv_files)
  on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)

  plain <- .fast_read_post_warmup_draws(prepared, "alpha[1]")
  flipped <- .fast_read_post_warmup_draws(prepared, "alpha[1]", flip = TRUE)
  expect_equal(as.numeric(flipped), -as.numeric(plain))
})

test_that(".prepare_fast_csv_read()/.cleanup_fast_csv_read() clean up their scratch files", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  prepared <- .prepare_fast_csv_read(fx$csv_files)

  expect_true(all(file.exists(prepared$clean_files)))
  expect_false(any(prepared$clean_files %in% fx$csv_files))

  .cleanup_fast_csv_read(prepared)
  expect_false(any(file.exists(prepared$clean_files)))
})

test_that(".prepare_fast_csv_read() respects an explicit scratch_dir", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  scratch <- tempfile()
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)

  prepared <- .prepare_fast_csv_read(fx$csv_files, scratch_dir = scratch)
  on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)

  expect_true(all(dirname(prepared$clean_files) == scratch))
  expect_true(all(file.exists(prepared$clean_files)))
})

test_that(".fast_read_post_warmup_draws() errors clearly on an unknown variable", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  prepared <- .prepare_fast_csv_read(fx$csv_files)
  on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)

  expect_error(
    .fast_read_post_warmup_draws(prepared, "not_a_real_variable"),
    "not found in the CmdStan CSV header"
  )
})
