# --- synthetic CmdStan-shaped CSV, no CmdStan/cmdstan_path() needed -------
# --- (cmdstanr:::read_csv_metadata() is a pure text parser; confirmed it ---
# --- accepts a hand-written file with no real Stan run behind it) ---------
# --- .make_synthetic_stan_csv() itself now lives in helper_fixtures.R, ---
# --- shared with test_diagnose_convergence.R's peak-RSS regression guard -

test_that(".prepare_fast_csv_read()/.fast_read_post_warmup_draws() handle a wide synthetic CSV (no CmdStan)", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  n_cols <- 20000L
  n_draws <- 50L
  .make_synthetic_stan_csv(f, n_cols = n_cols, n_draws = n_draws, n_warmup = 3L, save_warmup = FALSE)

  prepared <- .prepare_fast_csv_read(f)
  expect_equal(prepared$num_post_warmup_draws, n_draws)
  expect_equal(prepared$n_chains, 1L)
  expect_length(prepared$variables, n_cols + 1L) # + lp__
  expect_true(all(c("lp__", "x[1]", paste0("x[", n_cols, "]")) %in% prepared$variables))
  expect_gt(prepared$file_mb, 0)

  draws <- .fast_read_post_warmup_draws(prepared, c("x[1]", "x[500]", paste0("x[", n_cols, "]")))
  expect_equal(dim(draws), c(n_draws, 1L, 3L))
  expect_equal(posterior::variables(draws), c("x[1]", "x[500]", paste0("x[", n_cols, "]")))
  # draw d, column c holds d * 1e6 + c -- distinct per (draw, column)
  # pair, so this actually exercises column identity, not just row count
  expect_equal(as.numeric(posterior::extract_variable(draws, "x[1]")), seq_len(n_draws) * 1e6 + 1)
  expect_equal(as.numeric(posterior::extract_variable(draws, "x[500]")), seq_len(n_draws) * 1e6 + 500)
  expect_equal(
    as.numeric(posterior::extract_variable(draws, paste0("x[", n_cols, "]"))),
    seq_len(n_draws) * 1e6 + n_cols
  )

  # a request whose variable order deliberately differs from file order --
  # the only CmdStan-free coverage of fread(select = )'s column-order
  # guarantee the whole col_index design rests on (see
  # .fast_read_post_warmup_draws()'s docs): if select() ever silently
  # returned columns in FILE order instead of REQUEST order, this would
  # catch it, since "x[3]"'s values would show up where "lp__"'s value is
  # expected and vice versa
  reordered <- .fast_read_post_warmup_draws(prepared, c(paste0("x[", n_cols, "]"), "x[3]", "lp__"))
  expect_equal(posterior::variables(reordered), c(paste0("x[", n_cols, "]"), "x[3]", "lp__"))
  expect_equal(
    as.numeric(posterior::extract_variable(reordered, paste0("x[", n_cols, "]"))),
    seq_len(n_draws) * 1e6 + n_cols
  )
  expect_equal(as.numeric(posterior::extract_variable(reordered, "x[3]")), seq_len(n_draws) * 1e6 + 3)
  expect_equal(as.numeric(posterior::extract_variable(reordered, "lp__")), rep(-1, n_draws))
})

test_that(".prepare_fast_csv_read() computes data_skip correctly with save_warmup = 1 (warmup rows before the adaptation block)", {
  f <- tempfile(fileext = ".csv")
  on.exit(unlink(f), add = TRUE)
  n_cols <- 200L
  n_draws <- 10L
  n_warmup <- 4L
  .make_synthetic_stan_csv(f, n_cols = n_cols, n_draws = n_draws, n_warmup = n_warmup, save_warmup = TRUE)

  prepared <- .prepare_fast_csv_read(f)
  expect_equal(prepared$num_post_warmup_draws, n_draws)

  draws <- .fast_read_post_warmup_draws(prepared, "x[1]")
  # exactly the post-warmup values (d * 1e6 + 1), never the -99 warmup rows
  expect_equal(as.numeric(draws), seq_len(n_draws) * 1e6 + 1)
})

# --- .fast_read_post_warmup_draws() matches cmdstanr::read_cmdstan_csv() --
# --- exactly, reading straight from the raw CSVs, no scratch copy ---------

test_that(".fast_read_post_warmup_draws() matches cmdstanr::read_cmdstan_csv() for exact and base-name variables", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()

  prepared <- .prepare_fast_csv_read(fx$csv_files)

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

test_that(".prepare_fast_csv_read()/.fast_read_post_warmup_draws() never create a file next to the CSVs", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  outdir <- dirname(fx$csv_files[1])
  before <- list.files(outdir, all.files = TRUE, no.. = TRUE)

  prepared <- .prepare_fast_csv_read(fx$csv_files)
  .fast_read_post_warmup_draws(prepared, "theta")

  after <- list.files(outdir, all.files = TRUE, no.. = TRUE)
  expect_equal(sort(before), sort(after))
})

test_that(".prepare_fast_csv_read() errors clearly on a chain with too few post-warmup rows (B7)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  truncated <- tempfile(fileext = ".csv")
  lines <- readLines(fx$csv_files[1])
  # drop the last two lines of the file: on a fixture with a short (non-#)
  # timing-footer-free tail this removes real post-warmup data rows,
  # simulating a chain that ended mid-sampling
  data_line_idx <- which(!startsWith(lines, "#"))
  writeLines(lines[-utils::tail(data_line_idx, 2)], truncated)

  # fread() itself also warns ("Stopped early on line ...") when its
  # single-column probe runs into the footer comments a few rows short
  # of num_post_warmup_draws -- expected collateral from the short read
  # this test deliberately creates, not something to assert on.
  expect_error(
    suppressWarnings(.prepare_fast_csv_read(c(truncated, fx$csv_files[-1]))),
    "post-warmup row"
  )
})

test_that(".fast_read_post_warmup_draws() errors clearly on an unknown variable", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture()
  prepared <- .prepare_fast_csv_read(fx$csv_files)

  expect_error(
    .fast_read_post_warmup_draws(prepared, "not_a_real_variable"),
    "not found in the CmdStan CSV header"
  )
})
