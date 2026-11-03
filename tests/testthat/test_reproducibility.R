# Per-chain, seeded initial values (a fit with a `seed` is reproducible),
# and the fitted model recorded on the fit and inferred downstream.

repro_stan_data <- function(...) {
  events <- suppressMessages(recode_cameo(make_fake_events(), code_col = "EventCode"))
  agg <- suppressMessages(aggregate_events(events, grouping_var = "PentaClass"))
  suppressMessages(prepare_stan_data(agg, ...))
}

repro_one_dyad <- function() {
  events <- suppressMessages(recode_cameo(make_fake_events(), code_col = "EventCode"))
  agg <- suppressMessages(aggregate_events(events, grouping_var = "PentaClass"))
  dplyr::filter(agg, dyad == "USA_CHN")
}

# --- initial values, no CmdStan ---------------------------------------------

test_that("bilatr_init_fn() returns a function that draws fresh values on every call", {
  sd <- repro_stan_data()
  random_parts <- list(
    stable = c("alpha_raw", "z_theta0"),
    stable_gamma = c("alpha_raw", "z_theta0", "gamma_z"),
    ou = c("alpha_raw", "mu_dyad_raw")
  )
  for (model in names(random_parts)) {
    init_fn <- bilatr_init_fn(sd, model)
    a <- init_fn()
    b <- init_fn()
    for (par in random_parts[[model]]) {
      expect_false(identical(a[[par]], b[[par]]), label = paste(model, par))
    }
    # everything else is fixed
    fixed <- setdiff(names(a), random_parts[[model]])
    expect_identical(a[fixed], b[fixed])
    # alpha_raw sums to zero with a positive first element
    expect_equal(sum(a$alpha_raw), 0)
    expect_gt(a$alpha_raw[1], 0)
  }
})

test_that(".bilatr_draw_inits() gives each chain its own values, reproducibly with a seed", {
  sd <- repro_stan_data()
  inits <- .bilatr_draw_inits(sd, "stable", chains = 3, seed = 1)
  expect_length(inits, 3)
  expect_false(identical(inits[[1]]$alpha_raw, inits[[2]]$alpha_raw))
  expect_false(identical(inits[[2]]$z_theta0, inits[[3]]$z_theta0))
  expect_identical(.bilatr_draw_inits(sd, "stable", chains = 3, seed = 1), inits)
  # chain k is drawn under seed + k
  expect_identical(inits[[2]], .with_seed(1 + 2, bilatr_init_fn(sd, "stable")()))
  # with chain_ids, under seed + chain_ids[k]
  by_id <- .bilatr_draw_inits(sd, "stable", chains = 2, seed = 1, chain_ids = c(2L, 3L))
  expect_identical(by_id, inits[2:3])
  expect_error(.bilatr_draw_inits(sd, "stable", chains = 3, seed = 1, chain_ids = 1:2), "`chain_ids` has length 2 but `chains` is 3")
})

test_that("seeded inits leave the caller's random-number stream alone; seed = NULL draws from it", {
  sd <- repro_stan_data()
  set.seed(99)
  before <- .Random.seed
  .bilatr_draw_inits(sd, "stable", chains = 2, seed = 1)
  expect_identical(.Random.seed, before)

  set.seed(5)
  a <- .bilatr_draw_inits(sd, "stable", chains = 2, seed = NULL)
  set.seed(5)
  b <- .bilatr_draw_inits(sd, "stable", chains = 2, seed = NULL)
  expect_identical(a, b)
  expect_false(identical(a[[1]], a[[2]]))
})

# A fake compile that keeps what $sample() was given.
local_capture_sample <- function(env = parent.frame()) {
  local_mocked_bindings(
    .compile_stan_model = function(...) {
      list(sample = function(...) {
        fit <- new.env()
        fit$args <- list(...)
        fit
      })
    },
    .package = "bilatr",
    .env = env
  )
}

test_that("fit_*() pass one init list per chain, and a user-supplied init as given", {
  local_capture_sample()
  sd <- repro_stan_data()
  fit <- fit_panel(sd, chains = 3, seed = 4)
  expect_identical(fit$args$init, .bilatr_draw_inits(sd, "stable", chains = 3, seed = 4))
  expect_identical(fit$args$seed, 4)

  mine <- list(list(alpha_raw = c(1, -1, 0, 0, 0)))
  fit <- fit_panel(sd, chains = 1, init = mine)
  expect_identical(fit$args$init, mine)
  fit <- fit_panel(sd, chains = 2, init = 0.5)
  expect_identical(fit$args$init, 0.5)

  # chain_ids through `...`
  fit <- fit_panel(sd, chains = 2, seed = 4, chain_ids = c(5L, 6L))
  expect_identical(fit$args$init, .bilatr_draw_inits(sd, "stable", chains = 2, seed = 4, chain_ids = c(5L, 6L)))
  expect_identical(fit$args$chain_ids, c(5L, 6L))
})

# --- real fits ----------------------------------------------------------------

tiny_fit <- function(data, ..., fn = fit_dyad_ts) {
  suppressWarnings(suppressMessages(fn(
    data, ...,
    iter_warmup = 20, iter_sampling = 20, opt_level = 1,
    refresh = 0, show_messages = FALSE
  )))
}

test_that("two fits with the same seed give identical draws, from different per-chain inits", {
  skip_on_cran()
  skip_if_no_cmdstan()
  one <- repro_one_dyad()
  fit1 <- tiny_fit(one, seed = 1, chains = 2)
  fit2 <- tiny_fit(one, seed = 1, chains = 2)
  expect_identical(fit1$draws(), fit2$draws())

  # the inits CmdStan was given: different across chains, the same across fits
  inits <- fit1$init()
  expect_length(inits, 2)
  expect_false(isTRUE(all.equal(inits[[1]]$alpha_raw, inits[[2]]$alpha_raw)))
  expect_false(isTRUE(all.equal(inits[[1]]$z_theta0, inits[[2]]$z_theta0)))
  expect_equal(fit2$init(), inits)

  # a different seed gives different draws
  fit3 <- tiny_fit(one, seed = 2, chains = 2)
  expect_false(identical(fit3$draws(), fit1$draws()))
})

test_that("seed = NULL still fits, and a user-supplied init is used as given", {
  skip_on_cran()
  skip_if_no_cmdstan()
  one <- repro_one_dyad()
  fit <- tiny_fit(one, seed = NULL, chains = 1)
  expect_s3_class(fit, "CmdStanMCMC")

  sd <- attr(fit, "stan_data")
  mine <- bilatr_init_fn(sd, "stable")()
  mine$alpha_raw <- c(0.4, -0.1, -0.1, -0.1, -0.1)[seq_len(sd$A)]
  mine$alpha_raw <- mine$alpha_raw - mean(mine$alpha_raw)
  fit_mine <- tiny_fit(sd, seed = 3, chains = 1, init = list(mine))
  expect_equal(fit_mine$init()[[1]]$alpha_raw, mine$alpha_raw, tolerance = 1e-6)
})

test_that("a stable_gamma fit is checked as stable_gamma without stan_model, from the fit, the CSVs, or the sidecar", {
  skip_on_cran()
  skip_if_no_cmdstan()
  events <- suppressMessages(recode_cameo(make_fake_events(), code_col = "EventCode"))
  agg <- suppressMessages(aggregate_events(events, grouping_var = "PentaClass"))
  od <- withr::local_tempdir()
  fit <- tiny_fit(agg, reference_category = "0", stan_model = "stable_gamma", seed = 1, chains = 2, output_dir = od, fn = fit_panel)
  expect_identical(attr(fit, "stan_model"), "stable_gamma")
  expect_identical(attr(attr(fit, "stan_data"), "stan_model"), "stable_gamma")
  expect_identical(attr(readRDS(file.path(od, "bilatr_stan_data.rds")), "stan_model"), "stable_gamma")

  # inferred: the residual check subtracts g_d exactly as with an explicit stan_model
  inferred <- check_compositional_residuals(fit, n_dyads = 4, n_draws = 20)
  explicit <- check_compositional_residuals(fit, stan_model = "stable_gamma", n_dyads = 4, n_draws = 20)
  expect_identical(inferred, explicit)
  expect_identical(inferred$settings$stan_model, "stable_gamma")
  # ... and that matters: checked as "stable" (here a bare draws list with
  # nothing recorded), g_d is left out and the statistics differ
  bare_fit <- list(draws = function(variables) fit$draws(variables))
  bare_sd <- attr(fit, "stan_data")
  attr(bare_sd, "stan_model") <- NULL
  as_stable <- check_compositional_residuals(bare_fit, bare_sd, stan_model = "stable", n_dyads = 4, n_draws = 20)
  expect_false(identical(as_stable$global, inferred$global))

  # the diagnostics split gamma off, and icc_curves() accepts `country`
  diag <- suppressWarnings(diagnose_convergence(fit, tiers = 1))
  expect_false(is.null(diag$gamma))
  expect_false(any(startsWith(diag$tier1$variable, "gamma[")))
  expect_no_warning(icc_curves(fit, theta_range = c(-1, 1), n_grid = 5, country = "USA"))

  # a contradiction errors, naming both
  expect_error(
    check_compositional_residuals(fit, stan_model = "stable", n_dyads = 4, n_draws = 20),
    "`stan_model = \"stable\"` contradicts the model recorded in the fit \\(\"stable_gamma\"\\)"
  )
  expect_error(diagnose_convergence(fit, stan_model = "ou", tiers = 1), "contradicts")

  # CSV paths alone, without the sidecar: read from the CSV header
  csv_files <- fit$output_files()
  unlink(file.path(od, "bilatr_stan_data.rds"))
  expect_identical(.stan_model_from_csv_header(csv_files[1]), "stable_gamma")
  expect_identical(.resolve_stan_model_for_fit(csv_files, NULL, NULL, "f()"), "stable_gamma")
  diag_csv <- suppressWarnings(suppressMessages(diagnose_convergence(csv_files, tiers = 1)))
  expect_false(is.null(diag_csv$gamma))
  expect_error(.resolve_stan_model_for_fit(csv_files, NULL, "stable", "f()"), "recorded in the CmdStan CSV header \\(\"stable_gamma\"\\)")
})

test_that("with nothing recorded, the model falls back to \"stable\" with a message", {
  truth <- .make_fake_residual_data(D = 6, Tn = 4, A = 4, seed = 1)
  fake_fit <- .make_fake_residual_fit(truth, n_pseudo_draws = 30)
  expect_message(
    res <- check_compositional_residuals(fake_fit, truth$stan_data, n_dyads = 6, n_draws = 10),
    "no model is recorded on the fit, so `stan_model = \"stable\"` is assumed"
  )
  expect_identical(res$settings$stan_model, "stable")
  expect_message(out <- .resolve_stan_model_for_fit(list(), NULL, NULL, "f()"), "f\\(\\): no model is recorded")
  expect_identical(out, "stable")
  # an explicit stan_model needs nothing recorded, and gives no message
  expect_no_message(.resolve_stan_model_for_fit(list(), NULL, "ou", "f()"))
  # recorded on the Stan data only
  sd <- truth$stan_data
  attr(sd, "stan_model") <- "ou"
  expect_identical(.resolve_stan_model_for_fit(fake_fit, sd, NULL, "f()"), "ou")
  expect_error(.resolve_stan_model_for_fit(fake_fit, sd, "stable", "f()"), "recorded in `stan_data` \\(\"ou\"\\)")
})

test_that(".stan_model_from_csv_header() maps CmdStan model names through the registry", {
  path <- withr::local_tempfile(fileext = ".csv")
  .make_synthetic_stan_csv(path, n_cols = 2, n_draws = 3)
  # the synthetic CSV's model is "synthetic_model": not registered
  expect_null(.stan_model_from_csv_header(path))
  expect_identical(.scan_stan_csv_header_and_skip(path)$model, "synthetic_model")
  lines <- readLines(path)
  for (model in names(.bilatr_stan_models)) {
    cmdstan_name <- paste0(tools::file_path_sans_ext(.bilatr_stan_models[[model]]$file), "_model")
    writeLines(sub("^# model = .*$", paste("# model =", cmdstan_name), lines), path)
    expect_identical(.stan_model_from_csv_header(path), model)
  }
  # no model line at all
  writeLines(lines[!startsWith(lines, "# model =")], path)
  expect_null(.stan_model_from_csv_header(path))
  expect_true(is.na(.scan_stan_csv_header_and_skip(path)$model))
})
