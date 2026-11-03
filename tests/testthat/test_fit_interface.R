# Tests for the fit interface: fit_panel()/fit_dyad_ts() take a
# dyad-period table or a prepared stan_data, run every data check before
# compiling, attach stan_data (recording the model on it) to the fit and
# write the sidecar .rds; the downstream functions default to the
# attached/sidecar stan_data.
#
# Cases that only need to reach the pre-compile checks (or the call into
# the sampler) replace .compile_stan_model() with a fake, so they need no
# CmdStan. Cases that need real draws use short single-chain fits.

# --- fixtures -------------------------------------------------------------

fit_events <- function() {
  suppressMessages(recode_cameo(make_fake_events(), code_col = "EventCode"))
}

fit_agg <- function(...) {
  suppressMessages(aggregate_events(fit_events(), grouping_var = "PentaClass", ...))
}

# A user table with no actor columns: arbitrary dyad labels, a week_id
# index with a gap, three action columns.
user_table <- function() {
  set.seed(3)
  tibble::tibble(
    dyad = rep(c("A-B", "A-C", "B-C", "C-D"), each = 6),
    week_id = rep(c(1L, 2L, 3L, 5L, 6L, 7L), times = 4),
    meetings = stats::rpois(24, 4),
    protests = stats::rpois(24, 2),
    sanctions = stats::rpois(24, 1)
  )
}

# The same, undirected, with actor columns but no w_send.
user_table_actors <- function() {
  ut <- user_table()
  ut$actor1 <- sub("-.*", "", ut$dyad)
  ut$actor2 <- sub(".*-", "", ut$dyad)
  ut$dyad <- NULL
  ut
}

# Replaces .compile_stan_model() for the rest of the calling test. The
# fake model's $sample() returns an environment (like a CmdStanMCMC)
# holding the data it was given. `calls$n` counts compiles.
local_fake_compile <- function(env = parent.frame()) {
  calls <- new.env()
  calls$n <- 0L
  local_mocked_bindings(
    .compile_stan_model = function(stan_model, opt_level, ...) {
      calls$n <- calls$n + 1L
      list(sample = function(data, ...) {
        fit <- new.env()
        fit$data <- data
        fit$args <- list(...)
        fit
      })
    },
    .package = "bilatr",
    .env = env
  )
  calls
}

# Fails the test if anything gets compiled.
local_no_compile <- function(env = parent.frame()) {
  local_mocked_bindings(
    .compile_stan_model = function(...) stop("compiled: the check did not run before compile"),
    .package = "bilatr",
    .env = env
  )
}

real_fit <- function(data, ..., fn = fit_panel) {
  suppressWarnings(suppressMessages(fn(
    data, ...,
    chains = 1, parallel_chains = 1, threads_per_chain = 1,
    iter_warmup = 50, iter_sampling = 50, seed = 1, opt_level = 1,
    refresh = 0, show_messages = FALSE
  )))
}

# What fit_*() attach: the prepared stan_data, plus the model it was fit
# with.
with_model <- function(stan_data, stan_model = "stable") {
  attr(stan_data, "stan_model") <- stan_model
  stan_data
}

# --- dispatch on `data` (fake compile) --------------------------------------

test_that("a data frame goes through prepare_stan_data(), and the result is attached", {
  calls <- local_fake_compile()
  agg <- fit_agg()
  fit <- fit_panel(agg, reference_category = "0", chunk_size = 3)
  expect_identical(calls$n, 1L)
  expected <- with_model(suppressMessages(prepare_stan_data(agg, reference_category = "0", chunk_size = 3)))
  expect_identical(attr(fit, "stan_data"), expected)
  expect_identical(fit$data, expected)
  expect_identical(attr(fit, "stan_model"), "stable")
})

test_that("fit_*() pass on prepare_stan_data()'s message when reference_category is NULL", {
  local_fake_compile()
  agg <- fit_agg()
  msgs <- testthat::capture_messages(fit_panel(agg, chunk_size = 3))
  expect_match(
    grep("reference_category not set", msgs, value = TRUE),
    "^prepare_stan_data\\(\\): reference_category not set; anchoring on the first action column, \"0\""
  )
  msgs <- testthat::capture_messages(fit_panel(agg, reference_category = "1", chunk_size = 3))
  expect_false(any(grepl("reference_category not set", msgs)))
})

test_that("tidyselect and character action_cols both reach prepare_stan_data() through fit_panel()", {
  local_fake_compile()
  ut <- user_table()
  expected <- with_model(prepare_stan_data(ut, time = "week_id", action_cols = c("meetings", "protests", "sanctions")))
  by_range <- fit_panel(ut, time = "week_id", action_cols = c(meetings:sanctions))
  by_helper <- fit_panel(ut, time = "week_id", action_cols = tidyselect::starts_with(c("m", "p", "s")))
  by_chr <- fit_panel(ut, time = "week_id", action_cols = c("meetings", "protests", "sanctions"))
  cols <- c("meetings", "protests", "sanctions")
  by_var <- fit_panel(ut, time = "week_id", action_cols = cols)
  for (fit in list(by_range, by_helper, by_chr, by_var)) {
    expect_identical(attr(fit, "stan_data"), expected)
  }
})

test_that("a tidyselect action_cols works through a wrapper that forwards `...`", {
  local_fake_compile()
  wrapper <- function(data, ...) fit_panel(data, time = "week_id", ...)
  fit <- wrapper(user_table(), action_cols = c(protests, meetings))
  expect_identical(attr(attr(fit, "stan_data"), "event_classes"), c("protests", "meetings"))
})

test_that("a prepared stan_data is used as-is, and sampler arguments still pass positionally by name", {
  calls <- local_fake_compile()
  sd <- suppressMessages(prepare_stan_data(fit_agg()))
  fit <- fit_panel(sd, chains = 1, iter_warmup = 10, refresh = 0)
  expect_identical(calls$n, 1L)
  expect_identical(attr(fit, "stan_data"), with_model(sd))
  expect_identical(fit$args$chains, 1)
  expect_identical(fit$args$iter_warmup, 10)
  expect_identical(fit$args$refresh, 0)
})

test_that("(iv) a prepared stan_data plus a data-preparation argument errors, naming it, before compile", {
  local_no_compile()
  sd <- suppressMessages(prepare_stan_data(fit_agg()))
  expect_error(
    fit_panel(sd, reference_category = "0"),
    "already prepared stan_data list, so the data-preparation argument\\(s\\) `reference_category` cannot be applied"
  )
  # several, and an argument passed with its default value still counts
  expect_error(
    fit_panel(sd, min_n_events = 1, chunk_size = 600, chains = 1),
    "`min_n_events`, `chunk_size` cannot be applied"
  )
  expect_error(fit_dyad_ts(sd, time = "year"), "`time` cannot be applied")
})

test_that("anything other than a data frame or stan_data list errors, saying what it got", {
  local_no_compile()
  expect_error(fit_panel(1:3), "got an object of class <integer>")
  expect_error(fit_panel(list(a = 1)), "got an object of class <list>")
  expect_error(fit_panel("events.csv"), "got an object of class <character>")
})

test_that("a stan_data-like list missing required fields errors before compile", {
  local_no_compile()
  expect_error(
    fit_panel(list(D = 3, T = 3, A = 4)),
    "missing field\\(s\\) `Y`, `is_obs`"
  )
})

# --- checks before compile (no compile) ------------------------------------

test_that("(v) fit_dyad_ts() needs D == 1 and fit_panel() needs D >= 2, before compile", {
  local_no_compile()
  agg <- fit_agg()
  expect_error(
    suppressMessages(fit_dyad_ts(agg)),
    "fit_dyad_ts\\(\\) expects a single-dyad stan_data \\(D == 1\\), got D = 6\\. Use fit_panel\\(\\)"
  )
  one <- dplyr::filter(agg, dyad == "USA_CHN")
  expect_error(
    suppressMessages(fit_panel(one)),
    "fit_panel\\(\\) expects multiple dyads \\(D >= 2\\), got D = 1\\. Use fit_dyad_ts\\(\\)"
  )
})

test_that("the D error reports the dyad count before min_n_events filtering", {
  local_no_compile()
  ut <- user_table()
  # C-D has events only in week 1
  ut[ut$dyad == "C-D", c("meetings", "protests", "sanctions")] <- 0L
  ut[ut$dyad == "C-D" & ut$week_id == 1L, "meetings"] <- 1L
  ut <- ut[ut$dyad %in% c("A-B", "C-D"), ]
  expect_error(
    suppressMessages(fit_panel(ut, time = "week_id", action_cols = c(meetings, protests, sanctions), min_n_events = 5)),
    "got D = 1 after `min_n_events`/`time_range` filtering \\(`data` has 2 dyad\\(s\\); `min_n_events` = 5\\)"
  )
})

test_that("(ii) a user table without actors errors under stable_gamma before compile", {
  local_no_compile()
  expect_error(
    suppressMessages(fit_panel(
      user_table(),
      time = "week_id", action_cols = c(meetings, protests, sanctions),
      stan_model = "stable_gamma"
    )),
    "stable_gamma\" needs the country fields `n_countries`, `ctry_a`, `ctry_b`, `w_send`.*Supply the `actor1`/`actor2` columns"
  )
})

test_that("an unknown stan_model errors before anything else", {
  local_no_compile()
  expect_error(fit_panel(1:3, stan_model = "nope"), "Unknown stan_model")
})

test_that("(iii) an undirected table without w_send warns under stable_gamma only", {
  local_fake_compile()
  ut <- user_table_actors()
  args <- list(ut, time = "week_id", action_cols = c("meetings", "protests", "sanctions"), directed = FALSE)
  expect_warning(
    suppressMessages(do.call(fit_panel, c(args, stan_model = "stable_gamma"))),
    "`w_send` was not supplied, so 0.5 is assumed for every undirected dyad.*aggregate_events\\(\\) computes"
  )
  expect_no_warning(suppressMessages(do.call(fit_panel, c(args, stan_model = "stable"))))
  expect_no_warning(suppressMessages(do.call(fit_panel, c(args, stan_model = "ou"))))
  # directed, or undirected with w_send: no warning
  expect_no_warning(suppressMessages(do.call(
    fit_panel, c(list(ut, time = "week_id", action_cols = c("meetings", "protests", "sanctions")), stan_model = "stable_gamma")
  )))
  ut$w_send <- 0.4
  expect_no_warning(suppressMessages(do.call(fit_panel, c(list(ut), args[-1], stan_model = "stable_gamma"))))
})

# --- sidecar -----------------------------------------------------------------

test_that("the sidecar is written atomically, named by output_basename or 'bilatr'", {
  local_fake_compile()
  sd <- suppressMessages(prepare_stan_data(fit_agg()))
  od <- withr::local_tempfile()
  fit_panel(sd, output_dir = od)
  expect_identical(list.files(od, all.files = TRUE, no.. = TRUE), "bilatr_stan_data.rds")
  expect_identical(readRDS(file.path(od, "bilatr_stan_data.rds")), with_model(sd))
  fit_panel(sd, output_dir = od, output_basename = "m16")
  expect_setequal(list.files(od, all.files = TRUE, no.. = TRUE), c("bilatr_stan_data.rds", "m16_stan_data.rds"))
  # no output_dir, no sidecar
  fit <- fit_panel(sd)
  expect_null(fit$args$output_dir)
})

test_that("the sidecar is written before sampling, so a crashed run leaves it", {
  local_mocked_bindings(
    .compile_stan_model = function(...) list(sample = function(...) stop("sampler crashed")),
    .package = "bilatr"
  )
  sd <- suppressMessages(prepare_stan_data(fit_agg()))
  od <- withr::local_tempfile()
  expect_error(fit_panel(sd, output_dir = od), "sampler crashed")
  expect_identical(readRDS(file.path(od, "bilatr_stan_data.rds")), with_model(sd))
})

# --- real fits ------------------------------------------------------------------

test_that("(i) fit_panel(aggregate_events(...)) samples; the attribute and the sidecar match prepare_stan_data()", {
  skip_on_cran()
  skip_if_no_cmdstan()
  agg <- fit_agg()
  od <- withr::local_tempfile()
  fit <- real_fit(agg, reference_category = "0", output_dir = od)
  expect_s3_class(fit, "CmdStanMCMC")
  expect_gt(posterior::ndraws(fit$draws("alpha")), 0)
  expected <- with_model(suppressMessages(prepare_stan_data(agg, reference_category = "0")))
  expect_identical(attr(fit, "stan_data"), expected)
  expect_identical(readRDS(file.path(od, "bilatr_stan_data.rds")), expected)

  # (vi) the attribute survives save_object() and saveRDS() round trips
  f1 <- withr::local_tempfile(fileext = ".rds")
  fit$save_object(f1)
  expect_identical(attr(readRDS(f1), "stan_data"), expected)
  expect_identical(attr(readRDS(f1), "stan_model"), "stable")
  f2 <- withr::local_tempfile(fileext = ".rds")
  saveRDS(fit, f2)
  expect_identical(attr(readRDS(f2), "stan_data"), expected)
})

test_that("(ii) a user table without actors fits with stable", {
  skip_on_cran()
  skip_if_no_cmdstan()
  fit <- real_fit(user_table(), time = "week_id", action_cols = c(meetings, protests, sanctions))
  expect_s3_class(fit, "CmdStanMCMC")
  theta <- extract_theta(fit)
  expect_identical(nrow(theta), 4L * 7L)
  expect_true(all(is.na(theta$actor1)))
})

test_that("fit_dyad_ts() fits a single-dyad table", {
  skip_on_cran()
  skip_if_no_cmdstan()
  one <- dplyr::filter(fit_agg(), dyad == "USA_CHN")
  fit <- real_fit(one, fn = fit_dyad_ts)
  expect_identical(attr(fit, "stan_data")$D, 1L)
  expect_identical(unique(extract_theta(fit)$dyad), "USA_CHN")
})

# --- downstream defaults ----------------------------------------------------------

test_that("downstream functions use the attached stan_data, and the sidecar next to CSVs", {
  skip_on_cran()
  skip_if_no_cmdstan()
  agg <- fit_agg()
  od <- withr::local_tempfile()
  fit <- real_fit(agg, reference_category = "0", output_dir = od)
  sd <- attr(fit, "stan_data")

  expect_identical(extract_theta(fit), extract_theta(fit, sd))
  expect_identical(extract_alpha(fit)$event_class, attr(sd, "event_classes"))
  expect_identical(extract_mu_intercept(fit)$event_class, attr(sd, "event_classes"))
  # an explicit event_classes still wins
  expect_identical(extract_alpha(fit, event_classes = letters[1:5])$event_class, letters[1:5])
  expect_identical(
    diagnose_category_merges(fit)$categories,
    diagnose_category_merges(fit, sd)$categories
  )
  expect_identical(
    check_compositional_residuals(fit, n_dyads = 4, n_draws = 20)$global,
    check_compositional_residuals(fit, sd, n_dyads = 4, n_draws = 20)$global
  )
  icc <- icc_curves(fit, theta_range = c(-1, 1), n_grid = 5)
  expect_identical(unique(icc$event_class), attr(sd, "event_classes"))

  # CSV paths: the sidecar is found, with a message naming it
  csv_files <- fit$output_files()
  expect_message(
    theta_csv <- muffle_short_fit_warnings(extract_theta(csv_files, max_memory_mb = 1024)),
    "bilatr_stan_data.rds"
  )
  expect_identical(theta_csv$dyad, extract_theta(fit)$dyad)
  n_dt <- tibble::tibble(dyad_id = seq_len(sd$D), n_dt = apply(sd$Y, 1, sum))
  # One chain of 50 draws, so posterior caps the ESS and warns.
  expect_message(
    res <- muffle_short_fit_warnings(diagnose_and_extract_bilatr(csv_files, n_dt = n_dt, max_memory_mb = 1024)),
    "bilatr_stan_data.rds"
  )
  expect_identical(res$alpha$event_class, attr(sd, "event_classes"))
  expect_identical(res$theta$dyad, theta_csv$dyad)

  # n_dt defaults to the per-dyad event totals of the sidecar's Stan data
  res_default <- muffle_short_fit_warnings(suppressMessages(diagnose_and_extract_bilatr(csv_files, max_memory_mb = 1024)))
  expect_identical(res_default$diagnostics, res$diagnostics)

  # two sidecars: ambiguous, so an error that says to pass stan_data
  file.copy(file.path(od, "bilatr_stan_data.rds"), file.path(od, "other_stan_data.rds"))
  expect_error(
    muffle_short_fit_warnings(suppressMessages(diagnose_and_extract_bilatr(csv_files, n_dt = n_dt, max_memory_mb = 1024))),
    "there are 2 \\*_stan_data.rds files.*Pass `stan_data`"
  )
  expect_no_error(muffle_short_fit_warnings(suppressMessages(
    diagnose_and_extract_bilatr(csv_files, sd, n_dt = n_dt, max_memory_mb = 1024)
  )))
  # icc_curves() stays soft: no error, the raw action indices
  expect_no_error(icc_curves(csv_files, theta_range = c(-1, 1), n_grid = 5))
})

test_that("downstream functions error helpfully when no stan_data can be found", {
  skip_on_cran()
  skip_if_no_cmdstan()
  # Its own output_dir, with the sidecar removed, rather than cmdstanr's
  # default tempdir(): no test writes a sidecar there any more, but
  # tempdir() is shared by the whole R session (every test file under
  # R CMD check, or an interactive session that ran fit_*(output_dir =
  # tempdir())), so "no sidecar next to these CSVs" is only a fact in a
  # directory the test owns.
  od <- withr::local_tempdir()
  fit <- real_fit(fit_agg(), output_dir = od)
  unlink(file.path(od, "bilatr_stan_data.rds"))
  attr(fit, "stan_data") <- NULL
  for (f in list(extract_theta, extract_gamma, diagnose_category_merges, check_compositional_residuals)) {
    err <- expect_error(f(fit), "`stan_data` was not supplied and could not be found")
    expect_match(conditionMessage(err), "prepare_stan_data()", fixed = TRUE)
    expect_match(conditionMessage(err), "attr(fit, \"stan_data\")", fixed = TRUE)
    expect_no_match(conditionMessage(err), "assemble_stan_data")
  }
  expect_false("event_class" %in% names(extract_alpha(fit)))
  # no sidecar next to these CSVs
  expect_error(extract_theta(fit$output_files()), "no \\*_stan_data.rds file next to the CSV files")
})

test_that("extract_gamma() and the gamma-aware checks use the attached stan_data of a stable_gamma fit", {
  skip_on_cran()
  skip_if_no_cmdstan()
  fit <- real_fit(fit_agg(), stan_model = "stable_gamma")
  sd <- attr(fit, "stan_data")
  gamma <- extract_gamma(fit)
  expect_s3_class(gamma, "bilatr_gamma")
  expect_identical(gamma, extract_gamma(fit, sd))
  expect_setequal(gamma$country_code, attr(sd, "country_codes"))
  expect_identical(
    check_compositional_residuals(fit, stan_model = "stable_gamma", n_dyads = 4, n_draws = 20)$global,
    check_compositional_residuals(fit, sd, stan_model = "stable_gamma", n_dyads = 4, n_draws = 20)$global
  )
  expect_no_error(icc_curves(fit, theta_range = c(-1, 1), n_grid = 5, stan_model = "stable_gamma", country = "USA"))
})

test_that("a prepared stan_data that records another model gets the fitted model, with a message", {
  local_fake_compile()
  sd <- with_model(suppressMessages(prepare_stan_data(fit_agg())), "ou")
  expect_message(fit <- fit_panel(sd, chains = 1), "records stan_model = \"ou\"; replacing it with \"stable\"")
  expect_identical(attr(attr(fit, "stan_data"), "stan_model"), "stable")
  # the same model: no message
  expect_no_message(fit_panel(attr(fit, "stan_data"), chains = 1))
})
