# Compiled models live in a user cache directory, never in the installed
# package: .bilatr_model_dir(), the .stan copy, compile_bilatr_model() and
# bilatr_clear_model_cache().

test_that(".bilatr_model_dir() defaults to a per-version directory under R_user_dir()", {
  withr::local_options(bilatr.model_dir = NULL)
  dir <- .bilatr_model_dir(create = FALSE)
  expect_identical(
    dir,
    file.path(tools::R_user_dir("bilatr", "cache"), "models", as.character(utils::packageVersion("bilatr")))
  )
})

test_that(".bilatr_model_dir() uses options(bilatr.model_dir = ) and creates it", {
  dir <- file.path(withr::local_tempdir(), "a dir with spaces", "models")
  withr::local_options(bilatr.model_dir = dir)
  expect_false(dir.exists(dir))
  expect_identical(.bilatr_model_dir(), dir)
  expect_true(dir.exists(dir))

  withr::local_options(bilatr.model_dir = c("a", "b"))
  expect_error(.bilatr_model_dir(), "must be a single directory path")
})

test_that("the .stan file is copied only when missing or changed", {
  dir <- withr::local_tempdir(pattern = "model dir ")
  src <- .resolve_stan_model("stable")
  target <- .bilatr_stage_stan_file("stable", dir)
  expect_identical(target, file.path(dir, "bilatr_alphanorm.stan"))
  expect_identical(unname(tools::md5sum(target)), unname(tools::md5sum(src)))

  # unchanged: not rewritten, so its modification time stays put
  old_time <- as.POSIXct("2020-01-01", tz = "UTC")
  Sys.setFileTime(target, old_time)
  .bilatr_stage_stan_file("stable", dir)
  expect_equal(file.mtime(target), old_time, ignore_attr = TRUE)

  # changed (or stale): copied again, and the in-session cache entry for
  # the old file is dropped
  cat("// edited\n", file = target, append = TRUE)
  key <- paste(target, 3, TRUE, sep = "::")
  assign(key, "stale model", envir = .bilatr_model_cache)
  withr::defer(rm(list = intersect(key, ls(.bilatr_model_cache)), envir = .bilatr_model_cache))
  .bilatr_stage_stan_file("stable", dir)
  expect_identical(unname(tools::md5sum(target)), unname(tools::md5sum(src)))
  expect_false(exists(key, envir = .bilatr_model_cache, inherits = FALSE))
})

test_that("compile_bilatr_model() validates stan_model before anything else", {
  expect_error(compile_bilatr_model("nope"), "Unknown stan_model 'nope'")
  expect_error(compile_bilatr_model(c("stable", "ou")), "single string")
})

test_that(".require_cmdstan() requires the cmdstanr version DESCRIPTION names", {
  suggests <- read.dcf(system.file("DESCRIPTION", package = "bilatr"), fields = "Suggests")[1, 1]
  expect_match(suggests, paste0("cmdstanr \\(>= ", gsub(".", "\\.", .BILATR_MIN_CMDSTANR, fixed = TRUE), "\\)"))

  skip_if_not_installed("cmdstanr")
  local_mocked_bindings(.BILATR_MIN_CMDSTANR = "999.0.0")
  expect_error(.require_cmdstan("compile_bilatr_model()"), "cmdstanr.*999\\.0\\.0")
})

test_that(".require_cmdstan() names install_cmdstan() when CmdStan is missing", {
  skip_if_not_installed("cmdstanr")
  local_mocked_bindings(cmdstan_version = function(...) NULL, .package = "cmdstanr")
  expect_error(.require_cmdstan("compile_bilatr_model()"), "compile_bilatr_model\\(\\): CmdStan is not installed.*cmdstanr::install_cmdstan\\(\\)")
  # and the compile path checks it before doing anything
  expect_error(compile_bilatr_model(), "cmdstanr::install_cmdstan\\(\\)")
})

test_that("bilatr_cmdstan_available() is FALSE, without erroring, whatever is missing", {
  # cmdstanr too old (or, here also without cmdstanr, absent)
  local({
    local_mocked_bindings(.BILATR_MIN_CMDSTANR = "999.0.0")
    expect_false(bilatr_cmdstan_available())
    expect_identical(.cmdstan_problem(), "cmdstanr")
  })

  skip_if_not_installed("cmdstanr")
  local({
    local_mocked_bindings(cmdstan_version = function(...) NULL, .package = "cmdstanr")
    expect_false(bilatr_cmdstan_available())
    expect_identical(.cmdstan_problem(), "cmdstan")
  })
  local({
    local_mocked_bindings(cmdstan_version = function(...) stop("boom"), .package = "cmdstanr")
    expect_false(bilatr_cmdstan_available())
  })
  local({
    local_mocked_bindings(
      cmdstan_version = function(...) {
        warning("no path")
        NULL
      },
      .package = "cmdstanr"
    )
    expect_no_warning(expect_false(bilatr_cmdstan_available()))
  })
})

test_that("bilatr_cmdstan_available() agrees with whether CmdStan runs", {
  has_cmdstan <- rlang::is_installed("cmdstanr", version = .BILATR_MIN_CMDSTANR) &&
    !is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL))
  expect_identical(bilatr_cmdstan_available(), has_cmdstan)
  if (has_cmdstan) {
    expect_true(.require_cmdstan("test"))
  }
})

test_that("bilatr_clear_model_cache() deletes only bilatr's files in a bilatr.model_dir, and empties the session cache", {
  dir <- withr::local_tempdir(pattern = "model dir ")
  withr::local_options(bilatr.model_dir = dir)
  ours <- file.path(dir, c("bilatr_alphanorm.stan", "bilatr_alphanorm", "bilatr_alphanorm_gamma.stan", "bilatr_alphanorm_ou.exe"))
  theirs <- file.path(dir, c("notes.txt", "bilatr_alphanorm_old.stan"))
  file.create(c(ours, theirs))
  assign("some::key", "model", envir = .bilatr_model_cache)

  deleted <- withVisible(bilatr_clear_model_cache())
  expect_false(deleted$visible)
  expect_setequal(deleted$value, ours)
  expect_false(any(file.exists(ours)))
  expect_true(all(file.exists(theirs)))
  expect_length(ls(.bilatr_model_cache), 0)

  # nothing left to delete
  expect_identical(bilatr_clear_model_cache(), character(0))
  expect_error(bilatr_clear_model_cache(all = NA), "`all` must be TRUE or FALSE")
})

test_that("bilatr_clear_model_cache() deletes the per-version directory, or all versions with all = TRUE", {
  cache_root <- withr::local_tempdir(pattern = "user cache ")
  withr::local_options(bilatr.model_dir = NULL)
  local_mocked_bindings(.bilatr_cache_root = function() file.path(cache_root, "bilatr"))
  this_version <- .bilatr_model_dir()
  other_version <- file.path(dirname(this_version), "0.0.1")
  dir.create(other_version)
  file.create(file.path(c(this_version, other_version), "bilatr_alphanorm"))

  expect_identical(bilatr_clear_model_cache(), this_version)
  expect_false(dir.exists(this_version))
  expect_true(dir.exists(other_version))

  expect_identical(bilatr_clear_model_cache(all = TRUE), dirname(this_version))
  expect_false(dir.exists(dirname(this_version)))
})

test_that("compiling writes nothing into the installed package, and works from a directory with a space in its path", {
  skip_on_cran()
  skip_if_no_cmdstan()
  dir <- file.path(withr::local_tempdir(pattern = "bilatr cache "), "compiled models")
  withr::local_options(bilatr.model_dir = dir)

  installed <- system.file("stan", package = "bilatr")
  snapshot <- function() {
    files <- list.files(installed, recursive = TRUE, all.files = TRUE, full.names = TRUE)
    info <- file.info(files)[c("size", "mtime")]
    info[order(rownames(info)), ]
  }
  before <- snapshot()
  mod <- compile_bilatr_model(opt_level = 1)
  expect_identical(snapshot(), before)

  expect_identical(normalizePath(mod$stan_file()), normalizePath(file.path(dir, "bilatr_alphanorm.stan")))
  expect_true(file.exists(mod$exe_file()))
  expect_identical(normalizePath(dirname(mod$exe_file())), normalizePath(dir))

  # the same object again, from the session cache
  expect_identical(compile_bilatr_model(opt_level = 1), mod)

  # sampling from the spaced path works end to end
  agg <- suppressMessages(aggregate_events(
    recode_cameo(make_fake_events(), code_col = "EventCode"),
    grouping_var = "PentaClass"
  ))
  one <- dplyr::filter(agg, dyad == "USA_CHN")
  fit <- suppressWarnings(suppressMessages(fit_dyad_ts(
    one, chains = 1, iter_warmup = 20, iter_sampling = 20, seed = 1, opt_level = 1,
    refresh = 0, show_messages = FALSE
  )))
  expect_s3_class(fit, "CmdStanMCMC")
  expect_identical(snapshot(), before)

  bilatr_clear_model_cache()
  expect_false(file.exists(file.path(dir, "bilatr_alphanorm")))
})
