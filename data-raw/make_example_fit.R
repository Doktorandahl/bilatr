# Builds the tiny example fits shipped in inst/extdata/example_fit/, which
# bilatr_example_fit() returns. They exist so that every post-fit function
# has an example (and CRAN-visible tests) that run without CmdStan. They
# are far too short to be a valid inference.
#
#   inst/extdata/example_fit/stable/        stable model, with
#                                           compute_theta_filtered = 1 and
#                                           compute_log_lik = 1
#   inst/extdata/example_fit/stable_gamma/  stable_gamma model
#
# Each directory holds the CmdStan CSVs (example-1.csv, example-2.csv) and
# the Stan data sidecar (example_stan_data.rds) that fit_panel() writes.
# The CSV comment headers are scrubbed: every absolute path is replaced by
# its basename, and every other byte is left as CmdStan wrote it.
#
# Needs cmdstanr and CmdStan. Run from the package root:
#   Rscript data-raw/make_example_fit.R

library(dplyr)
devtools::load_all(".", quiet = TRUE)

out_root <- file.path("inst", "extdata", "example_fit")
# Size reductions, applied in this order until inst/extdata/ is below
# 400 KB (see the end of the script): 100 -> 60 draws per chain, then
# 4 -> 3 countries, then 4 significant figures in the CSVs.
iter_sampling <- 60
countries <- c("USA", "CHN", "RUS")
sig_figs <- 4

# -- 1. Simulated events --------------------------------------------------
#
# 4 countries, directed, 2015-2020. Each dyad has a latent cooperation
# level that follows a random walk over the years; an event's CAMEO code
# is drawn with probability proportional to exp(level * Goldstein / 4), so
# more cooperative dyad-years produce more cooperative codes and the fit
# has a real signal to find. Dyads differ in how active they are.

set.seed(20261103)
years <- 2015:2020
dyads <- expand.grid(a1 = countries, a2 = countries, stringsAsFactors = FALSE) %>%
  filter(a1 != a2) %>%
  mutate(rate = round(exp(runif(n(), log(8), log(60)))))

codes <- cameo_lookup$CAMEOEVENTCODE
goldstein <- cameo_lookup$GoldsteinScore

events <- purrr::pmap_dfr(dyads, function(a1, a2, rate) {
  level <- cumsum(c(rnorm(1, 0, 1), rnorm(length(years) - 1, 0, 0.4)))
  purrr::map2_dfr(years, level, function(year, lev) {
    n <- rpois(1, rate)
    if (n == 0) {
      return(NULL)
    }
    days <- seq(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")), by = "day")
    tibble::tibble(
      Actor1CountryCode = a1,
      Actor2CountryCode = a2,
      SQLDATE = as.integer(format(sample(days, n, replace = TRUE), "%Y%m%d")),
      EventCode = sample(codes, n, replace = TRUE, prob = exp(lev * goldstein / 4))
    )
  })
})
events <- recode_cameo(events, code_col = "EventCode")

agg <- aggregate_events(
  events,
  grouping_var = "PentaClass",
  resolution = "yearly",
  start = "2015-01-01",
  end = "2020-12-31"
)
message(
  nrow(events), " events, ", length(unique(agg$dyad)), " dyads, ",
  nrow(agg), " dyad-years"
)

# -- 2. Fits ----------------------------------------------------------------

fit_example <- function(stan_model, ...) {
  dir <- file.path(out_root, stan_model)
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  fit <- fit_panel(
    agg,
    reference_category = "1", # verbal cooperation
    stan_model = stan_model,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 150,
    iter_sampling = iter_sampling,
    seed = 1,
    output_dir = dir,
    output_basename = "example",
    refresh = 0,
    sig_figs = sig_figs,
    ...
  )
  # CmdStan 2.33+ also writes *_config.json and *_metric.json, which
  # bilatr does not read and which record absolute paths
  unlink(list.files(dir, pattern = "_(config|metric)\\.json$", full.names = TRUE))
  sort(fit$output_files())
}

csv <- list(
  stable = fit_example("stable", compute_theta_filtered = 1, compute_log_lik = 1),
  stable_gamma = fit_example("stable_gamma")
)

# -- 3. Scrub the CSV comment headers ---------------------------------------
#
# CmdStan records absolute paths in the comment header (`stancflags =
# --include-paths=...`, `data file = ...`, `output file = ...`,
# `init = ...`, and the `profile_file`/`diagnostic_file` lines when set).
# Each path value is replaced by its basename. Every other line, and every
# other character of the edited lines, is left byte-identical.

scrub_header_line <- function(line) {
  if (!startsWith(line, "#")) {
    return(line)
  }
  # `#   <key> = <value>` with an absolute path as the value
  m <- regmatches(line, regexec("^(#\\s+[A-Za-z_ ]+ = )(/.*)$", line))[[1]]
  if (length(m) == 3) {
    return(paste0(m[2], basename(m[3])))
  }
  # `stancflags = --include-paths=/abs/dir --name=...`
  sub("--include-paths=/[^ ]*/([^/ ]+)", "--include-paths=\\1", line)
}

scrub_csv <- function(path) {
  bytes <- readBin(path, "raw", file.size(path))
  text <- rawToChar(bytes)
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  # only `#` lines (the header, and the timing block at the end) carry
  # paths
  hashed <- startsWith(lines, "#")
  lines[hashed] <- vapply(lines[hashed], scrub_header_line, "", USE.NAMES = FALSE)
  out <- paste(lines, collapse = "\n")
  if (endsWith(text, "\n")) {
    out <- paste0(out, "\n")
  }
  writeBin(charToRaw(out), path)
  invisible(path)
}

vars_to_check <- c("lp__", "accept_stat__", "alpha", "mu_intercept", "theta", "phi")

for (model in names(csv)) {
  files <- csv[[model]]
  before_dir <- file.path(tempdir(), paste0("unscrubbed_", model))
  dir.create(before_dir, showWarnings = FALSE)
  file.copy(files, before_dir, overwrite = TRUE)
  before <- lapply(files, readLines)
  draws_before <- read_bilatr_draws(files, vars_to_check)

  for (f in files) scrub_csv(f)

  after <- lapply(files, readLines)
  for (i in seq_along(files)) {
    changed <- which(before[[i]] != after[[i]])
    message("\n", basename(files[i]), ": ", length(changed), " header line(s) changed")
    for (j in changed) {
      message("  - ", before[[i]][j], "\n  + ", after[[i]][j])
    }
    stopifnot(
      length(before[[i]]) == length(after[[i]]),
      all(startsWith(before[[i]][changed], "#")),
      !any(grepl(normalizePath("."), after[[i]], fixed = TRUE)),
      !any(grepl("^#.* = /", after[[i]]))
    )
    hdr <- .scan_stan_csv_header_and_skip(files[i])
    stopifnot(identical(hdr$model, paste0(tools::file_path_sans_ext(
      .bilatr_stan_models[[model]]$file
    ), "_model")))
    stopifnot(identical(.stan_model_from_csv_header(files[i]), model))
  }
  draws_after <- read_bilatr_draws(files, vars_to_check)
  stopifnot(identical(draws_before, draws_after))
  message(model, ": header parsed, model line intact, draws identical after scrubbing")
}

# -- 4. Sizes -----------------------------------------------------------------

sizes <- file.info(list.files(out_root, recursive = TRUE, full.names = TRUE))["size"]
print(sizes)
message("total: ", round(sum(sizes$size) / 1024), " KB")
