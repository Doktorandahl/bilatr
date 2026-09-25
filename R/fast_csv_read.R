#' Extract one `key = value` config field from a CmdStan CSV comment line
#'
#' Config lines look like `#     num_samples = 1000 (Default)` or
#' `#     save_warmup = false (Default)`; keys are indented by nesting
#' level, so the same key name can in principle appear more than once at
#' different levels (not observed in practice for the five keys this
#' package parses, but guarded against anyway -- see
#' [.scan_stan_csv_header_and_skip()]'s "first match wins" handling).
#'
#' @param line A single `#`-prefixed line.
#' @param key The config key to match, e.g. `"num_samples"`.
#' @return The matched value as a string (e.g. `"1000"`, `"false"`), or
#'   `NA_character_` if `line` doesn't match `key`.
#' @keywords internal
.match_stan_csv_config_field <- function(line, key) {
  m <- regmatches(line, regexpr(paste0("^#\\s*", key, "\\s*=\\s*([^ (]+)"), line, perl = TRUE))
  if (length(m) == 0) {
    return(NA_character_)
  }
  sub(paste0("^#\\s*", key, "\\s*=\\s*"), "", m, perl = TRUE)
}

#' Scan a raw CmdStan CSV for its header, config, and the offset of its
#' first post-warmup data row -- replaces `cmdstanr:::read_csv_metadata()`
#'
#' A CmdStan `sample` CSV is laid out, in order: a `#`-prefixed config
#' comment block; the header line (not a comment); `num_warmup_rows` data
#' rows if `save_warmup` (0 otherwise); an "Adaptation terminated"
#' comment block (`# Adaptation terminated`, `# Step size = ...`, `#
#' Diagonal elements of inverse mass matrix:`, and one more `#` line with
#' the values themselves -- 4 lines for the default `diag_e` metric
#' confirmed empirically to stay 4 lines even at 500 parameters, since
#' the diagonal is written as one long comma-separated line; a
#' non-default `dense_e` metric would print more lines, one per row of
#' the covariance matrix, which is why this scans for however many
#' consecutive `#` lines actually follow rather than hard-coding 4) --
#' this block is present even when `adapt_engaged = 0` (confirmed
#' empirically: CmdStan still reports the metric/step size used, just
#' the un-adapted defaults); post-warmup data rows; a `#`-prefixed timing
#' footer.
#'
#' Parses `method`/`num_samples`/`num_warmup`/`save_warmup`/`thin` out of
#' the config block in the same single pass that finds the header and the
#' post-warmup data offset (0.10.0; previously a separate
#' `cmdstanr:::read_csv_metadata()` call, which does a full-file `grep`
#' pass and is a `:::` call into cmdstanr's private API). Each key's
#' FIRST match wins, since keys are indented by nesting level and could in
#' principle repeat at a different level. `save_warmup` is normalized to
#' logical, accepting both `0`/`1` (older CmdStan) and `false`/`true`
#' (CmdStan >= 2.33, confirmed against a real 2.38 CSV). `method` must be
#' `"sample"`; anything else errors, naming the file -- this package only
#' ever reads `sample` output.
#'
#' Reads one line at a time via a text connection and discards each as
#' it's read (never keeps a warmup row), so peak memory here is one
#' line's text, not the file. In this package's own usage `save_warmup`
#' is always left at its `FALSE` default (grep confirms no caller sets
#' it), so this scan touches only the config block, the header, and the
#' (typically ~4-line) adaptation block -- tens of lines regardless of the
#' file's total size or column count.
#'
#' @param csv_file A single raw CmdStan CSV file path.
#' @return A list with `skip` (integer: number of lines to pass so the
#'   next line is the first post-warmup draw -- i.e. `fread(skip =
#'   this)`'s first row), `header_line` (the raw, unparsed header line),
#'   `num_samples`/`num_warmup`/`thin` (integers) and `save_warmup`
#'   (logical), parsed from the config block, and `num_post_warmup_draws`
#'   (`ceiling(num_samples / thin)`, the same arithmetic
#'   `cmdstanr:::read_csv_metadata()` used).
#' @keywords internal
.scan_stan_csv_header_and_skip <- function(csv_file) {
  con <- file(csv_file, open = "rt")
  on.exit(close(con), add = TRUE)

  config_keys <- c("method", "num_samples", "num_warmup", "save_warmup", "thin")
  config <- stats::setNames(as.list(rep(NA_character_, length(config_keys))), config_keys)

  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop("'", csv_file, "': reached end of file before finding the CSV header.", call. = FALSE)
    }
    skip <- skip + 1L
    if (!startsWith(line, "#")) break
    for (key in config_keys) {
      if (is.na(config[[key]])) {
        val <- .match_stan_csv_config_field(line, key)
        if (!is.na(val)) config[[key]] <- val
      }
    }
  }
  header_line <- line

  missing_keys <- config_keys[vapply(config, is.na, logical(1))]
  if (length(missing_keys) > 0) {
    stop(
      "'", csv_file, "': could not find config field(s) ",
      paste(missing_keys, collapse = ", "), " in the CSV's comment header.",
      call. = FALSE
    )
  }
  if (!identical(config$method, "sample")) {
    stop(
      "'", csv_file, "': method = '", config$method, "', expected 'sample' ",
      "-- this package only reads CmdStan `sample` output.",
      call. = FALSE
    )
  }
  save_warmup <- config$save_warmup %in% c("1", "true")
  num_samples <- as.integer(config$num_samples)
  num_warmup <- as.integer(config$num_warmup)
  thin <- as.integer(config$thin)
  num_warmup_rows <- if (save_warmup) ceiling(num_warmup / thin) else 0L

  if (num_warmup_rows > 0L) {
    for (i in seq_len(num_warmup_rows)) {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (length(line) == 0L) {
        stop(
          "'", csv_file, "': reached end of file while passing over ",
          num_warmup_rows, " warmup row(s) -- fewer than expected are present.",
          call. = FALSE
        )
      }
      skip <- skip + 1L
    }
  }

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop(
        "'", csv_file, "': reached end of file in the post-warmup ",
        "adaptation comment block, before any post-warmup draws.",
        call. = FALSE
      )
    }
    if (startsWith(line, "#")) {
      skip <- skip + 1L
    } else {
      break
    }
  }

  list(
    skip = skip, header_line = header_line,
    num_samples = num_samples, num_warmup = num_warmup,
    save_warmup = save_warmup, thin = thin,
    num_post_warmup_draws = ceiling(num_samples / thin)
  )
}

#' Prepare a set of raw CmdStan CSV files for repeated, memory-bounded reads
#'
#' The one-time setup shared by every CSV-file-path entry point
#' ([diagnose_convergence()], [diagnose_and_extract_bilatr()],
#' [extract_theta()], [extract_alpha()]/[extract_mu_intercept()] via
#' [.get_draws]). Scans every one of `csv_files` (see
#' [.scan_stan_csv_header_and_skip], which replaces
#' `cmdstanr:::read_csv_metadata()`, 0.10.0) for its config and the line
#' offset of its first post-warmup draw, and checks the chains agree
#' (unlike the `cmdstanr:::read_csv_metadata()`-based path this replaced,
#' which only ever looked at `csv_files[1]` and silently assumed the rest
#' matched). Everything downstream -- [.fast_read_post_warmup_draws] for
#' the Tier 1/2 read and every Tier 3 chunk -- reads `variables` directly
#' out of the raw CSV at that offset via `data.table::fread(file = ...,
#' skip = ..., nrows = ...)`, never a comment-stripped copy: no scratch
#' file is written anywhere, and no file is read more than once per
#' variable subset requested (see individual entry points for exactly
#' how many read passes each makes).
#'
#' Also validates, per file, that at least `num_post_warmup_draws` rows
#' actually follow the computed offset (a chain killed mid-sampling
#' produces fewer), via a single-column `fread(skip=, nrows=)` probe --
#' `fread` already has to touch the whole file once per call regardless
#' of how few columns are selected (see [.fast_read_post_warmup_draws]'s
#' docs), so this probe's real cost is that one file-sized touch, same
#' as any other read of the file, not something a slower line-by-line
#' scan would avoid.
#'
#' @param csv_files Character vector of raw CmdStan CSV file paths (one
#'   per chain).
#' @return A list with `csv_files` (as given), `variables` (all
#'   monitored quantities, bracket form, e.g. `"theta[3,12]"`),
#'   `dot_name_map` (named character vector: bracket name -> raw
#'   dot-form CSV column name, e.g. `"theta.3.12"`), `col_index` (named
#'   integer vector: bracket name -> 1-based column position in the raw
#'   CSV header -- sampler-diagnostic columns like `accept_stat__`
#'   occupy earlier positions and are counted in this numbering even
#'   though they're never named in `variables`), `num_post_warmup_draws`,
#'   `n_chains`, `data_skip` (integer vector, one per `csv_files` entry,
#'   for `fread(skip = )`), and `file_mb` (the LARGEST `csv_files` entry's
#'   size in MB -- `fread()` reads one file at a time, so this is the
#'   file-sized memory floor any single read pays regardless of
#'   `chunk_size`; see [.compute_chunk_size()]).
#' @keywords internal
.prepare_fast_csv_read <- function(csv_files) {
  scans <- lapply(csv_files, .scan_stan_csv_header_and_skip)

  ref <- scans[[1]]
  for (k in seq_along(scans)[-1]) {
    if (!identical(scans[[k]]$num_samples, ref$num_samples) ||
      !identical(scans[[k]]$thin, ref$thin) ||
      !identical(scans[[k]]$save_warmup, ref$save_warmup)) {
      stop(
        "'", csv_files[k], "' disagrees with '", csv_files[1], "' on ",
        "num_samples/thin/save_warmup -- chains from the same run must ",
        "share these.",
        call. = FALSE
      )
    }
    if (!identical(scans[[k]]$header_line, ref$header_line)) {
      stop(
        "'", csv_files[k], "' has a different CSV header than '", csv_files[1],
        "' -- chains from the same run must share the same variable layout.",
        call. = FALSE
      )
    }
  }

  num_post_warmup_draws <- ref$num_post_warmup_draws
  data_skip <- vapply(scans, `[[`, integer(1), "skip")

  header_fields <- strsplit(ref$header_line, ",", fixed = TRUE)[[1]]
  col_index_dot <- stats::setNames(seq_along(header_fields), header_fields)

  # Every registered CmdStan model's monitored quantities: every column
  # NOT ending in "__" (the sampler-diagnostic columns, e.g.
  # accept_stat__/treedepth__/energy__), plus "lp__" specifically --
  # confirmed against cmdstanr:::read_csv_metadata()$variables on a real
  # CmdStan 2.38 CSV, which keeps lp__ despite its trailing "__".
  variables_dot <- header_fields[!endsWith(header_fields, "__") | header_fields == "lp__"]

  bracket_vars <- vapply(variables_dot, .dot_name_to_bracket, character(1), USE.NAMES = FALSE)
  dot_name_map <- stats::setNames(variables_dot, bracket_vars)
  col_index <- stats::setNames(col_index_dot[variables_dot], bracket_vars)

  n_chains <- length(csv_files)

  for (k in seq_len(n_chains)) {
    probe <- data.table::fread(
      file = csv_files[k], skip = data_skip[k], nrows = num_post_warmup_draws,
      header = FALSE, select = 1L, data.table = FALSE, showProgress = FALSE
    )
    if (nrow(probe) < num_post_warmup_draws) {
      stop(
        "'", csv_files[k], "' has only ", nrow(probe), " post-warmup row(s), ",
        "expected ", num_post_warmup_draws, " (a chain that ended early, e.g. ",
        "killed mid-sampling, produces this). Re-run or exclude this file.",
        call. = FALSE
      )
    }
  }

  list(
    csv_files = csv_files,
    variables = bracket_vars,
    dot_name_map = dot_name_map,
    col_index = col_index,
    num_post_warmup_draws = num_post_warmup_draws,
    n_chains = n_chains,
    data_skip = data_skip,
    file_mb = max(file.size(csv_files)) / 1e6
  )
}

#' Expand variable filters to concrete variable names, cmdstanr-style
#'
#' Replicates `cmdstanr:::matching_variables()`'s two-step lookup so
#' [.fast_read_post_warmup_draws] accepts the same `variables` spellings
#' [cmdstanr::read_cmdstan_csv()] does: an exact match (e.g.
#' `"theta[3,12]"`, `"alpha[1]"`) is used as-is, and anything that
#' doesn't match exactly (e.g. a bare base name like `"alpha"` or
#' `"mu_intercept"`, as [.get_draws] passes for [extract_alpha()]/
#' [extract_mu_intercept()]) is expanded to every available variable
#' whose name starts with `"<filter>["` -- i.e. every indexed element of
#' that parameter.
#'
#' @param available Character vector of all bracket-form variable names
#'   in the CSV header (`prepared$variables`).
#' @param filters Character vector of requested variable names/base names.
#' @return Character vector of concrete, available variable names, in
#'   the order their filters were given (a base-name filter expands to
#'   its matches in `available`'s order), with duplicates dropped (B5): a
#'   request like `c("alpha", "alpha[1]")` would otherwise return
#'   `"alpha[1]"` twice, and `data.table::fread(select = )` errors opaquely
#'   ("Column number ... has been selected twice") on a duplicated column.
#' @keywords internal
.expand_bilatr_variable_filters <- function(available, filters) {
  matched <- as.list(match(filters, available))
  for (id in which(is.na(matched))) {
    matched[[id]] <- which(startsWith(available, paste0(filters[id], "[")))
  }
  not_found <- filters[vapply(matched, length, 0L) == 0]
  if (length(not_found) > 0) {
    stop(
      "Variable(s) not found in the CmdStan CSV header: ",
      paste(not_found, collapse = ", "),
      call. = FALSE
    )
  }
  unique(available[unlist(matched)])
}

#' Read named variables' post-warmup draws directly from raw CmdStan CSVs
#'
#' The direct replacement for
#' `cmdstanr::read_cmdstan_csv(csv_files, variables = variables)$post_warmup_draws`
#' used throughout the package's CSV-file-path branches: reads
#' `variables` straight out of each of `prepared$csv_files` with
#' `data.table::fread(file = ..., skip = prepared$data_skip[k], nrows =
#' prepared$num_post_warmup_draws, header = FALSE, select =
#' <column positions>)` -- the comment lines before the block are
#' skipped by line count, the timing footer is never reached because
#' `nrows` stops the parser there, and column identity comes from
#' `prepared$col_index` rather than a header row (there isn't one in the
#' slice being read). No intermediate comment-stripped copy of the file
#' is ever created.
#'
#' `fread` on a real file path (as opposed to `cmd =`/piped input) still
#' touches the whole file once per call -- it samples rows spread across
#' the full byte range for column-type detection before honoring
#' `skip`/`nrows` for the actual parse -- so this is not free at
#' production scale (the caller's chunking exists to bound how much is
#' actually *parsed and kept*, not this per-call file-sized touch; see
#' [.compute_chunk_size()]'s docs for how that's now counted). What it
#' avoids is `cmdstanr::read_cmdstan_csv()`'s `cmd = grep` path, which
#' additionally has to buffer/copy the file's full non-comment content
#' before parsing anything from it.
#'
#' Builds the result directly as a `posterior::draws_array` rather than
#' via a list of per-chain matrices: preallocates `array(NA_real_, dim =
#' c(n_draws, n_chains, n_vars))` and fills it one chain's `as.matrix()`
#' at a time (dropping that chain's `data.table` immediately after), so
#' only one chain's matrix is ever alive alongside the array, rather than
#' `n_chains` matrices plus whatever copies `as_draws_array()` would make
#' assembling a list of them.
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param variables Character vector of posterior-style variable names to
#'   read -- either exact (e.g. `"theta[3,12]"`, `"alpha[1]"`) or a bare
#'   base name (e.g. `"alpha"`, `"mu_intercept"`), expanded to every
#'   indexed element of that name via [.expand_bilatr_variable_filters],
#'   matching [cmdstanr::read_cmdstan_csv()]'s own `variables` matching.
#' @return A `posterior::draws_array`.
#' @keywords internal
.fast_read_post_warmup_draws <- function(prepared, variables) {
  resolved_vars <- .expand_bilatr_variable_filters(prepared$variables, variables)
  col_idx <- unname(prepared$col_index[resolved_vars])

  n_draws <- prepared$num_post_warmup_draws
  n_chains <- prepared$n_chains
  n_vars <- length(resolved_vars)

  arr <- array(
    NA_real_, dim = c(n_draws, n_chains, n_vars),
    dimnames = list(NULL, NULL, resolved_vars)
  )

  for (k in seq_len(n_chains)) {
    dt <- data.table::fread(
      file = prepared$csv_files[k], skip = prepared$data_skip[k],
      nrows = n_draws, header = FALSE, select = col_idx,
      data.table = TRUE, showProgress = FALSE
    )
    m <- as.matrix(dt)
    rm(dt)
    arr[, k, ] <- m
    rm(m)
  }

  posterior::as_draws_array(arr)
}
