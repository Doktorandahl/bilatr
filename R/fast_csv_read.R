#' Scan a raw CmdStan CSV for its header and the offset of its first
#' post-warmup data row
#'
#' A CmdStan `sample` CSV is laid out, in order: a `#`-prefixed config
#' comment block; the header line (not a comment); `num_warmup_rows` data
#' rows if `save_warmup = 1` (0 otherwise); an "Adaptation terminated"
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
#' Reads one line at a time via a text connection and discards each as
#' it's read (never keeps a warmup row), so peak memory here is one
#' line's text, not the file. In this package's own usage `save_warmup`
#' is always left at its `FALSE` default (grep confirms no caller sets
#' it), so `num_warmup_rows` is 0 and this scan touches only the config
#' block, the header, and the (typically ~4-line) adaptation block --
#' tens of lines regardless of the file's total size or column count.
#'
#' @param csv_file A single raw CmdStan CSV file path.
#' @param num_warmup_rows Number of warmup data rows to pass over (0 if
#'   `save_warmup = FALSE`), from metadata.
#' @return A list with `skip` (integer: number of lines to pass so the
#'   next line is the first post-warmup draw -- i.e. `fread(skip =
#'   this)`'s first row) and `header_line` (the raw, unparsed header
#'   line).
#' @keywords internal
.scan_stan_csv_header_and_skip <- function(csv_file, num_warmup_rows) {
  con <- file(csv_file, open = "rt")
  on.exit(close(con), add = TRUE)

  skip <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      stop("'", csv_file, "': reached end of file before finding the CSV header.", call. = FALSE)
    }
    skip <- skip + 1L
    if (!startsWith(line, "#")) break
  }
  header_line <- line

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

  list(skip = skip, header_line = header_line)
}

#' Prepare a set of raw CmdStan CSV files for repeated, memory-bounded reads
#'
#' The one-time setup shared by every CSV-file-path entry point
#' ([diagnose_convergence()], [diagnose_and_extract_bilatr()],
#' [extract_theta()], [extract_alpha()]/[extract_mu_intercept()] via
#' [.get_draws]). Calls `cmdstanr:::read_csv_metadata()` exactly once
#' (on `csv_files[1]`; every chain of the same run shares the same
#' variable/column layout, matching what
#' `cmdstanr::read_cmdstan_csv()` itself assumes) and scans each file
#' once (see [.scan_stan_csv_header_and_skip]) for the line offset of
#' its first post-warmup draw. Everything downstream --
#' [.fast_read_post_warmup_draws] for the Tier 1/2 read, the `alpha[1]`
#' flip check, and every Tier 3 chunk -- reads `variables` directly out
#' of the raw CSV at that offset via `data.table::fread(file = ...,
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
  meta <- cmdstanr:::read_csv_metadata(csv_files[1])
  num_warmup_rows <- if (isTRUE(meta$save_warmup == 1)) ceiling(meta$iter_warmup / meta$thin) else 0L
  num_post_warmup_draws <- ceiling(meta$iter_sampling / meta$thin)

  scans <- lapply(csv_files, .scan_stan_csv_header_and_skip, num_warmup_rows = num_warmup_rows)
  data_skip <- vapply(scans, `[[`, integer(1), "skip")

  header_fields <- strsplit(scans[[1]]$header_line, ",", fixed = TRUE)[[1]]
  col_index_dot <- stats::setNames(seq_along(header_fields), header_fields)

  bracket_vars <- vapply(meta$variables, .dot_name_to_bracket, character(1), USE.NAMES = FALSE)
  dot_name_map <- stats::setNames(meta$variables, bracket_vars)
  col_index <- stats::setNames(col_index_dot[meta$variables], bracket_vars)

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
#'   its matches in `available`'s order).
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
  available[unlist(matched)]
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
#' @param flip_vars Character vector of variable base names
#'   ([.bilatr_flip_variables()]'s output) to negate -- matched against
#'   `resolved_vars` via [.bilatr_match_draws_columns()] and negated on
#'   each chain's temporary matrix `m`, before it is ever assigned into
#'   `arr`. `m` is already about to be discarded at that point, so this
#'   costs nothing extra: no copy, no second pass, no data-frame
#'   round-trip. `character(0)` (the default) negates nothing. Column-
#'   selective by design (B5): a request mixing `theta`/`theta_raw` with
#'   `log_lik[d,t]` flips only the columns `flip_vars` actually lists.
#' @return A `posterior::draws_array`.
#' @keywords internal
.fast_read_post_warmup_draws <- function(prepared, variables, flip_vars = character(0)) {
  resolved_vars <- .expand_bilatr_variable_filters(prepared$variables, variables)
  col_idx <- unname(prepared$col_index[resolved_vars])

  n_draws <- prepared$num_post_warmup_draws
  n_chains <- prepared$n_chains
  n_vars <- length(resolved_vars)

  flip_cols <- .bilatr_match_draws_columns(resolved_vars, flip_vars)
  flip_pos <- match(flip_cols, resolved_vars)

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
    if (length(flip_pos) > 0) m[, flip_pos] <- -m[, flip_pos]
    arr[, k, ] <- m
    rm(m)
  }

  posterior::as_draws_array(arr)
}
