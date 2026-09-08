#' Strip comment lines from a raw CmdStan CSV, streaming, to a real file
#'
#' CmdStan's own CSV format intersperses `#`-prefixed comment lines (a
#' config block up front, an "Adaptation terminated" marker between
#' warmup and post-warmup rows, and a timing footer) with the actual
#' header/data rows. [cmdstanr::read_cmdstan_csv()] strips these via
#' `data.table::fread(cmd = "grep -v '^#' <file>")`; per `?data.table::fread`,
#' a `cmd=`/piped input is first run to completion with its *entire*
#' output written to a fresh file in `tempdir()`, and only then read "as
#' normal" -- so every single call, however few `variables` it requests,
#' pays for a full near-complete copy of that one chain's raw CSV (often
#' several GB at production scale) landing somewhere under `tempdir()`.
#' [.chunked_summarise_csv] calls this once per chunk, so a Tier 3 sweep
#' with hundreds of chunks repeats that copy hundreds of times per chain
#' -- and if `tempdir()`/`$TMPDIR` happens to be a RAM-backed `tmpfs` (a
#' common per-node HPC/SLURM default), each copy is a direct hit against
#' the job's memory allocation, independent of `max_memory_mb`/
#' `chunk_size`, which is what actually exhausts memory on jobs like
#' this one however tightly `max_memory_mb` is set.
#'
#' This does the same comment-stripping exactly once per chain file (see
#' [.prepare_fast_csv_read], which every CSV-path entry point calls a
#' single time up front, before any chunking), streaming line-by-line via
#' [readr::read_lines_chunked()] rather than shelling out to `grep`, so
#' peak memory here is bounded by `chunk_size` lines rather than the
#' whole file, and it depends on no external binary (portable to
#' Windows/any system with just R installed -- no `grep`/`grep.exe` PATH
#' assumptions). The output is written directly to `out_file`, a real
#' path the caller controls (see `scratch_dir` in
#' [.prepare_fast_csv_read]), never through `tempdir()`.
#'
#' @param csv_file A single raw CmdStan CSV file path.
#' @param out_file Destination path for the comment-stripped copy.
#' @return `out_file`, invisibly.
#' @keywords internal
.strip_stan_csv_comments <- function(csv_file, out_file, chunk_size = 200000L) {
  con_out <- file(out_file, open = "wt")
  on.exit(close(con_out), add = TRUE)

  readr::read_lines_chunked(
    csv_file,
    callback = readr::SideEffectChunkCallback$new(function(lines, pos) {
      keep <- lines[!startsWith(lines, "#")]
      if (length(keep) > 0) writeLines(keep, con_out)
    }),
    chunk_size = chunk_size,
    progress = FALSE
  )

  if (file.size(out_file) == 0) {
    stop(
      "Stripping comment lines from '", csv_file, "' produced an empty ",
      "file -- is this a valid CmdStan CSV?",
      call. = FALSE
    )
  }
  invisible(out_file)
}

#' Prepare a set of raw CmdStan CSV files for repeated, memory-bounded reads
#'
#' The one-time setup shared by every CSV-file-path entry point
#' ([diagnose_convergence()], [diagnose_and_extract_bilatr()],
#' [extract_theta()], [extract_alpha()]/[extract_mu_intercept()] via
#' [.get_draws]): comment-strips each chain file exactly once (see
#' [.strip_stan_csv_comments]) and reads the header metadata needed to
#' slice out post-warmup rows and map posterior-style variable names
#' (`"theta[3,12]"`) back to the raw CSV's dot-form column names
#' (`"theta.3.12"`). Everything downstream -- the Tier 1/2 read, the
#' `alpha[1]` flip check, every chunk of a Tier 3 sweep -- reads from the
#' resulting cleaned files via [.fast_read_post_warmup_draws] instead of
#' [cmdstanr::read_cmdstan_csv()], so the expensive comment-stripping
#' pass happens once per chain file per call, not once per chunk.
#'
#' @param csv_files Character vector of raw CmdStan CSV file paths (one
#'   per chain).
#' @param scratch_dir Directory to write the comment-stripped copies
#'   into. `NULL` (the default) writes each copy alongside its source
#'   file (`dirname(csv_file)`) -- deliberately, since that's already
#'   known-good storage for these (often many-GB) files, unlike
#'   `tempdir()`, which on some HPC systems is a RAM-backed `tmpfs` (see
#'   [.strip_stan_csv_comments]). Pass an explicit path (e.g. fast local
#'   scratch) if the output directory's filesystem is quota-constrained,
#'   read-only, or otherwise unsuitable for a temporary same-size copy of
#'   each chain file.
#' @return A list with `clean_files` (comment-stripped copies, one per
#'   `csv_files` entry, same order), `dot_name_map` (a named character
#'   vector mapping every posterior-style variable name to its raw
#'   dot-form CSV column name), and `num_post_warmup_draws`.
#' @keywords internal
.prepare_fast_csv_read <- function(csv_files, scratch_dir = NULL) {
  meta <- cmdstanr:::read_csv_metadata(csv_files[1])
  bracket_vars <- vapply(meta$variables, .dot_name_to_bracket, character(1), USE.NAMES = FALSE)
  dot_name_map <- stats::setNames(meta$variables, bracket_vars)
  num_post_warmup_draws <- ceiling(meta$iter_sampling / meta$thin)

  clean_files <- vapply(csv_files, function(f) {
    dir <- scratch_dir %||% dirname(f)
    out <- file.path(dir, paste0(basename(f), ".nocomments"))
    .strip_stan_csv_comments(f, out)
    out
  }, character(1), USE.NAMES = FALSE)

  list(
    clean_files = clean_files,
    dot_name_map = dot_name_map,
    num_post_warmup_draws = num_post_warmup_draws
  )
}

#' Delete the comment-stripped scratch copies made by [.prepare_fast_csv_read]
#' @param prepared Output of [.prepare_fast_csv_read].
#' @keywords internal
.cleanup_fast_csv_read <- function(prepared) {
  unlink(prepared$clean_files)
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
#'   in the CSV header (`names(prepared$dot_name_map)`).
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

#' Read named variables' post-warmup draws from prepared CmdStan CSVs
#'
#' The direct replacement for
#' `cmdstanr::read_cmdstan_csv(csv_files, variables = variables)$post_warmup_draws`
#' used throughout the package's CSV-file-path branches (see
#' [.prepare_fast_csv_read]'s docs for why). Reads `variables` from each
#' of `prepared`'s comment-stripped files with `data.table::fread(file =
#' ...)` -- a real file path, so `fread` can `mmap()` and parse only the
#' `select`ed columns, unlike the `cmd=`/piped path this replaces -- then
#' keeps only the last `num_post_warmup_draws` rows (warmup rows, if
#' present, always precede post-warmup rows in a CmdStan CSV, so this is
#' exact regardless of whether warmup was saved) and assembles the
#' per-chain matrices into a `posterior::draws_array` exactly as
#' [cmdstanr::read_cmdstan_csv()] does internally
#' (`posterior::as_draws_array()` on a list of same-column-name
#' matrices).
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param variables Character vector of posterior-style variable names to
#'   read -- either exact (e.g. `"theta[3,12]"`, `"alpha[1]"`) or a bare
#'   base name (e.g. `"alpha"`, `"mu_intercept"`), expanded to every
#'   indexed element of that name via [.expand_bilatr_variable_filters],
#'   matching [cmdstanr::read_cmdstan_csv()]'s own `variables` matching.
#' @param flip If `TRUE`, negate every chain's draws for `variables`
#'   before assembling -- applied pre-summary to the raw values, matching
#'   [bilatr_orient()]'s sign correction (see [.chunked_summarise_csv]).
#' @return A `posterior::draws_array`.
#' @keywords internal
.fast_read_post_warmup_draws <- function(prepared, variables, flip = FALSE) {
  resolved_vars <- .expand_bilatr_variable_filters(names(prepared$dot_name_map), variables)
  dot_vars <- prepared$dot_name_map[resolved_vars]

  chain_mats <- lapply(prepared$clean_files, function(f) {
    df <- data.table::fread(
      file = f, select = unname(dot_vars), header = TRUE,
      data.table = FALSE, showProgress = FALSE
    )
    df <- utils::tail(df, prepared$num_post_warmup_draws)
    m <- as.matrix(df)
    colnames(m) <- resolved_vars
    if (flip) m <- -m
    m
  })

  posterior::as_draws_array(chain_mats)
}
