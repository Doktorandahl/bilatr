#' Names of Stan parameters treated as global/shared (Tier 1)
#'
#' Matched against a monitored quantity's *base name* (its variable name
#' with any `[...]` index stripped), so this correctly matches vector
#' parameters like `alpha[1]`, `alpha[2]`, ... via their shared base name
#' `alpha`.
#' @keywords internal
.bilatr_tier1_names <- c(
  "alpha", "alpha_raw", "mu_intercept", "mu_intercept_raw",
  "mu_theta0", "sigma_theta0",
  "mu_log_phi", "sigma_log_phi", "mu_log_noise", "sigma_log_noise",
  "lp__"
)

#' Split a `posterior::summarise_draws()` variable name into base name and
#' index components
#'
#' `"phi[3]"` -> `base = "phi"`, `indices = list(3L)`. `"theta[3,12]"` ->
#' `base = "theta"`, `indices = list(3L, 12L)`. Scalars (no `[...]`) get
#' `indices = list()`.
#'
#' @param variable Character vector of `summarise_draws()` variable names.
#' @return A tibble with columns `variable`, `base_name`, `n_index`, `index_1`,
#'   `index_2` (the latter two `NA_integer_` when not applicable).
#' @keywords internal
.parse_variable_indices <- function(variable) {
  base_name <- stringr::str_remove(variable, "\\[.*\\]$")
  index_str <- stringr::str_match(variable, "\\[(.*)\\]$")[, 2]
  index_parts <- stringr::str_split(index_str, ",")

  tibble::tibble(
    variable = variable,
    base_name = base_name,
    n_index = purrr::map_int(index_parts, ~ if (all(is.na(.x))) 0L else length(.x)),
    index_1 = purrr::map_int(index_parts, ~ if (all(is.na(.x))) NA_integer_ else as.integer(.x[1])),
    index_2 = purrr::map_int(index_parts, ~ if (all(is.na(.x))) NA_integer_ else if (length(.x) >= 2) as.integer(.x[2]) else NA_integer_)
  )
}

#' Classify monitored quantities into diagnostic tiers by name/index shape
#'
#' Tier 1 (global/shared) is matched first, by base name, against
#' [.bilatr_tier1_names]. Everything else is classified structurally by
#' how many `[...]` indices it carries: a single index (`name[d]`) is
#' assumed to be a per-dyad hierarchical parameter (Tier 2, joined on
#' `d`); two indices (`name[d, t]`) is assumed to be a per-dyad-period
#' latent state (Tier 3, joined on `d`). This is deliberately structural
#' rather than a fixed per-parameter name list, so a future model variant
#' that changes a parameter's shape (e.g. makes `phi` per-dyad-period,
#' `phi[d, t]`, instead of the per-dyad `phi[d]` of the `stable` model) is
#' still classified consistently without special-casing.
#' Anything with no brackets that isn't in the Tier 1 name list (should
#' not occur for the package's own models, but could for a hand-edited
#' Stan file) is folded into Tier 1 rather than dropped, since its
#' sparsity profile is unknown and it should never be silently hidden.
#'
#' @param variable Character vector of `summarise_draws()` variable names.
#' @return A tibble with columns `variable`, `tier` (`1L`, `2L`, or `3L`),
#'   `dyad_id` (the first index, `NA` for Tier 1), and `time_index` (the
#'   second index, `NA` outside Tier 3).
#' @keywords internal
.classify_bilatr_tier <- function(variable) {
  parsed <- .parse_variable_indices(variable)

  dplyr::mutate(
    parsed,
    tier = dplyr::case_when(
      base_name %in% .bilatr_tier1_names ~ 1L,
      n_index == 0 ~ 1L,
      n_index == 1 ~ 2L,
      n_index >= 2 ~ 3L
    ),
    dyad_id = dplyr::if_else(tier %in% c(2L, 3L), index_1, NA_integer_),
    time_index = dplyr::if_else(tier == 3L, index_2, NA_integer_)
  ) %>%
    dplyr::select(variable, tier, dyad_id, time_index)
}

#' Normalize the `n_dt` argument to a two-column tibble
#'
#' @param n_dt A data frame with `dyad`/`n_dt`-like columns, or a named
#'   numeric vector keyed by dyad id.
#' @return A tibble with columns `dyad_id` (integer) and `n_dt` (numeric).
#' @keywords internal
.normalize_n_dt <- function(n_dt) {
  if (is.data.frame(n_dt)) {
    nm <- names(n_dt)
    dyad_col <- nm[stringr::str_detect(tolower(nm), "^dyad")][1]
    n_dt_col <- nm[stringr::str_detect(tolower(nm), "n_dt|n_obs|n_events")][1]
    if (is.na(dyad_col) || is.na(n_dt_col)) {
      stop(
        "`n_dt` must have a dyad-id-like column (matching \"dyad...\") and ",
        "an observation-count-like column (matching \"n_dt\"/\"n_obs\"/\"n_events\").",
        call. = FALSE
      )
    }
    out <- tibble::tibble(
      dyad_id = as.integer(n_dt[[dyad_col]]),
      n_dt = as.numeric(n_dt[[n_dt_col]])
    )
  } else if (is.numeric(n_dt) && !is.null(names(n_dt))) {
    out <- tibble::tibble(
      dyad_id = as.integer(names(n_dt)),
      n_dt = as.numeric(n_dt)
    )
  } else {
    stop(
      "`n_dt` must be a data frame with dyad-id and n_dt columns, or a ",
      "named numeric vector keyed by dyad id.",
      call. = FALSE
    )
  }

  if (anyNA(out$dyad_id)) {
    stop("`n_dt` has dyad ids that could not be coerced to integers.", call. = FALSE)
  }
  out
}

#' Flag diagnostics that are worse than a smooth sparsity trend predicts
#'
#' Bins dyads into `n_dt` quantile buckets and flags a dyad's diagnostic
#' value as "worse than expected for its sparsity bracket" when it falls
#' below `median - 1.5 * IQR` of its own bucket. This is deliberately a
#' within-bracket outlier rule rather than a smoother (loess/rank
#' regression): it needs no bandwidth/shape assumptions, degrades
#' gracefully with few dyads (bucket count shrinks automatically), and
#' the "1.5 * IQR" rule is the standard Tukey outlier convention, so a
#' flagged dyad is an outlier *relative to equally-sparse peers*, not
#' merely a sparse dyad.
#'
#' @param n_dt Numeric vector of per-dyad observation counts.
#' @param value Numeric vector (same length/order as `n_dt`) of the
#'   diagnostic to flag (e.g. ESS_bulk); lower is assumed worse.
#' @param n_bins Target number of quantile bins; automatically reduced
#'   for small `n_dt`.
#' @return Logical vector, same length as `n_dt`/`value`: `TRUE` where the
#'   value is worse than its bracket's expectation. `NA` where `value` is
#'   `NA`.
#' @keywords internal
.flag_worse_than_expected <- function(n_dt, value, n_bins = 6L) {
  n <- length(n_dt)
  bins <- dplyr::ntile(n_dt, min(n_bins, max(1L, n)))

  tibble::tibble(.row = seq_len(n), bin = bins, value = value) %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(
      .threshold = stats::median(value, na.rm = TRUE) - 1.5 * stats::IQR(value, na.rm = TRUE),
      .flag = value < .threshold
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.row) %>%
    dplyr::pull(.flag)
}

#' Flag a diagnostic row as breaching Rhat/ESS thresholds
#'
#' Shared threshold logic used for Tier 1 (and, per-entry, Tier 3)
#' flagging, factored out so it has a single, directly testable
#' definition.
#'
#' @param rhat,ess_bulk,ess_tail Numeric vectors (recycled against each
#'   other as usual).
#' @param rhat_threshold Rhat values strictly above this are flagged.
#' @param ess_threshold ESS_bulk/ESS_tail values strictly below this are
#'   flagged.
#' @return Logical vector: `TRUE` where any of the three diagnostics
#'   breaches its threshold.
#' @keywords internal
.flag_diagnostic <- function(rhat, ess_bulk, ess_tail, rhat_threshold, ess_threshold) {
  (rhat > rhat_threshold) | (ess_bulk < ess_threshold) | (ess_tail < ess_threshold)
}

#' Validate the `tiers` argument of [diagnose_convergence()]
#'
#' @param tiers Value passed as `diagnose_convergence()`'s `tiers`
#'   argument.
#' @return Sorted, deduplicated integer vector, a subset of `1:3`.
#' @keywords internal
.validate_tiers <- function(tiers) {
  tiers_int <- suppressWarnings(as.integer(tiers))
  if (length(tiers_int) == 0 || anyNA(tiers_int) || !all(tiers_int %in% 1:3)) {
    stop(
      "`tiers` must be a subset of 1:3, identifying which of Tier 1 ",
      "(global/shared), Tier 2 (per-dyad), and Tier 3 (per-dyad-period) ",
      "to compute. Got: ", paste(tiers, collapse = ", "), ".",
      call. = FALSE
    )
  }
  sort(unique(tiers_int))
}

#' Compute the Tier 1 (global/shared) diagnostics tibble
#' @keywords internal
.compute_tier1 <- function(summ, rhat_threshold, ess_threshold) {
  summ %>%
    dplyr::filter(tier == 1L) %>%
    dplyr::mutate(
      flagged = .flag_diagnostic(rhat, ess_bulk, ess_tail, rhat_threshold, ess_threshold)
    ) %>%
    dplyr::select(variable, rhat, ess_bulk, ess_tail, flagged) %>%
    dplyr::arrange(dplyr::desc(flagged))
}

#' Compute the Tier 2 (per-dyad hierarchical parameter) diagnostics tibble
#' @keywords internal
.compute_tier2 <- function(summ, n_dt_tbl) {
  tier2_raw <- summ %>%
    dplyr::filter(tier == 2L) %>%
    dplyr::mutate(base_name = stringr::str_remove(variable, "\\[.*\\]$"))

  dyads_in_draws <- unique(tier2_raw$dyad_id)
  dyads_in_n_dt <- n_dt_tbl$dyad_id

  tier2_joined <- dplyr::full_join(tier2_raw, n_dt_tbl, by = "dyad_id")

  tier2 <- tier2_joined %>%
    dplyr::select(dyad_id, n_dt, base_name, rhat, ess_bulk, ess_tail) %>%
    tidyr::pivot_wider(
      id_cols = c(dyad_id, n_dt),
      names_from = base_name,
      values_from = c(rhat, ess_bulk, ess_tail),
      names_glue = "{base_name}_{.value}"
    ) %>%
    dplyr::select(-dplyr::any_of("NA_rhat"), -dplyr::any_of("NA_ess_bulk"), -dplyr::any_of("NA_ess_tail"))

  worse_cols <- grep("_ess_bulk$", names(tier2), value = TRUE)
  worse_flags <- purrr::map(worse_cols, ~ .flag_worse_than_expected(tier2$n_dt, tier2[[.x]]))
  tier2$worse_than_expected <- if (length(worse_flags) > 0) {
    purrr::reduce(worse_flags, `|`, .init = rep(FALSE, nrow(tier2))) & !is.na(tier2$n_dt)
  } else {
    rep(NA, nrow(tier2))
  }
  tier2 <- dplyr::arrange(tier2, dplyr::desc(worse_than_expected))

  list(
    tier2 = tier2,
    n_dyads_missing_from_n_dt = length(setdiff(dyads_in_draws, dyads_in_n_dt)),
    n_dyads_missing_from_draws = length(setdiff(dyads_in_n_dt, dyads_in_draws))
  )
}

#' Compute the Tier 3 (per-dyad-period latent state) diagnostics tibble
#' @keywords internal
.compute_tier3 <- function(summ, n_dt_tbl, rhat_threshold, ess_threshold) {
  tier3_raw <- summ %>%
    dplyr::filter(tier == 3L) %>%
    dplyr::mutate(
      rhat_flag = rhat > rhat_threshold,
      ess_flag = (ess_bulk < ess_threshold) | (ess_tail < ess_threshold)
    )

  tier3_joined <- dplyr::left_join(tier3_raw, n_dt_tbl, by = "dyad_id")

  tier3_joined %>%
    dplyr::group_by(dyad_id, n_dt) %>%
    dplyr::summarise(
      n_theta = dplyr::n(),
      min_ess_bulk = min(ess_bulk, na.rm = TRUE),
      min_ess_tail = min(ess_tail, na.rm = TRUE),
      max_rhat = max(rhat, na.rm = TRUE),
      share_ess_below_threshold = mean(ess_flag, na.rm = TRUE),
      share_rhat_above_threshold = mean(rhat_flag, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dyad_id)
}

#' Assemble a `bilatr_diagnostics` object from a tier-classified summary
#'
#' Shared by both branches of [diagnose_convergence()] (in-memory fit and
#' raw CmdStan CSV files): everything past "we now have a
#' `summarise_draws()`-shaped tibble with `tier`/`dyad_id`/`time_index`
#' columns joined on" is identical regardless of how that tibble was
#' produced, so it lives here once rather than being duplicated per
#' branch.
#'
#' @param summ A tibble as produced by [posterior::summarise_draws()],
#'   left-joined with [.classify_bilatr_tier]'s `tier`/`dyad_id`/
#'   `time_index` columns.
#' @param n_dt_tbl Output of [.normalize_n_dt], or `NULL` if `tiers`
#'   excludes both 2 and 3.
#' @param tiers Validated (via [.validate_tiers]) tiers to compute.
#' @param rhat_threshold,ess_threshold See [diagnose_convergence()].
#' @return A list of class `bilatr_diagnostics`; see
#'   [diagnose_convergence()]'s `@return` for the element-by-element
#'   description.
#' @keywords internal
.assemble_bilatr_diagnostics <- function(summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold) {
  tier1 <- if (1L %in% tiers) .compute_tier1(summ, rhat_threshold, ess_threshold) else NULL
  tier2_result <- if (2L %in% tiers) .compute_tier2(summ, n_dt_tbl) else NULL
  tier2 <- tier2_result$tier2
  tier3 <- if (3L %in% tiers) .compute_tier3(summ, n_dt_tbl, rhat_threshold, ess_threshold) else NULL

  summary_info <- list(
    tiers_computed = tiers,
    n_tier1_flagged = if (!is.null(tier1)) sum(tier1$flagged, na.rm = TRUE) else NA_integer_,
    n_tier1_total = if (!is.null(tier1)) nrow(tier1) else NA_integer_,
    n_dyads_tier2 = if (!is.null(tier2)) nrow(tier2) else NA_integer_,
    n_dyads_tier2_worse_than_expected = if (!is.null(tier2)) sum(tier2$worse_than_expected, na.rm = TRUE) else NA_integer_,
    n_dyads_missing_from_n_dt = if (!is.null(tier2_result)) tier2_result$n_dyads_missing_from_n_dt else NA_integer_,
    n_dyads_missing_from_draws = if (!is.null(tier2_result)) tier2_result$n_dyads_missing_from_draws else NA_integer_,
    n_dyads_tier3 = if (!is.null(tier3)) nrow(tier3) else NA_integer_,
    n_dyads_tier3_min_ess_below_threshold = if (!is.null(tier3)) {
      sum(tier3$min_ess_bulk < ess_threshold | tier3$min_ess_tail < ess_threshold, na.rm = TRUE)
    } else {
      NA_integer_
    },
    rhat_threshold = rhat_threshold,
    ess_threshold = ess_threshold
  )

  structure(
    list(tier1 = tier1, tier2 = tier2, tier3 = tier3, summary = summary_info),
    class = "bilatr_diagnostics"
  )
}

#' Convert a CmdStan raw CSV dot-index name to posterior bracket-index form
#'
#' `"theta.3.12"` -> `"theta[3,12]"`; `"lp__"` (no indices) -> `"lp__"`
#' unchanged. Safe because Stan identifiers cannot themselves contain a
#' dot, so the first dot-separated segment is always the base name and
#' any remaining segments are always indices.
#'
#' @param name A single CmdStan raw CSV column name.
#' @return The equivalent posterior/bracket-style name.
#' @keywords internal
.dot_name_to_bracket <- function(name) {
  parts <- strsplit(name, ".", fixed = TRUE)[[1]]
  if (length(parts) == 1) {
    return(name)
  }
  paste0(parts[1], "[", paste(parts[-1], collapse = ","), "]")
}

#' Get posterior-style variable names from a CmdStan CSV without reading
#' any draws
#'
#' Uses `cmdstanr:::read_csv_metadata()` (unexported -- if a future
#' cmdstanr release changes or removes it, this needs revisiting), which
#' scans a CSV's header/comment lines only. This bounds R's memory
#' regardless of file size (draws are never touched), though not
#' necessarily wall-clock time on a very large file, since the scan is
#' still sequential over the whole file on disk.
#'
#' @param csv_file A single CmdStan CSV file path (any one chain's file;
#'   variable structure is identical across chains of the same run).
#' @return Character vector of variable names in posterior/bracket form.
#' @keywords internal
.stan_csv_variable_names <- function(csv_file) {
  meta <- cmdstanr:::read_csv_metadata(csv_file)
  vapply(meta$variables, .dot_name_to_bracket, character(1), USE.NAMES = FALSE)
}

#' Cheaply get per-chain draw count and chain count from CmdStan CSVs
#'
#' `n_draws` comes from the same header-only metadata read as
#' [.stan_csv_variable_names]; `n_chains` is simply `length(csv_files)`
#' (one file per chain, per this package's own SLURM submission
#' convention -- see `runscripts/submit_bilatr_runs.R`), not the
#' `num_chains` metadata field, which records only the single chain each
#' individual CSV file's own run was configured for.
#'
#' @param csv_files Character vector of CmdStan CSV file paths.
#' @return A list with `n_draws` and `n_chains`.
#' @keywords internal
.stan_csv_dims <- function(csv_files) {
  meta <- cmdstanr:::read_csv_metadata(csv_files[1])
  list(n_draws = meta$iter_sampling, n_chains = length(csv_files))
}

#' Derive a Tier 3 chunk size (variables per chunk) from a memory budget
#'
#' Solves, for `chunk_size`, the same memory model
#' [.estimate_diagnostics_memory_mb] reports in the other direction:
#' peak memory (MB) ~= `n_draws * n_chains * chunk_size * 8 bytes *
#' overhead_factor * (n_workers if parallel else 1) / 1e6`.
#' `overhead_factor` (default `2`) is a fixed safety margin for
#' `read_cmdstan_csv()`'s intermediate structures and
#' `summarise_draws()`'s own working memory, which a raw
#' 8-bytes-per-double count understates -- this is meant as a
#' sanity-check number, not a guarantee (see `max_memory_mb` in
#' [diagnose_convergence()]'s documentation).
#'
#' @param n_draws,n_chains From [.stan_csv_dims].
#' @param max_memory_mb See [diagnose_convergence()].
#' @param n_workers,parallel See [diagnose_convergence()].
#' @param overhead_factor Fixed safety multiplier; not user-facing.
#' @return Integer chunk size, at least `1`.
#' @keywords internal
.compute_chunk_size <- function(n_draws, n_chains, max_memory_mb, n_workers, parallel, overhead_factor = 2) {
  denom <- n_draws * n_chains * 8 * overhead_factor * (if (parallel) n_workers else 1)
  chunk_size <- floor((max_memory_mb * 1e6) / denom)
  if (chunk_size < 1) {
    warning(
      "max_memory_mb (", max_memory_mb, ") is too tight to fit even 1 ",
      "variable per chunk at this n_draws/n_chains/n_workers combination; ",
      "using chunk_size = 1. Per-chunk read overhead will dominate runtime.",
      call. = FALSE
    )
    chunk_size <- 1L
  }
  as.integer(chunk_size)
}

#' Estimate peak memory (MB) for a given Tier 3 chunk size
#'
#' The same memory model as [.compute_chunk_size], used in the other
#' direction (chunk size known, want the resulting estimate) for the
#' pre-flight `message()` in [diagnose_convergence()].
#'
#' @inheritParams .compute_chunk_size
#' @param chunk_size Variables per chunk.
#' @return Estimated peak memory in MB.
#' @keywords internal
.estimate_diagnostics_memory_mb <- function(n_draws, n_chains, chunk_size, n_workers, parallel, overhead_factor = 2) {
  bytes <- n_draws * n_chains * chunk_size * 8 * overhead_factor * (if (parallel) n_workers else 1)
  bytes / 1e6
}

#' Read and summarise CSV variables in memory-bounded chunks
#'
#' Shared by [diagnose_convergence()]'s CSV-path branch and
#' [extract_theta()]'s: reads and summarises `variables` in groups of
#' `chunk_size`, discarding each chunk's draws before moving to the next
#' -- this is what keeps peak memory bounded regardless of how many
#' variables there are in total. `parallel = TRUE` processes chunks
#' concurrently via `furrr::future_map_dfr()`, trading the sequential
#' path's memory bound (now multiplied by `n_workers`, since that many
#' chunks are in memory at once) for wall-clock speed. The active
#' `future::plan()` is saved and restored on exit, so this never
#' permanently changes the caller's parallel backend.
#'
#' Reads via [.fast_read_post_warmup_draws] against `prepared`'s
#' comment-stripped files (see [.prepare_fast_csv_read]), not
#' [cmdstanr::read_cmdstan_csv()] directly: the latter's `cmd = grep`
#' read re-materializes a full near-complete copy of each chain's raw
#' CSV on *every* call (see [.prepare_fast_csv_read]'s docs), which
#' would otherwise happen once per chunk here -- exactly the redundant,
#' memory-model-breaking cost this function's `max_memory_mb`/
#' `chunk_size` accounting is meant to avoid.
#'
#' Benchmarked against a real production-scale single-chain CSV (~1.9M
#' Tier 3 columns): reading is strongly I/O-bound, not parsing-bound --
#' per-call wall-time barely depends on how many `variables` are
#' requested (a 1000x range in chunk size changed per-call time by under
#' 10%), because extracting even one column from a row-oriented CSV
#' requires scanning the full row width regardless of how many fields
#' are kept. Consequence: total sweep time scales with the NUMBER OF
#' CHUNKS, not with memory saved -- a small `chunk_size` does not make
#' this cheaper, it makes it much slower for a given `variables` list,
#' potentially by orders of magnitude. See `max_memory_mb` in
#' [diagnose_convergence()]'s documentation, which this benchmark
#' motivated.
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param variables Character vector of variable names to summarise.
#' @param chunk_size Variables per chunk.
#' @param parallel,n_workers See [diagnose_convergence()].
#' @param flip If `TRUE`, negate each chunk's raw draws before
#'   summarising (used by [extract_theta()]'s CSV path to apply
#'   [bilatr_orient()]'s sign correction; always `FALSE` for
#'   [diagnose_convergence()], since Rhat/ESS are invariant to a
#'   deterministic sign flip and it would be pointless work there).
#'   Applied pre-summary, per chunk, matching how [bilatr_orient()]
#'   flips raw draws for the in-memory path -- not a post-hoc
#'   transformation of the summary columns, which would need to swap
#'   the quantile columns (`quantile(-X, p) == -quantile(X, 1 - p)`),
#'   not just negate them.
#' @return A tibble, the row-bound [posterior::summarise_draws()] output
#'   across all chunks.
#' @keywords internal
.chunked_summarise_csv <- function(prepared, variables, chunk_size, parallel, n_workers, flip = FALSE) {
  chunks <- split(variables, ceiling(seq_along(variables) / chunk_size))

  summarise_one_chunk <- function(chunk_vars) {
    draws <- .fast_read_post_warmup_draws(prepared, chunk_vars, flip = flip)
    posterior::summarise_draws(draws)
  }

  if (!parallel) {
    return(purrr::map_dfr(chunks, summarise_one_chunk))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (.Platform$OS.type == "unix") {
    future::plan(future::multicore, workers = n_workers)
  } else {
    warning(
      "parallel = TRUE on Windows falls back to future::multisession, ",
      "which copies data to each worker rather than sharing it via ",
      "copy-on-write (unlike future::multicore on Unix-like systems); ",
      "peak memory will run higher than the max_memory_mb estimate ",
      "assumes. Consider parallel = FALSE if memory is tight.",
      call. = FALSE
    )
    future::plan(future::multisession, workers = n_workers)
  }

  furrr::future_map_dfr(chunks, summarise_one_chunk)
}

#' Resolve a chunk size for a CSV-file-path chunked read, and report it
#'
#' Shared by [diagnose_convergence()]'s and [extract_theta()]'s CSV-path
#' branches. If `max_memory_mb` was left at its default, emits a one-time
#' `message()` naming the I/O-bound wall-time/memory tradeoff a real
#' benchmark against a production-scale CSV found (see
#' [.chunked_summarise_csv]): a smaller `chunk_size` does NOT make
#' reading cheaper -- it multiplies wall-time roughly by the number of
#' chunks, since each chunk pays nearly the same full-row-scan cost
#' regardless of how many columns it keeps. A second `message()`, always
#' emitted, reports the resolved chunk count/size and estimated peak
#' memory, so a caller can see the actual tradeoff being made before a
#' long run commits to it.
#'
#' @param n_vars Number of variables the chunked sweep will cover (for
#'   the reporting message only; does not affect the chunk_size
#'   calculation itself).
#' @param csv_files,max_memory_mb,chunk_size,parallel,n_workers See
#'   [diagnose_convergence()].
#' @param max_memory_mb_missing Whether the caller left `max_memory_mb`
#'   at its default (via `missing()` in the calling function).
#' @return The resolved integer chunk size.
#' @keywords internal
.resolve_chunk_size_and_report <- function(
  n_vars, csv_files, max_memory_mb, chunk_size, parallel, n_workers,
  max_memory_mb_missing
) {
  dims <- .stan_csv_dims(csv_files)

  if (max_memory_mb_missing) {
    message(
      "Using the default max_memory_mb = ", max_memory_mb, " (",
      round(max_memory_mb / 1024, 1), " GB). Reading is I/O-bound, not ",
      "parsing-bound: a benchmark against a production-scale CSV found ",
      "per-chunk wall-time barely depends on how many variables are ",
      "requested (a 1000x range in chunk size changed per-call time by ",
      "under 10%), because extracting even one column from a row-",
      "oriented CSV means scanning the full row regardless of how much ",
      "of it is kept. So chunk COUNT, not chunk size, drives total ",
      "wall-time: a smaller max_memory_mb produces more chunks and can ",
      "multiply total time by orders of magnitude for a modest memory ",
      "saving. Prefer the LARGEST max_memory_mb your job's memory ",
      "allocation can afford; only lower it if memory, not time, is the ",
      "binding constraint."
    )
  }

  chunk_size_used <- chunk_size %||% .compute_chunk_size(
    n_draws = dims$n_draws, n_chains = dims$n_chains,
    max_memory_mb = max_memory_mb, n_workers = n_workers, parallel = parallel
  )

  est_mb <- .estimate_diagnostics_memory_mb(
    n_draws = dims$n_draws, n_chains = dims$n_chains,
    chunk_size = chunk_size_used, n_workers = n_workers, parallel = parallel
  )
  message(
    n_vars, " variable(s) in ",
    ceiling(n_vars / chunk_size_used), " chunk(s) of ",
    chunk_size_used, " variable(s) each; estimated peak memory ~",
    round(est_mb), " MB",
    if (parallel) paste0(" across ", n_workers, " worker(s)") else "", "."
  )

  chunk_size_used
}

#' Build the tier-classified summary tibble directly from raw CmdStan CSVs
#'
#' The CSV-path counterpart of the in-memory branch in
#' [diagnose_convergence()]: reads only Tier 1/2 variables in one small
#' read (cheap, as today), and Tier 3 variables (typically, by far, the
#' most numerous of the three tiers) in memory-bounded chunks via
#' [.chunked_summarise_csv]. Never materializes the full multi-chain
#' draws array in memory, unlike the in-memory branch, which necessarily
#' receives an already-fully-read `fit`.
#'
#' @param csv_files Character vector of CmdStan CSV file paths.
#' @param tiers,max_memory_mb,chunk_size,parallel,n_workers See
#'   [diagnose_convergence()].
#' @param max_memory_mb_missing Whether the caller left `max_memory_mb`
#'   at its default (via `missing()` in [diagnose_convergence()]) --
#'   gates the one-time "this is a default, not a calibrated value"
#'   message.
#' @param scratch_dir See [diagnose_convergence()].
#' @return A tibble in the same shape [diagnose_convergence()]'s
#'   in-memory branch produces: [posterior::summarise_draws()] columns
#'   left-joined with [.classify_bilatr_tier]'s `tier`/`dyad_id`/
#'   `time_index`.
#' @keywords internal
.read_diagnostics_summary_from_csv <- function(
  csv_files, tiers, max_memory_mb, chunk_size, parallel, n_workers,
  max_memory_mb_missing, scratch_dir = NULL
) {
  prepared <- .prepare_fast_csv_read(csv_files, scratch_dir)
  on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)

  all_vars <- .stan_csv_variable_names(csv_files[1])
  var_tiers <- .classify_bilatr_tier(all_vars)
  keep_tiers <- var_tiers[var_tiers$tier %in% tiers, ]

  tier12_vars <- keep_tiers$variable[keep_tiers$tier %in% c(1L, 2L)]
  tier3_vars <- keep_tiers$variable[keep_tiers$tier == 3L]

  tier12_summ <- if (length(tier12_vars) > 0) {
    draws <- .fast_read_post_warmup_draws(prepared, tier12_vars)
    posterior::summarise_draws(draws)
  } else {
    NULL
  }

  tier3_summ <- if (length(tier3_vars) > 0) {
    chunk_size_used <- .resolve_chunk_size_and_report(
      length(tier3_vars), csv_files, max_memory_mb, chunk_size, parallel, n_workers,
      max_memory_mb_missing
    )
    .chunked_summarise_csv(prepared, tier3_vars, chunk_size_used, parallel, n_workers)
  } else {
    NULL
  }

  summ_raw <- dplyr::bind_rows(tier12_summ, tier3_summ)
  dplyr::left_join(summ_raw, var_tiers, by = "variable")
}

#' Triage MCMC convergence diagnostics by parameter tier
#'
#' Runs [posterior::summarise_draws()] over a fitted `bilatr` model and
#' splits the result into three diagnostic tiers before applying Rhat/ESS
#' thresholds, so that pathologies in global/shared parameters are never
#' masked by the wide, poorly-identified posteriors expected for
#' sparsely-observed dyads. See `vignette("diagnostics")` for a worked
#' example and the rationale behind the tiering.
#'
#' Tier 1 (global/shared parameters plus `lp__`) is always reported in
#' full and never summarized away. Tier 2 (per-dyad hierarchical
#' parameters, e.g. `phi`, `process_noise`, `theta0`) and Tier 3
#' (per-dyad-period latent states, e.g. `theta`, `theta_raw`) are joined
#' against `n_dt` and screened for dyads whose diagnostics are worse than
#' a smooth degradation-with-sparsity trend predicts (see
#' [.flag_worse_than_expected]), so sparse-but-unremarkable dyads don't
#' flood the report. See [.classify_bilatr_tier] for how parameter names
#' are assigned to tiers.
#'
#' `tiers` limits *which* of these are computed at all: quantities not
#' assigned to a requested tier are never read into memory in the first
#' place (see `fit` below), rather than being read and then hidden. This
#' matters in practice because Tier 3 (`theta`/`theta_raw`) typically
#' has, by far, the most monitored quantities of the three tiers (one per
#' dyad-period); requesting only `tiers = 1` or `tiers = 1:2` skips
#' reading/computing Rhat/ESS for all of them.
#'
#' @param fit Either (a) a `CmdStanMCMC`/`CmdStanFit`-like fit object
#'   (anything with a `$draws()` method) or a `posterior::draws_array`/
#'   `draws_df` -- the whole object is already in memory, so `tiers`
#'   controls what gets summarised but not what gets read, and
#'   `max_memory_mb`/`chunk_size`/`parallel`/`n_workers` are unused; or
#'   (b) a character vector of raw CmdStan CSV file paths (one per
#'   chain, e.g. from a completed SLURM run never loaded into this R
#'   session) -- in this case Tier 1/2 variables are read in one small
#'   call, and Tier 3 variables (typically the overwhelming majority for
#'   production-sized panels) are read and summarised in memory-bounded
#'   chunks via [cmdstanr::read_cmdstan_csv()]'s `variables` argument,
#'   discarding each chunk's draws before moving to the next, so the
#'   full draws array is never materialized at once. This is the path
#'   that matters for production Stan output too large to read normally
#'   (tested against panels with millions of Tier 3 columns).
#' @param n_dt A data frame with a dyad-id column (matching `"dyad..."`)
#'   and an observation-count column (matching `"n_dt"`/`"n_obs"`/
#'   `"n_events"`), or a named numeric vector of counts keyed by dyad id.
#'   Only required when `tiers` includes `2` and/or `3` (Tier 1 has no
#'   per-dyad structure to join against `n_dt`).
#' @param rhat_threshold Rhat values strictly above this are flagged.
#'   Defaults to `1.01` (Vehtari et al. 2021).
#' @param ess_threshold ESS_bulk/ESS_tail values strictly below this are
#'   flagged. Defaults to `400` (Vehtari et al. 2021).
#' @param tiers Integer vector, a subset of `1:3`, naming which tier(s) to
#'   compute: `1` (global/shared), `2` (per-dyad hierarchical parameters),
#'   `3` (per-dyad-period latent states). Defaults to `1:3` (all tiers).
#'   A tier not requested is left as `NULL` in the returned object rather
#'   than an empty tibble, so `is.null(diag$tier3)` distinguishes "not
#'   computed" from "computed, nothing to report".
#' @param max_memory_mb Only used when `fit` is CSV file paths and `tiers`
#'   includes `3`. Target ceiling, in MB, for Tier 3's peak memory;
#'   drives the automatically-derived `chunk_size` (see `chunk_size`
#'   below) so you don't have to guess a variable count yourself.
#'   Defaults to `8192` (8 GB) -- a guess, not a calibration against your
#'   hardware, and this function says so via `message()` the first time
#'   you rely on that default rather than setting it explicitly. This is
#'   a sanity-check number, not a guarantee: actual peak memory depends
#'   on `read_cmdstan_csv()`/`summarise_draws()` internals this function
#'   doesn't control. IMPORTANT, from a real benchmark against a
#'   production-scale CSV (see [.chunked_summarise_csv]): reading is
#'   strongly I/O-bound, so a SMALLER `max_memory_mb` (more, smaller
#'   chunks) does not make this cheaper -- it can make it dramatically
#'   *slower*, since each chunk pays nearly the same full-row-scan cost
#'   regardless of how many columns it keeps, and total wall-time scales
#'   with chunk count. Prefer the LARGEST `max_memory_mb` your job's
#'   memory allocation can afford; only lower it if memory, not time, is
#'   the binding constraint.
#' @param chunk_size Only used when `fit` is CSV file paths and `tiers`
#'   includes `3`. Explicit override: number of Tier 3 variables read per
#'   chunk. `NULL` (the default) derives this from `max_memory_mb`
#'   instead; set this directly only if you want to bypass that
#'   calculation (e.g. you've measured actual memory use and want to
#'   tune it by hand).
#' @param parallel Only used when `fit` is CSV file paths and `tiers`
#'   includes `3`. `FALSE` (default) processes Tier 3 chunks
#'   sequentially, one at a time -- this is what makes `max_memory_mb`'s
#'   bound hold regardless of `n_workers`. `TRUE` processes chunks
#'   concurrently via `furrr::future_map_dfr()`, trading that memory
#'   bound (now effectively `max_memory_mb` times up to `n_workers`,
#'   since that many chunks are in memory at once) for wall-clock speed.
#'   Prefer `parallel = FALSE` when memory is already tight (e.g. a
#'   memory-constrained HPC allocation) and `parallel = TRUE` when you
#'   have memory headroom to spend on speed instead. On Windows this
#'   falls back from `future::multicore` to `future::multisession` with
#'   a loud `warning()`, since `multisession` copies data to each worker
#'   rather than sharing it via copy-on-write -- the memory math changes
#'   substantially, and `max_memory_mb` is less trustworthy there.
#' @param n_workers Only used when `parallel = TRUE`. Defaults to
#'   `max(1, parallel::detectCores() - 1)`.
#' @param scratch_dir Only used when `fit` is CSV file paths. Directory
#'   to write a one-time, comment-stripped copy of each chain file into
#'   before reading (see [.prepare_fast_csv_read]). `NULL` (the default)
#'   writes each copy alongside its source file, which is deliberate:
#'   that's already known-good storage for these (often many-GB) files,
#'   unlike `tempdir()`/`$TMPDIR`, which on some HPC systems is a
#'   RAM-backed `tmpfs` -- silently turning a routine disk read into a
#'   direct hit against the job's memory allocation. Point this at an
#'   explicit path only if the CSVs' own directory is unsuitable (e.g.
#'   quota-constrained, read-only, or a slow network filesystem when
#'   faster local scratch is available).
#' @return A list of class `bilatr_diagnostics` with elements:
#'   \describe{
#'     \item{tier1}{Tibble of global/shared diagnostics, one row per
#'       monitored quantity, with a `flagged` column; `NULL` if `1` was
#'       not in `tiers`.}
#'     \item{tier2}{Tibble with one row per dyad found in `n_dt` (or in
#'       the draws, if unmatched), per-dyad-parameter Rhat/ESS columns,
#'       and a `worse_than_expected` column; `NULL` if `2` was not in
#'       `tiers`.}
#'     \item{tier3}{Tibble with one row per dyad summarizing its
#'       per-dyad-period latent-state diagnostics (min ESS, share of
#'       entries breaching thresholds); `NULL` if `3` was not in `tiers`.}
#'     \item{summary}{A short named list of headline counts (see
#'       [print.bilatr_diagnostics]).}
#'   }
#' @examples
#' \dontrun{
#' diag <- diagnose_convergence(fit, n_dt = dplyr::count(events, dyad, wt = 1))
#' diag
#' diag$tier1
#'
#' # only the cheap, always-important global/shared parameters, and no
#' # need to supply n_dt at all:
#' diagnose_convergence(fit, tiers = 1)
#'
#' # a completed SLURM run, never read into this R session: Tier 3 is
#' # read and summarised in memory-bounded chunks rather than all at once
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' diagnose_convergence(csv_files, n_dt = n_dt, max_memory_mb = 4096)
#' }
#' @export
diagnose_convergence <- function(
  fit,
  n_dt = NULL,
  rhat_threshold = 1.01,
  ess_threshold = 400,
  tiers = 1:3,
  max_memory_mb = 8192,
  chunk_size = NULL,
  parallel = FALSE,
  n_workers = max(1L, parallel::detectCores() - 1L),
  scratch_dir = NULL
) {
  max_memory_mb_missing <- missing(max_memory_mb)

  tiers <- .validate_tiers(tiers)
  if (any(c(2L, 3L) %in% tiers) && is.null(n_dt)) {
    stop(
      "`n_dt` is required when `tiers` includes 2 and/or 3 (Tier 1 alone ",
      "needs no per-dyad join).",
      call. = FALSE
    )
  }
  n_dt_tbl <- if (any(c(2L, 3L) %in% tiers)) .normalize_n_dt(n_dt) else NULL

  if (is.character(fit)) {
    summ <- .read_diagnostics_summary_from_csv(
      fit, tiers, max_memory_mb, chunk_size, parallel, n_workers, max_memory_mb_missing,
      scratch_dir
    )
  } else {
    draws <- if (posterior::is_draws(fit)) fit else fit$draws()
    var_tiers <- .classify_bilatr_tier(posterior::variables(draws))
    keep_vars <- var_tiers$variable[var_tiers$tier %in% tiers]
    summ <- posterior::summarise_draws(posterior::subset_draws(draws, variable = keep_vars))
    summ <- dplyr::left_join(summ, var_tiers, by = "variable")
  }

  .assemble_bilatr_diagnostics(summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold)
}

#' Print a `bilatr_diagnostics` object
#'
#' Tier 1 (global/shared parameters) is always printed in full, since it
#' should never be silently summarized away. Tier 2 (per-dyad
#' hierarchical parameters) is printed as a compact table sorted with
#' dyads flagged as "worse than expected for their sparsity" first. Tier
#' 3 (per-dyad-period latent states) is expected to be noisy for sparse
#' dyads, so it is reported only as aggregate one-line statistics rather
#' than flooding the console with per-dyad-period rows.
#'
#' @param x A `bilatr_diagnostics` object, as returned by
#'   [diagnose_convergence()].
#' @param n_tier2 Maximum number of Tier 2 rows to print.
#' @param ... Ignored; present for S3 consistency.
#' @return `x`, invisibly.
#' @export
print.bilatr_diagnostics <- function(x, n_tier2 = 20, ...) {
  cat("<bilatr_diagnostics>\n\n")

  if (1L %in% x$summary$tiers_computed) {
    cat(sprintf(
      "== Tier 1: global/shared parameters (%d/%d flagged) ==\n",
      x$summary$n_tier1_flagged, x$summary$n_tier1_total
    ))
    if (x$summary$n_tier1_flagged > 0) {
      print(dplyr::filter(x$tier1, flagged), n = Inf)
    } else {
      cat("No Tier 1 issues: all global/shared parameters (and lp__) meet threshold.\n")
    }
    cat("\n")
  }

  if (2L %in% x$summary$tiers_computed) {
    cat(sprintf(
      "== Tier 2: per-dyad hierarchical parameters (%d dyads, %d worse than expected for their n_dt) ==\n",
      x$summary$n_dyads_tier2, x$summary$n_dyads_tier2_worse_than_expected
    ))
    if (x$summary$n_dyads_missing_from_n_dt > 0) {
      cat(sprintf(
        "Note: %d dyad(s) in the draws have no matching entry in `n_dt`.\n",
        x$summary$n_dyads_missing_from_n_dt
      ))
    }
    if (x$summary$n_dyads_missing_from_draws > 0) {
      cat(sprintf(
        "Note: %d dyad(s) in `n_dt` have no matching Tier 2 parameters in the draws.\n",
        x$summary$n_dyads_missing_from_draws
      ))
    }
    print(utils::head(x$tier2, n_tier2), n = Inf)
    cat("\n")
  }

  if (3L %in% x$summary$tiers_computed) {
    share_below <- x$summary$n_dyads_tier3_min_ess_below_threshold / max(x$summary$n_dyads_tier3, 1)
    low_n_dt_share <- tryCatch({
      flagged <- x$tier3$min_ess_bulk < x$summary$ess_threshold | x$tier3$min_ess_tail < x$summary$ess_threshold
      if (any(flagged, na.rm = TRUE)) {
        stats::median(x$tier3$n_dt[flagged], na.rm = TRUE)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)

    cat(sprintf(
      "== Tier 3: per-dyad-period latent states (aggregated; %d dyads) ==\n",
      x$summary$n_dyads_tier3
    ))
    cat(sprintf(
      "%.0f%% of dyads have min theta ESS below threshold (%s)%s.\n",
      100 * share_below,
      x$summary$ess_threshold,
      if (!is.na(low_n_dt_share)) {
        sprintf(", concentrated around dyads with n_dt <= %.0f", low_n_dt_share)
      } else {
        ""
      }
    ))
  }

  invisible(x)
}
