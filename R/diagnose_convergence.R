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

#' Memory-model constant: a fixed floor for R and its loaded packages
#'
#' Unlike every other `.BILATR_CHUNK_*` constant, this is not a
#' multiplier of `raw_mb` -- it is a fixed number of MB, added once
#' alongside `file_mb` in [.compute_chunk_size()]/
#' [.estimate_diagnostics_memory_mb()], because it does not shrink with
#' chunk size any more than `file_mb` does: it is the cost of `library(
#' bilatr)` and its dependencies (`cmdstanr`, `data.table`, `posterior`,
#' `dplyr`, etc.) simply being loaded, before a single byte of any CSV
#' is read. Measured via `dev/bench_memory.R` (2026-09; `devtools::
#' load_all()` alone, R 4.5.2, macOS): ~160 MB peak RSS; set to 200 here
#' for margin against a heavier dependency set on a different R/package
#' version. If your own session or job already has substantially more
#' loaded before calling into this package, this floor will
#' under-count it -- it models a fresh `Rscript` invocation, which is
#' what every SLURM job in this package's own `runscripts/` does.
#' @keywords internal
.BILATR_CHUNK_BASELINE_MB <- 200

#' Memory-model constant: the preallocated result array, the fread parse
#' buffer, and the one chain matrix in flight, combined
#'
#' See [.bilatr_chunk_overhead_multiplier()] for how this and the other
#' `.BILATR_CHUNK_*_FACTOR` constants combine into the `k` multiplier
#' [.compute_chunk_size()]/[.estimate_diagnostics_memory_mb()] apply to a
#' chunk's raw byte count (`raw_mb = n_draws * n_chains * chunk_size * 8
#' / 1e6`). Re-derived (2026-09, `dev/bench_memory.R`, R 4.5.2, macOS)
#' against the array-native flip (no more `as_draws_df()` round-trip --
#' see `NEWS.md`), from two measurements at `n_cores = 1`, no flip: a
#' single-chain, 60,000-Tier-3-column, 300-draw synthetic CSV
#' (`raw_mb = 144`, `file_mb` ~= 202) measured 1002 MB peak RSS, and the
#' same shape at 2 chains (`raw_mb = 288`) measured 1409 MB. Subtracting
#' [.BILATR_CHUNK_BASELINE_MB] (measured separately) and `file_mb` from
#' each and solving the two-chain-count system for `ARRAY +
#' PER_CHAIN / n_chains` gives `ARRAY ~= 2.8`, `PER_CHAIN ~= 1.6`;
#' rounded up here for margin. Read alone (before Step 1's flip fix)
#' measured close to the old estimate (the read chain itself is
#' unchanged by that fix), so this jump from the old `ARRAY = 1` is
#' `summarise_draws()`'s own working set (rank-normalisation, split-Rhat,
#' bulk/tail ESS all build same-shaped intermediate arrays), which the
#' pre-refactor model never isolated from the read.
#' @keywords internal
.BILATR_CHUNK_ARRAY_FACTOR <- 3

#' Memory-model constant: `fread()`'s parse buffer and the one chain
#' matrix in flight, together sized for a single chain
#'
#' Applied as `.BILATR_CHUNK_PER_CHAIN_FACTOR / n_chains` (each term is
#' `n_draws * chunk_size * 8`, i.e. `raw_mb / n_chains`, not `raw_mb`
#' itself) -- see [.bilatr_chunk_overhead_multiplier()] and
#' [.BILATR_CHUNK_ARRAY_FACTOR]'s docs for the joint derivation.
#' @keywords internal
.BILATR_CHUNK_PER_CHAIN_FACTOR <- 2

#' Memory-model constant: the post-Step-1 array-native sign flip's
#' residual cost
#'
#' Applied unconditionally (not just when a flip actually happens):
#' whether `alpha[1]`'s posterior median is negative isn't known until
#' after the first chunk is read, so [.compute_chunk_size()] must budget
#' for the possibility on every call, not just the calls that end up
#' needing it. Since [.fast_read_post_warmup_draws()]'s `flip_vars`
#' negates `flip_pos` columns of `m` in place before it's assigned into
#' `arr` (no extra copy there), and [.chunked_summarise_csv_with_orientation()]/
#' [.read_and_orient_draws_summary()]'s array-native
#' `draws[, , flip_cols] <- -draws[, , flip_cols, drop = FALSE]` negates
#' only the matched subset of an already-materialized array, the
#' remaining cost is R's own copy-on-modify for that replacement.
#' Measured (2026-09, `dev/bench_memory.R`): flip added ~91-142 MB over
#' otherwise-identical no-flip runs (~0.6-1.0x `raw_mb`) -- far below
#' the ~3.9x `raw_mb` the pre-Step-1 `as_draws_df()` round-trip cost
#' (see `dev/refactor_verification_2026-09-10.md` section 1), but not
#' zero, so (per that document's own instruction) a term stays here
#' rather than assuming it away.
#' @keywords internal
.BILATR_CHUNK_FLIP_FACTOR <- 1

#' Memory-model constant: the one-time cost of `posterior::
#' summarise_draws()` forking at all (`.cores > 1`)
#'
#' Applied once, whenever `n_cores > 1`, together with
#' [.BILATR_CHUNK_CORES_PER_CORE_FACTOR] -- see
#' [.bilatr_chunk_overhead_multiplier()]. Measured, not assumed: see
#' that constant's docs for why this and the per-core term replace the
#' old single flat `CORES_FACTOR` (which this package's own pinned
#' tests, before this measurement, asserted was identical from
#' `n_cores = 2` through `24` -- that assumption was never measured
#' against a real fork and was wrong).
#' @keywords internal
.BILATR_CHUNK_CORES_STEP_FACTOR <- 4.5

#' Memory-model constant: the additional cost of each forked worker
#' beyond the first, when `.cores > 1`
#'
#' `posterior::summarise_draws(.cores = k)` forks `k` worker processes
#' via `parallel::mclapply()`, each computing rank-normalised Rhat/bulk
#' ESS/tail ESS (each of which builds several same-shaped intermediate
#' arrays) over its own slice of variables. This is NOT a flat cost
#' regardless of `n_cores`, as the pre-measurement model assumed:
#' measured (2026-09, `dev/bench_memory.R`, single chain, 60,000
#' Tier-3 columns, 300 draws, no flip, peak RSS summed across the whole
#' process tree via the `ps` package -- `/usr/bin/time`'s own child-RSS
#' accounting does not see a further-forked grandchild's memory, and
#' was confirmed empirically to under-count this exact case) at
#' `n_cores` = 1/2/4/8: 1002 / 1822 / 2378 / 3571 MB. Fitting `extra
#' beyond n_cores = 1 baseline = STEP + PER_CORE * (n_cores - 1)` in
#' units of `raw_mb` (144 MB here) across the 2->4 and 4->8 intervals
#' gives `PER_CORE ~= 1.9-2.1`, consistently; `STEP ~= 3.0-3.3`. Both
#' rounded up here for margin. Summing RSS across a fork tree likely
#' over-counts pages the OS still shares copy-on-write between workers
#' (a SLURM cgroup would charge shared pages once, not once per
#' process), so this constant is a conservative upper bound, not a
#' precise physical-memory figure -- an acceptable direction of error
#' for a job-sizing floor, given the alternative is under-counting and
#' risking exactly the OOM this package's memory model exists to
#' prevent. **Only measured up to `n_cores = 8`**: a caller sizing a job
#' with `n_workers` well beyond that (the audit's own per-chain example
#' used 24) is extrapolating this linear term 3x past its measured
#' range -- validate with `dev/bench_memory.R` at your actual `n_cores`
#' before trusting the estimate for a real SLURM allocation that large.
#' @keywords internal
.BILATR_CHUNK_CORES_PER_CORE_FACTOR <- 2.5

#' The `k` memory-model multiplier [.compute_chunk_size()]/
#' [.estimate_diagnostics_memory_mb()] apply to a chunk's raw byte count
#'
#' The post-[.fast_read_post_warmup_draws] pipeline holds, per chunk of
#' `v` variables: the preallocated `array(n_draws, n_chains, v)` plus
#' `summarise_draws()`'s own working set, both sized for the full
#' `n_chains` array (`.BILATR_CHUNK_ARRAY_FACTOR`, ~3x `raw_mb`);
#' `fread()`'s parse buffer and one chain's resulting matrix, each sized
#' for ONE chain (`raw_mb / n_chains`, i.e.
#' `.BILATR_CHUNK_PER_CHAIN_FACTOR / n_chains`); the array-native sign
#' flip's residual cost, budgeted unconditionally
#' (`.BILATR_CHUNK_FLIP_FACTOR`, see its own docs for why); and, only
#' when `n_cores > 1`, `posterior::summarise_draws()`'s forked-worker
#' cost -- a one-time step
#' (`.BILATR_CHUNK_CORES_STEP_FACTOR`) plus a per-additional-worker term
#' (`.BILATR_CHUNK_CORES_PER_CORE_FACTOR * (n_cores - 1)`), NOT a flat
#' penalty regardless of the exact core count (see that constant's docs
#' -- this replaced an assumption that was never measured and was
#' wrong).
#'
#' @param n_chains From [.prepare_fast_csv_read].
#' @param n_cores See [diagnose_convergence()]'s `parallel`/`n_workers`.
#' @return The multiplier `k`, in units of `raw_mb`.
#' @keywords internal
.bilatr_chunk_overhead_multiplier <- function(n_chains, n_cores) {
  .BILATR_CHUNK_ARRAY_FACTOR +
    .BILATR_CHUNK_PER_CHAIN_FACTOR / n_chains +
    .BILATR_CHUNK_FLIP_FACTOR +
    (if (n_cores > 1) {
      .BILATR_CHUNK_CORES_STEP_FACTOR + .BILATR_CHUNK_CORES_PER_CORE_FACTOR * (n_cores - 1)
    } else {
      0
    })
}

#' Derive a Tier 3 chunk size (variables per chunk) from a memory budget
#'
#' Solves, for `chunk_size`, the same memory model
#' [.estimate_diagnostics_memory_mb] reports in the other direction:
#' peak memory (MB) ~= `baseline_mb + file_mb + raw_mb * k`, `raw_mb =
#' n_draws * n_chains * chunk_size * 8 / 1e6`, `k` from
#' [.bilatr_chunk_overhead_multiplier()]. `baseline_mb`
#' ([.BILATR_CHUNK_BASELINE_MB]) and `file_mb` are both fixed costs
#' `chunk_size` cannot shrink away -- `fread()` touches the largest
#' chain file's full byte range on every call, however few columns are
#' selected (see [.fast_read_post_warmup_draws]), and R plus its loaded
#' packages cost the same regardless of how much of any file is read --
#' so `max_memory_mb` must cover both before any budget is left for
#' variables at all; if it doesn't, this `stop()`s rather than silently
#' returning a `chunk_size` that can never keep the job under budget.
#'
#' @param n_draws,n_chains From [.prepare_fast_csv_read]'s
#'   `num_post_warmup_draws`/`n_chains`.
#' @param file_mb From [.prepare_fast_csv_read].
#' @param max_memory_mb See [diagnose_convergence()].
#' @param n_cores See [diagnose_convergence()]'s `parallel`/`n_workers`.
#' @return Integer chunk size, at least `1`.
#' @keywords internal
.compute_chunk_size <- function(n_draws, n_chains, file_mb, max_memory_mb, n_cores) {
  k <- .bilatr_chunk_overhead_multiplier(n_chains, n_cores)
  per_var_mb <- n_draws * n_chains * 8 * k / 1e6
  budget_mb <- max_memory_mb - file_mb - .BILATR_CHUNK_BASELINE_MB

  if (budget_mb < per_var_mb) {
    stop(
      "max_memory_mb (", round(max_memory_mb), " MB) leaves no room for ",
      "even one Tier 3 variable per chunk: the largest chain file (~",
      round(file_mb), " MB) is mapped in full during every read, and ~",
      .BILATR_CHUNK_BASELINE_MB, " MB is a fixed floor for R and its ",
      "loaded packages -- both before any variables are selected, ",
      "regardless of chunk_size. Allocate at least ",
      ceiling(file_mb + .BILATR_CHUNK_BASELINE_MB + per_var_mb),
      " MB for this job.",
      call. = FALSE
    )
  }

  as.integer(floor(budget_mb / per_var_mb))
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
.estimate_diagnostics_memory_mb <- function(n_draws, n_chains, chunk_size, n_cores, file_mb) {
  k <- .bilatr_chunk_overhead_multiplier(n_chains, n_cores)
  raw_mb <- n_draws * n_chains * chunk_size * 8 / 1e6
  .BILATR_CHUNK_BASELINE_MB + file_mb + raw_mb * k
}

#' Summarise draws with the fixed measure set the CSV-path functions use
#'
#' `mean`, [posterior::quantile2()] (default `q5`/`q95`), `median`,
#' `rhat`, `ess_bulk`, `ess_tail` -- every column `.compute_tier1()`/
#' `.compute_tier2()`/`.compute_tier3()` (which each select their own
#' fixed subset regardless of what else is present) and
#' `diagnose_and_extract_bilatr()`'s `extract_from_summ()`/
#' `action_extract()` (which need `mean`/`q5`/`median`/`q95` by name)
#' actually use. `posterior::summarise_draws()`'s own defaults also
#' compute `sd` and `mad`, silently discarded by every caller today --
#' dropping them isn't free at the scale this package targets
#' (~1.9M Tier 3 variables per chain), so this is used everywhere a CSV
#' chunk gets summarised instead of relying on the defaults.
#'
#' @param draws A `posterior::draws_array`/`draws_df` to summarise.
#' @param n_cores Passed through to `summarise_draws()`'s `.cores`:
#'   `posterior::summarise_draws.draws()` splits `draws`' variables into
#'   `n_cores` slices and forks one process per slice
#'   (`parallel::mclapply()`, copy-on-write on Unix), so this parallelises
#'   the CPU-bound Rhat/ESS computation, not the (sequential, disk-bound)
#'   read that produced `draws`.
#' @return A tibble, one row per variable in `draws`.
#' @keywords internal
.summarise_bilatr_draws <- function(draws, n_cores = 1L) {
  posterior::summarise_draws(
    draws, "mean", posterior::quantile2, "median", "rhat", "ess_bulk", "ess_tail",
    .cores = n_cores
  )
}

#' Read and summarise CSV variables in memory-bounded chunks
#'
#' Shared by [diagnose_convergence()]'s CSV-path branch and
#' [extract_theta()]'s: reads and summarises `variables` in groups of
#' `chunk_size`, discarding each chunk's draws before moving to the next
#' -- this is what keeps peak memory bounded regardless of how many
#' variables there are in total. Chunks are always read strictly
#' sequentially (reading is disk-bound and single-threaded regardless --
#' see [.fast_read_post_warmup_draws] -- so there is nothing to gain, and
#' several chunks' worth of memory to lose, by reading more than one at a
#' time); `n_cores` instead parallelises the CPU-bound part, each chunk's
#' own [.summarise_bilatr_draws()] call, via `posterior::summarise_draws()`'s
#' `.cores` argument.
#'
#' Reads each chunk via [.fast_read_post_warmup_draws] (see its docs for
#' what that read costs per call), discarding the chunk's draws once
#' summarised so peak memory reflects one chunk at a time, per
#' [.compute_chunk_size()]'s accounting.
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param variables Character vector of variable names to summarise.
#' @param chunk_size Variables per chunk.
#' @param n_cores Cores for [.summarise_bilatr_draws()]'s `.cores`; `1L`
#'   for a fully sequential run. See [diagnose_convergence()]'s
#'   `parallel`/`n_workers` for how a caller arrives at this number.
#' @param flip_vars Character vector of variable base names to negate
#'   before summarising (used by [extract_theta()]'s CSV path to apply
#'   [bilatr_orient()]'s sign correction; `character(0)`, the default,
#'   for [diagnose_convergence()], since Rhat/ESS are invariant to a
#'   deterministic sign flip and it would be pointless work there).
#'   Passed straight through to [.fast_read_post_warmup_draws()], which
#'   negates matching columns of each chain's temporary matrix before it
#'   is ever assigned into the read's result array -- not a post-hoc
#'   transformation of the summary columns (which would need to swap the
#'   quantile columns, `quantile(-X, p) == -quantile(X, 1 - p)`, not just
#'   negate them), and not a `draws_df` round-trip (which would force
#'   `summarise_draws()` to rebuild a `draws_array` internally at the
#'   cost of a second full copy of the chunk -- see
#'   [.fast_read_post_warmup_draws()]'s own docs). Column-selective by
#'   construction (B5): `chunk_vars` can legitimately mix a flip-needing
#'   variable (`theta`/`theta_raw`) with one that must NOT flip
#'   (`log_lik[d,t]`, also Tier 3, when `compute_log_lik = 1`), and only
#'   the columns `flip_vars` actually lists are ever negated.
#' @return A tibble, the row-bound [.summarise_bilatr_draws()] output
#'   across all chunks.
#' @keywords internal
.chunked_summarise_csv <- function(prepared, variables, chunk_size, n_cores = 1L, flip_vars = character(0)) {
  chunks <- split(variables, ceiling(seq_along(variables) / chunk_size))

  summarise_one_chunk <- function(chunk_vars) {
    draws <- .fast_read_post_warmup_draws(prepared, chunk_vars, flip_vars = flip_vars)
    .summarise_bilatr_draws(draws, n_cores = n_cores)
  }

  purrr::map_dfr(chunks, summarise_one_chunk)
}

#' Chunked-summarise a Tier 3 (or theta-only) variable set, determining
#' sign orientation from the first chunk instead of a separate `alpha[1]`
#' read
#'
#' Shared by [diagnose_and_extract_bilatr()] and [extract_theta()]'s CSV
#' branches for the case where no Tier 1/2 read is already happening to
#' piggyback `alpha[1]` onto (i.e. `tiers` excludes both 1 and 2, or
#' `extract_theta()`, which never reads anything but `theta`). If
#' `flip_vars` is empty (the model has no reflection symmetry) or
#' `variables` is empty, this is exactly [.chunked_summarise_csv()] with
#' `flip = FALSE` -- no `alpha[1]` read at all. Otherwise, `"alpha[1]"`
#' is prepended to the FIRST chunk only (unless it's already in it),
#' read once unflipped; its posterior median decides `flip`; that
#' chunk's `flip_vars`-matching columns (see
#' [.bilatr_match_draws_columns()]) are negated directly on the
#' `draws_array` (mirroring [.read_and_orient_draws_summary()]'s
#' array-native flip -- no `draws_df` round-trip, which would force
#' `summarise_draws()` to rebuild a `draws_array` internally at the cost
#' of a second full copy) and, if `"alpha[1]"` was only added for this
#' orientation check, its row is dropped from the first chunk's summary
#' afterward, the same way [.read_and_orient_draws_summary()] does it;
#' every remaining chunk is then read via [.chunked_summarise_csv()]
#' with `flip` already known. A full sweep still costs exactly one read
#' per chunk, with no extra pass just to check `alpha[1]`'s sign.
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param variables Character vector of Tier 3 (or theta-only) variable
#'   names to summarise.
#' @param chunk_size Variables per chunk.
#' @param n_cores See [.chunked_summarise_csv()].
#' @param flip_vars Output of [.bilatr_flip_variables()] for the model
#'   being read.
#' @return A list with `summ` (the row-bound summary tibble, same shape
#'   [.chunked_summarise_csv()] returns) and `flip` (logical, the
#'   orientation decision -- for a caller that also needs to flip a
#'   Tier 1/2 read using the same decision).
#' @keywords internal
.chunked_summarise_csv_with_orientation <- function(prepared, variables, chunk_size, n_cores, flip_vars) {
  if (length(flip_vars) == 0 || length(variables) == 0) {
    return(list(
      summ = .chunked_summarise_csv(prepared, variables, chunk_size, n_cores),
      flip = FALSE
    ))
  }

  first_chunk_vars <- utils::head(variables, chunk_size)
  rest_vars <- utils::tail(variables, -length(first_chunk_vars))

  alpha1_injected <- !("alpha[1]" %in% first_chunk_vars)
  read_vars <- if (alpha1_injected) c("alpha[1]", first_chunk_vars) else first_chunk_vars
  first_draws <- .fast_read_post_warmup_draws(prepared, read_vars)
  flip <- stats::median(posterior::extract_variable(first_draws, "alpha[1]")) < 0

  if (flip) {
    flip_cols <- .bilatr_match_draws_columns(posterior::variables(first_draws), flip_vars)
    if (length(flip_cols) > 0) {
      first_draws[, , flip_cols] <- -first_draws[, , flip_cols, drop = FALSE]
    }
  }
  first_summ <- .summarise_bilatr_draws(first_draws, n_cores = n_cores)
  if (alpha1_injected) {
    first_summ <- dplyr::filter(first_summ, variable != "alpha[1]")
  }

  rest_summ <- if (length(rest_vars) > 0) {
    .chunked_summarise_csv(prepared, rest_vars, chunk_size, n_cores, flip_vars = if (flip) flip_vars else character(0))
  } else {
    NULL
  }

  list(summ = dplyr::bind_rows(first_summ, rest_summ), flip = flip)
}

#' Resolve a chunk size for a CSV-file-path chunked read, and report it
#'
#' Shared by [diagnose_convergence()]'s and [extract_theta()]'s CSV-path
#' branches. If `max_memory_mb` was left at its default, emits a one-time
#' `message()` explaining the new cost model: reading is single-threaded
#' and each chunk's cost is roughly a full parse of the largest chain
#' file regardless of how many variables are kept (see
#' [.fast_read_post_warmup_draws]), so fewer, larger chunks are cheaper
#' in wall-time, and the CPU-bound summary step's wall-time (not memory)
#' is what `n_workers`/`parallel` trade against. A second `message()`,
#' always emitted, reports the resolved chunk count/size and the
#' estimated peak broken into its file-sized and per-chunk terms
#' separately (see [.compute_chunk_size()]), so a caller can see the
#' actual budget being spent before a long run commits to it.
#'
#' @param n_vars Number of variables the chunked sweep will cover (for
#'   the reporting message only; does not affect the chunk_size
#'   calculation itself).
#' @param prepared Output of [.prepare_fast_csv_read] (for `n_draws`,
#'   `n_chains`, `file_mb`).
#' @param max_memory_mb,chunk_size See [diagnose_convergence()].
#' @param n_cores See [diagnose_convergence()]'s `parallel`/`n_workers`.
#' @param max_memory_mb_missing Whether the caller left `max_memory_mb`
#'   at its default (via `missing()` in the calling function).
#' @return The resolved integer chunk size.
#' @keywords internal
.resolve_chunk_size_and_report <- function(
  n_vars, prepared, max_memory_mb, chunk_size, n_cores,
  max_memory_mb_missing
) {
  n_draws <- prepared$num_post_warmup_draws
  n_chains <- prepared$n_chains
  file_mb <- prepared$file_mb

  if (max_memory_mb_missing) {
    message(
      "Using the default max_memory_mb = ", max_memory_mb, " (",
      round(max_memory_mb / 1024, 1), " GB). Reading is single-threaded: ",
      "each chunk's cost is roughly a full parse of the largest chain ",
      "file regardless of how many variables are kept (see ",
      "?.fast_read_post_warmup_draws), so fewer, larger chunks are ",
      "cheaper in wall-time -- prefer the LARGEST max_memory_mb your ",
      "job's memory allocation can afford. `parallel`/`n_workers` trade ",
      "wall-time in the summary step (Rhat/rank-normalised ESS) against ",
      "memory, NOT independent of chunk_size: forking n_workers costs ",
      "memory roughly proportional to n_workers (see ",
      "?.bilatr_chunk_overhead_multiplier), which comes out of this same ",
      "budget and can force smaller/more chunks -- check the resolved ",
      "chunk count below before committing a large n_workers to a long run."
    )
  }

  chunk_size_used <- chunk_size %||% .compute_chunk_size(
    n_draws = n_draws, n_chains = n_chains, file_mb = file_mb,
    max_memory_mb = max_memory_mb, n_cores = n_cores
  )

  # The resolved chunk_size can exceed n_vars (the whole sweep fits in
  # one chunk with room to spare); report against the variable count
  # actually read in that case, not the theoretical, unused capacity.
  k <- .bilatr_chunk_overhead_multiplier(n_chains, n_cores)
  effective_chunk <- min(chunk_size_used, n_vars)
  raw_mb <- n_draws * n_chains * effective_chunk * 8 / 1e6
  est_mb <- .BILATR_CHUNK_BASELINE_MB + file_mb + raw_mb * k
  message(
    n_vars, " variable(s) in ",
    ceiling(n_vars / chunk_size_used), " chunk(s) of ",
    chunk_size_used, " variable(s) each; estimated peak ~", round(est_mb),
    " MB (", .BILATR_CHUNK_BASELINE_MB, " MB R/package baseline + ",
    round(file_mb), " MB file-sized term + ", round(raw_mb),
    " MB raw chunk x ", signif(k, 3), " multiplier)",
    if (n_cores > 1) paste0(" [", n_cores, " core(s) for summarising]") else "",
    "."
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
#' @return A tibble in the same shape [diagnose_convergence()]'s
#'   in-memory branch produces: [posterior::summarise_draws()] columns
#'   left-joined with [.classify_bilatr_tier]'s `tier`/`dyad_id`/
#'   `time_index`.
#' @keywords internal
.read_diagnostics_summary_from_csv <- function(
  csv_files, tiers, max_memory_mb, chunk_size, parallel, n_workers,
  max_memory_mb_missing
) {
  prepared <- .prepare_fast_csv_read(csv_files)

  var_tiers <- .classify_bilatr_tier(prepared$variables)
  keep_tiers <- var_tiers[var_tiers$tier %in% tiers, ]

  tier12_vars <- keep_tiers$variable[keep_tiers$tier %in% c(1L, 2L)]
  tier3_vars <- keep_tiers$variable[keep_tiers$tier == 3L]

  n_cores <- if (parallel) n_workers else 1L

  chunk_size_used <- if (length(tier3_vars) > 0) {
    .resolve_chunk_size_and_report(
      length(tier3_vars), prepared, max_memory_mb, chunk_size, n_cores,
      max_memory_mb_missing
    )
  } else {
    NULL
  }
  # No sign orientation to worry about here (Rhat/ESS are flip-invariant,
  # unlike diagnose_and_extract_bilatr()'s equivalent), so folding
  # Tier 1/2 and Tier 3 together when Tier 3 fits one chunk anyway is
  # just a single unconditional read.
  tier3_fits_one_chunk <- length(tier3_vars) > 0 && length(tier3_vars) <= chunk_size_used

  tier12_summ <- NULL
  tier3_summ <- NULL

  if (length(tier12_vars) > 0 && tier3_fits_one_chunk) {
    draws <- .fast_read_post_warmup_draws(prepared, c(tier12_vars, tier3_vars))
    tier12_summ <- .summarise_bilatr_draws(draws, n_cores = n_cores)
  } else {
    if (length(tier12_vars) > 0) {
      draws <- .fast_read_post_warmup_draws(prepared, tier12_vars)
      tier12_summ <- .summarise_bilatr_draws(draws, n_cores = n_cores)
    }
    if (length(tier3_vars) > 0) {
      tier3_summ <- .chunked_summarise_csv(prepared, tier3_vars, chunk_size_used, n_cores)
    }
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
#' Every `fit` case (a)-(c) below summarises via the same fixed measure
#' set ([.summarise_bilatr_draws()]: `mean`, `quantile2`, `median`,
#' `rhat`, `ess_bulk`, `ess_tail`), not `posterior::summarise_draws()`'s
#' own defaults (which also compute `sd`/`mad`, unused by every consumer
#' here) -- so the tibble you get back has the same columns regardless
#' of which kind of `fit` produced it.
#'
#' @param fit One of three things. (a) A `posterior::draws_array`/
#'   `draws_df` -- already fully materialized in memory, so `tiers`
#'   controls what gets summarised but not what gets read (there is
#'   nothing left to avoid reading), and `max_memory_mb`/`chunk_size`/
#'   `parallel`/`n_workers` are unused. (b) A `CmdStanMCMC`/`CmdStanFit`-
#'   like fit object (anything with `$metadata()` and
#'   `$draws(variables = ...)` methods) -- variable names are read via
#'   `$metadata()$variables` without touching a single draw, classified
#'   into tiers, and only the tiers actually requested are read via
#'   `$draws(variables = keep_vars)`; a quantity outside `tiers` is
#'   therefore never read into memory at all, same as case (c) below.
#'   `max_memory_mb`/`chunk_size`/`parallel`/`n_workers` are still unused
#'   here -- chunking only applies to the raw-CSV path, since `$draws()`
#'   already holds whatever it returns as one in-memory array. (c) A
#'   character vector of raw CmdStan CSV file paths (one per chain, e.g.
#'   from a completed SLURM run never loaded into this R session) -- in
#'   this case Tier 1/2 variables are read in one small call, and Tier 3
#'   variables (typically the overwhelming majority for production-sized
#'   panels) are read and summarised directly from the raw CSVs in
#'   memory-bounded chunks (see [.fast_read_post_warmup_draws]),
#'   discarding each chunk's draws before moving to the next, so the
#'   full draws array is never materialized at once. Total file touches
#'   per chain file: one metadata scan plus one validation probe (both in
#'   [.prepare_fast_csv_read]), one Tier 1/2 read, and one per Tier 3
#'   chunk -- or, when Tier 3 fits in a single chunk, Tier 1/2 is folded
#'   into that one read instead of being separate. This is the path that
#'   matters for production Stan output too large to read normally
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
#'   includes `3`. Target ceiling, in MB, for Tier 3's peak memory per
#'   chunk; drives the automatically-derived `chunk_size` (see
#'   `chunk_size` below) so you don't have to guess a variable count
#'   yourself. Defaults to `8192` (8 GB) -- a guess, not a calibration
#'   against your hardware, and this function says so via `message()`
#'   the first time you rely on that default rather than setting it
#'   explicitly. This is a sanity-check number, not a guarantee: actual
#'   peak memory depends on `data.table::fread()`/`summarise_draws()`
#'   internals this function doesn't control (see
#'   [.compute_chunk_size()] for the model and its derivation). Reading
#'   is single-threaded, and each chunk's cost is roughly a full parse
#'   of the largest chain file regardless of how many variables are
#'   kept -- `max_memory_mb` must be large enough to cover that file
#'   size, plus a fixed R/package-loading floor, at minimum (this
#'   `stop()`s with a clear message naming the shortfall if it isn't).
#'   Beyond that, fewer/larger chunks are cheaper in read wall-time, but
#'   -- see `parallel` below -- a larger `n_workers` shrinks the
#'   largest chunk size the SAME `max_memory_mb` budget can afford, so
#'   the two are not independent: prefer the LARGEST `max_memory_mb`
#'   your job's allocation can afford, and check the resolved chunk
#'   count in the pre-flight `message()` before committing to a large
#'   `n_workers` on a tight budget.
#' @param chunk_size Only used when `fit` is CSV file paths and `tiers`
#'   includes `3`. Explicit override: number of Tier 3 variables read per
#'   chunk. `NULL` (the default) derives this from `max_memory_mb`
#'   instead; set this directly only if you want to bypass that
#'   calculation (e.g. you've measured actual memory use and want to
#'   tune it by hand).
#' @param parallel Only used when `fit` is CSV file paths and `tiers`
#'   includes `3`. Chunks are always read strictly sequentially regardless
#'   of this argument -- reading is disk-bound and single-threaded no
#'   matter how many cores are available (see
#'   [.fast_read_post_warmup_draws]), so there is nothing to gain, and
#'   `n_workers` chunks' worth of memory to lose, by reading more than one
#'   chunk at a time. `parallel` instead controls whether the CPU-bound
#'   part -- each chunk's Rhat/rank-normalised-ESS computation -- uses
#'   `n_workers` cores (via `posterior::summarise_draws()`'s `.cores`
#'   argument, which forks `n_workers` worker processes) or just one.
#'   Measured (`dev/bench_memory.R`, see [.BILATR_CHUNK_CORES_PER_CORE_FACTOR]):
#'   this genuinely costs memory roughly proportional to `n_workers`, not
#'   a fixed amount regardless of it, and that cost comes out of the SAME
#'   `max_memory_mb` budget the read uses -- so a larger `n_workers`
#'   indirectly means MORE, not fewer, passes over each chain file at a
#'   fixed `max_memory_mb`, trading read wall-time for summary wall-time
#'   rather than being free on top of it. Worth checking the resolved
#'   chunk count (in the pre-flight `message()`) at your intended
#'   `n_workers` before committing a long run to it -- a very large
#'   `n_workers` on a large `n_vars` sweep can end up costing MORE total
#'   wall-time than a smaller one, if it forces enough extra chunks.
#' @param n_workers Only used when `parallel = TRUE`. Defaults to
#'   `parallelly::availableCores()`, which -- unlike
#'   `parallel::detectCores()` -- respects a SLURM allocation's
#'   `SLURM_CPUS_PER_TASK` (among other cluster/container schedulers)
#'   rather than reporting the whole node's core count. See `parallel`
#'   above: this is a genuine memory/wall-time trade-off now, not a
#'   free choice, so consider an explicit, moderate value (rather than
#'   the default, which can be large on a big allocation) for a
#'   memory-constrained job.
#' @param scratch_dir Deprecated and ignored since 0.4.1: reads are made
#'   directly against the raw CSV files (see [.prepare_fast_csv_read]),
#'   so no scratch copy is ever made any more. Passing a non-`NULL`
#'   value emits a warning.
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
  n_workers = parallelly::availableCores(),
  scratch_dir = NULL
) {
  max_memory_mb_missing <- missing(max_memory_mb)
  if (!is.null(scratch_dir)) {
    warning(
      "`scratch_dir` is deprecated and ignored since 0.4.1: no scratch ",
      "copy is made any more.",
      call. = FALSE
    )
  }

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
      fit, tiers, max_memory_mb, chunk_size, parallel, n_workers, max_memory_mb_missing
    )
  } else if (posterior::is_draws(fit)) {
    var_tiers <- .classify_bilatr_tier(posterior::variables(fit))
    keep_vars <- var_tiers$variable[var_tiers$tier %in% tiers]
    summ <- .summarise_bilatr_draws(posterior::subset_draws(fit, variable = keep_vars))
    summ <- dplyr::left_join(summ, var_tiers, by = "variable")
  } else {
    # A CmdStanMCMC-like fit object still backed by its own CSV files:
    # $metadata()$variables lists every posterior-style variable name
    # without touching a single draw, so tiers not requested can be
    # excluded from $draws()'s own read via `variables =`, instead of
    # reading everything and subsetting afterward.
    var_tiers <- .classify_bilatr_tier(fit$metadata()$variables)
    keep_vars <- var_tiers$variable[var_tiers$tier %in% tiers]
    summ <- .summarise_bilatr_draws(fit$draws(variables = keep_vars))
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
