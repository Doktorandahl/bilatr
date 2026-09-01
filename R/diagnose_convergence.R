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
  "beta_logn", "lp__"
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
#' rather than a fixed per-parameter name list, so it classifies
#' consistently across Stan model variants that change a parameter's
#' shape (e.g. `phi` is per-dyad in the `stable` model but per-dyad-period
#' in `phi_logn`; see `R/model_registry.R`) without special-casing either.
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
#' assigned to a requested tier are dropped before
#' [posterior::summarise_draws()] ever runs on them (via
#' [posterior::subset_draws()]), rather than just being hidden afterward.
#' This matters in practice because Tier 3 (`theta`/`theta_raw`) typically
#' has, by far, the most monitored quantities of the three tiers (one per
#' dyad-period); requesting only `tiers = 1` or `tiers = 1:2` skips
#' computing Rhat/ESS for all of them.
#'
#' @param fit A `CmdStanMCMC`/`CmdStanFit`-like fit object (anything with
#'   a `$draws()` method), or a `posterior::draws_array`/`draws_df`.
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
#' }
#' @export
diagnose_convergence <- function(
  fit,
  n_dt = NULL,
  rhat_threshold = 1.01,
  ess_threshold = 400,
  tiers = 1:3
) {
  tiers <- .validate_tiers(tiers)
  if (any(c(2L, 3L) %in% tiers) && is.null(n_dt)) {
    stop(
      "`n_dt` is required when `tiers` includes 2 and/or 3 (Tier 1 alone ",
      "needs no per-dyad join).",
      call. = FALSE
    )
  }

  draws <- if (posterior::is_draws(fit)) fit else fit$draws()
  var_tiers <- .classify_bilatr_tier(posterior::variables(draws))
  keep_vars <- var_tiers$variable[var_tiers$tier %in% tiers]
  summ <- posterior::summarise_draws(posterior::subset_draws(draws, variable = keep_vars))
  summ <- dplyr::left_join(summ, var_tiers, by = "variable")

  n_dt_tbl <- if (any(c(2L, 3L) %in% tiers)) .normalize_n_dt(n_dt) else NULL

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
