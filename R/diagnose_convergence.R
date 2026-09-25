#' Names of Stan parameters treated as global/shared (Tier 1)
#'
#' Matched against a monitored quantity's *base name* (its variable name
#' with any `[...]` index stripped), so this correctly matches vector
#' parameters like `alpha[1]`, `alpha[2]`, ... via their shared base name
#' `alpha`. `alpha_raw` is NOT listed here (despite also being indexed by
#' action type, not dyad): it is sign-ambiguous under the alpha/theta
#' reflection symmetry (see [.bilatr_sign_ambiguous_raw_names()]) and
#' caught by that check, earlier in [.classify_bilatr_tier()]'s
#' `case_when()`, before this list is ever consulted.
#'
#' `mu_intercept_raw`/`mu_theta0` (0.5.0 tidy-up, see NEWS.md) were
#' dropped from this list: neither exists in any currently registered
#' model (`mu_intercept` is a plain `sum_to_zero_vector`, not a
#' non-centered `_raw` parameter, and `mu_theta0` was removed when
#' `theta0` became `sigma_theta0 * z_theta0` with no separate location
#' parameter) -- listing them here was dead code that could only ever
#' silently match nothing.
#'
#' `theta_filtered`/`theta_filtered_sd` (0.5.0) need no entry here or
#' anywhere else in [.classify_bilatr_tier()]: Tier 3 is matched
#' structurally, by bracket-index count (`n_index >= 2`, see that
#' function), not by name, so a `theta_filtered[d, t]`/
#' `theta_filtered_sd[d, t]` variable lands there automatically, the same
#' way `log_lik[d, t]` already does. One caveat: when
#' `assemble_stan_data()`'s `filter_dyads` narrows the dyad set, the
#' `dyad_id`/`index_1` these two variables get in Tier 3 output is the
#' *position within the filtered subset* (`1..n_filter_dyads`), not the
#' true `D`-space `dyad_id` -- there is no translation layer back to the
#' true `dyad_id` here (out of scope for 0.5.0); joining these two
#' variables' Tier 3 rows to `dyad_ids`/other per-dyad metadata by
#' `dyad_id` is only valid when `filter_dyads` was `NULL` (i.e. every
#' dyad was filtered, in `Y`'s original order).
#' `gamma` (0.7.0, `stable_gamma` only: `A x n_countries`) is listed here
#' too, for the same by-name reason as `alpha`: it has two `[...]`
#' indices (action, country), which the structural
#' `n_index >= 2 -> Tier 3` rule below would otherwise misclassify as a
#' dyad-indexed quantity, reading the country index as a `dyad_id` and
#' sweeping it into the expensive Tier 3 per-dyad-period pass (see
#' inst/stan/bilatr_alphanorm_gamma.stan's header and
#' dev/claude_code_prompt_0.7.0_country_offsets.md, "Part 4"). `gamma` is
#' fully identified (unlike `gamma_z`, which is prior-only pinned in its
#' projected-out directions -- see [.bilatr_sign_tied_names()] in
#' `R/sign_ambiguity.R` for where THAT is excluded instead, via the
#' `tier = NA` path below, for a different reason than sign ambiguity).
#' @keywords internal
.bilatr_tier1_names <- c(
  "alpha", "mu_intercept", "sigma_theta0",
  "mu_log_phi", "sigma_log_phi", "mu_log_noise", "sigma_log_noise",
  "gamma", "sigma_gamma",
  "lp__"
)

#' Names of raw, sampled parameters that are sign-ambiguous under the
#' alpha/theta reflection symmetry, across every registered model
#'
#' `stable`/`ou`'s orientation fold (see each `.stan` file's header,
#' "IDENTIFICATION: ORIENTATION FOLD") corrects the reported
#' `alpha`/`theta`/etc., but the RAW parameters it's built from
#' (`alpha_raw` itself, plus whichever of `z_theta0`/`mu_dyad_raw` and
#' `theta_raw` a given model has) remain genuinely sign-ambiguous: if two
#' chains land in opposite raw-space basins, those variables' own
#' cross-chain Rhat is meaningless even though everything reported is
#' fine. [.classify_bilatr_tier()] uses this to keep them out of the
#' tiered diagnostics tables entirely (`tier = NA`), rather than let a
#' meaningless Rhat surface as a false Tier 1 alarm.
#'
#' Derived from [.bilatr_sign_tied_names()] (`R/sign_ambiguity.R`), the
#' single source of truth for which raw names are sign-tied per model --
#' unioning its `raw` component across every currently-registered model.
#' Computed inside the function body, not as a top-level constant, so it
#' doesn't depend on `R/model_registry.R` having been sourced first.
#'
#' This function's name and output must not change across versions: the
#' runscripts call `bilatr:::.bilatr_sign_ambiguous_raw_names()` directly.
#'
#' @return Character vector of variable base names.
#' @keywords internal
.bilatr_sign_ambiguous_raw_names <- function() {
  unique(unlist(lapply(
    names(.bilatr_stan_models),
    function(m) .bilatr_sign_tied_names(m)$raw
  )))
}

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
#' Sign-ambiguous raw parameters (see [.bilatr_sign_ambiguous_raw_names()]:
#' `alpha_raw`, `z_theta0`/`mu_dyad_raw`, `theta_raw`) are excluded first,
#' by base name -- `tier = NA`, dropped from every downstream tier table
#' by the `tier %in% c(1L, 2L, 3L)` filters already used throughout this
#' file, rather than surfacing their meaningless cross-chain Rhat as a
#' false convergence alarm (see the alpha/theta reflection symmetry
#' discussion in each `.stan` file's header). Tier 1 (global/shared) is
#' matched next, by base name, against [.bilatr_tier1_names]. Everything
#' else is classified structurally by how many `[...]` indices it
#' carries: a single index (`name[d]`) is assumed to be a per-dyad
#' hierarchical parameter (Tier 2, joined on `d`); two indices (`name[d,
#' t]`) is assumed to be a per-dyad-period latent state (Tier 3, joined
#' on `d`). This is deliberately structural rather than a fixed
#' per-parameter name list, so a future model variant that changes a
#' parameter's shape (e.g. makes `phi` per-dyad-period, `phi[d, t]`,
#' instead of the per-dyad `phi[d]` of the `stable` model) is still
#' classified consistently without special-casing.
#' Anything with no brackets that isn't in the Tier 1 name list (should
#' not occur for the package's own models, but could for a hand-edited
#' Stan file) is folded into Tier 1 rather than dropped, since its
#' sparsity profile is unknown and it should never be silently hidden.
#'
#' @param variable Character vector of `summarise_draws()` variable names.
#' @return A tibble with columns `variable`, `tier` (`1L`, `2L`, `3L`, or
#'   `NA` for excluded sign-ambiguous raw parameters), `dyad_id` (the
#'   first index, `NA` outside Tier 2/3), and `time_index` (the second
#'   index, `NA` outside Tier 3).
#' @keywords internal
.classify_bilatr_tier <- function(variable) {
  parsed <- .parse_variable_indices(variable)
  sign_ambiguous <- .bilatr_sign_ambiguous_raw_names()

  dplyr::mutate(
    parsed,
    tier = dplyr::case_when(
      base_name %in% sign_ambiguous ~ NA_integer_,
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
#'
#' @param exclude_pattern Optional regex; `variable`s matching it are
#'   excluded before `flagged` is computed. Used (0.7.1) to pull `gamma`
#'   out of this table and into its own (see [.compute_gamma_tier()] and
#'   [.assemble_bilatr_diagnostics()]) -- `gamma` is correctly Tier 1 for
#'   *classification* (see `.bilatr_tier1_names`), but at production
#'   scale it is `A x n_countries` (~thousands of parameters), and Tier 1
#'   is also a *reporting* category designed around a handful of global
#'   scalars: `sigma_gamma[...]` (a per-category scalar, not
#'   country-indexed) is NOT excluded by the `stable_gamma` caller's
#'   pattern and stays in this table.
#' @keywords internal
.compute_tier1 <- function(summ, rhat_threshold, ess_threshold, exclude_pattern = NULL) {
  t1 <- dplyr::filter(summ, tier == 1L)
  if (!is.null(exclude_pattern)) {
    t1 <- dplyr::filter(t1, !grepl(exclude_pattern, variable))
  }
  t1 %>%
    dplyr::mutate(
      flagged = .flag_diagnostic(rhat, ess_bulk, ess_tail, rhat_threshold, ess_threshold)
    ) %>%
    dplyr::select(variable, rhat, ess_bulk, ess_tail, flagged) %>%
    dplyr::arrange(dplyr::desc(flagged))
}

#' Compute the `gamma` (country-level offset) diagnostics tibble
#'
#' Split out of Tier 1 (0.7.1) so that a routine number of flagged
#' `gamma` elements -- dozens to hundreds at production scale, just from
#' `ess_threshold` against ~4,000 draws -- doesn't flood
#' [print.bilatr_diagnostics()]'s Tier 1 block or dominate
#' `n_tier1_flagged`. `country_code`/`action_index` are parsed directly
#' from the `gamma[k,c]` variable name (self-contained; does not depend
#' on [.classify_bilatr_tier()] exposing raw bracket indices), the same
#' pattern [extract_gamma()] uses.
#'
#' @param summ A tibble as produced by [posterior::summarise_draws()],
#'   left-joined with [.classify_bilatr_tier]'s `tier` column.
#' @param country_codes Optional character vector (in `country_index`
#'   order, e.g. `stan_data`'s `"country_codes"` attribute) for labelling
#'   `country_index`. `NULL` (e.g. plain [diagnose_convergence()], which
#'   has no `stan_data`) leaves the table with `country_index` only.
#' @return A tibble: `variable`, `action_index`, `country_index` (and
#'   `country_code` if `country_codes` is supplied), `rhat`, `ess_bulk`,
#'   `ess_tail`, `flagged`.
#' @keywords internal
.compute_gamma_tier <- function(summ, rhat_threshold, ess_threshold, country_codes = NULL) {
  g <- summ %>%
    dplyr::filter(tier == 1L, startsWith(variable, "gamma[")) %>%
    dplyr::mutate(
      idx = stringr::str_match(variable, "\\[(\\d+),(\\d+)\\]"),
      action_index = as.integer(idx[, 2]),
      country_index = as.integer(idx[, 3]),
      flagged = .flag_diagnostic(rhat, ess_bulk, ess_tail, rhat_threshold, ess_threshold)
    ) %>%
    dplyr::select(-idx)

  if (!is.null(country_codes)) {
    g <- dplyr::mutate(g, country_code = country_codes[country_index])
  }

  dplyr::select(
    g,
    dplyr::any_of(c("variable", "action_index", "country_index", "country_code")),
    rhat, ess_bulk, ess_tail, flagged
  )
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
#' @param has_gamma (0.7.1) Whether the model being diagnosed has a
#'   `gamma` (country-level offset) parameter -- from
#'   [.bilatr_model_has_gamma()], never inferred by checking `summ` for a
#'   `gamma` name. When `TRUE` and Tier 1 was requested, `gamma` is
#'   pulled out of `tier1` into its own `gamma` element (see
#'   [.compute_gamma_tier()]); when `FALSE`, `gamma` is `NULL` and
#'   `n_gamma_*` are `NA_integer_`, so a model without `gamma` gets no
#'   trace of the section at all.
#' @param country_codes (0.7.1) Optional character vector for labelling
#'   `gamma`'s `country_index`; see [.compute_gamma_tier()].
#' @return A list of class `bilatr_diagnostics`; see
#'   [diagnose_convergence()]'s `@return` for the element-by-element
#'   description.
#' @keywords internal
.assemble_bilatr_diagnostics <- function(summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold,
                                          has_gamma = FALSE, country_codes = NULL) {
  tier1 <- if (1L %in% tiers) {
    .compute_tier1(summ, rhat_threshold, ess_threshold, exclude_pattern = if (has_gamma) "^gamma\\[" else NULL)
  } else {
    NULL
  }
  gamma_diag <- if (1L %in% tiers && has_gamma) {
    .compute_gamma_tier(summ, rhat_threshold, ess_threshold, country_codes = country_codes)
  } else {
    NULL
  }
  tier2_result <- if (2L %in% tiers) .compute_tier2(summ, n_dt_tbl) else NULL
  tier2 <- tier2_result$tier2
  tier3 <- if (3L %in% tiers) .compute_tier3(summ, n_dt_tbl, rhat_threshold, ess_threshold) else NULL

  summary_info <- list(
    tiers_computed = tiers,
    n_tier1_flagged = if (!is.null(tier1)) sum(tier1$flagged, na.rm = TRUE) else NA_integer_,
    n_tier1_total = if (!is.null(tier1)) nrow(tier1) else NA_integer_,
    n_gamma_flagged = if (!is.null(gamma_diag)) sum(gamma_diag$flagged, na.rm = TRUE) else NA_integer_,
    n_gamma_total = if (!is.null(gamma_diag)) nrow(gamma_diag) else NA_integer_,
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
    list(tier1 = tier1, gamma = gamma_diag, tier2 = tier2, tier3 = tier3, summary = summary_info),
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

#' Memory-model constant: a documented FLOOR for R and its loaded packages,
#' used only where the runtime measurement in
#' [.bilatr_estimate_chunk_baseline_mb()] is unavailable
#'
#' Unlike every other `.BILATR_CHUNK_*` constant, this is not a
#' multiplier of `raw_mb` -- it is a fixed number of MB, added once
#' alongside `file_mb` in [.compute_chunk_size()]/
#' [.estimate_diagnostics_memory_mb()], because it does not shrink with
#' chunk size any more than `file_mb` does: it is the cost of `library(
#' bilatr)` and its dependencies (`cmdstanr`, `data.table`, `posterior`,
#' `dplyr`, etc.) simply being loaded, before a single byte of any CSV
#' is read. 200 was measured via `dev/bench_memory.R` on exactly ONE
#' machine (2026-09; `devtools::load_all()` alone, R 4.5.2, macOS: ~160
#' MB peak RSS, rounded up for margin) -- an artefact of that machine's R
#' build and loaded namespaces, not a universal constant, and the reason
#' [.bilatr_estimate_chunk_baseline_mb()] measures the CALLER's own process
#' instead of assuming this number applies. This value is retained only
#' as the floor that measurement falls back to when it fails (a platform
#' without `/proc`, or `gc()` itself erroring), so the model never
#' silently under-counts to zero.
#' @keywords internal
.BILATR_CHUNK_BASELINE_MB <- 200

#' Measure the caller's own process memory footprint, for use as the
#' memory model's baseline term
#'
#' Replaces a hardcoded assumption ([.BILATR_CHUNK_BASELINE_MB]) with an
#' actual reading of the current process: on Linux, `VmRSS` from
#' `/proc/self/status` (the same per-process resident-memory figure the
#' rest of this package's benchmarking uses); everywhere else, `sum(gc()[
#' , "(Mb)"])` (Ncells + Vcells currently used) as a portable fallback --
#' cruder (R's own view of its heap, not the OS's view of the whole
#' process, so it can miss non-R allocations made by loaded C libraries)
#' but available anywhere. Whichever succeeds is floored at
#' [.BILATR_CHUNK_BASELINE_MB] so a session that happens to measure
#' lighter than that (e.g. right after a `gc()`) doesn't understate the
#' margin the constant was chosen to provide, and a measurement that
#' fails outright (unreadable `/proc`, `gc()` erroring) falls back to it
#' entirely rather than propagating `NA` into the memory model.
#'
#' @return Estimated baseline memory in MB.
#' @keywords internal
.bilatr_estimate_chunk_baseline_mb <- function() {
  measured <- tryCatch(
    {
      if (identical(Sys.info()[["sysname"]], "Linux") && file.exists("/proc/self/status")) {
        status <- readLines("/proc/self/status")
        vmrss_line <- grep("^VmRSS:", status, value = TRUE)
        if (length(vmrss_line) == 1) {
          kb <- as.numeric(sub("[^0-9]+([0-9]+).*", "\\1", vmrss_line))
          kb / 1024
        } else {
          NA_real_
        }
      } else {
        NA_real_
      }
    },
    error = function(e) NA_real_
  )

  if (is.na(measured)) {
    measured <- tryCatch(sum(gc()[, "(Mb)"]), error = function(e) NA_real_)
  }

  if (is.na(measured)) {
    .BILATR_CHUNK_BASELINE_MB
  } else {
    max(measured, .BILATR_CHUNK_BASELINE_MB)
  }
}

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
#'
#' **PENDING RE-FIT (0.4.2), NOT YET DONE**: the measurement below sums
#' RSS across the fork tree, which double-counts copy-on-write shared
#' pages a SLURM cgroup only charges once -- `dev/bench_memory.R` was
#' corrected (0.4.2) to poll summed Pss and the cgroup's own peak-usage
#' counter instead (see that file's header), and a small Linux
#' measurement using the corrected poller found roughly 3.8x `raw_mb` at
#' the first extra worker and ~2.6x `raw_mb` per worker after that --
#' both LOWER than the 4.5/2.5 below, consistent with this constant
#' being ~1.5-2x conservative in the `n_cores > 1` range. This value has
#' NOT yet been updated to reflect that: doing so needs a full
#' Linux-sourced run of the corrected benchmark (this constant's own
#' derivation below is macOS sum-RSS, and the jobs it sizes run on
#' Linux/SLURM), which had not happened as of 0.4.2's release. Treat the
#' current 4.5 as a documented-conservative placeholder, not a
#' recalibrated figure, until that run happens.
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
#'
#' **PENDING RE-FIT (0.4.2), NOT YET DONE**: see
#' [.BILATR_CHUNK_CORES_STEP_FACTOR]'s docs -- the corrected (Pss/cgroup)
#' benchmark suggests ~2.6 here, consistent with this constant's 2.5
#' within measurement noise (unlike the step term, which moved more).
#' Still pending a full Linux-sourced re-fit before treating either as
#' recalibrated rather than the pre-correction placeholder.
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
#' `.BILATR_CHUNK_PER_CHAIN_FACTOR / n_chains`); and, only when `n_cores >
#' 1`, `posterior::summarise_draws()`'s forked-worker cost -- a one-time
#' step (`.BILATR_CHUNK_CORES_STEP_FACTOR`) plus a per-additional-worker
#' term (`.BILATR_CHUNK_CORES_PER_CORE_FACTOR * (n_cores - 1)`), NOT a
#' flat penalty regardless of the exact core count (see that constant's
#' docs -- this replaced an assumption that was never measured and was
#' wrong).
#'
#' 0.10.0 dropped the `.BILATR_CHUNK_FLIP_FACTOR` term that used to sit
#' here: the post-hoc sign-flip machinery it budgeted for
#' (`bilatr_orient()` and its array-native equivalents) was retired along
#' with the soft-anchor stack that needed it (see NEWS.md), so `k` is
#' exactly 1 lower than before for every `(n_chains, n_cores)` pair.
#'
#' @param n_chains From [.prepare_fast_csv_read].
#' @param n_cores See [diagnose_convergence()]'s `parallel`/`n_workers`.
#' @return The multiplier `k`, in units of `raw_mb`.
#' @keywords internal
.bilatr_chunk_overhead_multiplier <- function(n_chains, n_cores) {
  .BILATR_CHUNK_ARRAY_FACTOR +
    .BILATR_CHUNK_PER_CHAIN_FACTOR / n_chains +
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
  baseline_mb <- .bilatr_estimate_chunk_baseline_mb()
  budget_mb <- max_memory_mb - file_mb - baseline_mb

  if (budget_mb < per_var_mb) {
    stop(
      "max_memory_mb (", round(max_memory_mb), " MB) leaves no room for ",
      "even one Tier 3 variable per chunk: the largest chain file (~",
      round(file_mb), " MB) is mapped in full during every read, and ~",
      round(baseline_mb), " MB is this process's own measured baseline ",
      "(R and its loaded packages) -- both before any variables are ",
      "selected, regardless of chunk_size. Allocate at least ",
      ceiling(file_mb + baseline_mb + per_var_mb),
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
  .bilatr_estimate_chunk_baseline_mb() + file_mb + raw_mb * k
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
#' @return A tibble, the row-bound [.summarise_bilatr_draws()] output
#'   across all chunks.
#' @keywords internal
.chunked_summarise_csv <- function(prepared, variables, chunk_size, n_cores = 1L) {
  chunks <- split(variables, ceiling(seq_along(variables) / chunk_size))

  summarise_one_chunk <- function(chunk_vars) {
    draws <- .fast_read_post_warmup_draws(prepared, chunk_vars)
    .summarise_bilatr_draws(draws, n_cores = n_cores)
  }

  purrr::map_dfr(chunks, summarise_one_chunk)
}

#' Report the wall-time/core-hour trade-off across candidate `n_workers`
#'
#' Chunk size shrinks as `n_cores` rises (forking `n_workers` costs memory
#' roughly proportional to `n_workers`, out of the same `max_memory_mb`
#' budget -- see [.bilatr_chunk_overhead_multiplier()]), and each chunk's
#' read is a full parse of the largest chain file (see
#' [.fast_read_post_warmup_draws]) -- so raising `n_workers` doesn't just
#' speed up the CPU-bound summary step, it also forces MORE chunks, i.e.
#' more full-file re-parses. Past some point this makes wall time worse,
#' not better, and core-seconds (`n_workers * wall time`) rise
#' monotonically with `n_workers` even where wall time itself is still
#' falling -- nothing before this helper made either consequence visible
#' to a caller sizing a job.
#'
#' For each candidate in `worker_levels`, resolves `chunk_size` ([.compute_chunk_size()])
#' and the resulting `n_chunks = ceiling(n_vars / chunk_size)`, estimated
#' peak memory ([.estimate_diagnostics_memory_mb()]), estimated wall time
#' `n_chunks * read_seconds + n_vars * summarise_seconds_per_variable /
#' n_workers`, and estimated core-seconds (`n_workers * wall time`). A
#' `n_workers` level whose budget can't fit even one variable per chunk
#' (see [.compute_chunk_size()]'s own `stop()`) is reported as infeasible
#' (`NA` `chunk_size`/`n_chunks`/wall/core-seconds) rather than erroring
#' the whole table.
#'
#' @param n_vars Number of Tier 3 variables the sweep will cover.
#' @param n_draws,n_chains,file_mb From [.prepare_fast_csv_read].
#' @param max_memory_mb See [diagnose_convergence()].
#' @param worker_levels Candidate `n_workers` values to evaluate. Defaults
#'   to `c(1, 2, 4, 8, 16, 24)`; the
#'   [.BILATR_CHUNK_CORES_PER_CORE_FACTOR]/[.BILATR_CHUNK_CORES_STEP_FACTOR]
#'   memory-model constants this helper's peak-MB estimate depends on are
#'   only measured up to `n_cores = 8`, so treat levels above that as an
#'   extrapolation (see those constants' docs).
#' @param read_seconds Wall-time (seconds) for ONE chunk's read -- a full
#'   parse of the largest chain file, roughly independent of chunk size
#'   (see [.fast_read_post_warmup_draws]). No default: this is the one
#'   input this helper cannot estimate on its own, and the whole
#'   wall-time column is only as good as it is. Take it from your own
#'   job's logs (a `diagnose_convergence()`/`extract_theta()` run against
#'   the same files reports its own read time), not a guess.
#' @param summarise_seconds_per_variable Wall-time (seconds) for the
#'   CPU-bound summary step ([posterior::summarise_draws()]'s
#'   rank-normalised Rhat/bulk-tail-ESS), per Tier 3 variable, on ONE
#'   core. Defaults to `0.0026` (~2.6 ms/variable), measured against this
#'   package's own Tier 3 sweep (`dev/bench_memory.R`, 2026-09); scales
#'   with however expensive computing those statistics is for your
#'   posterior (number of draws, mostly), so treat the default as a
#'   starting point, not a calibration for your own data.
#' @return A tibble, one row per `worker_levels` entry: `n_workers`,
#'   `chunk_size`, `n_chunks`, `peak_mb`, `wall_seconds`, `core_seconds`
#'   (the last four `NA` where infeasible).
#' @keywords internal
.bilatr_worker_tradeoff <- function(
  n_vars, n_draws, n_chains, file_mb, max_memory_mb,
  worker_levels = c(1, 2, 4, 8, 16, 24),
  read_seconds,
  summarise_seconds_per_variable = 0.0026
) {
  rows <- lapply(worker_levels, function(n_workers) {
    chunk_size <- tryCatch(
      .compute_chunk_size(
        n_draws = n_draws, n_chains = n_chains, file_mb = file_mb,
        max_memory_mb = max_memory_mb, n_cores = n_workers
      ),
      error = function(e) NA_integer_
    )
    if (is.na(chunk_size)) {
      return(tibble::tibble(
        n_workers = n_workers, chunk_size = NA_integer_, n_chunks = NA_integer_,
        peak_mb = NA_real_, wall_seconds = NA_real_, core_seconds = NA_real_
      ))
    }

    n_chunks <- ceiling(n_vars / chunk_size)
    effective_chunk <- min(chunk_size, n_vars)
    peak_mb <- .estimate_diagnostics_memory_mb(
      n_draws = n_draws, n_chains = n_chains, chunk_size = effective_chunk,
      n_cores = n_workers, file_mb = file_mb
    )
    wall_seconds <- n_chunks * read_seconds + n_vars * summarise_seconds_per_variable / n_workers

    tibble::tibble(
      n_workers = n_workers, chunk_size = chunk_size, n_chunks = n_chunks,
      peak_mb = peak_mb, wall_seconds = wall_seconds, core_seconds = n_workers * wall_seconds
    )
  })

  dplyr::bind_rows(rows)
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
#' actual budget being spent before a long run commits to it. A third
#' `message()`, only when `read_seconds` is supplied, names a
#' `n_workers` level (via [.bilatr_worker_tradeoff()]) that would give
#' both lower estimated wall time AND lower estimated core-seconds than
#' the chosen `n_cores`, if one exists in a small grid around it --
#' silently skipped without `read_seconds`, since wall-time can't be
#' estimated at all without it (see that helper's docs).
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
#' @param read_seconds Optional; see [.bilatr_worker_tradeoff()]. `NULL`
#'   (the default) skips the trade-off message entirely.
#' @return The resolved integer chunk size.
#' @keywords internal
.resolve_chunk_size_and_report <- function(
  n_vars, prepared, max_memory_mb, chunk_size, n_cores,
  max_memory_mb_missing, read_seconds = NULL
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
  baseline_mb <- .bilatr_estimate_chunk_baseline_mb()
  effective_chunk <- min(chunk_size_used, n_vars)
  raw_mb <- n_draws * n_chains * effective_chunk * 8 / 1e6
  est_mb <- baseline_mb + file_mb + raw_mb * k
  message(
    n_vars, " variable(s) in ",
    ceiling(n_vars / chunk_size_used), " chunk(s) of ",
    chunk_size_used, " variable(s) each; estimated peak ~", round(est_mb),
    " MB (", round(baseline_mb), " MB R/package baseline (measured) + ",
    round(file_mb), " MB file-sized term + ", round(raw_mb),
    " MB raw chunk x ", signif(k, 3), " multiplier)",
    if (n_cores > 1) paste0(" [", n_cores, " core(s) for summarising]") else "",
    "."
  )

  if (!is.null(read_seconds)) {
    # A small neighborhood around the chosen n_cores, not just {1,
    # n_cores}: dominance (lower wall time AND lower core-seconds) is
    # common between ADJACENT levels on the high side of the wall-time-
    # minimizing point (e.g. n_workers = 16 dominating 24), but rare
    # between 1 and a large n_cores directly, since low n_workers
    # typically trades better core-seconds for worse wall time rather
    # than dominating outright.
    grid <- sort(unique(pmax(1, c(1, n_cores, round(n_cores / 2), n_cores * 2))))
    tradeoff <- .bilatr_worker_tradeoff(
      n_vars = n_vars, n_draws = n_draws, n_chains = n_chains, file_mb = file_mb,
      max_memory_mb = max_memory_mb, worker_levels = grid, read_seconds = read_seconds
    )
    current <- tradeoff[tradeoff$n_workers == n_cores, ]
    if (nrow(current) == 1 && !is.na(current$wall_seconds)) {
      better <- tradeoff[
        !is.na(tradeoff$wall_seconds) & tradeoff$n_workers != n_cores &
          tradeoff$wall_seconds < current$wall_seconds &
          tradeoff$core_seconds < current$core_seconds,
      ]
      if (nrow(better) > 0) {
        best <- better[which.min(better$wall_seconds), ]
        message(
          "n_workers = ", best$n_workers, " is estimated to give BOTH lower ",
          "wall time (", round(best$wall_seconds), "s vs ", round(current$wall_seconds),
          "s) and lower core-seconds (", round(best$core_seconds), " vs ",
          round(current$core_seconds), ") than the current n_workers = ", n_cores,
          " -- see .bilatr_worker_tradeoff() for the full grid."
        )
      }
    }
  }

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
  max_memory_mb_missing, read_seconds = NULL
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
      max_memory_mb_missing, read_seconds = read_seconds
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
#'   wall-time than a smaller one, if it forces enough extra chunks. More
#'   workers is NOT monotonically better on either axis: pass
#'   `read_seconds` (below) to see the actual trade-off
#'   ([.bilatr_worker_tradeoff()]) for your job rather than guessing at
#'   it from this description alone.
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
#' @param read_seconds Only used when `fit` is CSV file paths and `tiers`
#'   includes `3`. Optional wall-time (seconds) for ONE chunk's read,
#'   from your own job's logs -- see [.bilatr_worker_tradeoff()], which
#'   this is passed straight through to. `NULL` (the default) skips the
#'   trade-off entirely: there is no way to estimate wall time at all
#'   without it, so nothing is reported rather than guessed. When
#'   supplied, an extra `message()` names a `n_workers` level that would
#'   give both lower estimated wall time and lower estimated
#'   core-seconds than your current `n_workers`, if the small grid
#'   checked (`n_workers` itself, `1`, half, and double) finds one.
#' @param stan_model (0.7.1) Name registered in `.bilatr_stan_models`; see
#'   [.canonical_stan_model()]. Used
#'   only to decide whether `gamma` (the experimental `stable_gamma`
#'   variant's country-level offset) gets its own report element instead
#'   of flooding Tier 1 -- see `@return`'s `gamma` element. Rhat/ESS
#'   themselves need no orientation info regardless of `stan_model`.
#' @return A list of class `bilatr_diagnostics` with elements:
#'   \describe{
#'     \item{tier1}{Tibble of global/shared diagnostics, one row per
#'       monitored quantity, with a `flagged` column; `NULL` if `1` was
#'       not in `tiers`. Excludes `gamma` when `stan_model` has one (see
#'       `gamma` below) -- `n_tier1_flagged`/`n_tier1_total` in `summary`
#'       count only this table, same meaning as pre-0.7.1.}
#'     \item{gamma}{(0.7.1) Tibble of `stable_gamma`'s country-level
#'       offset diagnostics -- `variable`, `action_index`,
#'       `country_index`, `rhat`, `ess_bulk`, `ess_tail`, `flagged`; one
#'       row per `gamma[k,c]` element. `NULL` whenever `1` was not in
#'       `tiers` or `stan_model` has no `gamma` (silently absent, not an
#'       empty tibble) -- see `.bilatr_model_has_gamma()`. `summary`'s
#'       `n_gamma_flagged`/`n_gamma_total` count this table.}
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
  scratch_dir = NULL,
  read_seconds = NULL,
  stan_model = .BILATR_DEFAULT_MODEL
) {
  # stan_model (0.7.1) is used ONLY to decide has_gamma below (via
  # .bilatr_model_has_gamma()) -- Rhat/ESS are already invariant to the
  # alpha/theta reflection symmetry's sign flip (see
  # .chunked_summarise_csv()'s own docs), so this adds no orientation
  # logic here, unlike extract_theta()/extract_alpha()/etc.
  has_gamma <- .bilatr_model_has_gamma(stan_model)
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
      fit, tiers, max_memory_mb, chunk_size, parallel, n_workers, max_memory_mb_missing,
      read_seconds = read_seconds
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

  .assemble_bilatr_diagnostics(summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold, has_gamma = has_gamma)
}

#' Print a `bilatr_diagnostics` object
#'
#' Tier 1 (global/shared parameters) is always printed in full, since it
#' should never be silently summarized away. `stable_gamma`'s `gamma`
#' (0.7.1) gets its own section instead of being part of that full-print
#' promise -- at production scale it is thousands of parameters, and a
#' routine number will flag as a matter of course, so only the worst few
#' by Rhat and by ESS are shown (silently absent for a model without
#' `gamma`; see [diagnose_convergence()]'s `stan_model` argument). Tier 2
#' (per-dyad hierarchical parameters) is printed as a compact table
#' sorted with dyads flagged as "worse than expected for their sparsity"
#' first. Tier 3 (per-dyad-period latent states) is expected to be noisy
#' for sparse dyads, so it is reported only as aggregate one-line
#' statistics rather than flooding the console with per-dyad-period rows.
#'
#' @param x A `bilatr_diagnostics` object, as returned by
#'   [diagnose_convergence()].
#' @param n_tier2 Maximum number of Tier 2 rows to print.
#' @param n_gamma (0.7.1) Maximum number of `gamma` rows to print, PER
#'   ranking (worst by `rhat`, worst by `ess_bulk` -- so up to `2 *
#'   n_gamma` rows total, fewer if the two rankings overlap or there are
#'   fewer flagged elements than that).
#' @param ... Ignored; present for S3 consistency.
#' @return `x`, invisibly.
#' @export
print.bilatr_diagnostics <- function(x, n_tier2 = 20, n_gamma = 10, ...) {
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

  if (!is.null(x$gamma)) {
    cat(sprintf(
      "== gamma: country-level offsets (%d/%d flagged) ==\n",
      x$summary$n_gamma_flagged, x$summary$n_gamma_total
    ))
    if (x$summary$n_gamma_flagged > 0) {
      flagged_gamma <- dplyr::filter(x$gamma, flagged)
      cat("worst by rhat:\n")
      print(utils::head(dplyr::arrange(flagged_gamma, dplyr::desc(rhat)), n_gamma), n = Inf)
      cat("worst by ess_bulk:\n")
      print(utils::head(dplyr::arrange(flagged_gamma, ess_bulk), n_gamma), n = Inf)
    } else {
      cat("No gamma issues: all country-level offsets meet threshold.\n")
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
