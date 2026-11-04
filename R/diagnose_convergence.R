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
#' `theta_filtered`/`theta_filtered_sd` need no entry here: Tier 3 is
#' matched structurally, by bracket-index count (`n_index >= 2`, see
#' [.classify_bilatr_tier()]), so they land there automatically, as
#' `log_lik[d, t]` does. Their first index is the position within
#' `filter_dyads`, not the `D`-space `dyad_id`, so [.compute_tier3()]
#' gives every Tier 3 base name other than `theta` one row in a separate
#' `tier3_other` table, with no per-dyad breakdown.
#'
#' `gamma` (`stable_gamma` only: `A x n_countries`) is listed here for the
#' same by-name reason as `alpha`: it has two `[...]` indices (action,
#' country), which the structural `n_index >= 2 -> Tier 3` rule would
#' otherwise read as a dyad-indexed quantity, taking the country index for
#' a `dyad_id`. `gamma` is fully identified, unlike `gamma_z`, which is
#' informed only by its prior in its projected-out directions and is
#' excluded instead, via [.bilatr_sign_tied_names()] and the `tier = NA`
#' path.
#' @noRd
.bilatr_tier1_names <- c(
  "alpha", "mu_intercept", "sigma_theta0",
  "mu_log_phi", "sigma_log_phi", "mu_log_noise", "sigma_log_noise",
  "gamma", "sigma_gamma",
  "lp__"
)

#' Names of raw, sampled parameters that are sign-ambiguous under the
#' alpha/theta reflection symmetry, across every registered model
#'
#' Every model's orientation fold (see each `.stan` file's header,
#' "Sign: the ORIENTATION FOLD") corrects the reported
#' `alpha`/`theta`/etc., but the RAW parameters it's built from
#' (`alpha_raw` itself, plus whichever of `z_theta0`/`mu_dyad_raw` and
#' `theta_raw` a given model has) remain genuinely sign-ambiguous: if two
#' chains land in opposite raw-space basins, those variables' own
#' cross-chain Rhat is meaningless even though everything reported is
#' fine. [.classify_bilatr_tier()] uses this to keep them out of the
#' tiered diagnostics tables entirely (`tier = NA`), rather than let a
#' meaningless Rhat surface as a false Tier 1 alarm.
#'
#' Derived from [.bilatr_sign_tied_names()], the single source of truth
#' for which raw names are sign-tied per model, by taking the union of its
#' `raw` component across every registered model. Computed inside the
#' function body, not as a top-level constant, so it does not depend on
#' the order in which the package's files are collated.
#'
#' @return Character vector of variable base names.
#' @noRd
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
#' @noRd
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
#' @return A tibble with columns `variable`, `base_name` (kept so
#'   [.compute_tier3()] can split by base name without re-deriving it),
#'   `tier` (`1L`, `2L`, `3L`, or `NA` for
#'   excluded sign-ambiguous raw parameters), `dyad_id` (the first index,
#'   `NA` outside Tier 2/3), and `time_index` (the second index, `NA`
#'   outside Tier 3).
#' @noRd
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
    dplyr::select(variable, base_name, tier, dyad_id, time_index)
}

#' Normalize the `n_dt` argument to a two-column tibble
#'
#' Accepts a data frame with a dyad column and a count column, or a named
#' numeric vector. The dyads may be given by integer `dyad_id` or by
#' dyad name (as in the `dyad` column of `attr(stan_data, "dyad_ids")`);
#' names are mapped to `dyad_id` through `stan_data`, and names not found
#' there are an error.
#'
#' @param n_dt A data frame with a dyad column (`dyad_id`, else the first
#'   column whose name starts with `"dyad"`) and a count column (matching
#'   `"n_dt"`/`"n_obs"`/`"n_events"`, or named `n`, as from
#'   [dplyr::count()]), or a named numeric vector keyed by dyad id or
#'   name.
#' @param stan_data The Stan data list, or `NULL`. Needed only for
#'   name-keyed input.
#' @return A tibble with columns `dyad_id` (integer) and `n_dt` (numeric).
#' @noRd
.normalize_n_dt <- function(n_dt, stan_data = NULL) {
  if (is.data.frame(n_dt)) {
    nm <- names(n_dt)
    dyad_col <- if ("dyad_id" %in% nm) "dyad_id" else nm[stringr::str_detect(tolower(nm), "^dyad")][1]
    n_dt_col <- nm[stringr::str_detect(tolower(nm), "n_dt|n_obs|n_events")][1]
    if (is.na(n_dt_col) && "n" %in% nm) {
      n_dt_col <- "n"
    }
    if (is.na(dyad_col) || is.na(n_dt_col)) {
      stop(
        "`n_dt` must have a dyad column (`dyad_id`, or a name starting with \"dyad\") and ",
        "a count column (matching \"n_dt\"/\"n_obs\"/\"n_events\", or named \"n\").",
        call. = FALSE
      )
    }
    key <- n_dt[[dyad_col]]
    if (is.factor(key)) key <- as.character(key)
    counts <- n_dt[[n_dt_col]]
  } else if (is.numeric(n_dt) && !is.null(names(n_dt))) {
    key <- names(n_dt)
    counts <- n_dt
  } else {
    stop(
      "`n_dt` must be a data frame with dyad and count columns, or a ",
      "numeric vector named by dyad id or dyad name.",
      call. = FALSE
    )
  }

  dyad_id <- suppressWarnings(as.integer(key))
  if (anyNA(dyad_id[!is.na(key)]) && is.character(key)) {
    dyad_id <- .dyad_names_to_ids(key, stan_data)
  }
  if (anyNA(dyad_id)) {
    stop("`n_dt` has dyad ids that could not be coerced to integers.", call. = FALSE)
  }
  tibble::tibble(dyad_id = dyad_id, n_dt = as.numeric(counts))
}

#' Map dyad names to `dyad_id` through a Stan data list's `dyad_ids`
#'
#' @param dyads Character vector of dyad names.
#' @param stan_data The Stan data list, or `NULL`.
#' @return An integer vector of `dyad_id`s, aligned with `dyads`.
#' @noRd
.dyad_names_to_ids <- function(dyads, stan_data) {
  ids <- attr(stan_data, "dyad_ids", exact = TRUE)
  if (is.null(ids)) {
    stop(
      "`n_dt` is keyed by dyad name, which needs the Stan data (its `dyad_ids` ",
      "attribute) to map names to `dyad_id`. Key `n_dt` by `dyad_id` instead, ",
      "or use a fit that carries its Stan data.",
      call. = FALSE
    )
  }
  lookup <- dplyr::distinct(ids, .data$dyad, .data$dyad_id)
  out <- lookup$dyad_id[match(dyads, lookup$dyad)]
  not_found <- unique(dyads[is.na(out)])
  if (length(not_found) > 0) {
    stop(
      "`n_dt` names dyad(s) not found in the Stan data: ",
      paste(utils::head(not_found, 10), collapse = ", "),
      if (length(not_found) > 10) sprintf(" (and %d more)", length(not_found) - 10),
      ".",
      call. = FALSE
    )
  }
  as.integer(out)
}

#' Per-dyad event totals from a Stan data list
#'
#' The default `n_dt` of [diagnose_convergence()] and
#' [diagnose_and_extract_bilatr()]: `apply(stan_data$Y, 1, sum)`, keyed by
#' `dyad_id`.
#'
#' @param stan_data A Stan data list with `D` and `Y`.
#' @return A tibble with columns `dyad_id` and `n_dt`.
#' @noRd
.n_dt_from_stan_data <- function(stan_data) {
  tibble::tibble(
    dyad_id = seq_len(stan_data$D),
    n_dt = as.numeric(apply(stan_data$Y, 1, sum))
  )
}

#' Resolve the `n_dt` table for tiers 2 and 3
#'
#' @param n_dt The caller's `n_dt`, or `NULL`.
#' @param stan_data The resolved Stan data list, or `NULL`.
#' @param tiers Validated tiers.
#' @param caller The calling function's name, for the error.
#' @return A tibble from [.normalize_n_dt()], or `NULL` when neither tier
#'   2 nor 3 is requested.
#' @noRd
.resolve_n_dt <- function(n_dt, stan_data, tiers, caller) {
  if (!any(c(2L, 3L) %in% tiers)) {
    return(NULL)
  }
  if (!is.null(n_dt)) {
    return(.normalize_n_dt(n_dt, stan_data))
  }
  if (is.null(stan_data) || is.null(stan_data$Y)) {
    stop(
      caller, ": `n_dt` (each dyad's number of events) is needed when `tiers` ",
      "includes 2 or 3, and there is no Stan data to compute it from. Pass `n_dt`, ",
      "e.g. tibble::tibble(dyad_id = 1:D, n_dt = apply(stan_data$Y, 1, sum)); pass ",
      "a fit from fit_panel()/fit_dyad_ts(), or CSV files with their ",
      "*_stan_data.rds, which carry the Stan data; or use `tiers = 1`.",
      call. = FALSE
    )
  }
  .n_dt_from_stan_data(stan_data)
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
#' @noRd
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
#' @noRd
.flag_diagnostic <- function(rhat, ess_bulk, ess_tail, rhat_threshold, ess_threshold) {
  (rhat > rhat_threshold) | (ess_bulk < ess_threshold) | (ess_tail < ess_threshold)
}

#' Validate the `tiers` argument of [diagnose_convergence()]
#'
#' @param tiers Value passed as `diagnose_convergence()`'s `tiers`
#'   argument.
#' @return Sorted, deduplicated integer vector, a subset of `1:3`.
#' @noRd
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
#'   excluded before `flagged` is computed. Used to pull `gamma`
#'   out of this table and into its own (see [.compute_gamma_tier()] and
#'   [.assemble_bilatr_diagnostics()]) -- `gamma` is correctly Tier 1 for
#'   *classification* (see `.bilatr_tier1_names`), but in a large panel
#'   it is `A x n_countries` (thousands of parameters), and Tier 1
#'   is also a *reporting* category designed around a handful of global
#'   scalars: `sigma_gamma[...]` (a per-category scalar, not
#'   country-indexed) is NOT excluded by the `stable_gamma` caller's
#'   pattern and stays in this table.
#' @noRd
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
#' Split out of Tier 1 so that a routine number of flagged `gamma`
#' elements -- dozens to hundreds in a large panel, just from
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
#'   order, e.g. `stan_data`'s `"country_codes"` attribute) for labeling
#'   `country_index`. `NULL` (e.g. plain [diagnose_convergence()], which
#'   has no `stan_data`) leaves the table with `country_index` only.
#' @return A tibble: `variable`, `action_index`, `country_index` (and
#'   `country_code` if `country_codes` is supplied), `rhat`, `ess_bulk`,
#'   `ess_tail`, `flagged`.
#' @noRd
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
#' @noRd
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

#' Compute the Tier 3 (per-dyad-period latent state) diagnostics tibbles
#'
#' Tier 3 is matched structurally by index count (see
#' [.classify_bilatr_tier()]), so it can hold more than `theta[d, t]`:
#' `theta_filtered[k, t]`/`theta_filtered_sd[k, t]` (`k` is the
#' *position within* `filter_dyads`, not `dyad_id`) and `log_lik[d, t]`
#' (`compute_log_lik = 1`) also land here. Grouping every Tier 3 base
#' name by its first index as `dyad_id` would mis-join `theta_filtered*`
#' onto the wrong dyads whenever `filter_dyads` narrows the dyad set,
#' fold `log_lik` into `theta`'s own `min_ess`/`max_rhat`, and inflate
#' `n_theta` to count all of them. So Tier 3 is split in two: `tier3`
#' aggregates `base_name == "theta"` only (matching its name and the
#' per-dyad `dyad_id`/`time_index` semantics the rest of the package
#' assumes), and every other Tier 3 base name gets one row in
#' `tier3_other` -- no per-dyad breakdown, so no `filter_dyads`
#' index-translation is needed there either.
#' @noRd
.compute_tier3 <- function(summ, n_dt_tbl, rhat_threshold, ess_threshold) {
  tier3_raw <- summ %>%
    dplyr::filter(tier == 3L) %>%
    dplyr::mutate(
      rhat_flag = rhat > rhat_threshold,
      ess_flag = (ess_bulk < ess_threshold) | (ess_tail < ess_threshold)
    )

  theta_joined <- tier3_raw %>%
    dplyr::filter(.data$base_name == "theta") %>%
    dplyr::left_join(n_dt_tbl, by = "dyad_id")

  tier3 <- theta_joined %>%
    dplyr::group_by(dyad_id, n_dt) %>%
    dplyr::summarise(
      n_theta = dplyr::n(),
      min_ess_bulk = .min_or_na(ess_bulk),
      min_ess_tail = .min_or_na(ess_tail),
      max_rhat = .max_or_na(rhat),
      share_ess_below_threshold = mean(ess_flag, na.rm = TRUE),
      share_rhat_above_threshold = mean(rhat_flag, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dyad_id)

  tier3_other <- tier3_raw %>%
    dplyr::filter(.data$base_name != "theta") %>%
    dplyr::group_by(.data$base_name) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_rhat_above = sum(.data$rhat_flag, na.rm = TRUE),
      n_ess_below = sum(.data$ess_flag, na.rm = TRUE),
      max_rhat = .max_or_na(.data$rhat),
      min_ess_bulk = .min_or_na(.data$ess_bulk),
      min_ess_tail = .min_or_na(.data$ess_tail),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$base_name)

  list(tier3 = tier3, tier3_other = tier3_other)
}

#' `max()`/`min()` that return `NA` for an all-`NA` input
#'
#' For per-group summaries: a group whose values are all `NA` (e.g.
#' R-hat of constant draws) gives `NA`, not -Inf/Inf with a warning.
#'
#' @param x A numeric vector.
#' @return A single number, or `NA_real_`.
#' @noRd
.max_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

#' @noRd
.min_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
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
#' @param has_gamma Whether the model being diagnosed has a
#'   `gamma` (country-level offset) parameter -- from
#'   [.bilatr_model_has_gamma()], never inferred by checking `summ` for a
#'   `gamma` name. When `TRUE` and Tier 1 was requested, `gamma` is
#'   pulled out of `tier1` into its own `gamma` element (see
#'   [.compute_gamma_tier()]); when `FALSE`, `gamma` is `NULL` and
#'   `n_gamma_*` are `NA_integer_`, so a model without `gamma` gets no
#'   trace of the section at all.
#' @param country_codes Optional character vector for labeling
#'   `gamma`'s `country_index`; see [.compute_gamma_tier()].
#' @return A list of class `bilatr_diagnostics`; see
#'   [diagnose_convergence()]'s `@return` for the element-by-element
#'   description.
#' @noRd
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
  tier3_result <- if (3L %in% tiers) .compute_tier3(summ, n_dt_tbl, rhat_threshold, ess_threshold) else NULL
  tier3 <- tier3_result$tier3
  tier3_other <- tier3_result$tier3_other

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
    n_tier3_other_groups = if (!is.null(tier3_other)) nrow(tier3_other) else NA_integer_,
    n_tier3_other_rhat_above = if (!is.null(tier3_other)) sum(tier3_other$n_rhat_above, na.rm = TRUE) else NA_integer_,
    n_tier3_other_ess_below = if (!is.null(tier3_other)) sum(tier3_other$n_ess_below, na.rm = TRUE) else NA_integer_,
    rhat_threshold = rhat_threshold,
    ess_threshold = ess_threshold
  )

  structure(
    list(
      tier1 = tier1, gamma = gamma_diag, tier2 = tier2, tier3 = tier3,
      tier3_other = tier3_other, summary = summary_info
    ),
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
#' @noRd
.dot_name_to_bracket <- function(name) {
  parts <- strsplit(name, ".", fixed = TRUE)[[1]]
  if (length(parts) == 1) {
    return(name)
  }
  paste0(parts[1], "[", paste(parts[-1], collapse = ","), "]")
}

#' Memory-model constant: a FLOOR for R and its loaded packages, used
#' only where the runtime measurement in
#' [.bilatr_estimate_chunk_baseline_mb()] is unavailable
#'
#' Unlike every other `.BILATR_CHUNK_*` constant, this is not a
#' multiplier of `raw_mb` -- it is a fixed number of MB, added once
#' alongside `file_mb` in [.compute_chunk_size()]/
#' [.estimate_diagnostics_memory_mb()], because it does not shrink with
#' chunk size any more than `file_mb` does: it is the cost of bilatr and
#' its dependencies (`data.table`, `posterior`, `dplyr`, etc.) being
#' loaded, before a single byte of any CSV is read. A session with only
#' bilatr loaded measured about 160 MB peak RSS (R 4.5, macOS); 200
#' rounds that up for margin. It depends on the R build and the loaded
#' namespaces, which is why [.bilatr_estimate_chunk_baseline_mb()]
#' measures the caller's own process and uses this value only as the
#' floor it falls back to (a platform without `/proc`, or `gc()` itself
#' erroring), so the model never silently under-counts to zero.
#' @noRd
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
#' @noRd
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

#' Memory-model constant: the preallocated result array plus
#' `summarise_draws()`'s working set
#'
#' See [.bilatr_chunk_overhead_multiplier()] for how this and the other
#' `.BILATR_CHUNK_*_FACTOR` constants combine into the `k` multiplier
#' [.compute_chunk_size()]/[.estimate_diagnostics_memory_mb()] apply to a
#' chunk's raw byte count (`raw_mb = n_draws * n_chains * chunk_size * 8
#' / 1e6`). Derived from two measurements at `n_cores = 1` (R 4.5,
#' macOS): a single-chain, 60,000-Tier-3-column, 300-draw synthetic CSV
#' (`raw_mb = 144`, `file_mb` ~= 202) peaked at 1002 MB RSS, and the
#' same shape with 2 chains (`raw_mb = 288`) at 1409 MB. Subtracting
#' [.BILATR_CHUNK_BASELINE_MB] and `file_mb` from each and solving the
#' two-chain-count system for `ARRAY + PER_CHAIN / n_chains` gives
#' `ARRAY ~= 2.8`, `PER_CHAIN ~= 1.6`, rounded up here for margin. Most
#' of `ARRAY` is `summarise_draws()`'s own working set
#' (rank-normalization, split-Rhat, bulk/tail ESS all build same-shaped
#' intermediate arrays), not the read.
#' @noRd
.BILATR_CHUNK_ARRAY_FACTOR <- 3

#' Memory-model constant: `fread()`'s parse buffer and the one chain
#' matrix in flight, together sized for a single chain
#'
#' Applied as `.BILATR_CHUNK_PER_CHAIN_FACTOR / n_chains` (each term is
#' `n_draws * chunk_size * 8`, i.e. `raw_mb / n_chains`, not `raw_mb`
#' itself) -- see [.bilatr_chunk_overhead_multiplier()] and
#' [.BILATR_CHUNK_ARRAY_FACTOR]'s docs for the joint derivation.
#' @noRd
.BILATR_CHUNK_PER_CHAIN_FACTOR <- 2

#' Memory-model constant: the one-time cost of `posterior::
#' summarise_draws()` forking at all (`.cores > 1`)
#'
#' Applied once, whenever `n_cores > 1`, together with
#' [.BILATR_CHUNK_CORES_PER_CORE_FACTOR] -- see
#' [.bilatr_chunk_overhead_multiplier()] and that constant's docs for the
#' measurement.
#'
#' **Conservative by design:** the measurement sums RSS across the fork
#' tree on macOS, which double-counts copy-on-write pages that a Linux
#' cgroup charges only once. A smaller Linux measurement (summed Pss and
#' the cgroup's peak-usage counter) found about 3.8x `raw_mb` at the
#' first extra worker and 2.6x per worker after that, so 4.5 is roughly
#' 1.5-2x conservative for `n_cores > 1`. The larger value is kept: an
#' estimate used to size jobs is more useful too large than too small.
#' @noRd
.BILATR_CHUNK_CORES_STEP_FACTOR <- 4.5

#' Memory-model constant: the additional cost of each forked worker
#' beyond the first, when `.cores > 1`
#'
#' `posterior::summarise_draws(.cores = k)` forks `k` worker processes
#' via `parallel::mclapply()`, each computing rank-normalized Rhat/bulk
#' ESS/tail ESS (each of which builds several same-shaped intermediate
#' arrays) over its own slice of variables, so the cost grows with
#' `n_cores`. Measured on a single chain, 60,000 Tier-3 columns and 300
#' draws, with peak RSS summed across the whole process tree (R 4.5,
#' macOS), at `n_cores` = 1/2/4/8: 1002 / 1822 / 2378 / 3571 MB. Fitting
#' `extra beyond n_cores = 1 = STEP + PER_CORE * (n_cores - 1)` in units
#' of `raw_mb` (144 MB here) over the 2->4 and 4->8 intervals gives
#' `PER_CORE ~= 1.9-2.1` and `STEP ~= 3.0-3.3`; both are rounded up here
#' for margin. Summing RSS across a fork tree over-counts shared pages
#' (see [.BILATR_CHUNK_CORES_STEP_FACTOR]), so this is a conservative
#' upper bound; the Linux measurement gave about 2.6, consistent within
#' noise. **Only measured up to `n_cores = 8`**: well beyond that the
#' linear term is extrapolated, so check actual memory use before
#' relying on the estimate for a much larger `n_cores`.
#' @noRd
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
#' term (`.BILATR_CHUNK_CORES_PER_CORE_FACTOR * (n_cores - 1)`).
#'
#' @param n_chains From [.prepare_fast_csv_read].
#' @param n_cores See [diagnose_convergence()]'s `n_cores`.
#' @return The multiplier `k`, in units of `raw_mb`.
#' @noRd
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
#' @param n_cores See [diagnose_convergence()]'s `n_cores`.
#' @return Integer chunk size, at least `1`.
#' @noRd
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
#' @noRd
.estimate_diagnostics_memory_mb <- function(n_draws, n_chains, chunk_size, n_cores, file_mb) {
  k <- .bilatr_chunk_overhead_multiplier(n_chains, n_cores)
  raw_mb <- n_draws * n_chains * chunk_size * 8 / 1e6
  .bilatr_estimate_chunk_baseline_mb() + file_mb + raw_mb * k
}

#' Summarize draws with the fixed measure set the CSV-path functions use
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
#' chunk gets summarized instead of relying on the defaults.
#'
#' @param draws A `posterior::draws_array`/`draws_df` to summarize.
#' @param n_cores Passed through to `summarise_draws()`'s `.cores`:
#'   `posterior::summarise_draws.draws()` splits `draws`' variables into
#'   `n_cores` slices and forks one process per slice
#'   (`parallel::mclapply()`, copy-on-write on Unix), so this parallelizes
#'   the CPU-bound Rhat/ESS computation, not the (sequential, disk-bound)
#'   read that produced `draws`.
#' @return A tibble, one row per variable in `draws`.
#' @noRd
.summarise_bilatr_draws <- function(draws, n_cores = 1L) {
  posterior::summarise_draws(
    draws, "mean", posterior::quantile2, "median", "rhat", "ess_bulk", "ess_tail",
    .cores = n_cores
  )
}

#' Read and summarize CSV variables in memory-bounded chunks
#'
#' Shared by [diagnose_convergence()]'s CSV-path branch and
#' [extract_theta()]'s: reads and summarizes `variables` in groups of
#' `chunk_size`, discarding each chunk's draws before moving to the next
#' -- this is what keeps peak memory bounded regardless of how many
#' variables there are in total. Chunks are always read strictly
#' sequentially (reading is disk-bound and single-threaded regardless --
#' see [.fast_read_post_warmup_draws] -- so there is nothing to gain, and
#' several chunks' worth of memory to lose, by reading more than one at a
#' time); `n_cores` instead parallelizes the CPU-bound part, each chunk's
#' own [.summarise_bilatr_draws()] call, via `posterior::summarise_draws()`'s
#' `.cores` argument.
#'
#' Reads each chunk via [.fast_read_post_warmup_draws] (see its docs for
#' what that read costs per call), discarding the chunk's draws once
#' summarized so peak memory reflects one chunk at a time, per
#' [.compute_chunk_size()]'s accounting.
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param variables Character vector of variable names to summarize.
#' @param chunk_size Variables per chunk.
#' @param n_cores Cores for [.summarise_bilatr_draws()]'s `.cores`; `1L`
#'   for a fully sequential run. See [diagnose_convergence()]'s `n_cores`.
#' @return A tibble, the row-bound [.summarise_bilatr_draws()] output
#'   across all chunks.
#' @noRd
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
#'   rank-normalized Rhat/bulk-tail-ESS), per Tier 3 variable, on ONE
#'   core. Defaults to `0.0026` (~2.6 ms/variable), measured on a
#'   60,000-variable, 300-draw Tier 3 sweep (R 4.5, macOS); scales
#'   with however expensive computing those statistics is for your
#'   posterior (number of draws, mostly), so treat the default as a
#'   starting point, not a calibration for your own data.
#' @return A tibble, one row per `worker_levels` entry: `n_workers`,
#'   `chunk_size`, `n_chunks`, `peak_mb`, `wall_seconds`, `core_seconds`
#'   (the last four `NA` where infeasible).
#' @noRd
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
#' branches. If `max_memory_mb` was left at its default, emits a
#' `message()` summarizing the cost model: reading is single-threaded
#' and each chunk's cost is roughly a full parse of the largest chain
#' file regardless of how many variables are kept (see
#' [.fast_read_post_warmup_draws]), so fewer, larger chunks are cheaper
#' in wall-time, and the CPU-bound summary step's wall-time (not memory)
#' is what `n_cores` trades against. A second `message()`,
#' always emitted, reports the resolved chunk count/size and the
#' estimated peak broken into its file-sized and per-chunk terms
#' separately (see [.compute_chunk_size()]), so a caller can see the
#' actual budget being spent before a long run starts. For the
#' wall-time/core-hour trade-off across candidate `n_cores` values, see
#' [.bilatr_worker_tradeoff()], which takes a wall-time reading from your
#' own job's logs (there is no way to estimate it without one, so this
#' function does not attempt to).
#'
#' @param n_vars Number of variables the chunked sweep will cover (for
#'   the reporting message only; does not affect the chunk_size
#'   calculation itself).
#' @param prepared Output of [.prepare_fast_csv_read] (for `n_draws`,
#'   `n_chains`, `file_mb`).
#' @param max_memory_mb,chunk_size See [diagnose_convergence()].
#' @param n_cores See [diagnose_convergence()]'s `n_cores`.
#' @param max_memory_mb_missing Whether the caller left `max_memory_mb`
#'   at its default (via `missing()` in the calling function).
#' @return The resolved integer chunk size.
#' @noRd
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
      round(max_memory_mb / 1024, 1), " GB); each chunk re-reads the largest ",
      "chain file, so set the largest budget you can afford. More n_cores ",
      "speeds up the summaries but takes memory from the same budget, so ",
      "check the chunk count below."
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
    if (n_cores > 1) paste0(" [", n_cores, " core(s) for summarizing]") else "",
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
#' @param tiers,max_memory_mb,chunk_size,n_cores See
#'   [diagnose_convergence()].
#' @param max_memory_mb_missing Whether the caller left `max_memory_mb`
#'   at its default (via `missing()` in [diagnose_convergence()]) --
#'   gates the one-time "this is a default, not a calibrated value"
#'   message.
#' @return A tibble in the same shape [diagnose_convergence()]'s
#'   in-memory branch produces: [posterior::summarise_draws()] columns
#'   left-joined with [.classify_bilatr_tier]'s `tier`/`dyad_id`/
#'   `time_index`.
#' @noRd
.read_diagnostics_summary_from_csv <- function(
  csv_files, tiers, max_memory_mb, chunk_size, n_cores,
  max_memory_mb_missing
) {
  prepared <- .prepare_fast_csv_read(csv_files)

  var_tiers <- .classify_bilatr_tier(prepared$variables)
  keep_tiers <- var_tiers[var_tiers$tier %in% tiers, ]

  tier12_vars <- keep_tiers$variable[keep_tiers$tier %in% c(1L, 2L)]
  tier3_vars <- keep_tiers$variable[keep_tiers$tier == 3L]

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
#' Computes R-hat and effective sample sizes for a fitted bilatr model
#' and splits them into three tiers before applying thresholds, so that
#' problems in the global parameters are never masked by the wide,
#' weakly identified posteriors expected for sparsely observed dyads. See
#' `vignette("diagnostics", package = "bilatr")` for a worked example and
#' the rationale.
#'
#' * **Tier 1**, global parameters (`alpha`, `mu_intercept`, the
#'   hierarchical scales and `lp__`), is always reported in full.
#' * **Tier 2**, per-dyad parameters (e.g. `phi`, `process_noise`,
#'   `theta0`), and **Tier 3**, the per-dyad-period latent states
#'   (`theta`), are joined to each dyad's event count `n_dt`. A dyad is
#'   flagged only when its diagnostics are worse than the trend of
#'   diagnostics against sparsity predicts, so sparse but unremarkable
#'   dyads do not flood the report.
#'
#' Quantities are assigned to tiers by name and by their number of
#' indices: Tier 1 names are listed explicitly, a single index means a
#' per-dyad parameter (Tier 2), and two indices a per-dyad-period
#' quantity (Tier 3). The raw parameters the models' orientation fold is
#' built from (`alpha_raw`, `z_theta0` or `mu_dyad_raw`, `theta_raw`, and
#' `gamma_z` in `stable_gamma`) are sign-ambiguous or informed only by
#' their prior, so their R-hat is meaningless; they are left out of every
#' tier.
#'
#' `tiers` limits which tiers are computed: quantities outside the
#' requested tiers are not read at all, rather than read and hidden. Tier
#' 3 usually has by far the most quantities (one per dyad-period), so
#' `tiers = 1` or `tiers = 1:2` is much faster on a large fit.
#'
#' Every kind of `fit` is summarized with the same measures (`mean`, the
#' 5% and 95% quantiles, `median`, `rhat`, `ess_bulk`, `ess_tail`), so the
#' returned tables have the same columns whichever kind was given.
#'
#' @section Large fits from CSV files:
#' With CmdStan CSV files, Tier 1 and 2 are read in one call and Tier 3 is
#' read and summarized in chunks, discarding each chunk before reading the
#' next, so the full set of draws is never in memory at once. Reading is
#' single-threaded and each chunk costs about one pass over the largest
#' chain file, however few columns it keeps, so fewer, larger chunks are
#' faster.
#'
#' * `max_memory_mb` sets the chunk size. Use the largest budget your
#'   machine or job allows. The default, 8 GB, is a guess, and the
#'   function says so when you rely on it. Peak memory is estimated as
#'   the size of the largest chain file, plus about 200 MB for R and the
#'   loaded packages, plus a multiple of the chunk's raw size (draws x
#'   chains x variables x 8 bytes): 3.5 times it for four chains on one
#'   core.
#' * `n_cores` > 1 computes each chunk's R-hat and ESS in that many forked
#'   processes, which is faster but costs memory from the same budget:
#'   roughly 4.5x the chunk's raw size for the first extra process and
#'   2.5x for each one after that (estimates that are conservative and
#'   were measured up to 8 cores). At a fixed `max_memory_mb`, more cores
#'   therefore mean smaller chunks and more passes over the files, and a
#'   large `n_cores` can be slower overall. A message before reading
#'   reports the chunk size, the number of chunks and the estimated peak
#'   memory; check it before a long run.
#'
#' @param fit One of: a `posterior` draws object (`draws_array` or
#'   `draws_df`), already in memory; a `CmdStanMCMC` fit (or any object
#'   with `$metadata()` and `$draws(variables = )` methods), whose
#'   requested tiers are read with `$draws()`; or a character vector of
#'   CmdStan CSV files, one per chain, for example from a run in another
#'   session, which are read in memory-bounded chunks (see "Large fits
#'   from CSV files"). `max_memory_mb`, `chunk_size` and `n_cores` apply
#'   only to CSV files.
#' @param n_dt Each dyad's total number of events, which tiers 2 and 3
#'   screen diagnostics against. `NULL` (default) computes it from the
#'   fit's Stan data (`apply(stan_data$Y, 1, sum)`), found as for
#'   [extract_theta()]: `attr(fit, "stan_data")`, or the
#'   `*_stan_data.rds` file next to CSV paths. To pass your own: a data
#'   frame with a dyad column (`dyad_id`, or `dyad` with dyad names as in
#'   `attr(stan_data, "dyad_ids")`) and a count column (`n_dt`, `n_obs`,
#'   `n_events` or `n`), or a numeric vector named by dyad id or name.
#'   Not used when `tiers = 1`. With a `draws` object, which carries no
#'   Stan data, `n_dt` must be supplied for tiers 2 and 3.
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
#' @param max_memory_mb CSV files with Tier 3 only. The memory budget, in
#'   MB, that sets how many Tier 3 variables are read per chunk. Defaults
#'   to `8192` (8 GB). It must cover at least the largest chain file plus
#'   about 200 MB; otherwise the function stops and says how much is
#'   needed. The estimate is a guide, not a guarantee. See "Large fits
#'   from CSV files".
#' @param chunk_size CSV files with Tier 3 only. The number of Tier 3
#'   variables per chunk, to override the one derived from
#'   `max_memory_mb`. `NULL` (default) derives it.
#' @param n_cores CSV files with Tier 3 only. The number of processes that
#'   compute R-hat and ESS for each chunk (chunks are always read one at a
#'   time). Defaults to `1`. More cores speed up the summaries but use
#'   memory from the same `max_memory_mb` budget; see "Large fits from CSV
#'   files". Under a job scheduler or in a container,
#'   `parallelly::availableCores()` gives the cores actually allocated,
#'   unlike `parallel::detectCores()`.
#'   On Unix-alikes, a session that fits a model between two `n_cores > 1`
#'   calls can print "Error while shutting down parallel" when R exits;
#'   this is harmless.
#' @param stan_model The model that was fit: `"stable"`, `"ou"` or
#'   `"stable_gamma"`. `NULL` (default) uses the model recorded on the fit
#'   (by [fit_panel()]/[fit_dyad_ts()]), on its Stan data, or in the
#'   CmdStan CSV header, else assumes `"stable"` with a message. Giving a
#'   model that contradicts the recorded one is an error. Used only to
#'   decide whether `gamma` (the `stable_gamma` model's country-level
#'   offset) gets its own report element instead of flooding Tier 1 (see
#'   the `gamma` element below).
#' @return A list of class `bilatr_diagnostics` with elements:
#'   \describe{
#'     \item{tier1}{Tibble of global/shared diagnostics, one row per
#'       monitored quantity, with a `flagged` column; `NULL` if `1` was
#'       not in `tiers`. Excludes `gamma` when the model has one (see
#'       `gamma` below); `n_tier1_flagged`/`n_tier1_total` in `summary`
#'       count only this table.}
#'     \item{gamma}{Tibble of `stable_gamma`'s country-level offset
#'       diagnostics -- `variable`, `action_index`, `country_index` (and
#'       `country_code` when the Stan data is available), `rhat`,
#'       `ess_bulk`, `ess_tail`, `flagged`; one row per `gamma[k,c]`
#'       element. `NULL` when `1` was not in `tiers` or the model has no
#'       `gamma`. `summary`'s `n_gamma_flagged`/`n_gamma_total` count this
#'       table.}
#'     \item{tier2}{Tibble with one row per dyad found in `n_dt` (or in
#'       the draws, if unmatched), per-dyad-parameter Rhat/ESS columns,
#'       and a `worse_than_expected` column; `NULL` if `2` was not in
#'       `tiers`.}
#'     \item{tier3}{Tibble with one row per dyad summarizing `theta`'s
#'       per-dyad-period diagnostics (minimum ESS, maximum R-hat, share of
#'       entries breaching the thresholds); `NULL` if `3` was not in
#'       `tiers`. Holds `theta` only; other per-dyad-period quantities are
#'       in `tier3_other`.}
#'     \item{tier3_other}{Tibble with one row per Tier 3 quantity other
#'       than `theta` (`theta_filtered`, `theta_filtered_sd`, `log_lik`,
#'       whichever are present): `n`, `n_rhat_above`, `n_ess_below`,
#'       `max_rhat`, `min_ess_bulk`, `min_ess_tail`, with no per-dyad
#'       breakdown (`theta_filtered`'s first index is the position in
#'       `filter_dyads`, not a dyad id). `NULL` if `3` was not in `tiers`;
#'       zero rows if the fit has no Tier 3 quantity besides `theta`.}
#'     \item{summary}{A short named list of headline counts (see
#'       [print.bilatr_diagnostics()]), including `n_tier3_other_groups`,
#'       `n_tier3_other_rhat_above` and `n_tier3_other_ess_below`.}
#'   }
#' @references Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and
#'   Burkner, P.-C. (2021). Rank-normalization, folding, and localization:
#'   an improved R-hat for assessing convergence of MCMC. *Bayesian
#'   Analysis* 16(2), 667-718. \doi{10.1214/20-BA1221}
#' @examples
#' # CmdStan CSV files from an earlier run, with the *_stan_data.rds file
#' # that fit_panel(output_dir = ) wrote next to them (here, the example
#' # fit shipped with bilatr): Tier 3 is read and summarized in
#' # memory-bounded chunks rather than all at once. The example fit is far
#' # too short to converge, so most of Tier 1 is flagged.
#' csv_files <- bilatr_example_fit()
#' diag <- diagnose_convergence(csv_files, max_memory_mb = 4096)
#' diag
#' diag$tier1
#'
#' # only the cheap, always-important global/shared parameters:
#' diagnose_convergence(csv_files, tiers = 1)
#'
#' # a stable_gamma fit gets its own `gamma` report element
#' diagnose_convergence(bilatr_example_fit("stable_gamma"), tiers = 1)$gamma
#' @export
diagnose_convergence <- function(
  fit,
  n_dt = NULL,
  rhat_threshold = 1.01,
  ess_threshold = 400,
  tiers = 1:3,
  max_memory_mb = 8192,
  chunk_size = NULL,
  n_cores = 1L,
  stan_model = NULL
) {
  caller <- "diagnose_convergence()"
  max_memory_mb_missing <- missing(max_memory_mb)
  tiers <- .validate_tiers(tiers)

  # The Stan data is looked up only when it is needed: for the model
  # (when not given) and for the default n_dt.
  needs_n_dt <- any(c(2L, 3L) %in% tiers) && is.null(n_dt)
  stan_data <- if (needs_n_dt || is.null(stan_model) || !is.null(n_dt)) {
    .stan_data_from_fit(fit, NULL, caller, required = FALSE)
  }
  # stan_model only decides has_gamma below: Rhat/ESS are invariant to
  # the alpha/theta reflection's sign flip, so no orientation logic here.
  stan_model <- .resolve_stan_model_for_fit(fit, stan_data, stan_model, caller)
  has_gamma <- .bilatr_model_has_gamma(stan_model)
  n_dt_tbl <- .resolve_n_dt(n_dt, stan_data, tiers, caller)

  if (is.character(fit)) {
    summ <- .read_diagnostics_summary_from_csv(
      fit, tiers, max_memory_mb, chunk_size, n_cores, max_memory_mb_missing
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

  .assemble_bilatr_diagnostics(
    summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold,
    has_gamma = has_gamma, country_codes = attr(stan_data, "country_codes")
  )
}

#' Print a `bilatr_diagnostics` object
#'
#' Tier 1 (global/shared parameters) is always printed in full, since it
#' should never be silently summarized away. `stable_gamma`'s `gamma`
#' gets its own section instead: in a large panel it has thousands of
#' elements, and a routine number will be flagged, so only the worst few
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
#' @param n_gamma Maximum number of `gamma` rows to print, PER
#'   ranking (worst by `rhat`, worst by `ess_bulk` -- so up to `2 *
#'   n_gamma` rows total, fewer if the two rankings overlap or there are
#'   fewer flagged elements than that).
#' @param ... Ignored; present for S3 consistency.
#' @return `x`, invisibly.
#' @examples
#' diag <- diagnose_convergence(bilatr_example_fit(), tiers = 1:2)
#' print(diag, n_tier2 = 3)
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
      "== Tier 3: per-dyad-period latent states (theta, aggregated; %d dyads) ==\n",
      x$summary$n_dyads_tier3
    ))
    cat(sprintf(
      "%.0f%% of dyads have min theta ESS below threshold (%s)%s.\n",
      100 * share_below,
      x$summary$ess_threshold,
      if (!is.na(low_n_dt_share)) {
        sprintf(", median n_dt among them: %.0f", low_n_dt_share)
      } else {
        ""
      }
    ))

    if (!is.null(x$tier3_other) && nrow(x$tier3_other) > 0) {
      for (i in seq_len(nrow(x$tier3_other))) {
        row <- x$tier3_other[i, ]
        cat(sprintf(
          "  %s: %d entries, %d above Rhat threshold, %d below ESS threshold (max Rhat %.3f, min ESS bulk/tail %.0f/%.0f).\n",
          row$base_name, row$n, row$n_rhat_above, row$n_ess_below,
          row$max_rhat, row$min_ess_bulk, row$min_ess_tail
        ))
      }
    }
  }

  invisible(x)
}
