#' Run convergence diagnostics and extract theta/alpha/mu_intercept from raw
#' CmdStan CSVs in a single shared read
#'
#' A fused counterpart to calling [diagnose_convergence()], [extract_theta()],
#' [extract_alpha()], and [extract_mu_intercept()] separately against the
#' same set of CSV files. Calling them separately re-reads overlapping data
#' from disk: [diagnose_convergence()]'s Tier 3 sweep already covers every
#' `theta`/`theta_raw` column [extract_theta()] needs, and its Tier 1/2
#' read already covers every `alpha`/`mu_intercept` column
#' [extract_alpha()]/[extract_mu_intercept()] need. `data.table::fread()`
#' touches each file's full byte range once per call regardless of how
#' few columns are requested (see [.fast_read_post_warmup_draws]), so
#' each of those extra calls pays close to a full-file touch for data
#' already summarised once. This function reads Tier 1/2 once and Tier 3
#' once (chunked, same as [diagnose_convergence()]) via
#' [.read_diagnostics_summary_from_csv()] -- the same helper
#' [diagnose_convergence()]'s own CSV-path branch uses -- and derives
#' every output from that one summary rather than re-reading anything.
#' Sharing that helper (rather than a second, orientation-aware copy of
#' its tier/chunk logic, as before 0.10.0) means the two paths cannot
#' drift apart.
#'
#' Total file touches per chain file, at production scale (Tier 3 too
#' large for one chunk): one header/config scan and one row-count
#' validation probe, both in [.prepare_fast_csv_read]; one read for Tier
#' 1/2; and one per Tier 3 chunk. When Tier 3 fits in a single chunk
#' (small panels, or a large `max_memory_mb`), Tier 1/2 is folded into
#' that same read instead of being a separate one -- see
#' [.read_diagnostics_summary_from_csv()].
#'
#' Only meaningful for the CSV-file-paths case: an in-memory `fit` has no
#' repeated-disk-read cost to avoid, so there is nothing to fuse there --
#' call the four functions separately instead.
#'
#' Unlike [extract_alpha()]/[extract_mu_intercept()]'s standalone `probs`
#' argument, the `alpha`/`mu_intercept` outputs here always carry the same
#' fixed `mean`/`5%`/`50%`/`95%` columns as [extract_theta()]'s CSV path
#' (it reuses the same underlying [posterior::summarise_draws()] default
#' summary, not a second call with different quantiles); call
#' [extract_alpha()]/[extract_mu_intercept()] directly if you need other
#' quantiles.
#'
#' @inheritParams diagnose_convergence
#' @param csv_files Character vector of raw CmdStan CSV file paths (one per
#'   chain). Required; this function does not accept an in-memory fit
#'   (see Details).
#' @param stan_data The Stan data list used to produce the fit, from
#'   [prepare_stan_data()] (must carry its `dyad_ids` attribute), for
#'   [extract_theta()]'s dyad-identifier join. `NULL` (default) reads the
#'   single `*_stan_data.rds` file next to `csv_files` (which
#'   [fit_panel()]/[fit_dyad_ts()] write to their `output_dir`), with a
#'   `message()`; errors if there is none or more than one.
#' @param n_dt Each dyad's total number of events, for tiers 2 and 3.
#'   `NULL` (default) computes it from `stan_data`; see
#'   [diagnose_convergence()] for passing your own.
#' @param stan_model The model that was fit. `NULL` (default) uses the
#'   model recorded on `stan_data` or in the CmdStan CSV header; see
#'   [diagnose_convergence()].
#' @param event_classes Optional character vector of event-class labels,
#'   passed through to the `alpha`/`mu_intercept` outputs exactly as in
#'   [extract_alpha()]. `NULL` (default) uses `stan_data`'s
#'   `"event_classes"` attribute, if present; an explicit value wins. When supplied, `alpha`/`mu_intercept` also get
#'   a `class_label` column (0.10.1), from `event_class_labels(stan_data)`
#'   (see [event_class_labels()]) unless overridden by `class_labels`.
#' @param class_labels Optional; see [.resolve_class_labels()]. Ignored
#'   when `event_classes` is `NULL` (nothing to key labels by). Defaults
#'   to `NULL`, which uses `event_class_labels(stan_data)` (see
#'   [event_class_labels()]).
#' @return A list with elements `diagnostics` (a `bilatr_diagnostics`
#'   object, as from [diagnose_convergence()]), `theta`, `alpha`, and
#'   `mu_intercept` (tibbles, in the same shape [extract_theta()]/
#'   [extract_alpha()]/[extract_mu_intercept()] return, plus a
#'   `class_label` column on `alpha`/`mu_intercept` when `event_classes`
#'   is supplied), plus
#'   `theta_filtered`/`theta_filtered_sd` (empty tibbles unless
#'   `compute_theta_filtered = 1` was set when the data was prepared) --
#'   same `dyad_id`/`time_index`/`dyad_ids`-joined shape as `theta`, with
#'   the position-within-`filter_dyads` index already translated back to
#'   the true `dyad_id` via `stan_data$filter_dyads`.
#'
#'   `theta`, `theta_filtered` and `theta_filtered_sd` are built from the
#'   Tier 3 summary, so they are empty tibbles unless `tiers` includes 3.
#'   `alpha` and `mu_intercept` come from Tier 1 and are empty unless
#'   `tiers` includes 1.
#' @examples
#' # CmdStan CSVs with the *_stan_data.rds file that fit_panel(output_dir = )
#' # wrote next to them; here, the example fit shipped with bilatr
#' result <- diagnose_and_extract_bilatr(bilatr_example_fit())
#' result$diagnostics
#' result$theta
#' result$alpha
#' @export
diagnose_and_extract_bilatr <- function(
  csv_files, stan_data = NULL, n_dt = NULL,
  stan_model = NULL,
  rhat_threshold = 1.01, ess_threshold = 400, tiers = 1:3,
  event_classes = NULL,
  class_labels = NULL,
  max_memory_mb = 8192, chunk_size = NULL, n_cores = 1L
) {
  caller <- "diagnose_and_extract_bilatr()"
  if (!is.character(csv_files)) {
    stop(
      "diagnose_and_extract_bilatr() only supports raw CmdStan CSV file ",
      "paths -- its entire purpose is sharing a single read across ",
      "diagnostics and extraction, and an in-memory fit has no repeated-",
      "disk-read cost to avoid. Call diagnose_convergence()/",
      "extract_theta()/extract_alpha()/extract_mu_intercept() separately ",
      "for an in-memory fit.",
      call. = FALSE
    )
  }
  max_memory_mb_missing <- missing(max_memory_mb)

  tiers <- .validate_tiers(tiers)

  stan_data <- .stan_data_from_fit(csv_files, stan_data, caller)
  stan_model <- .resolve_stan_model_for_fit(csv_files, stan_data, stan_model, caller)
  n_dt_tbl <- .resolve_n_dt(n_dt, stan_data, tiers, caller)
  event_classes <- event_classes %||% attr(stan_data, "event_classes")
  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be a Stan data list from prepare_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }

  # has_gamma/country_codes for .assemble_bilatr_diagnostics()'s
  # gamma-vs-Tier-1 split, with country_code labels on the gamma table.
  has_gamma <- .bilatr_model_has_gamma(stan_model)
  country_codes <- attr(stan_data, "country_codes")

  summ <- .read_diagnostics_summary_from_csv(
    csv_files, tiers, max_memory_mb, chunk_size, n_cores,
    max_memory_mb_missing
  )

  diagnostics <- .assemble_bilatr_diagnostics(
    summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold,
    has_gamma = has_gamma, country_codes = country_codes
  )

  extract_from_summ <- function(prefix) {
    summ %>%
      dplyr::filter(startsWith(variable, prefix)) %>%
      dplyr::select(variable, mean, `5%` = q5, `50%` = median, `95%` = q95)
  }

  theta <- extract_from_summ("theta[") %>%
    dplyr::mutate(variable = stringr::str_remove_all(variable, "theta\\[|\\]")) %>%
    tidyr::separate(variable, into = c("dyad_id", "time_index"), sep = ",", convert = TRUE) %>%
    dplyr::left_join(dyad_ids, by = c("dyad_id", "time_index"))

  # theta_filtered/theta_filtered_sd (0.5.0+, present only if
  # compute_theta_filtered = 1 was set): their first index is the
  # POSITION within filter_dyads, not the true D-space dyad_id, whenever
  # filter_dyads narrowed the dyad set -- see R/diagnose_convergence.R's
  # .bilatr_tier1_names docs for the caveat this resolves. stan_data$
  # filter_dyads (prepare_stan_data()'s own resolved integer vector, 1:D
  # when filter_dyads = NULL was used) maps that position back to the
  # true dyad_id before joining dyad_ids, so this is correct in both the
  # "all dyads filtered" and "an explicit subset" cases -- not just the
  # former.
  filtered_extract <- function(prefix) {
    extract_from_summ(prefix) %>%
      dplyr::mutate(idx = stringr::str_match(variable, "\\[(\\d+),(\\d+)\\]")) %>%
      dplyr::mutate(
        filter_index = as.integer(idx[, 2]),
        time_index = as.integer(idx[, 3])
      ) %>%
      dplyr::select(-idx) %>%
      dplyr::mutate(dyad_id = stan_data$filter_dyads[filter_index]) %>%
      dplyr::select(-filter_index) %>%
      dplyr::left_join(dyad_ids, by = c("dyad_id", "time_index"))
  }

  theta_filtered <- filtered_extract("theta_filtered[")
  theta_filtered_sd <- filtered_extract("theta_filtered_sd[")

  action_extract <- function(prefix) {
    extract_from_summ(prefix) %>%
      dplyr::mutate(action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])")))
  }

  alpha <- action_extract("alpha[")
  mu_intercept <- action_extract("mu_intercept[")

  if (!is.null(event_classes)) {
    labels <- .resolve_class_labels(class_labels, event_classes, stan_data = stan_data)
    alpha <- dplyr::mutate(
      alpha,
      event_class = event_classes[action_index],
      class_label = labels[action_index]
    )
    mu_intercept <- dplyr::mutate(
      mu_intercept,
      event_class = event_classes[action_index],
      class_label = labels[action_index]
    )
  }

  list(
    diagnostics = diagnostics, theta = theta,
    theta_filtered = theta_filtered, theta_filtered_sd = theta_filtered_sd,
    alpha = alpha, mu_intercept = mu_intercept
  )
}
