#' Run convergence diagnostics and extract theta/alpha/mu_intercept from raw
#' CmdStan CSVs in a single shared read
#'
#' A fused counterpart to calling [diagnose_convergence()], [extract_theta()],
#' [extract_alpha()], and [extract_mu_intercept()] separately against the
#' same set of CSV files. Calling them separately re-reads overlapping data
#' from disk: [diagnose_convergence()]'s Tier 3 sweep already covers every
#' `theta`/`theta_raw` column [extract_theta()] needs, and its Tier 1/2
#' read already covers every `alpha`/`mu_intercept` column
#' [extract_alpha()]/[extract_mu_intercept()] need. Since
#' [cmdstanr::read_cmdstan_csv()] is I/O-bound (see [.chunked_summarise_csv]),
#' each of those extra calls pays close to a full-file-scan cost for data
#' already summarised once. This function reads Tier 1/2 once and Tier 3
#' once (chunked, same as [diagnose_convergence()]), reorients whichever
#' columns [bilatr_orient()] would flip for `stan_model` as part of that
#' single read, and derives all four outputs from those two summaries
#' rather than re-reading anything.
#'
#' Only meaningful for the CSV-file-paths case: an in-memory `fit` has no
#' repeated-disk-read cost to avoid, so there is nothing to fuse there --
#' call the four functions separately instead.
#'
#' Rhat/ESS are invariant to a deterministic sign flip, so reorienting
#' before summarising never changes anything [diagnose_convergence()]
#' reports; it only affects the mean/quantile columns, which is exactly
#' what [extract_theta()]/[extract_alpha()]/[extract_mu_intercept()] need
#' reoriented anyway. For every currently-registered model, Tier 3 is
#' exactly `{theta, theta_raw}`, and [.bilatr_flip_variables()] lists both
#' whenever either needs flipping -- so the whole Tier 3 chunked sweep can
#' be flipped uniformly, the same way [extract_theta()]'s own CSV path
#' already does for `theta` alone.
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
#' @param stan_data The Stan data list used to produce the fit, as returned
#'   by [assemble_stan_data()] (must still carry its `dyad_ids` attribute),
#'   for [extract_theta()]'s dyad-identifier join.
#' @param stan_model Name registered in `.bilatr_stan_models`; see
#'   [bilatr_orient()]. Defaults to `.BILATR_DEFAULT_MODEL`.
#' @param event_classes Optional character vector of event-class labels,
#'   passed through to the `alpha`/`mu_intercept` outputs exactly as in
#'   [extract_alpha()]. Defaults to `stan_data`'s `"event_classes"`
#'   attribute, if present.
#' @return A list with elements `diagnostics` (a `bilatr_diagnostics`
#'   object, as from [diagnose_convergence()]), `theta`, `alpha`, and
#'   `mu_intercept` (tibbles, in the same shape [extract_theta()]/
#'   [extract_alpha()]/[extract_mu_intercept()] return).
#' @examples
#' \dontrun{
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' result <- diagnose_and_extract_bilatr(
#'   csv_files, stan_data, n_dt = n_dt, stan_model = "alphanorm"
#' )
#' result$diagnostics
#' result$theta
#' }
#' @export
diagnose_and_extract_bilatr <- function(
  csv_files, stan_data, n_dt,
  stan_model = .BILATR_DEFAULT_MODEL,
  rhat_threshold = 1.01, ess_threshold = 400, tiers = 1:3,
  event_classes = attr(stan_data, "event_classes"),
  max_memory_mb = 8192, chunk_size = NULL, parallel = FALSE,
  n_workers = max(1L, parallel::detectCores() - 1L)
) {
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
  if (any(c(2L, 3L) %in% tiers) && missing(n_dt)) {
    stop(
      "`n_dt` is required when `tiers` includes 2 and/or 3 (Tier 1 alone ",
      "needs no per-dyad join).",
      call. = FALSE
    )
  }
  n_dt_tbl <- if (any(c(2L, 3L) %in% tiers)) .normalize_n_dt(n_dt) else NULL

  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be the output of assemble_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }

  all_vars <- .stan_csv_variable_names(csv_files[1])
  var_tiers <- .classify_bilatr_tier(all_vars)
  keep_tiers <- var_tiers[var_tiers$tier %in% tiers, ]

  tier12_vars <- keep_tiers$variable[keep_tiers$tier %in% c(1L, 2L)]
  tier3_vars <- keep_tiers$variable[keep_tiers$tier == 3L]

  # Orientation is determined once, from a single dedicated alpha[1] read,
  # rather than relying on alpha[1] happening to already be among
  # tier12_vars (it won't be if `tiers` excludes 1) or reading it twice
  # (once per tier group). Both branches below apply the same flip
  # decision to whichever of their own columns .bilatr_flip_variables()
  # lists, mirroring bilatr_orient()'s exact column-selective negation.
  flip_vars <- .bilatr_flip_variables(stan_model)
  flip <- FALSE
  if (length(flip_vars) > 0) {
    alpha1_draws <- cmdstanr::read_cmdstan_csv(csv_files, variables = "alpha[1]")$post_warmup_draws
    flip <- stats::median(posterior::extract_variable(alpha1_draws, "alpha[1]")) < 0
  }

  # Tier 1/2: rhat/ess are unaffected by the flip, so this single read/
  # summary serves both diagnose_convergence()'s Tier 1/2 tables below
  # and extract_alpha()/extract_mu_intercept()'s already-correctly-
  # oriented mean/quantile columns, with no second read.
  tier12_summ <- NULL
  if (length(tier12_vars) > 0) {
    tier12_draws <- cmdstanr::read_cmdstan_csv(csv_files, variables = tier12_vars)$post_warmup_draws
    if (flip) {
      df <- as.data.frame(posterior::as_draws_df(tier12_draws))
      flip_cols <- .bilatr_match_draws_columns(names(df), flip_vars)
      if (length(flip_cols) > 0) {
        df[flip_cols] <- -df[flip_cols]
        tier12_draws <- posterior::as_draws_df(df)
      }
    }
    tier12_summ <- posterior::summarise_draws(tier12_draws)
  }

  # Tier 3 (theta + theta_raw together), flipped as a whole chunked sweep
  # -- see Details above for why this is exact, not an approximation,
  # for every currently-registered model. Serves both
  # diagnose_convergence()'s Tier 3 table and extract_theta()'s already-
  # correctly-oriented "theta[" subset, with no second sweep over the
  # same columns.
  tier3_summ <- NULL
  if (length(tier3_vars) > 0) {
    chunk_size_used <- .resolve_chunk_size_and_report(
      length(tier3_vars), csv_files, max_memory_mb, chunk_size, parallel, n_workers,
      max_memory_mb_missing
    )
    tier3_summ <- .chunked_summarise_csv(
      csv_files, tier3_vars, chunk_size_used, parallel, n_workers, flip = flip
    )
  }

  summ <- dplyr::left_join(dplyr::bind_rows(tier12_summ, tier3_summ), var_tiers, by = "variable")
  diagnostics <- .assemble_bilatr_diagnostics(summ, n_dt_tbl, tiers, rhat_threshold, ess_threshold)

  extract_from_summ <- function(prefix) {
    summ %>%
      dplyr::filter(startsWith(variable, prefix)) %>%
      dplyr::select(variable, mean, `5%` = q5, `50%` = median, `95%` = q95)
  }

  theta <- extract_from_summ("theta[") %>%
    dplyr::mutate(variable = stringr::str_remove_all(variable, "theta\\[|\\]")) %>%
    tidyr::separate(variable, into = c("dyad_id", "time_index"), sep = ",", convert = TRUE) %>%
    dplyr::left_join(dyad_ids, by = c("dyad_id", "time_index"))

  action_extract <- function(prefix) {
    extract_from_summ(prefix) %>%
      dplyr::mutate(action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])")))
  }

  alpha <- action_extract("alpha[")
  mu_intercept <- action_extract("mu_intercept[")

  if (!is.null(event_classes)) {
    alpha <- dplyr::mutate(alpha, event_class = event_classes[action_index])
    mu_intercept <- dplyr::mutate(mu_intercept, event_class = event_classes[action_index])
  }

  list(diagnostics = diagnostics, theta = theta, alpha = alpha, mu_intercept = mu_intercept)
}
