#' Extract theta trajectories with dyad identifiers reattached
#'
#' Pulls posterior summaries of the latent conflict trajectory `theta`
#' out of a fitted model, indexed by dyad and time period, and reattaches
#' the human-readable dyad identifiers (`dyad`, the undirected `dyad2`,
#' `year`, and `month` if applicable) via the `dyad_ids` attribute that
#' [assemble_stan_data()] attaches to its output. Works the same way for
#' single-dyad ([fit_dyad_ts()]) and panel ([fit_panel()]) fits.
#'
#' Draws are passed through [bilatr_orient()] before summarizing, so
#' `theta`'s sign always reflects the canonical orientation (higher theta
#' = better relations) regardless of which reflection-symmetry basin the
#' sampler actually landed in. This is a no-op for `stan_model`s other
#' than `"alphanorm"`/`"alphanorm_ou"`, so the default preserves existing
#' behavior for `stable`-model fits (its `alpha[1] = 1` is fixed, never
#' negative).
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths
#' (one per chain), matching [diagnose_convergence()]'s CSV-path mode --
#' reads and summarises `theta` in memory-bounded chunks via
#' [.chunked_summarise_csv] rather than materializing the full draws
#' array, for the same production-scale-panel reason `diagnose_convergence()`
#' needed it. **This trades memory for wall-time, and the exchange rate
#' can be steep**: a real benchmark against a production-scale CSV found
#' `read_cmdstan_csv()` strongly I/O-bound (per-chunk cost barely depends
#' on how many columns are requested, so total sweep time scales with
#' chunk COUNT, not memory saved) -- see `max_memory_mb` below and in
#' [diagnose_convergence()]'s documentation. Prefer the in-memory path
#' when you already have `fit` in memory (no benefit to re-reading), and
#' the CSV path's default `max_memory_mb` as large as your job's memory
#' allocation can afford, not as small as "safely" possible. Sign
#' orientation in this mode requires one small extra read of `alpha[1]`
#' up front (to determine the flip, applied per chunk before
#' summarising, matching the in-memory path's raw-draws flip rather than
#' a post-hoc adjustment of already-computed quantiles); the flip only
#' happens for `stan_model`s with a reflection symmetry, same as the
#' in-memory path. This mode only supports the default `probs`, since it
#' reuses [diagnose_convergence()]'s chunked-read helper (which always
#' computes a fixed `q5`/`median`/`q95` summary, not user-configurable
#' quantiles) rather than re-reading data already summarised once; pass
#' an in-memory `fit` if you need other quantiles.
#'
#' @param fit A `CmdStanMCMC` fit object from [fit_dyad_ts()] or
#'   [fit_panel()], or a character vector of raw CmdStan CSV file paths
#'   (see Details).
#' @param stan_data The Stan data list used to produce `fit`, as returned
#'   by [assemble_stan_data()] (must still carry its `dyad_ids`
#'   attribute).
#' @param probs Posterior quantiles to report alongside the mean. Only
#'   the default is supported when `fit` is CSV file paths (see Details).
#' @param stan_model Name registered in `.bilatr_stan_models` identifying
#'   which model produced `fit`; see [bilatr_orient()]. Defaults to
#'   `.BILATR_DEFAULT_MODEL` (`"stable"`), matching what [fit_dyad_ts()]/
#'   [fit_panel()] always fit.
#' @param max_memory_mb,chunk_size,parallel,n_workers,scratch_dir Only
#'   used when `fit` is CSV file paths; identical in meaning to
#'   [diagnose_convergence()]'s arguments of the same name (including
#'   the one-time default-`max_memory_mb` `message()`), applied here to
#'   `theta` alone rather than all of Tier 3.
#' @return A tibble with one row per dyad-period: `dyad_id`,
#'   `time_index`, `dyad`, `dyad2`, `year` (and `month`, if applicable),
#'   the posterior `mean` of theta, and one column per requested quantile.
#' @examples
#' \dontrun{
#' theta <- extract_theta(fit, stan_data)
#'
#' # a completed SLURM run, never read into this R session
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' theta <- extract_theta(csv_files, stan_data, stan_model = "alphanorm", max_memory_mb = 16384)
#' }
#' @export
extract_theta <- function(
  fit, stan_data, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL,
  max_memory_mb = 8192, chunk_size = NULL, parallel = FALSE,
  n_workers = max(1L, parallel::detectCores() - 1L), scratch_dir = NULL
) {
  max_memory_mb_missing <- missing(max_memory_mb)

  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be the output of assemble_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }

  if (is.character(fit)) {
    if (!identical(probs, c(0.05, 0.5, 0.95))) {
      stop(
        "extract_theta()'s CSV-file-path mode only supports the default ",
        "`probs` (0.05, 0.5, 0.95): it reuses diagnose_convergence()'s ",
        "chunked-read helper, which always computes the fixed q5/median/",
        "q95 summary rather than user-configurable quantiles, to avoid ",
        "re-reading data already summarised once. Pass an in-memory ",
        "`fit` object if you need other quantiles.",
        call. = FALSE
      )
    }

    csv_files <- fit
    prepared <- .prepare_fast_csv_read(csv_files, scratch_dir)
    on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)

    all_vars <- .stan_csv_variable_names(csv_files[1])
    var_tiers <- .classify_bilatr_tier(all_vars)
    theta_vars <- var_tiers$variable[var_tiers$tier == 3L & startsWith(var_tiers$variable, "theta[")]

    flip <- FALSE
    if (length(.bilatr_flip_variables(stan_model)) > 0) {
      alpha1_draws <- .fast_read_post_warmup_draws(prepared, "alpha[1]")
      flip <- stats::median(posterior::extract_variable(alpha1_draws, "alpha[1]")) < 0
    }

    chunk_size_used <- .resolve_chunk_size_and_report(
      length(theta_vars), csv_files, max_memory_mb, chunk_size, parallel, n_workers,
      max_memory_mb_missing
    )
    theta_summ <- .chunked_summarise_csv(
      prepared, theta_vars, chunk_size_used, parallel, n_workers, flip = flip
    ) %>%
      dplyr::select(variable, mean, `5%` = q5, `50%` = median, `95%` = q95)
  } else {
    draws <- fit$draws(variables = c("alpha[1]", "theta"))
    draws <- bilatr_orient(draws, stan_model = stan_model, variables = "theta")
    theta_summ <- posterior::summarise_draws(
      draws,
      mean = mean,
      ~ stats::quantile(.x, probs = probs)
    )
  }

  theta_summ %>%
    dplyr::mutate(variable = stringr::str_remove_all(variable, "theta\\[|\\]")) %>%
    tidyr::separate(variable, into = c("dyad_id", "time_index"), sep = ",", convert = TRUE) %>%
    dplyr::left_join(dyad_ids, by = c("dyad_id", "time_index"))
}

#' Read named variables from either an in-memory fit or raw CmdStan CSVs
#'
#' Shared one-line dispatch behind [extract_alpha()]/
#' [extract_mu_intercept()]'s `fit`-as-CSV-paths support: unlike
#' [extract_theta()]/[diagnose_convergence()]'s CSV-path branches, these
#' two already know exactly which (small) set of variable names they
#' want, so there is no tier enumeration or chunking to do -- just one
#' plain, unchunked read across all chain files, via
#' [.prepare_fast_csv_read]/[.fast_read_post_warmup_draws] rather than
#' [cmdstanr::read_cmdstan_csv()] (see [.prepare_fast_csv_read]'s docs
#' for why: the latter's `cmd = grep` read re-materializes a full
#' near-complete copy of each chain's raw CSV on every call, however few
#' columns are requested). Not free in absolute terms, though: reading
#' still scans each file's full row width regardless of how few columns
#' are requested (see [.chunked_summarise_csv]'s docs for the benchmark
#' this is based on), so this is the minimum possible number of file
#' scans (one), not an instant lookup -- there is no way to make reading
#' a handful of columns out of a very wide CSV cheaper than that at the
#' file-format level.
#'
#' @param fit A `CmdStanMCMC`-like fit object, or a character vector of
#'   CmdStan CSV file paths.
#' @param variables Character vector of variable names to read.
#' @param scratch_dir See [diagnose_convergence()]. Only used when `fit`
#'   is CSV file paths.
#' @return A `posterior::draws_array`.
#' @keywords internal
.get_draws <- function(fit, variables, scratch_dir = NULL) {
  if (is.character(fit)) {
    prepared <- .prepare_fast_csv_read(fit, scratch_dir)
    on.exit(.cleanup_fast_csv_read(prepared), add = TRUE)
    .fast_read_post_warmup_draws(prepared, variables)
  } else {
    fit$draws(variables = variables)
  }
}

#' Extract discrimination parameters (alpha) with action-class labels
#'
#' Pulls posterior summaries of the action-type discrimination parameters
#' `alpha` out of a fitted model. `alpha[1]` is fixed at 1 (the model's
#' scale/sign reference) for `stable`/`ou`; every other element is freely
#' estimated. See the package's identification notes in
#' `vignette("dyad_time_series")`.
#'
#' Draws are passed through [bilatr_orient()] before summarizing, so
#' `alpha`'s sign always reflects the canonical orientation regardless of
#' which reflection-symmetry basin the sampler actually landed in
#' (`stan_model = "alphanorm"`/`"alphanorm_ou"` only; a no-op otherwise).
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths
#' (one per chain), for the case where there is no in-memory fit at all
#' -- e.g. a post-hoc, cross-chain roll-up script run after several
#' independently-submitted single-chain SLURM jobs have all completed,
#' with only their saved CSVs on disk. Unlike [extract_theta()]'s
#' CSV-path mode, `alpha` is small (Tier 1: a handful of values
#' regardless of dyad-set size), so there is no chunking/`parallel`
#' machinery here -- just one plain read (via [.get_draws]) across all
#' of `fit`, then the same extraction logic either way.
#'
#' @inheritParams extract_theta
#' @param event_classes Optional character vector of event-class labels,
#'   in the same order used to build `stan_data` (i.e. the
#'   `"event_classes"` attribute attached by [assemble_stan_data()]). If
#'   supplied, an `event_class` column is added alongside the raw action
#'   index.
#' @param scratch_dir See [diagnose_convergence()]. Only used when `fit`
#'   is CSV file paths.
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of alpha, and one column per requested quantile.
#' @examples
#' \dontrun{
#' alpha <- extract_alpha(fit, event_classes = attr(stan_data, "event_classes"))
#'
#' # cross-chain roll-up from saved CSVs, no in-memory fit
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' alpha <- extract_alpha(csv_files, event_classes = event_classes, stan_model = "alphanorm")
#' }
#' @export
extract_alpha <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL, scratch_dir = NULL) {
  draws <- bilatr_orient(.get_draws(fit, "alpha", scratch_dir), stan_model = stan_model, variables = "alpha")

  out <- posterior::summarise_draws(
    draws,
    mean = mean,
    ~ stats::quantile(.x, probs = probs)
  ) %>%
    dplyr::mutate(
      action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])"))
    )

  if (!is.null(event_classes)) {
    out <- dplyr::mutate(out, event_class = event_classes[action_index])
  }
  out
}

#' Extract global action-type intercepts (mu_intercept) with action-class
#' labels
#'
#' Pulls posterior summaries of the global action-type intercepts
#' `mu_intercept` out of a fitted model. `mu_intercept[1]` is fixed at 0
#' (the model's softmax level-shift reference) for `stable`/`ou`.
#'
#' `mu_intercept` is unaffected by alphanorm/alphanorm_ou's reflection
#' symmetry (`alpha .* theta` is invariant under the joint negation, so
#' `mu_intercept` never needs to flip; see [bilatr_orient()]) -- `fit` is
#' passed through it anyway, for a single consistent code path across the
#' `extract_*()` functions, but it is always a no-op here.
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths,
#' for the same post-hoc/no-in-memory-fit case described in
#' [extract_alpha()]'s documentation (also Tier 1/2: small regardless of
#' dyad-set size, so no chunking here either).
#'
#' @inheritParams extract_alpha
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of mu_intercept, and one column per requested quantile.
#' @examples
#' \dontrun{
#' mu_intercept <- extract_mu_intercept(fit, event_classes = attr(stan_data, "event_classes"))
#' }
#' @export
extract_mu_intercept <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL, scratch_dir = NULL) {
  draws <- .get_draws(fit, c("alpha[1]", "mu_intercept"), scratch_dir)
  draws <- bilatr_orient(draws, stan_model = stan_model, variables = "mu_intercept")

  out <- posterior::summarise_draws(
    draws,
    mean = mean,
    ~ stats::quantile(.x, probs = probs)
  ) %>%
    dplyr::mutate(
      action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])"))
    )

  if (!is.null(event_classes)) {
    out <- dplyr::mutate(out, event_class = event_classes[action_index])
  }
  out
}
