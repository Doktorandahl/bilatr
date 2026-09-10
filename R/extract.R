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
#' sampler actually landed in. This matters for the default `stan_model`
#' (`"stable"`) itself, not just the experimental `"ou"`: both normalize
#' `alpha` via a `sum_to_zero_vector` with only a *soft* sign anchor on
#' `alpha[1]` (see [assemble_stan_data()]'s `anchor_scale`), so a chain's
#' init can still land in the wrong-sign basin; this reorientation is
#' what makes the result correct regardless. A no-op only for a
#' hypothetical `stan_model` whose identification hard-fixes `alpha[1]`'s
#' sign instead (as `stable`/`ou` themselves did before 0.4.0's
#' promotion; see NEWS.md).
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths
#' (one per chain), matching [diagnose_convergence()]'s CSV-path mode --
#' reads and summarises `theta` directly from the raw CSVs in
#' memory-bounded chunks via [.chunked_summarise_csv_with_orientation]
#' rather than materializing the full draws array, for the same
#' production-scale-panel reason `diagnose_convergence()` needed it.
#' **This trades memory for wall-time, and the exchange rate can be
#' steep**: reading is single-threaded, and each chunk's cost is roughly
#' a full parse of the largest chain file regardless of how many columns
#' are requested (see `max_memory_mb` below and in
#' [diagnose_convergence()]'s documentation), so total sweep time scales
#' with chunk COUNT, not memory saved. Prefer the in-memory path when
#' you already have `fit` in memory (no benefit to re-reading), and the
#' CSV path's default `max_memory_mb` as large as your job's memory
#' allocation can afford, not as small as "safely" possible. Sign
#' orientation in this mode is decided from the first chunk itself
#' (`alpha[1]` prepended, read once, dropped before summarising --
#' see [.chunked_summarise_csv_with_orientation]) rather than a separate
#' pass, applied per chunk before summarising, matching the in-memory
#' path's raw-draws flip rather than a post-hoc adjustment of
#' already-computed quantiles; the flip only happens for `stan_model`s
#' with a reflection symmetry, same as the in-memory path. This mode
#' only supports the default `probs`, since it reuses
#' [diagnose_convergence()]'s chunked-read helper (which always computes
#' a fixed `q5`/`median`/`q95` summary, not user-configurable quantiles)
#' rather than re-reading data already summarised once; pass an
#' in-memory `fit` if you need other quantiles.
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
#'   which model produced `fit`, or a recognized pre-0.4.0 alias
#'   (`"alphanorm"`/`"alphanorm_ou"`, mapped to `"stable"`/`"ou"` with a
#'   message; see [.canonical_stan_model()]) -- an unrecognized name
#'   errors immediately rather than silently skipping sign orientation.
#'   See [bilatr_orient()]. Defaults to `.BILATR_DEFAULT_MODEL`
#'   (`"stable"`), matching what [fit_dyad_ts()]/[fit_panel()] always fit.
#' @param max_memory_mb,chunk_size,parallel,n_workers Only used when
#'   `fit` is CSV file paths; identical in meaning to
#'   [diagnose_convergence()]'s arguments of the same name (including
#'   the one-time default-`max_memory_mb` `message()`), applied here to
#'   `theta` alone rather than all of Tier 3.
#' @param scratch_dir Deprecated and ignored since 0.4.1; see
#'   [diagnose_convergence()].
#' @return A tibble with one row per dyad-period: `dyad_id`,
#'   `time_index`, `dyad`, `dyad2`, `year` (and `month`, if applicable),
#'   the posterior `mean` of theta, and one column per requested quantile.
#' @examples
#' \dontrun{
#' theta <- extract_theta(fit, stan_data)
#'
#' # a completed SLURM run, never read into this R session
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' theta <- extract_theta(csv_files, stan_data, stan_model = "ou", max_memory_mb = 16384)
#' }
#' @export
extract_theta <- function(
  fit, stan_data, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL,
  max_memory_mb = 8192, chunk_size = NULL, parallel = FALSE,
  n_workers = parallelly::availableCores(), scratch_dir = NULL
) {
  stan_model <- .canonical_stan_model(stan_model)
  max_memory_mb_missing <- missing(max_memory_mb)
  if (!is.null(scratch_dir)) {
    warning(
      "`scratch_dir` is deprecated and ignored since 0.4.1: no scratch ",
      "copy is made any more.",
      call. = FALSE
    )
  }

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
    prepared <- .prepare_fast_csv_read(csv_files)

    var_tiers <- .classify_bilatr_tier(prepared$variables)
    theta_vars <- var_tiers$variable[var_tiers$tier == 3L & startsWith(var_tiers$variable, "theta[")]

    n_cores <- if (parallel) n_workers else 1L
    chunk_size_used <- .resolve_chunk_size_and_report(
      length(theta_vars), prepared, max_memory_mb, chunk_size, n_cores,
      max_memory_mb_missing
    )
    # No standalone alpha[1] pass: orientation is decided from the first
    # Tier 3 chunk itself (alpha[1] prepended, read once, dropped before
    # summarising) -- see .chunked_summarise_csv_with_orientation().
    theta_summ <- .chunked_summarise_csv_with_orientation(
      prepared, theta_vars, chunk_size_used, n_cores,
      .bilatr_flip_variables(stan_model)
    )$summ %>%
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
#' [.prepare_fast_csv_read]/[.fast_read_post_warmup_draws]. Not free in
#' absolute terms, though: `data.table::fread()` touches each file's
#' full byte range once per call regardless of how few columns are
#' requested (see [.fast_read_post_warmup_draws]'s docs), so this is the
#' minimum possible number of file touches (one), not an instant lookup.
#'
#' @param fit A `CmdStanMCMC`-like fit object, or a character vector of
#'   CmdStan CSV file paths.
#' @param variables Character vector of variable names to read.
#' @return A `posterior::draws_array`.
#' @keywords internal
.get_draws <- function(fit, variables) {
  if (is.character(fit)) {
    prepared <- .prepare_fast_csv_read(fit)
    .fast_read_post_warmup_draws(prepared, variables)
  } else {
    fit$draws(variables = variables)
  }
}

#' Extract discrimination parameters (alpha) with action-class labels
#'
#' Pulls posterior summaries of the action-type discrimination parameters
#' `alpha` out of a fitted model. `alpha` sums to exactly 0 and has RMS 1
#' (a `sum_to_zero_vector`, not a single fixed-to-1 element) for
#' `stable`/`ou`; `alpha[1]`, the reference/neutral action class supplied
#' via `reference_category`, is only softly anchored positive (see
#' [assemble_stan_data()]'s `anchor_scale`). See the package's
#' identification notes in `vignette("dyad_time_series")`.
#'
#' Draws are passed through [bilatr_orient()] before summarizing, so
#' `alpha`'s sign always reflects the canonical orientation regardless of
#' which reflection-symmetry basin the sampler actually landed in --
#' `stable`/`ou`'s soft anchor above doesn't guarantee that on its own.
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
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of alpha, and one column per requested quantile.
#' @examples
#' \dontrun{
#' alpha <- extract_alpha(fit, event_classes = attr(stan_data, "event_classes"))
#'
#' # cross-chain roll-up from saved CSVs, no in-memory fit
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' alpha <- extract_alpha(csv_files, event_classes = event_classes, stan_model = "ou")
#' }
#' @export
extract_alpha <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL, scratch_dir = NULL) {
  stan_model <- .canonical_stan_model(stan_model)
  if (!is.null(scratch_dir)) {
    warning(
      "`scratch_dir` is deprecated and ignored since 0.4.1: no scratch ",
      "copy is made any more.",
      call. = FALSE
    )
  }
  draws <- bilatr_orient(.get_draws(fit, "alpha"), stan_model = stan_model, variables = "alpha")

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
#' `mu_intercept` out of a fitted model. `mu_intercept` sums to exactly 0
#' (a `sum_to_zero_vector`, not a single fixed-to-0 element) for
#' `stable`/`ou`, so no residual location degree of freedom hides in the
#' softmax level-shift.
#'
#' `mu_intercept` is unaffected by stable/ou's reflection symmetry
#' (`alpha .* theta` is invariant under the joint negation, so
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
  stan_model <- .canonical_stan_model(stan_model)
  if (!is.null(scratch_dir)) {
    warning(
      "`scratch_dir` is deprecated and ignored since 0.4.1: no scratch ",
      "copy is made any more.",
      call. = FALSE
    )
  }
  draws <- .get_draws(fit, c("alpha[1]", "mu_intercept"))
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
