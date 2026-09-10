#' Read a variable set in one pass, determining sign orientation from
#' `alpha[1]` if it's needed but wasn't already among `vars`
#'
#' Shared by [diagnose_and_extract_bilatr()]'s Tier 1/2-only and
#' Tier-1/2-plus-Tier-3-folded-together cases (see that function's
#' Details): reads `vars`, prepending `"alpha[1]"` first if `flip_vars`
#' is non-empty and `vars` doesn't already include it (tier 1 wasn't
#' requested), determines `flip` from its median sign, negates whichever
#' of `vars`' own columns [.bilatr_match_draws_columns()] matches
#' `flip_vars`, summarises, and drops the injected `alpha[1]` row
#' afterward so it never surfaces as a silently-incomplete alpha output.
#'
#' @param prepared Output of [.prepare_fast_csv_read].
#' @param vars Character vector of variables to read (Tier 1/2, or
#'   Tier 1/2 + Tier 3 together).
#' @param flip_vars Output of [.bilatr_flip_variables()] for the model
#'   being read.
#' @param n_cores See [.chunked_summarise_csv()].
#' @return A list with `summ` (summary tibble) and `flip` (logical).
#' @keywords internal
.read_and_orient_draws_summary <- function(prepared, vars, flip_vars, n_cores = 1L) {
  alpha1_injected <- length(flip_vars) > 0 && !("alpha[1]" %in% vars)
  read_vars <- if (alpha1_injected) c("alpha[1]", vars) else vars
  draws <- .fast_read_post_warmup_draws(prepared, read_vars)

  flip <- FALSE
  if (length(flip_vars) > 0) {
    flip <- stats::median(posterior::extract_variable(draws, "alpha[1]")) < 0
    if (flip) {
      flip_cols <- .bilatr_match_draws_columns(posterior::variables(draws), flip_vars)
      if (length(flip_cols) > 0) {
        draws[, , flip_cols] <- -draws[, , flip_cols, drop = FALSE]
      }
    }
  }

  summ <- .summarise_bilatr_draws(draws, n_cores = n_cores)
  if (alpha1_injected) {
    summ <- dplyr::filter(summ, variable != "alpha[1]")
  }
  list(summ = summ, flip = flip)
}

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
#' once (chunked, same as [diagnose_convergence()]), reorients whichever
#' columns [bilatr_orient()] would flip for `stan_model` as part of that
#' single read, and derives all four outputs from those two summaries
#' rather than re-reading anything.
#'
#' Total file touches per chain file, at production scale (Tier 3 too
#' large for one chunk): one `cmdstanr:::read_csv_metadata()` call and
#' one row-count validation probe, both in [.prepare_fast_csv_read]; one
#' read for Tier 1/2 (`alpha[1]` included at no extra cost whenever
#' tier 1 itself is requested, which is the default); and one per Tier 3
#' chunk. When Tier 3 fits in a single chunk (small panels, or a large
#' `max_memory_mb`), Tier 1/2 is folded into that same read instead of
#' being a separate one -- see [.read_and_orient_draws_summary()].
#'
#' Only meaningful for the CSV-file-paths case: an in-memory `fit` has no
#' repeated-disk-read cost to avoid, so there is nothing to fuse there --
#' call the four functions separately instead.
#'
#' Rhat/ESS are invariant to a deterministic sign flip, so reorienting
#' before summarising never changes anything [diagnose_convergence()]
#' reports; it only affects the mean/quantile columns, which is exactly
#' what [extract_theta()]/[extract_alpha()]/[extract_mu_intercept()] need
#' reoriented anyway. Tier 3 is NOT always exactly `{theta, theta_raw}`:
#' with `compute_log_lik = 1` it also includes `log_lik[d,t]` (also
#' two-indexed), which must never flip even when `theta`/`theta_raw` in
#' the same chunk do. Each Tier 3 chunk is therefore negated
#' column-selectively, matching exactly the variables
#' [.bilatr_flip_variables()] lists (see [.chunked_summarise_csv()]'s
#' `flip_vars`), not flipped as a whole uniformly.
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
#' @param stan_model Name registered in `.bilatr_stan_models`, or a
#'   recognized pre-0.4.0 alias (`"alphanorm"`/`"alphanorm_ou"`, mapped to
#'   `"stable"`/`"ou"` with a message; see [.canonical_stan_model()]) --
#'   an unrecognized name errors immediately rather than silently
#'   skipping sign orientation. See [bilatr_orient()]. Defaults to
#'   `.BILATR_DEFAULT_MODEL`.
#' @param event_classes Optional character vector of event-class labels,
#'   passed through to the `alpha`/`mu_intercept` outputs exactly as in
#'   [extract_alpha()]. Defaults to `stan_data`'s `"event_classes"`
#'   attribute, if present.
#' @param scratch_dir Deprecated and ignored since 0.4.1; see
#'   [diagnose_convergence()].
#' @return A list with elements `diagnostics` (a `bilatr_diagnostics`
#'   object, as from [diagnose_convergence()]), `theta`, `alpha`, and
#'   `mu_intercept` (tibbles, in the same shape [extract_theta()]/
#'   [extract_alpha()]/[extract_mu_intercept()] return).
#' @examples
#' \dontrun{
#' csv_files <- list.files("model_output/some_spec", pattern = "\\.csv$", full.names = TRUE)
#' result <- diagnose_and_extract_bilatr(
#'   csv_files, stan_data, n_dt = n_dt, stan_model = "ou"
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
  n_workers = parallelly::availableCores(), scratch_dir = NULL
) {
  stan_model <- .canonical_stan_model(stan_model)
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
  if (!is.null(scratch_dir)) {
    warning(
      "`scratch_dir` is deprecated and ignored since 0.4.1: no scratch ",
      "copy is made any more.",
      call. = FALSE
    )
  }

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

  prepared <- .prepare_fast_csv_read(csv_files)

  var_tiers <- .classify_bilatr_tier(prepared$variables)
  keep_tiers <- var_tiers[var_tiers$tier %in% tiers, ]

  tier12_vars <- keep_tiers$variable[keep_tiers$tier %in% c(1L, 2L)]
  tier3_vars <- keep_tiers$variable[keep_tiers$tier == 3L]

  # Orientation is folded into whichever read already includes alpha[1],
  # rather than a separate dedicated pass: if tier 1 is requested,
  # alpha[1] is already part of tier12_vars (alpha is Tier 1); otherwise
  # [.read_and_orient_draws_summary] prepends it to the Tier 1/2 read (if
  # any Tier 2 variables are also being read) or, if none are (tiers = 3
  # alone), [.chunked_summarise_csv_with_orientation] prepends it to the
  # first Tier 3 chunk instead -- the same position extract_theta()'s
  # CSV branch is always in.
  flip_vars <- .bilatr_flip_variables(stan_model)
  flip <- FALSE
  n_cores <- if (parallel) n_workers else 1L

  chunk_size_used <- if (length(tier3_vars) > 0) {
    .resolve_chunk_size_and_report(
      length(tier3_vars), prepared, max_memory_mb, chunk_size, n_cores,
      max_memory_mb_missing
    )
  } else {
    NULL
  }
  # Whether the whole Tier 3 sweep fits in a single chunk anyway -- if
  # so, and Tier 1/2 also has variables to read, fold them into ONE read
  # together rather than two (a separate Tier 1/2 pass plus a
  # single-chunk Tier 3 "sweep" would otherwise touch every chain file
  # twice for no reason).
  tier3_fits_one_chunk <- length(tier3_vars) > 0 && length(tier3_vars) <= chunk_size_used

  tier12_summ <- NULL
  tier3_summ <- NULL

  if (length(tier12_vars) > 0 && tier3_fits_one_chunk) {
    combined <- .read_and_orient_draws_summary(prepared, c(tier12_vars, tier3_vars), flip_vars, n_cores)
    # combined result carries both tiers' rows; tier3_summ stays NULL so
    # dplyr::bind_rows(tier12_summ, tier3_summ) below isn't duplicated
    tier12_summ <- combined$summ
    flip <- combined$flip
  } else {
    # Tier 1/2: rhat/ess are unaffected by the flip, so this single
    # read/summary serves both diagnose_convergence()'s Tier 1/2 tables
    # below and extract_alpha()/extract_mu_intercept()'s already-
    # correctly-oriented mean/quantile columns, with no second read.
    if (length(tier12_vars) > 0) {
      tier12 <- .read_and_orient_draws_summary(prepared, tier12_vars, flip_vars, n_cores)
      tier12_summ <- tier12$summ
      flip <- tier12$flip
    }

    # Tier 3 (theta/theta_raw, and possibly log_lik if compute_log_lik = 1)
    # -- only the columns .bilatr_flip_variables() actually lists get
    # negated within each chunk (B5), never the whole chunk uniformly,
    # since a chunk can mix theta/theta_raw (which must flip) with
    # log_lik[d,t] (also Tier 3, which must NOT). Serves both
    # diagnose_convergence()'s Tier 3 table and extract_theta()'s
    # already-correctly-oriented "theta[" subset, with no second sweep
    # over the same columns.
    if (length(tier3_vars) > 0) {
      if (length(tier12_vars) > 0) {
        # flip already determined above from the Tier 1/2 read
        tier3_summ <- .chunked_summarise_csv(
          prepared, tier3_vars, chunk_size_used, n_cores,
          flip_vars = if (flip) flip_vars else character(0)
        )
      } else {
        oriented <- .chunked_summarise_csv_with_orientation(
          prepared, tier3_vars, chunk_size_used, n_cores, flip_vars
        )
        tier3_summ <- oriented$summ
        flip <- oriented$flip
      }
    }
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
