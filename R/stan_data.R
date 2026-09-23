#' Assemble Stan-ready data for the bilatr dyadic IRT model
#'
#' Aggregates CAMEO-coded event data to dyad-period action-class counts
#' and packages it as the data list expected by the package's Stan model
#' (`inst/stan/bilatr_dirmult_irt.stan`).
#'
#' @param data A data frame of event-level records, as produced by
#'   [extract_all_relevant_gdelt()] or [ingest_icews()] and recoded via
#'   [recode_cameo()].
#' @param years Integer vector of years to cover.
#' @param resolution Either `"monthly"` or `"yearly"`.
#' @param grouping_var Name of the event-class column to aggregate on
#'   (e.g. `"QuadClass"`, `"PentaClass"`).
#' @param directed If `TRUE` (default), dyads are directed; if `FALSE`,
#'   actor order is ignored.
#' @param reference_category Value of `grouping_var` to anchor as the
#'   model's scale/sign reference: `stable`/`ou` fold `alpha[1]`'s sign
#'   into the reported `alpha`/`theta` (see each `.stan` file's header,
#'   "IDENTIFICATION: ORIENTATION FOLD"), and RMS-normalize the whole
#'   `alpha` vector to 1, so positive `alpha[1]` means better relations at
#'   this reference/neutral class. Should typically be a
#'   low-conflict/cooperative class. Every other action class's
#'   discrimination (`alpha[2:A]`) is freely estimated.
#' @param min_n_events Minimum total events for a dyad to be retained.
#' @param weighted Defunct. The `dyad_weight`/`period_weight`/
#'   `action_weight` likelihood-weighting scheme was removed in 0.4.6 (see
#'   NEWS.md) -- never used in production, and `action_weight` in
#'   particular made the Dirichlet-multinomial concentration depend on
#'   `theta`, complicating the hand-differentiated forward filter built on
#'   top of this likelihood. Kept only as a formal so old call sites that
#'   actually requested weighting fail loudly instead of silently fitting
#'   an unweighted model: must be `FALSE` or `"none"` (both accepted as
#'   the only meaningful values now, matching this project's own SLURM
#'   runscripts, which map their own `weighted = "none"` CLI argument to
#'   `FALSE` before calling this function); any other value errors.
#' @param chunk_size `reduce_sum` grainsize used to chunk the likelihood
#'   across dyads. 16 cores with `chunk_size = 600` was the
#'   Pareto-optimal setting found in this project's own benchmarking for
#'   panel-sized data; smaller panels or single-dyad fits should use a
#'   smaller value (down to 1) since chunking overhead outweighs the
#'   benefit below that scale. See [fit_panel()] for the corresponding
#'   `threads_per_chain` argument.
#' @param rho_prior_a,rho_prior_b Shape parameters of the
#'   `beta(rho_prior_a, rho_prior_b)` prior on the OU/AR(1) persistence
#'   parameter `rho`, consumed only by the experimental `ou` Stan model
#'   variant (see `R/model_registry.R`); ignored by `stable`, which keeps
#'   a random-walk `theta` with no persistence parameter. Default `8, 2`
#'   (weighted toward strong persistence).
#' @param compute_log_lik `0` (default) or `1`. Data flag consumed by
#'   both registered Stan model variants, gating a per-dyad-period
#'   `log_lik` in `generated quantities`; left off by default since it is
#'   `D x T x` draws and file size already scales with dyad count.
#' @param prior_only `0` (default) or `1`. Data flag consumed by both
#'   registered Stan model variants, gating the `reduce_sum` likelihood
#'   call in `model` exactly as `compute_log_lik` gates its generated
#'   quantity -- `1` fits the prior alone, useful for prior-predictive
#'   checks and for validating [alpha_prior_moments()] against an actual
#'   fit.
#' @param compute_theta_filtered `0` (default) or `1`. Data flag consumed
#'   by both registered Stan model variants, gating
#'   `theta_filtered`/`theta_filtered_sd` in `generated quantities`: the
#'   forward-filtered state (conditional on each draw's hyperparameters,
#'   a West-Harrison/Fisher-scoring linear-Bayes update; see each `.stan`
#'   file's `generated quantities` block) as opposed to the smoothed
#'   `theta` every fit already returns. Left off by default: computing it
#'   for every dyad is `D * T` iterations with `A` `digamma()` calls each,
#'   single-threaded (no `reduce_sum` in `generated quantities`), on the
#'   order of a minute or two per chain for a full production-sized dyad
#'   set, and adds two more `D x T`-scaled Tier 3 blocks to the output.
#' @param filter_dyads Optional character vector of dyad names (matched
#'   against the `dyad_ids` attribute's `dyad` column -- the directed
#'   `"AAA_BBB"` key, same as everywhere else in this package) naming
#'   which dyads to compute `theta_filtered` for. Ignored if
#'   `compute_theta_filtered` is `0`. `NULL` (default) means all dyads,
#'   once `compute_theta_filtered = 1` is set -- pass a small subset (e.g.
#'   fifty dyads) to compute the filter cheaply for just those. Errors,
#'   listing the offending names, if any requested dyad was dropped by
#'   `min_n_events` or never existed.
#' @param anchor_scale Scale of the soft sign anchor
#'   `target += log_inv_logit(alpha[1] * inv(anchor_scale))`. Since 0.4.2
#'   (see NEWS.md), this is consumed only by the LEGACY
#'   `stable_soft_anchor`/`ou_soft_anchor` Stan programs (see
#'   `R/model_registry.R`), which normalize `alpha` via a free
#'   `sum_to_zero_vector` with no fixed element -- leaving an exact
#'   reflection symmetry (`alpha`, `theta` -> `-alpha`, `-theta` is
#'   likelihood-invariant) the anchor only softly penalizes (`alpha[1] <
#'   0`), not reliably breaks across independently-initialized chains;
#'   see `bilatr_orient()` for the post-hoc fix those two programs still
#'   need. The current `stable`/`ou` programs fold `alpha[1]`'s sign into
#'   the reported `alpha`/`theta` instead (see each `.stan` file's
#'   header, "IDENTIFICATION: ORIENTATION FOLD") and don't declare
#'   `anchor_scale` at all; it is still supplied unconditionally
#'   here (CmdStan ignores data a program doesn't declare), so this
#'   function doesn't need to branch on `stan_model`. `alpha[1]` is
#'   already the reference/neutral action class supplied via
#'   `reference_category` -- no separate index is needed. Default `0.1`.
#' @return A named list suitable as the `data` argument to
#'   `cmdstanr::CmdStanModel$sample()` for the bilatr Stan model: `D`,
#'   `T`, `A`, `C`, `is_obs`, `Y`, `rho_prior_a`, `rho_prior_b`,
#'   `compute_log_lik`, `prior_only`, `compute_theta_filtered`,
#'   `n_filter_dyads`, `filter_dyads`, `anchor_scale`, and (0.7.0,
#'   unconditionally -- see the fields' own inline comments in this
#'   function's body) `n_countries`, `ctry_a`, `ctry_b`, `w_send` (only
#'   consumed by the experimental `stable_gamma` variant; see
#'   `R/model_registry.R`). `n_countries` counts only the countries
#'   actually present in the RETAINED dyads (after `min_n_events`
#'   filtering) -- a country appearing only in a dropped dyad consumes no
#'   index. `ctry_a`/`ctry_b` are 1-indexed into that country set, in
#'   `dyad_id` order (`ctry_a` is the sender for directed dyads, side A
#'   for undirected); `w_send` is per-dyad (not per dyad-period -- see
#'   [grouped_events_to_dyad_period()]'s `w_send` attribute), exactly `1`
#'   for every directed dyad. A `message()` at assembly time reports
#'   `n_countries` and the dyads-per-country distribution (min/median/
#'   max), so a country identified from only one or two dyads (whose
#'   `gamma_c` is then close to confounded with those dyads' own `theta`)
#'   is visible before fitting, not after. Also carries a `dyad_ids`
#'   attribute (the output of [make_dyad_ids()]) for reattaching
#'   identifiers to posterior draws (see [extract_theta()]), an
#'   `event_classes` attribute, and (0.7.0) a `country_codes` attribute
#'   (character, in `ctry_a`/`ctry_b` index order) for labelling `gamma`
#'   (see [extract_gamma()]).
#' @examples
#' \dontrun{
#' events <- extract_all_relevant_gdelt("data/gdelt_raw/20200101.zip")
#' events <- recode_cameo(events)
#' agg <- grouped_events_to_dyad_period(
#'   events,
#'   resolution = "yearly",
#'   grouping_var = "PentaClass",
#'   reference_category = 0
#' )
#' stan_data <- assemble_stan_data(
#'   events,
#'   years = 2015:2020,
#'   resolution = "yearly",
#'   grouping_var = "PentaClass",
#'   reference_category = 0
#' )
#' }
#' @export
assemble_stan_data <- function(
  data,
  years,
  resolution = c("monthly", "yearly"),
  grouping_var,
  directed = TRUE,
  reference_category = NULL,
  min_n_events = 1,
  weighted = FALSE,
  chunk_size = 100,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  anchor_scale = 0.1
) {
  resolution <- match.arg(resolution)

  if (!isFALSE(weighted) && !identical(weighted, "none")) {
    stop(
      "likelihood weighting was removed in 0.4.6; see NEWS. ",
      "`weighted` must be `FALSE` or \"none\".",
      call. = FALSE
    )
  }

  agg <- grouped_events_to_dyad_period(
    data,
    resolution = resolution,
    grouping_var = grouping_var,
    directed = directed,
    reference_category = reference_category
  )
  # Captured here, immediately, rather than relied on to survive the
  # dplyr pipeline below (fill_dyad_period_skeleton()'s right_join/
  # arrange/mutate chain, then the min_n_events dplyr::filter()): custom
  # attributes are not guaranteed to pass through dplyr verbs, so this is
  # kept as a plain local tibble and joined back onto the SURVIVING
  # `dyads` explicitly below (see @return / Part 2a of the build
  # prompt) -- achieving the same "dropped dyads drop their weights with
  # them" outcome without depending on attribute pass-through.
  w_send_tbl <- attr(agg, "w_send")
  agg <- fill_dyad_period_skeleton(agg, years, resolution)

  drop_dyads <- agg %>%
    dplyr::group_by(dyad) %>%
    dplyr::summarise(total_events = sum(total_events), .groups = "drop") %>%
    dplyr::filter(total_events < min_n_events) %>%
    dplyr::pull(dyad)
  agg <- dplyr::filter(agg, !(dyad %in% drop_dyads))

  dyads <- unique(agg$dyad)
  D <- length(dyads)
  if (D == 0) {
    stop(
      "No dyads have at least `min_n_events` (", min_n_events, ") total events; ",
      "nothing to assemble. Lower min_n_events or check the input data.",
      call. = FALSE
    )
  }
  Tn <- nrow(agg) / D

  # Country index and per-dyad ctry_a/ctry_b/w_send (0.7.0; consumed only
  # by the experimental stable_gamma Stan variant, but built and attached
  # UNCONDITIONALLY for every model -- see @param n_countries below for
  # why). Row i of `country_info` corresponds to `dyads[i]` (matched
  # explicitly via `match()`, not assumed pre-sorted the same way) --
  # `dyads`' own D-indexing is exactly the `dyad_id` order
  # [make_dyad_ids()] assigns below, so ctry_a/ctry_b/w_send end up
  # indexed 1:D the same way Y/is_obs/dyad_ids already are.
  country_info <- w_send_tbl[match(dyads, w_send_tbl$dyad), ]
  countries_present <- sort(unique(c(country_info$ctry_a_code, country_info$ctry_b_code)))
  n_countries <- length(countries_present)
  country_index <- stats::setNames(seq_len(n_countries), countries_present)
  ctry_a <- unname(country_index[country_info$ctry_a_code])
  ctry_b <- unname(country_index[country_info$ctry_b_code])
  w_send <- country_info$w_send

  dyads_per_country <- table(c(ctry_a, ctry_b))
  message(sprintf(
    "assemble_stan_data(): n_countries = %d (dyads per country -- min %d, median %s, max %d)",
    n_countries, min(dyads_per_country),
    format(stats::median(as.numeric(dyads_per_country)), nsmall = 1),
    max(dyads_per_country)
  ))

  event_classes <- stringr::str_remove(
    grep("^EventClass_", names(agg), value = TRUE), "^EventClass_"
  )
  Anum <- length(event_classes)

  obs_matrix <- matrix(agg$is_obs, nrow = D, byrow = TRUE)

  events_list <- agg %>%
    dplyr::group_by(dyad) %>%
    dplyr::group_split() %>%
    purrr::map(~ unname(as.matrix(dplyr::select(.x, dplyr::starts_with("EventClass_")))))
  events_array <- aperm(simplify2array(events_list), c(3, 1, 2))

  dyad_ids <- make_dyad_ids(
    agg,
    years = years,
    resolution = resolution,
    min_n_events = min_n_events
  )

  if (compute_theta_filtered == 1) {
    dyad_lookup <- dplyr::distinct(dyad_ids, dyad_id, dyad)
    if (is.null(filter_dyads)) {
      filter_dyad_ids <- dyad_lookup$dyad_id
    } else {
      filter_dyad_ids <- dyad_lookup$dyad_id[match(filter_dyads, dyad_lookup$dyad)]
      missing_dyads <- filter_dyads[is.na(filter_dyad_ids)]
      if (length(missing_dyads) > 0) {
        stop(
          "`filter_dyads` names not found (dropped by min_n_events, or never existed): ",
          paste(missing_dyads, collapse = ", "),
          call. = FALSE
        )
      }
    }
  } else {
    filter_dyad_ids <- integer(0)
  }

  stan_data <- list(
    D = D,
    T = Tn,
    A = Anum,
    C = chunk_size,
    is_obs = obs_matrix,
    Y = events_array,
    # Consumed only by the experimental ou variant (R/model_registry.R);
    # unused by stable, whose Stan program doesn't declare it -- CmdStan
    # ignores data fields a program doesn't declare, so this is passed
    # unconditionally rather than branching on stan_model.
    rho_prior_a = rho_prior_a,
    rho_prior_b = rho_prior_b,
    # Consumed by both registered variants; gates a per-dyad-period
    # `log_lik` in generated quantities (D x T x draws, so off by
    # default).
    compute_log_lik = compute_log_lik,
    # Consumed by both registered variants; gates the reduce_sum
    # likelihood call in `model` (off by default -- fits the prior alone
    # when on).
    prior_only = prior_only,
    # Consumed by both registered variants; gates theta_filtered/
    # theta_filtered_sd in generated quantities for filter_dyad_ids (off
    # by default -- see @param compute_theta_filtered/filter_dyads).
    compute_theta_filtered = compute_theta_filtered,
    n_filter_dyads = length(filter_dyad_ids),
    filter_dyads = filter_dyad_ids,
    # Only consumed by the legacy stable_soft_anchor/ou_soft_anchor
    # programs (soft sign anchor on alpha[1]); stable/ou identify
    # alpha[1]'s sign by construction and don't declare this field, but
    # it's passed unconditionally regardless (CmdStan ignores data a
    # program doesn't declare) -- see @param anchor_scale above.
    anchor_scale = anchor_scale,
    # Consumed only by the experimental stable_gamma variant
    # (R/model_registry.R); unused by every other registered program,
    # whose Stan code doesn't declare these fields -- CmdStan ignores
    # data a program doesn't declare, so (matching rho_prior_a/b and
    # anchor_scale above) these are passed unconditionally rather than
    # gated on which program will be fit, so a caller doesn't need to
    # know the model name at assembly time. BREAKING: a `stan_data.rds`
    # saved before 0.7.0 lacks these four fields and cannot be used with
    # `stable_gamma` -- see NEWS.md.
    n_countries = n_countries,
    ctry_a = ctry_a,
    ctry_b = ctry_b,
    w_send = w_send
  )

  attr(stan_data, "dyad_ids") <- dyad_ids
  attr(stan_data, "event_classes") <- event_classes
  attr(stan_data, "country_codes") <- countries_present

  stan_data
}
