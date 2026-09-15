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
#'   `compute_log_lik`, `anchor_scale`. Also carries a `dyad_ids` attribute
#'   (the output of [make_dyad_ids()]) for reattaching identifiers to
#'   posterior draws; see [extract_theta()].
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
    # Only consumed by the legacy stable_soft_anchor/ou_soft_anchor
    # programs (soft sign anchor on alpha[1]); stable/ou identify
    # alpha[1]'s sign by construction and don't declare this field, but
    # it's passed unconditionally regardless (CmdStan ignores data a
    # program doesn't declare) -- see @param anchor_scale above.
    anchor_scale = anchor_scale
  )

  attr(stan_data, "dyad_ids") <- make_dyad_ids(
    agg,
    years = years,
    resolution = resolution,
    min_n_events = min_n_events
  )
  attr(stan_data, "event_classes") <- event_classes

  stan_data
}
