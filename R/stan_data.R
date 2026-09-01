#' Compute inverse-frequency dyad/period/action reweighting vectors
#'
#' Replicates the informal reweighting scheme used in earlier iterations of
#' this model: dyads and periods with very high total event counts are
#' downweighted (via `1 / log1p(total)`, normalized to mean 1), and action
#' types that occur rarely are upweighted the same way, so a handful of
#' high-volume dyad-periods or common action types don't dominate the
#' Dirichlet-multinomial likelihood. Only the reweighting schemes named in
#' `components` are computed.
#'
#' @param events_array A `D x T x A` array of event counts.
#' @param components Character vector, any non-empty subset of `c("dyad",
#'   "period", "action")`, naming which reweighting scheme(s) to compute.
#'   Defaults to all three.
#' @return A named list containing `dyad_weight` (length D), `period_weight`
#'   (length T), and/or `action_weight` (length A) — whichever were
#'   requested in `components`. Components not requested are simply absent
#'   from the returned list, not present as `NULL`.
#' @keywords internal
compute_default_weights <- function(events_array, components = c("dyad", "period", "action")) {
  components <- match.arg(components, choices = c("dyad", "period", "action"), several.ok = TRUE)
  normalize <- function(w) w / mean(w)
  out <- list()

  if ("dyad" %in% components) {
    dyad_totals <- apply(events_array, 1, sum)
    out$dyad_weight <- normalize(1 / log1p(pmax(dyad_totals, 1)))
  }
  if ("period" %in% components) {
    period_totals <- apply(events_array, 2, sum)
    out$period_weight <- normalize(1 / log1p(pmax(period_totals, 1)))
  }
  if ("action" %in% components) {
    action_means <- apply(events_array, 3, function(x) {
      nonzero <- x[x > 0]
      if (length(nonzero) == 0) {
        return(1)
      }
      mean(nonzero)
    })
    out$action_weight <- normalize(1 / log1p(action_means))
  }

  out
}

#' Parse the `weighted` argument of [assemble_stan_data()]
#'
#' @param weighted `FALSE`, or a single string: `"all"`, or a
#'   `"-"`-separated combination of `"dyad"`, `"period"`, `"action"` (e.g.
#'   `"dyad"`, `"dyad-period"`, `"period-action"`).
#' @return A character vector, a subset of `c("dyad", "period", "action")`,
#'   suitable as `compute_default_weights()`'s `components` argument.
#'   Empty (`character(0)`) if `weighted` is `FALSE`.
#' @keywords internal
parse_weighted_arg <- function(weighted) {
  if (isFALSE(weighted)) {
    return(character(0))
  }
  if (isTRUE(weighted)) {
    stop(
      "`weighted = TRUE` is no longer supported; use `weighted = \"all\"` instead ",
      "(or a specific combination, e.g. `weighted = \"dyad-period\"`).",
      call. = FALSE
    )
  }
  if (!is.character(weighted) || length(weighted) != 1 || is.na(weighted)) {
    stop(
      "`weighted` must be `FALSE` or a single string (\"all\", or a \"-\"-separated ",
      "combination of \"dyad\", \"period\", \"action\").",
      call. = FALSE
    )
  }
  if (identical(weighted, "all")) {
    return(c("dyad", "period", "action"))
  }

  components <- unique(stringr::str_split(weighted, "-")[[1]])
  valid <- c("dyad", "period", "action")
  invalid <- setdiff(components, valid)
  if (length(invalid) > 0) {
    stop(
      "Unrecognized weighting component(s) in `weighted`: ", paste(invalid, collapse = ", "), ". ",
      "Valid components are \"dyad\", \"period\", \"action\" (combine with \"-\", e.g. ",
      "\"dyad-period\"), or \"all\".",
      call. = FALSE
    )
  }
  components
}

#' Assemble Stan-ready data for the bilatr dyadic IRT model
#'
#' Aggregates CAMEO-coded event data to dyad-period action-class counts
#' and packages it as the data list expected by the package's Stan model
#' (`inst/stan/bilatr_dirmult_irt.stan`). This single function covers all
#' three previously-separate model variants (base, event-weighted,
#' dyad+period+event-weighted): `dyad_weight`, `period_weight`, and
#' `action_weight` default to vectors of 1s (recovering the base,
#' unweighted model), can be supplied explicitly, or computed
#' automatically per-component via `weighted`.
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
#'   model's scale/sign reference (`alpha[1] = 1`). Should typically be
#'   a low-conflict/cooperative class. Every other action class's
#'   discrimination (`alpha[2:A]`) is freely estimated.
#' @param min_n_events Minimum total events for a dyad to be retained.
#' @param weighted `FALSE` (default), or a single string naming which
#'   default reweighting scheme(s) to compute via
#'   [compute_default_weights()] for any of `dyad_weight`/`period_weight`/
#'   `action_weight` left as `NULL`: `"all"` for all three, or a
#'   `"-"`-separated combination of `"dyad"`, `"period"`, `"action"` for
#'   just those (e.g. `"dyad"`, `"dyad-period"`, `"period-action"`).
#'   Components not named here, and not requested, still default to 1s.
#'   `weighted = TRUE` is no longer accepted; use `weighted = "all"`.
#' @param dyad_weight Optional numeric vector of length D (number of
#'   retained dyads) rescaling each dyad's contribution to the
#'   likelihood. Defaults to 1s (or computed, if `"dyad"` is requested via
#'   `weighted`).
#' @param period_weight Optional numeric vector of length T (number of
#'   time periods) rescaling each period's contribution. Defaults to 1s
#'   (or computed, if `"period"` is requested via `weighted`).
#' @param action_weight Optional numeric vector of length A (number of
#'   action classes) rescaling the concentration per action type.
#'   Defaults to 1s (or computed, if `"action"` is requested via
#'   `weighted`).
#' @param chunk_size `reduce_sum` grainsize used to chunk the likelihood
#'   across dyads. 16 cores with `chunk_size = 600` was the
#'   Pareto-optimal setting found in this project's own benchmarking for
#'   panel-sized data; smaller panels or single-dyad fits should use a
#'   smaller value (down to 1) since chunking overhead outweighs the
#'   benefit below that scale. See [fit_panel()] for the corresponding
#'   `threads_per_chain` argument.
#' @return A named list suitable as the `data` argument to
#'   `cmdstanr::CmdStanModel$sample()` for the bilatr Stan model: `D`,
#'   `T`, `A`, `C`, `is_obs`, `Y`, `dyad_weight`, `period_weight`,
#'   `action_weight`. Also carries a `dyad_ids` attribute (the output of
#'   [make_dyad_ids()]) for reattaching identifiers to posterior draws;
#'   see [extract_theta()].
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
#'
#' # compute default dyad- and period-level reweighting, but leave
#' # action_weight at 1s:
#' stan_data_dp <- assemble_stan_data(
#'   events,
#'   years = 2015:2020,
#'   resolution = "yearly",
#'   grouping_var = "PentaClass",
#'   reference_category = 0,
#'   weighted = "dyad-period"
#' )
#'
#' # compute all three default reweighting schemes:
#' stan_data_all <- assemble_stan_data(
#'   events,
#'   years = 2015:2020,
#'   resolution = "yearly",
#'   grouping_var = "PentaClass",
#'   reference_category = 0,
#'   weighted = "all"
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
  dyad_weight = NULL,
  period_weight = NULL,
  action_weight = NULL,
  chunk_size = 100
) {
  resolution <- match.arg(resolution)

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

  weight_components <- parse_weighted_arg(weighted)
  weights <- if (length(weight_components) > 0) {
    compute_default_weights(events_array, components = weight_components)
  } else {
    NULL
  }
  dyad_weight <- dyad_weight %||% weights$dyad_weight %||% rep(1, D)
  period_weight <- period_weight %||% weights$period_weight %||% rep(1, Tn)
  action_weight <- action_weight %||% weights$action_weight %||% rep(1, Anum)

  stan_data <- list(
    D = D,
    T = Tn,
    A = Anum,
    C = chunk_size,
    is_obs = obs_matrix,
    Y = events_array,
    dyad_weight = dyad_weight,
    period_weight = period_weight,
    action_weight = action_weight
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
