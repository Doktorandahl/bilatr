#' Assemble Stan-ready data for the bilatr dyadic IRT model
#'
#' Aggregates event data to dyad-period action-class counts and packages
#' it as the data list expected by the package's Stan model (the
#' registered `stable` program; see `R/model_registry.R`).
#'
#' Validates `data` against `?bilatr_event_data` **first, against every
#' row** -- including rows outside `years` -- so a data problem in an
#' out-of-window row (e.g. an unparseable date) is still caught (0.9.1,
#' audit 0c). Only after validation does it drop rows outside `years` and
#' determine the set of observed event classes (see `@param years`), then
#' aggregate the survivors itself (the worker also used by
#' [grouped_events_to_dyad_period()], which this function no longer calls
#' -- see NEWS 0.9.1, audit 0f).
#'
#' @param data A data frame of event-level records; see
#'   `?bilatr_event_data`. Two helpers produce data in this format:
#'   [read_gdelt()]/[download_gdelt()] and [recode_cameo()].
#' @param years Integer vector of years to cover. Rows outside `years` are
#'   dropped -- after validation (see above), and before the class order
#'   is built (0.9.0; previously an action class seen only outside `years`
#'   still got an all-zero column of `Y`, identified only by its prior).
#'   A `message()` reports how many rows were dropped, and names any class
#'   present in `data` that has zero in-window events (and so is absent
#'   from `event_classes`/`A`/`Y`). If no rows of `data` fall within
#'   `years`, errors and reports the data's actual year range (0.9.1,
#'   audit 0e); an empty `data` errors too.
#' @param resolution Either `"monthly"` or `"yearly"`.
#' @param grouping_var Name of the event-class column to aggregate on
#'   (e.g. `"QuadClass"`, `"PentaClass"`).
#' @param directed If `TRUE` (default), dyads are directed; if `FALSE`,
#'   actor order is ignored.
#' @param actor1,actor2 Names of the sender/target actor columns; see
#'   `?bilatr_event_data`.
#' @param date Name of the event-date column; see `?bilatr_event_data`.
#' @param reference_category Value of `grouping_var` to anchor as the
#'   model's scale/sign reference: `stable`/`ou` fold `alpha[1]`'s sign
#'   into the reported `alpha`/`theta` (see each `.stan` file's header,
#'   "IDENTIFICATION: ORIENTATION FOLD"), and RMS-normalize the whole
#'   `alpha` vector to 1, so positive `alpha[1]` means better relations at
#'   this reference/neutral class. Should typically be a
#'   low-conflict/cooperative class. Every other action class's
#'   discrimination (`alpha[2:A]`) is freely estimated. `NULL` (default)
#'   anchors on the first class in C-locale numeric sort instead; if
#'   supplied but absent from the (in-window) data, errors (0.9.0; see
#'   [validate_reference_class()]).
#' @param min_n_events Minimum total events for a dyad to be retained.
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
#' @return A named list suitable as the `data` argument to
#'   `cmdstanr::CmdStanModel$sample()` for the bilatr Stan model: `D`,
#'   `T`, `A`, `C`, `is_obs`, `Y`, `rho_prior_a`, `rho_prior_b`,
#'   `compute_log_lik`, `prior_only`, `compute_theta_filtered`,
#'   `n_filter_dyads`, `filter_dyads`, and (0.7.0,
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
#'   `event_classes` attribute, a `grouping_var` attribute (0.9.0; the
#'   `grouping_var` argument itself, for downstream class-label
#'   derivation), and (0.7.0) a `country_codes` attribute (character, in
#'   `ctry_a`/`ctry_b` index order) for labelling `gamma` (see
#'   [extract_gamma()]).
#' @examples
#' \dontrun{
#' events <- download_gdelt("2020-01-01")
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
  chunk_size = 100,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  actor1 = "Actor1CountryCode",
  actor2 = "Actor2CountryCode",
  date = "SQLDATE"
) {
  resolution <- match.arg(resolution)

  validate_bilatr_events(data, grouping_var, actor1, actor2, date)
  # Validated and slimmed exactly once per call (0.9.1, audit 0f): the
  # date column is parsed here and nowhere else in this function, and
  # the slim table -- not the original `data` -- is what's windowed and
  # handed to .dyad_period_from_slim() below (the same worker
  # grouped_events_to_dyad_period() calls), so that function no longer
  # re-validates or re-parses.
  slim <- .slim_event_table(data, grouping_var, actor1, actor2, date)

  # Window-first (0.9.0, audit D3): drop rows outside `years` before the
  # class order is built, so a class seen only outside the window never
  # gets an all-zero column of `Y` identified only by its prior. `
  # all_classes` is captured from the FULL slim table, before windowing,
  # so the dropped-class message below can compare against every class
  # in `data`, not just the in-window survivors.
  all_classes <- unique(slim$.class)
  n_total <- nrow(slim)
  full_year_range <- if (n_total > 0) range(slim$.year) else NULL
  in_window <- slim$.year %in% years
  n_dropped <- sum(!in_window)
  if (n_dropped > 0) {
    message(sprintf(
      "assemble_stan_data(): dropping %d row(s) outside `years`.", n_dropped
    ))
  }
  slim <- slim[in_window, , drop = FALSE]

  # Empty-window guard (0.9.1, audit 0e): before this, an empty `data` or
  # a `years` matching no rows fell through to the "No dyads have at
  # least min_n_events" error below, which points at the wrong knob
  # (min_n_events, not years/data).
  if (nrow(slim) == 0) {
    if (n_total == 0) {
      stop("assemble_stan_data(): `data` is empty; nothing to assemble.", call. = FALSE)
    }
    stop(sprintf(
      "assemble_stan_data(): no events fall within `years` (data's years span %d-%d).",
      full_year_range[1], full_year_range[2]
    ), call. = FALSE)
  }

  zero_window_classes <- setdiff(all_classes, unique(slim$.class))
  if (length(zero_window_classes) > 0) {
    message(
      "assemble_stan_data(): the following class(es) of `", grouping_var,
      "` have zero in-window events and will not appear in `event_classes`/`Y`: ",
      paste(sort(zero_window_classes), collapse = ", ")
    )
  }

  agg <- .dyad_period_from_slim(slim, resolution, directed, reference_category)
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
  # why). Row i of `country_info` corresponds to `dyads[i]` (the left
  # side of the join below, so join order -- not `w_send_tbl`'s own row
  # order -- determines it): `dyads`' own D-indexing is exactly the
  # `dyad_id` order [make_dyad_ids()] assigns below, so
  # ctry_a/ctry_b/w_send end up indexed 1:D the same way Y/is_obs/
  # dyad_ids already are. Built from the assembled table's own carried
  # actor_a/actor_b (0.9.0; no more parsing 3-letter codes out of `dyad`)
  # -- every retained dyad has both an actor pair (present since it
  # survived to `agg`) and a `w_send` entry (guaranteed since
  # validate_bilatr_events() bans the NA actors that used to cause
  # misses here), so the 0.7.1 "no matching entry" guard is no longer
  # reachable and has been removed.
  dyad_actors <- dplyr::distinct(agg, dyad, actor_a, actor_b)
  country_info <- tibble::tibble(dyad = dyads) %>%
    dplyr::left_join(dyad_actors, by = "dyad") %>%
    dplyr::left_join(dplyr::select(w_send_tbl, dyad, w_send), by = "dyad")

  # method = "radix" (0.7.1): base sort()'s default method is locale-
  # collation-dependent for character input, and this ordering *defines*
  # what gamma's columns mean -- a stan_data.rds re-derived on a machine
  # with a different locale must not silently relabel countries. radix
  # is C-locale/byte order, so this is deterministic across machines (see
  # order_event_classes()'s locale = "C" str_sort() for the same
  # treatment of the action-class index).
  countries_present <- sort(
    unique(c(country_info$actor_a, country_info$actor_b)),
    method = "radix"
  )
  n_countries <- length(countries_present)
  country_index <- stats::setNames(seq_len(n_countries), countries_present)
  ctry_a <- unname(country_index[country_info$actor_a])
  ctry_b <- unname(country_index[country_info$actor_b])
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

  dyad_ids <- make_dyad_ids(agg)

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
    # Consumed only by the experimental stable_gamma variant
    # (R/model_registry.R); unused by every other registered program,
    # whose Stan code doesn't declare these fields -- CmdStan ignores
    # data a program doesn't declare, so (matching rho_prior_a/b above)
    # these are passed unconditionally rather than gated on which program
    # will be fit, so a caller doesn't need to know the model name at
    # assembly time. BREAKING: a `stan_data.rds` saved before 0.7.0 lacks
    # these four fields and cannot be used with `stable_gamma` -- see
    # NEWS.md.
    n_countries = n_countries,
    ctry_a = ctry_a,
    ctry_b = ctry_b,
    w_send = w_send
  )

  attr(stan_data, "dyad_ids") <- dyad_ids
  attr(stan_data, "event_classes") <- event_classes
  attr(stan_data, "country_codes") <- countries_present
  attr(stan_data, "grouping_var") <- grouping_var

  stan_data
}

#' Prepare Stan-ready data from a dyad-period table
#'
#' Packages a dyad x time-period table of action counts (see
#' `?bilatr_dyad_period_data`), from [aggregate_events()] or any other
#' source, as the data list expected by the package's Stan models (see
#' `R/model_registry.R`).
#'
#' In order:
#'
#' 1. validates `data` with [validate_bilatr_dyad_period()], building the
#'    dyad column from the actor columns if it is absent;
#' 2. drops action columns that are zero in every row, with a `message()`
#'    naming them;
#' 3. orders the action columns: `reference_category` first, then the
#'    rest in `data`'s column order;
#' 4. fills every dyad over the whole `time_range` (by default
#'    `min(time):max(time)` of the table), missing periods becoming
#'    all-zero and unobserved (`is_obs = 0`; a period is observed when its
#'    counts sum to more than zero);
#' 5. drops dyads with fewer than `min_n_events` events in total;
#' 6. if both actor columns are present, builds the country index
#'    (`n_countries`, `ctry_a`, `ctry_b`) and `w_send`, with a `message()`
#'    reporting `n_countries` and the dyads-per-country distribution;
#' 7. resolves `filter_dyads` and builds the list and its attributes.
#'
#' @section Size:
#' `Y` is a dense `D x T x A` integer array, about `D * T * A * 4` bytes
#' (the filled grid has `D * T` rows). `T` counts every period in
#' `time_range`, observed or not, so finer time units grow it fast:
#' 16,750 dyads over 35 years with 16 action types is about 37 MB yearly,
#' 450 MB monthly and 13.7 GB daily. Daily (or similarly fine) resolution
#' is only practical for few dyads or short windows; use `time_range` to
#' pin `T`.
#'
#' @param data A dyad-period data frame; see `?bilatr_dyad_period_data`.
#' @param time Name of the time-period column. `NULL` (default) takes
#'   `attr(data, "time_col")` (set by [aggregate_events()]); failing that,
#'   whichever single one of `year`, `month_id`, `week_id` or `day_id` is
#'   present, with a `message()`; otherwise errors.
#' @param action_cols The action-count columns: a character vector or a
#'   tidyselect expression. `NULL` (default) means every column except
#'   the dyad/actor/time/`w_send` columns, with a warning; see
#'   `?bilatr_dyad_period_data`.
#' @param dyad,actor1,actor2,w_send Names of the dyad, actor and
#'   sender-share columns; see `?bilatr_dyad_period_data`. The actor and
#'   `w_send` columns are optional.
#' @param directed Whether the dyads are directed. Defaults to
#'   `attr(data, "directed")` if set (by [aggregate_events()]), else
#'   `TRUE`. Only matters for `w_send`: with `directed = TRUE` a `w_send`
#'   column must be all `1`, and a missing one is taken as `1`; with
#'   `directed = FALSE` a missing `w_send` is taken as `0.5` and the
#'   `w_send_imputed` attribute is set to `TRUE`.
#' @param reference_category Name of the action column to anchor as the
#'   model's scale/sign reference: `stable`/`ou` fold `alpha[1]`'s sign
#'   into the reported `alpha`/`theta` (see each `.stan` file's header,
#'   "IDENTIFICATION: ORIENTATION FOLD"), so positive `alpha[1]` means
#'   better relations at this class. Should typically be a
#'   low-conflict/cooperative class. `NULL` (default) anchors on the first
#'   action column. Errors if it is not one of the (non-all-zero) action
#'   columns.
#' @param min_n_events Minimum total events for a dyad to be retained.
#' @param time_range `NULL` (default), or two whole numbers `c(min, max)`
#'   on the time column's scale. Rows outside it are dropped (with a
#'   `message()` giving the count), and every dyad is filled over
#'   `min:max`, so `T = max - min + 1` whatever periods the data covers.
#'   `NULL` means `min:max` of the data. Errors if it is not two whole
#'   numbers, if `min > max`, or if no rows fall inside it.
#' @param time_labels A data frame keyed by the time column, whose other
#'   columns (e.g. `year`, `month`) are joined into the `dyad_ids`
#'   attribute. Defaults to `attr(data, "time_labels")` (set by
#'   [aggregate_events()]); `NULL` adds none.
#' @param grouping_var Stored as the `grouping_var` attribute, for
#'   [event_class_labels()]. Defaults to `attr(data, "grouping_var")`.
#' @inheritParams assemble_stan_data
#' @return A named list suitable as the `data` argument to
#'   `cmdstanr::CmdStanModel$sample()`: `D`, `T`, `A`, `C`, `is_obs`,
#'   `Y`, `rho_prior_a`, `rho_prior_b`, `compute_log_lik`, `prior_only`,
#'   `compute_theta_filtered`, `n_filter_dyads`, `filter_dyads`, and, only
#'   when both actor columns are present, `n_countries`, `ctry_a`,
#'   `ctry_b` and `w_send` (consumed only by the experimental
#'   `stable_gamma` variant). `Y` is an integer `D x T x A` array, in
#'   `dyad_id` order (dyads sorted in C-locale byte order) and time order.
#'   Attributes:
#'
#'   * `dyad_ids`: one row per dyad-period, with `dyad_id`, `time_index`,
#'     `dyad`, `dyad2` (the undirected dyad key), `actor1`, `actor2`
#'     (`dyad2`/`actor1`/`actor2` are `NA` without actor columns), the time
#'     column, and any `time_labels` columns; see [extract_theta()];
#'   * `event_classes`: the action-column names, in `Y`'s order;
#'   * `country_codes`: the countries in `ctry_a`/`ctry_b` index order
#'     (`NULL` without actor columns); see [extract_gamma()];
#'   * `grouping_var`, `time_col`, `time_origin` (from
#'     `attr(data, "time_origin")`, else `NA`), `directed`;
#'   * `w_send_imputed`: `TRUE` if `w_send` was absent and taken as `0.5`
#'     for undirected dyads, else `FALSE`.
#' @examples
#' # A non-event source: arbitrary dyad labels, no actor columns, a weekly
#' # period index with gaps, and three action types.
#' dp <- tibble::tibble(
#'   dyad = c("a-b", "a-b", "a-b", "a-c", "a-c", "b-c"),
#'   week = c(1L, 2L, 6L, 1L, 3L, 4L),
#'   talks = c(3, 0, 2, 1, 4, 2),
#'   threats = c(0, 1, 0, 2, 0, 1),
#'   aid = c(1, 1, 0, 0, 2, 0)
#' )
#' stan_data <- prepare_stan_data(
#'   dp,
#'   time = "week", action_cols = c(talks, threats, aid),
#'   reference_category = "aid"
#' )
#' dim(stan_data$Y)
#'
#' # From event data
#' events <- tibble::tibble(
#'   Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN"),
#'   Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA"),
#'   SQLDATE = c(20180105L, 20180610L, 20190302L, 20190401L, 20200101L),
#'   PentaClass = c(0, 1, 0, 1, 3)
#' )
#' agg <- aggregate_events(events, grouping_var = "PentaClass")
#' stan_data <- prepare_stan_data(agg, action_cols = c("0", "1", "3"))
#' @export
prepare_stan_data <- function(
  data,
  time = NULL,
  action_cols = NULL,
  dyad = "dyad",
  actor1 = "actor1",
  actor2 = "actor2",
  w_send = "w_send",
  directed = attr(data, "directed") %||% TRUE,
  reference_category = NULL,
  min_n_events = 1,
  time_range = NULL,
  chunk_size = 100,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  time_labels = attr(data, "time_labels"),
  grouping_var = attr(data, "grouping_var")
) {
  fn <- "prepare_stan_data()"
  if (!rlang::is_bool(directed)) {
    stop(fn, ": `directed` must be TRUE or FALSE.", call. = FALSE)
  }
  time <- .resolve_time_col(data, time)

  # 1. validate; resolves action_cols exactly once
  cols <- .validate_dyad_period(data, time, rlang::enquo(action_cols), dyad, actor1, actor2, w_send)
  has_actors <- !is.null(cols$actor1)
  dyad_vec <- if (dyad %in% names(data)) {
    as.character(data[[dyad]])
  } else {
    paste(data[[actor1]], data[[actor2]], sep = "_")
  }
  time_vec <- as.integer(data[[time]])
  w_send_vec <- if (!is.null(cols$w_send)) as.double(data[[w_send]])
  if (directed && !is.null(w_send_vec) && any(w_send_vec != 1)) {
    stop(
      fn, ": `directed = TRUE`, but `", w_send, "` is not 1 in every row (a directed ",
      "dyad's sender is always actor1). Set `directed = FALSE` for undirected data.",
      call. = FALSE
    )
  }
  counts <- purrr::map(data[cols$action_cols], as.integer)
  actor1_vec <- if (has_actors) as.character(data[[actor1]])
  actor2_vec <- if (has_actors) as.character(data[[actor2]])

  # 1b. time_range: drop rows outside it; the grid is filled over it below
  if (is.null(time_range)) {
    time_range <- range(time_vec)
  } else {
    if (!is.numeric(time_range) || length(time_range) != 2 || anyNA(time_range) ||
        !all(.is_whole_number(time_range))) {
      stop(
        fn, ": `time_range` must be two whole numbers, c(min, max), on the scale of `",
        time, "`.",
        call. = FALSE
      )
    }
    if (time_range[1] > time_range[2]) {
      stop(
        fn, ": `time_range` min (", time_range[1], ") is greater than its max (",
        time_range[2], ").",
        call. = FALSE
      )
    }
    in_range <- time_vec >= time_range[1] & time_vec <= time_range[2]
    if (!any(in_range)) {
      stop(
        fn, ": no rows fall within `time_range` (", time_range[1], " to ", time_range[2],
        "); `", time, "` spans ", min(time_vec), " to ", max(time_vec), ".",
        call. = FALSE
      )
    }
    if (!all(in_range)) {
      message(sprintf("%s: dropping %d row(s) outside `time_range`.", fn, sum(!in_range)))
      dyad_vec <- dyad_vec[in_range]
      time_vec <- time_vec[in_range]
      w_send_vec <- w_send_vec[in_range]
      counts <- purrr::map(counts, ~ .x[in_range])
      actor1_vec <- actor1_vec[in_range]
      actor2_vec <- actor2_vec[in_range]
    }
  }
  time_values <- seq.int(as.integer(time_range[1]), as.integer(time_range[2]))

  # 2. drop all-zero action columns, before the min_n_events filter
  all_zero <- purrr::map_lgl(counts, ~ all(.x == 0L))
  if (any(all_zero)) {
    message(
      fn, ": dropping action column(s) with zero events in every row: ",
      paste(names(counts)[all_zero], collapse = ", ")
    )
    counts <- counts[!all_zero]
  }
  if (length(counts) < 2) {
    stop(
      fn, ": fewer than two action columns have any events (",
      paste(names(counts), collapse = ", "), "); the model needs at least two.",
      call. = FALSE
    )
  }

  # 3. reference first, then the user's column order
  reference_category <- validate_reference_class(
    reference_category, names(counts), "reference_category",
    null_hint = "the first action column"
  )
  event_classes <- c(reference_category, setdiff(names(counts), reference_category))
  counts <- counts[event_classes]
  A <- length(event_classes)

  # 4. fill the dyad x time grid. Counts get positional internal names
  # (.k1, .k2, ...) so no action column can collide with the key columns.
  count_names <- paste0(".k", seq_len(A))
  long <- tibble::tibble(.dyad = dyad_vec, .time = time_vec, !!!rlang::set_names(counts, count_names))
  # C-locale (radix) sort: the same order dplyr::group_by() gives, and so
  # the same dyad_id order make_dyad_ids() assigns.
  dyads <- sort(unique(dyad_vec), method = "radix")
  grid <- tidyr::expand_grid(.dyad = dyads, .time = time_values) %>%
    dplyr::left_join(long, by = c(".dyad", ".time")) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(count_names), ~ tidyr::replace_na(.x, 0L)))
  total <- rowSums(as.matrix(grid[count_names]))

  # 5. min_n_events filter
  dyad_totals <- tapply(total, factor(grid$.dyad, levels = dyads), sum)
  keep <- dyad_totals >= min_n_events
  dyads <- dyads[keep]
  D <- length(dyads)
  if (D == 0) {
    stop(
      "No dyads have at least `min_n_events` (", min_n_events, ") total events; ",
      "nothing to assemble. Lower min_n_events or check the input data.",
      call. = FALSE
    )
  }
  keep_rows <- grid$.dyad %in% dyads
  grid <- grid[keep_rows, , drop = FALSE]
  total <- total[keep_rows]
  # nrow / D, as assemble_stan_data() computes it (a double).
  Tn <- nrow(grid) / D

  is_obs <- matrix(as.integer(total > 0), nrow = D, byrow = TRUE)
  Y <- array(
    unlist(purrr::map(count_names, ~ matrix(grid[[.x]], nrow = D, byrow = TRUE)), use.names = FALSE),
    dim = c(D, Tn, A)
  )

  # 6. country index, only with actor columns
  dyad_info <- tibble::tibble(
    dyad = dyad_vec,
    actor1 = actor1_vec %||% NA_character_,
    actor2 = actor2_vec %||% NA_character_,
    w_send = w_send_vec %||% NA_real_
  ) %>%
    dplyr::distinct(.data$dyad, .keep_all = TRUE)
  dyad_info <- dplyr::left_join(tibble::tibble(dyad = dyads), dyad_info, by = "dyad")

  # 7. w_send (1 when directed; 0.5, flagged, when undirected and absent)
  w_send_imputed <- is.null(w_send_vec) && !directed
  if (is.null(w_send_vec)) {
    dyad_info$w_send <- rep(if (directed) 1 else 0.5, D)
  }

  country_fields <- NULL
  countries_present <- NULL
  if (has_actors) {
    # radix = C-locale byte order, as in assemble_stan_data(): this order
    # defines what gamma's columns mean, so it must not depend on locale.
    countries_present <- sort(unique(c(dyad_info$actor1, dyad_info$actor2)), method = "radix")
    n_countries <- length(countries_present)
    country_index <- stats::setNames(seq_len(n_countries), countries_present)
    ctry_a <- unname(country_index[dyad_info$actor1])
    ctry_b <- unname(country_index[dyad_info$actor2])
    dyads_per_country <- table(c(ctry_a, ctry_b))
    message(sprintf(
      "%s: n_countries = %d (dyads per country -- min %d, median %s, max %d)",
      fn, n_countries, min(dyads_per_country),
      format(stats::median(as.numeric(dyads_per_country)), nsmall = 1),
      max(dyads_per_country)
    ))
    country_fields <- list(
      n_countries = n_countries,
      ctry_a = ctry_a,
      ctry_b = ctry_b,
      w_send = dyad_info$w_send
    )
  }

  dyad_ids <- .make_dyad_ids_generic(
    dyads,
    actor1 = if (has_actors) dyad_info$actor1,
    actor2 = if (has_actors) dyad_info$actor2,
    time_values = time_values,
    time_col = time,
    time_labels = time_labels
  )

  # 8. filter_dyads, the list, and attributes
  if (compute_theta_filtered == 1) {
    if (is.null(filter_dyads)) {
      filter_dyad_ids <- seq_len(D)
    } else {
      filter_dyad_ids <- match(filter_dyads, dyads)
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

  stan_data <- c(
    list(
      D = D,
      T = Tn,
      A = A,
      C = chunk_size,
      is_obs = is_obs,
      Y = Y,
      rho_prior_a = rho_prior_a,
      rho_prior_b = rho_prior_b,
      compute_log_lik = compute_log_lik,
      prior_only = prior_only,
      compute_theta_filtered = compute_theta_filtered,
      n_filter_dyads = length(filter_dyad_ids),
      filter_dyads = filter_dyad_ids
    ),
    # Only stable_gamma declares these; CmdStan ignores undeclared fields.
    country_fields
  )

  attr(stan_data, "dyad_ids") <- dyad_ids
  attr(stan_data, "event_classes") <- event_classes
  attr(stan_data, "country_codes") <- countries_present
  attr(stan_data, "grouping_var") <- grouping_var
  attr(stan_data, "time_col") <- time
  attr(stan_data, "time_origin") <- attr(data, "time_origin") %||% as.Date(NA)
  attr(stan_data, "directed") <- directed
  attr(stan_data, "w_send_imputed") <- w_send_imputed

  stan_data
}

#' Resolve [prepare_stan_data()]'s `time = NULL`
#'
#' @param data The dyad-period data frame.
#' @param time The `time` argument.
#' @return A single column name.
#' @keywords internal
.resolve_time_col <- function(data, time) {
  if (!is.null(time)) {
    return(time)
  }
  from_attr <- attr(data, "time_col")
  if (!is.null(from_attr)) {
    return(from_attr)
  }
  known <- intersect(c("year", "month_id", "week_id", "day_id"), names(data))
  if (length(known) == 1) {
    message("prepare_stan_data(): using `", known, "` as the time column.")
    return(known)
  }
  stop(
    "prepare_stan_data(): cannot tell which column is the time period",
    if (length(known) > 1) paste0(" (found ", paste(known, collapse = ", "), ")"),
    "; set `time`.",
    call. = FALSE
  )
}
