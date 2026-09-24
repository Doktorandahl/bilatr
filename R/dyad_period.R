#' Validate a reference-category argument against the observed classes
#'
#' @param value The candidate reference value, or `NULL`.
#' @param classes Vector of observed event-class values.
#' @param arg_name Name of the argument, used in the error message.
#' @return `value` (coerced to character) if it is present in `classes`.
#'   `NULL` is returned as-is: the only way to say "no preference", which
#'   then anchors on the first class in [order_event_classes()]'s
#'   C-locale numeric sort.
#' @keywords internal
validate_reference_class <- function(value, classes, arg_name) {
  if (is.null(value)) {
    return(NULL)
  }
  value <- as.character(value)
  observed <- unique(as.character(classes))
  if (!(value %in% observed)) {
    stop(
      arg_name, " = '", value, "' was requested, but is not present in the ",
      "(in-window) data. Classes present: ",
      paste(sort(observed), collapse = ", "),
      ". `", arg_name, " = NULL` anchors on the first class instead ",
      "(C-locale numeric sort).",
      call. = FALSE
    )
  }
  value
}

#' Order event classes with the reference category first
#'
#' Puts `reference_category` first (if supplied and present), with all other
#' classes sorted alphabetically after it. This ordering is what implements
#' the model's identification constraint on the R side: `stable`/`ou` fold
#' `alpha[1]`'s (the first column's) sign into the reported `alpha`/`theta`
#' (see each `.stan` file's header, "IDENTIFICATION: ORIENTATION FOLD"),
#' so the reference/neutral class anchors alpha's sign and scale. All
#' remaining `alpha[2:A]` are freely estimated.
#'
#' @param classes Vector of observed event-class values.
#' @param reference_category Value to place first, or `NULL`.
#' @return Character vector of unique classes in anchor order.
#' @keywords internal
order_event_classes <- function(classes, reference_category = NULL) {
  # locale = "C" (0.7.1): str_sort()'s default locale is the system's
  # ("en" here), which is locale-COLLATION-dependent -- and this ordering
  # *defines* what alpha/mu_intercept/gamma's action-index columns mean
  # (the same class of concern as assemble_stan_data()'s country index;
  # see R/stan_data.R). "C" gives byte/ASCII-order collation, so a
  # stan_data.rds re-derived on a machine with a different locale can't
  # silently relabel action classes; `numeric = TRUE` still gives the
  # intended natural/numeric-aware ordering ("2" before "10") under "C".
  classes <- stringr::str_sort(unique(as.character(classes)), numeric = TRUE, locale = "C")
  middle <- setdiff(classes, reference_category)
  c(reference_category, middle)
}

#' Aggregate event data to dyad-period class counts
#'
#' Collapses event-level data (in the format documented at
#' `?bilatr_event_data`; as produced by e.g. [extract_all_relevant_gdelt()]
#' and [recode_cameo()]) to dyad-by-time-period counts of a chosen
#' event-class column, one column per class plus a `total_events` column.
#' This is the shared aggregation step feeding [assemble_stan_data()].
#'
#' Validates `data` with [validate_bilatr_events()] first, then builds a
#' slim internal table (`actor1`/`actor2`/`date`/`grouping_var`, renamed
#' and re-typed, plus a derived `year`/`month`) that every later step in
#' this function works from -- a user's own `dyad`/`date`/`year`/`month`/
#' `event_type` columns, if present in `data`, are never read or
#' overwritten.
#'
#' @param data A data frame of event-level records; see
#'   `?bilatr_event_data`.
#' @param resolution Either `"monthly"` or `"yearly"`.
#' @param grouping_var Name of the event-class column to aggregate on
#'   (e.g. `"QuadClass"`, `"PentaClass"`).
#' @param directed If `TRUE` (default), dyads are directed
#'   (actor1 -> actor2); if `FALSE`, actor order is ignored and dyads are
#'   collapsed to an unordered pair, with side A/B assigned by
#'   [.order_pair_c_locale()] (deterministic, locale-independent).
#' @param reference_category Value of `grouping_var` to place first in
#'   the class ordering (the model's scale/sign-reference, neutral
#'   action: `stable`/`ou` build `alpha[1]` positive by construction --
#'   see [order_event_classes()]). `NULL` (default) anchors on the first
#'   class in C-locale numeric sort instead. If supplied but not present
#'   among the (in-window) classes, errors -- see
#'   [validate_reference_class()]. All other action classes'
#'   discrimination is freely estimated.
#' @param years Integer vector of years, or `NULL` (default, meaning
#'   "use every row"). When supplied, rows outside `years` are dropped
#'   before any aggregation -- for every computation this function does,
#'   including `w_send` below, not only (as before 0.9.0) the `w_send`
#'   computation.
#' @param actor1,actor2 Names of the sender/target actor columns; see
#'   `?bilatr_event_data`.
#' @param date Name of the event-date column; see `?bilatr_event_data`.
#' @return A data frame with columns `dyad`, `year` (and `month` if
#'   `resolution = "monthly"`), `actor_a`/`actor_b` (the directed dyad's
#'   sender/target, or the undirected pair in
#'   [.order_pair_c_locale()]'s deterministic order), one
#'   `EventClass_<value>` column per observed class (ordered per
#'   `reference_category`), and `total_events`. Also carries a `w_send`
#'   attribute (0.7.0; consumed by [assemble_stan_data()] for the
#'   experimental `stable_gamma` Stan variant -- see
#'   `R/model_registry.R`): a tibble with one row per dyad, `actor_a`/
#'   `actor_b`, and `w_send` (the share of that pair's events with side A
#'   as sender, over the same rows the main aggregation uses). For
#'   directed data this is exactly `1` by construction (side A IS the
#'   sender in a directed dyad key), asserted internally below as a free
#'   check on the dyad-key convention.
#' @examples
#' \dontrun{
#' events <- extract_all_relevant_gdelt("data/gdelt_raw/20200101.zip")
#' events <- recode_cameo(events)
#' grouped_events_to_dyad_period(
#'   events,
#'   resolution = "yearly",
#'   grouping_var = "PentaClass",
#'   reference_category = 0
#' )
#' }
#' @export
grouped_events_to_dyad_period <- function(
  data,
  resolution = c("monthly", "yearly"),
  grouping_var,
  directed = TRUE,
  reference_category = NULL,
  years = NULL,
  actor1 = "Actor1CountryCode",
  actor2 = "Actor2CountryCode",
  date = "SQLDATE"
) {
  resolution <- match.arg(resolution)

  validate_bilatr_events(data, grouping_var, actor1, actor2, date)

  slim <- tibble::tibble(
    .actor1 = as.character(data[[actor1]]),
    .actor2 = as.character(data[[actor2]]),
    .date = .parse_event_date(data[[date]]),
    .class = as.character(data[[grouping_var]])
  )
  slim$.year <- as.integer(format(slim$.date, "%Y"))
  slim$.month <- as.integer(format(slim$.date, "%m"))

  if (!is.null(years)) {
    slim <- dplyr::filter(slim, .year %in% years)
  }

  if (directed) {
    slim$actor_a <- slim$.actor1
    slim$actor_b <- slim$.actor2
  } else {
    pair <- .order_pair_c_locale(slim$.actor1, slim$.actor2)
    slim$actor_a <- pair$actor_a
    slim$actor_b <- pair$actor_b
  }
  slim$dyad <- paste(slim$actor_a, slim$actor_b, sep = "_")
  slim$year <- slim$.year
  slim$month <- slim$.month
  slim$event_type <- slim$.class

  reference_category <- validate_reference_class(reference_category, slim$event_type, "reference_category")

  class_order <- order_event_classes(slim$event_type, reference_category)
  group_cols <- c("dyad", "year", if (resolution == "monthly") "month")
  column_order <- c(group_cols, "actor_a", "actor_b", paste0("EventClass_", class_order), "total_events")

  dyad_actors <- dplyr::distinct(slim, dyad, actor_a, actor_b)

  counts <- slim %>%
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), event_type) %>%
    tidyr::pivot_wider(
      names_from = event_type,
      values_from = n,
      values_fill = 0,
      names_prefix = "EventClass_"
    )

  totals <- slim %>%
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), name = "total_events")

  result <- counts %>%
    dplyr::left_join(totals, by = group_cols) %>%
    dplyr::left_join(dyad_actors, by = "dyad") %>%
    dplyr::select(dplyr::any_of(column_order))

  # w_send (0.7.0; see @return): per-dyad sender share, over the same
  # (already year-windowed) rows the main aggregation uses -- computed
  # from the raw event rows rather than from `result`'s per-period
  # counts, since it needs the sender per event, not per action class.
  # For directed data this is exactly 1 by construction (actor_a IS
  # .actor1 above), so no invariant check is needed -- unlike before
  # 0.9.0, there is also no NA-actor case to guard against:
  # validate_bilatr_events() above already guarantees neither actor
  # column has NA, so every row contributes and no dyad can be dropped
  # from this computation.
  w_send <- slim %>%
    dplyr::mutate(is_side_a_sender = .actor1 == actor_a) %>%
    dplyr::group_by(dyad) %>%
    dplyr::summarise(
      actor_a = dplyr::first(actor_a),
      actor_b = dplyr::first(actor_b),
      w_send = mean(is_side_a_sender),
      .groups = "drop"
    )

  attr(result, "w_send") <- w_send
  result
}

#' Build a full dyad x period skeleton, filling in gaps as zero-count rows
#'
#' Right-joins `data` onto the cross-product of observed dyads and the
#' requested time range, so every dyad has an entry for every period even
#' when no events occurred, and flags which periods were actually observed.
#'
#' @param data Output of [grouped_events_to_dyad_period()].
#' @param years Integer vector of years to cover.
#' @param resolution Either `"monthly"` or `"yearly"`.
#' @return `data` expanded to the full dyad x period grid, with count
#'   columns and `total_events` set to 0 for unobserved periods and an
#'   `is_obs` indicator column (1 if `total_events > 0`, else 0).
#' @keywords internal
fill_dyad_period_skeleton <- function(data, years, resolution = c("monthly", "yearly")) {
  resolution <- match.arg(resolution)

  dyad_actors <- dplyr::distinct(data, dyad, actor_a, actor_b)
  skeleton <- if (resolution == "yearly") {
    tidyr::expand_grid(dyad_actors, year = years)
  } else {
    tidyr::expand_grid(dyad_actors, year = years, month = 1:12)
  }

  join_cols <- c("dyad", "actor_a", "actor_b", "year", if (resolution == "monthly") "month")
  count_cols <- c(grep("^EventClass_", names(data), value = TRUE), "total_events")

  data %>%
    dplyr::right_join(skeleton, by = join_cols) %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(join_cols))) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(count_cols), ~ tidyr::replace_na(.x, 0))) %>%
    dplyr::mutate(is_obs = dplyr::if_else(total_events == 0, 0L, 1L)) %>%
    dplyr::ungroup()
}

#' Reattach dyad string identifiers to a fitted model's integer dyad index
#'
#' [assemble_stan_data()] indexes dyads by an integer `dyad_id` (1..D) in
#' row order. This recovers the `dyad_id` <-> `dyad` (and undirected
#' `dyad2`) mapping, so posterior draws indexed by `dyad_id` can be joined
#' back to human-readable identifiers. See [extract_theta()].
#'
#' A pure function of the already-assembled table: `data` must already be
#' skeleton-filled ([fill_dyad_period_skeleton()]) and `min_n_events`-
#' filtered, exactly as [assemble_stan_data()] passes it -- this function
#' does not re-fill or re-filter (0.9.0; previously it duplicated
#' [assemble_stan_data()]'s own fill/filter on top of the "agg" table it
#' was already given, one of which was redundant).
#'
#' @param data Output of [fill_dyad_period_skeleton()] after
#'   `min_n_events` filtering (carrying `dyad`, `actor_a`, `actor_b`,
#'   `year`, and `month` if monthly).
#' @return A data frame with `dyad_id`, `time_index`, `dyad`, `dyad2`
#'   (undirected dyad key), `actor_a`, `actor_b`, `year`, and (if
#'   monthly) `month`.
#' @keywords internal
make_dyad_ids <- function(data) {
  pair <- .order_pair_c_locale(data$actor_a, data$actor_b)
  data$dyad2 <- paste(pair$actor_a, pair$actor_b, sep = "_")

  data %>%
    dplyr::group_by(dyad) %>%
    dplyr::mutate(dyad_id = dplyr::cur_group_id(), time_index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::any_of(c(
      "dyad_id", "time_index", "dyad", "dyad2", "actor_a", "actor_b", "year", "month"
    )))
}
