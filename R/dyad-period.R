#' Validate a reference-category argument against the observed classes
#'
#' @param value The candidate reference value, or `NULL`.
#' @param classes Vector of observed event-class values.
#' @param arg_name Name of the argument, used in the warning message.
#' @return `value` (coerced to character) if it is present in `classes`,
#'   otherwise `NULL` with a warning.
#' @keywords internal
validate_reference_class <- function(value, classes, arg_name) {
  if (is.null(value)) {
    return(NULL)
  }
  value <- as.character(value)
  if (!(value %in% unique(classes))) {
    warning(
      arg_name, " '", value, "' not found in the data; ",
      "proceeding without it.",
      call. = FALSE
    )
    return(NULL)
  }
  value
}

#' Order event classes with reference categories at the ends
#'
#' Puts `reference_category` first (if supplied and present) and
#' `reference_hostile` last (if supplied and present), with all other
#' classes sorted alphabetically in between. This ordering is what
#' implements the model's identification constraints on the R side:
#' `alpha[1] = 1` anchors on the first column and `alpha[A]` anchors
#' (with a forced-negative sign) on the last.
#'
#' @param classes Vector of observed event-class values.
#' @param reference_category Value to place first, or `NULL`.
#' @param reference_hostile Value to place last, or `NULL`.
#' @return Character vector of unique classes in anchor order.
#' @keywords internal
order_event_classes <- function(classes, reference_category = NULL, reference_hostile = NULL) {
  classes <- sort(unique(as.character(classes)))
  middle <- setdiff(classes, c(reference_category, reference_hostile))
  c(reference_category, middle, reference_hostile)
}

#' Aggregate CAMEO-coded events to dyad-period class counts
#'
#' Collapses event-level data (as produced by [extract_all_relevant_gdelt()]
#' or [ingest_icews()]) to dyad-by-time-period counts of a chosen
#' event-class column (e.g. `QuadClass`, `PentaClass`), one column per
#' class plus a `total_events` column. This is the shared aggregation step
#' feeding [assemble_stan_data()].
#'
#' @param data A data frame of event-level records with
#'   `Actor1CountryCode`, `Actor2CountryCode`, `SQLDATE`, and the column
#'   named by `grouping_var`.
#' @param resolution Either `"monthly"` or `"yearly"`.
#' @param grouping_var Name of the event-class column to aggregate on
#'   (e.g. `"QuadClass"`, `"PentaClass"`).
#' @param directed If `TRUE` (default), dyads are directed
#'   (actor1 -> actor2); if `FALSE`, actor order is ignored and dyads are
#'   collapsed to an unordered pair.
#' @param reference_category Value of `grouping_var` to place first in
#'   the class ordering (the model's scale-reference / neutral action).
#'   If `NULL` or not present in the data, ignored with a warning.
#' @param reference_hostile Value of `grouping_var` to place last in the
#'   class ordering (the model's known-hostile anchor). If `NULL` or not
#'   present in the data, ignored with a warning.
#' @return A data frame with columns `dyad`, `year` (and `month` if
#'   `resolution = "monthly"`), one `EventClass_<value>` column per
#'   observed class (ordered per `reference_category`/`reference_hostile`),
#'   and `total_events`.
#' @examples
#' \dontrun{
#' events <- extract_all_relevant_gdelt("data/gdelt_raw/20200101.zip")
#' events <- recode_cameo(events)
#' grouped_events_to_dyad_period(
#'   events,
#'   resolution = "yearly",
#'   grouping_var = "PentaClass",
#'   reference_category = 0,
#'   reference_hostile = 4
#' )
#' }
#' @export
grouped_events_to_dyad_period <- function(
  data,
  resolution = c("monthly", "yearly"),
  grouping_var,
  directed = TRUE,
  reference_category = NULL,
  reference_hostile = NULL
) {
  resolution <- match.arg(resolution)

  data <- data %>%
    dplyr::mutate(
      dyad = if (directed) {
        paste(Actor1CountryCode, Actor2CountryCode, sep = "_")
      } else {
        paste(
          pmin(Actor1CountryCode, Actor2CountryCode),
          pmax(Actor1CountryCode, Actor2CountryCode),
          sep = "_"
        )
      },
      date = lubridate::ymd(SQLDATE),
      year = lubridate::year(date),
      month = lubridate::month(date),
      event_type = as.character(.data[[grouping_var]])
    )

  reference_category <- validate_reference_class(reference_category, data$event_type, "reference_category")
  reference_hostile <- validate_reference_class(reference_hostile, data$event_type, "reference_hostile")

  class_order <- order_event_classes(data$event_type, reference_category, reference_hostile)
  group_cols <- c("dyad", "year", if (resolution == "monthly") "month")
  column_order <- c(group_cols, paste0("EventClass_", class_order), "total_events")

  counts <- data %>%
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), event_type) %>%
    tidyr::pivot_wider(
      names_from = event_type,
      values_from = n,
      values_fill = 0,
      names_prefix = "EventClass_"
    )

  totals <- data %>%
    dplyr::count(dplyr::across(dplyr::all_of(group_cols)), name = "total_events")

  counts %>%
    dplyr::left_join(totals, by = group_cols) %>%
    dplyr::select(dplyr::any_of(column_order))
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

  skeleton <- if (resolution == "yearly") {
    tidyr::expand_grid(dyad = unique(data$dyad), year = years)
  } else {
    tidyr::expand_grid(dyad = unique(data$dyad), year = years, month = 1:12)
  }

  join_cols <- c("dyad", "year", if (resolution == "monthly") "month")
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
#' row order. This reconstructs the same dyad x period skeleton and
#' recovers the `dyad_id` <-> `dyad` (and undirected `dyad2`) mapping, so
#' posterior draws indexed by `dyad_id` can be joined back to
#' human-readable identifiers. See [extract_theta()].
#'
#' @inheritParams grouped_events_to_dyad_period
#' @param years Integer vector of years covered by the fitted model.
#' @param min_n_events Minimum total events for a dyad to have been
#'   retained by [assemble_stan_data()]; must match the value used there.
#' @return A data frame with `dyad_id`, `time_index`, `dyad`, `dyad2`
#'   (undirected dyad key), `year`, and (if monthly) `month`.
#' @keywords internal
make_dyad_ids <- function(
  data,
  years,
  resolution = c("monthly", "yearly"),
  min_n_events = 1
) {
  resolution <- match.arg(resolution)

  dyad_secondid <- data %>%
    dplyr::group_by(dyad) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      a1 = stringr::str_sub(dyad, 1, 3),
      a2 = stringr::str_sub(dyad, 5, 7),
      dyad2 = dplyr::if_else(a1 < a2, paste(a1, a2, sep = "_"), paste(a2, a1, sep = "_"))
    ) %>%
    dplyr::select(dyad, dyad2)

  data <- fill_dyad_period_skeleton(data, years, resolution)

  drop_dyads <- data %>%
    dplyr::group_by(dyad) %>%
    dplyr::summarise(total_events = sum(total_events), .groups = "drop") %>%
    dplyr::filter(total_events < min_n_events) %>%
    dplyr::pull(dyad)

  data %>%
    dplyr::filter(!(dyad %in% drop_dyads)) %>%
    dplyr::left_join(dyad_secondid, by = "dyad") %>%
    dplyr::group_by(dyad) %>%
    dplyr::mutate(dyad_id = dplyr::cur_group_id(), time_index = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::any_of(c("dyad_id", "time_index", "dyad", "dyad2", "year", "month")))
}
