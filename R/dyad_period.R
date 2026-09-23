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
#'   the class ordering (the model's scale/sign-reference, neutral
#'   action: `stable`/`ou` build `alpha[1]` positive by construction --
#'   see [order_event_classes()]). If `NULL` or not present in the data,
#'   ignored with a warning. All other action classes' discrimination is
#'   freely estimated.
#' @param years Integer vector of years, or `NULL` (default). Used
#'   **only** for the `w_send` computation described below -- restricts
#'   the event rows pooled into `w_send` to `year %in% years`, matching
#'   the analysis window [assemble_stan_data()] restricts `Y`'s counts to
#'   (via [fill_dyad_period_skeleton()]'s right join, downstream). The
#'   main dyad-period aggregation returned by this function is
#'   unaffected -- it is still built from every row of `data`, exactly as
#'   before; only `w_send` narrows. `NULL` pools every year in `data`,
#'   preserving this function's behaviour for a direct caller that isn't
#'   going through [assemble_stan_data()] (0.7.0's original behaviour).
#' @return A data frame with columns `dyad`, `year` (and `month` if
#'   `resolution = "monthly"`), one `EventClass_<value>` column per
#'   observed class (ordered per `reference_category`), and `total_events`.
#'   Also carries a `w_send` attribute (0.7.0; consumed by
#'   [assemble_stan_data()] for the experimental `stable_gamma` Stan
#'   variant -- see `R/model_registry.R`): a tibble with one row per
#'   dyad, `ctry_a_code`/`ctry_b_code` (side A/B's 3-letter country code,
#'   `str_sub(dyad, 1, 3)`/`str_sub(dyad, 5, 7)`, the same convention
#'   [make_dyad_ids()] uses for its own `dyad2`) and `w_send` (the share
#'   of that pair's events, POOLED over `years` if supplied (else every
#'   period in `data`) -- not per dyad-period, since the country offset
#'   this feeds is dyad-constant and reused across `t`; a per-period
#'   weight would cost an `A x D x T` object in [assemble_stan_data()]'s
#'   output for nothing -- with side A as `Actor1CountryCode`). For
#'   directed data this is exactly `1` by construction (side A IS
#'   `Actor1CountryCode` in the directed dyad key), asserted internally
#'   below as a free check on the dyad-key convention (it would silently
#'   invert if a future change swapped the `paste()` argument order
#'   building `dyad` above). A row with `NA` `Actor1CountryCode` is
#'   dropped from this computation (with a `warning()` naming how many)
#'   rather than propagating `NA` into `w_send` -- an `NA` reaching
#'   `stan_data$w_send` would otherwise surface as an opaque CmdStan data
#'   error much later.
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
  years = NULL
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

  class_order <- order_event_classes(data$event_type, reference_category)
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

  result <- counts %>%
    dplyr::left_join(totals, by = group_cols) %>%
    dplyr::select(dplyr::any_of(column_order))

  # w_send (0.7.0; see @return): per-dyad sender share, pooled over
  # `years` (0.7.1 -- previously every period in `data` regardless of the
  # analysis window; see @param years) not per dyad-period, computed from
  # the raw event rows rather than from `result`'s per-period counts,
  # since it needs Actor1CountryCode per event, not per action class.
  # Side A is str_sub(dyad, 1, 3), matching make_dyad_ids()'s own
  # convention for dyad2 -- both assume 3-letter country codes joined by
  # a single "_", the same assumption `dyad` itself already relies on
  # above.
  w_send_data <- if (is.null(years)) data else dplyr::filter(data, year %in% years)

  w_send_data <- dplyr::mutate(
    w_send_data,
    ctry_a_code = stringr::str_sub(dyad, 1, 3),
    ctry_b_code = stringr::str_sub(dyad, 5, 7),
    is_side_a_sender = Actor1CountryCode == ctry_a_code
  )

  # NA Actor1CountryCode -> NA is_side_a_sender (comparing NA to anything
  # is NA in R, regardless of what paste() upstream turned a missing
  # actor into inside `dyad` itself) -> mean() would silently return NA
  # for that whole dyad, surfacing as an opaque CmdStan data error only
  # once `w_send` reaches stan_data (see @return). Dropped here instead,
  # with a warning naming how many rows/dyads were affected. If an
  # affected dyad has EVERY row dropped this way, it disappears from the
  # w_send tibble entirely -- assemble_stan_data()'s own guard (see
  # R/stan_data.R, 0.7.1's country_info NA check) then stop()s naming
  # that dyad, rather than silently shipping an NA into ctry_a/w_send.
  na_rows <- is.na(w_send_data$is_side_a_sender)
  if (any(na_rows)) {
    n_dyads_affected <- length(unique(w_send_data$dyad[na_rows]))
    warning(
      sum(na_rows), " event row(s) (across ", n_dyads_affected,
      " dyad(s)) had NA Actor1CountryCode and were dropped from the ",
      "w_send computation.",
      call. = FALSE
    )
    w_send_data <- w_send_data[!na_rows, , drop = FALSE]
  }

  w_send <- w_send_data %>%
    dplyr::group_by(dyad) %>%
    dplyr::summarise(
      ctry_a_code = dplyr::first(ctry_a_code),
      ctry_b_code = dplyr::first(ctry_b_code),
      w_send = mean(is_side_a_sender),
      .groups = "drop"
    )

  if (directed && !all(w_send$w_send == 1)) {
    stop(
      "Internal invariant violated: computed w_send != 1 for directed ",
      "data. Side A (str_sub(dyad, 1, 3)) should always equal ",
      "Actor1CountryCode for a directed dyad key -- this points to a bug ",
      "in how `dyad` is built above (e.g. the paste() argument order), ",
      "not a data issue.",
      call. = FALSE
    )
  }

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
