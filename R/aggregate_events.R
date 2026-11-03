#' Aggregate event data to a dyad-period table of class counts
#'
#' Collapses event-level data (see `?bilatr_event_data`; as produced by
#' e.g. [read_gdelt()]/[download_gdelt()] and [recode_cameo()]) to a
#' dyad x time-period table of event-class counts in the format
#' documented at `?bilatr_dyad_period_data`, ready for
#' [prepare_stan_data()].
#'
#' Validates `data` with [validate_bilatr_events()] first, against every
#' row, including rows outside `start`/`end`, so a data problem in an
#' out-of-window row (e.g. an unparseable date) is still caught. Only
#' after validation does it drop rows outside the window and determine
#' the set of observed event classes. A `message()` reports how many rows
#' were dropped, and names any class present in `data` with zero
#' in-window events (it gets no column). If no rows fall within the
#' window, it errors and reports the data's date range; an empty `data`
#' errors too.
#'
#' Every dyad with at least one in-window event gets a row for every
#' period from the period containing `start` to the period containing
#' `end`, with zero counts where it had no events. So `start = "1991-01-01",
#' end = "2025-12-31"` covers every year (or month) of 1991-2025.
#'
#' @section Time IDs:
#' The time column holds absolute period IDs, named by resolution:
#'
#' * `"yearly"`: `year`, the calendar year;
#' * `"monthly"`: `month_id = year * 12 + month - 1`;
#' * `"weekly"`: `week_id = (date - origin) %/% 7`, whole weeks since the
#'   origin date, so weeks run from the origin's weekday (1946-01-01, the
#'   default origin, is a Tuesday);
#' * `"daily"`: `day_id = date - origin`, whole days since the origin.
#'
#' The origin is 1946-01-01, so weekly and daily IDs are comparable across
#' datasets. If the window starts earlier (through `start` or the data's
#' own earliest date), the origin moves to 1 January of that year instead,
#' with a warning that the IDs are then not comparable with tables built
#' on the default origin. Yearly and monthly IDs need no origin.
#'
#' @section Attributes:
#' The returned tibble carries these attributes, which
#' [prepare_stan_data()] picks up by default:
#'
#' * `time_labels`: a tibble keyed by the time column, giving calendar
#'   labels for every period in the range (yearly: `year`; monthly:
#'   `month_id`, `year`, `month`; weekly: `week_id`, `week_start` (a
#'   `Date`); daily: `day_id`, `date`);
#' * `time_col`: the name of the time column;
#' * `time_origin`: the origin `Date` (weekly/daily), or `NA`;
#' * `resolution`, `grouping_var`, `directed`: the arguments;
#' * `action_cols`: the class columns, in output order. With
#'   `action_cols = NULL`, [prepare_stan_data()] uses them silently, as
#'   long as every one is still a column.
#'
#' These are conveniences. Some dplyr verbs drop them (e.g. `summarise()`
#' does; `filter()`, `select()`, `mutate()` and `arrange()` keep them).
#' [prepare_stan_data()] works without them, given its `time`,
#' `action_cols` (and, for undirected data, `directed`) arguments.
#'
#' @section Size:
#' The Stan data's `Y` is a dense `D x T x A` integer array (dyads x
#' periods x classes), about `D * T * A * 4` bytes, and the dyad-period
#' table itself has `D * T` rows. Finer resolutions multiply `T`: 16,750
#' dyads over 35 years with 16 classes is about 37 MB at yearly
#' resolution, 450 MB monthly, and about 13.7 GB daily. Daily resolution
#' is only practical for few dyads or short windows.
#'
#' @param data A data frame of event-level records; see
#'   `?bilatr_event_data`.
#' @param grouping_var Name of the event-class column to count (e.g.
#'   `"QuadClass"`, `"PentaClass"`, `"ModifiedRootCode"`).
#' @param resolution One of `"yearly"` (default), `"monthly"`, `"weekly"`
#'   or `"daily"`; see the Time IDs section.
#' @param directed If `TRUE` (default), dyads are directed (`actor1` is
#'   the sender, `actor2` the target). If `FALSE`, actor order is ignored:
#'   each pair becomes one dyad, with `actor1`/`actor2` in a deterministic,
#'   locale-independent (C-locale byte) order.
#' @param start,end Window bounds: a `Date`, or a string `as.Date()`
#'   parses (e.g. `"1991-01-01"`). Rows dated outside `[start, end]` are
#'   dropped. `NULL` (default) means the data's own earliest/latest date.
#'   Giving both with `start > end` is an error; a one-sided window that
#'   misses the data entirely gets the empty-window error.
#' @param class_order Optional order of the class columns: a permutation
#'   of the in-window classes (values of `grouping_var`, as character).
#'   `NULL` (default) sorts them in C-locale numeric order (`"2"` before
#'   `"10"`). The order matters: [prepare_stan_data()] keeps it, after
#'   moving its `reference_category` to the front.
#' @param actor1,actor2 Names of the sender/target actor columns; see
#'   `?bilatr_event_data`.
#' @param date Name of the event-date column; see `?bilatr_event_data`.
#' @return A tibble with exactly the columns `dyad`, `actor1`, `actor2`,
#'   the time column (see Time IDs), `w_send`, and one integer count
#'   column per event class, named by the raw class value. Sorted by
#'   `dyad` (C-locale) then time. `w_send` is per dyad, repeated on every
#'   row: exactly `1` when `directed = TRUE`; when undirected, the share of
#'   the pair's in-window events sent by `actor1`. See the Attributes
#'   section for the attributes.
#' @examples
#' events <- tibble::tibble(
#'   Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN"),
#'   Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA"),
#'   SQLDATE = c(20180105L, 20180610L, 20190302L, 20190401L, 20200101L),
#'   PentaClass = c(0, 1, 0, 1, 3)
#' )
#' aggregate_events(events, grouping_var = "PentaClass")
#' aggregate_events(
#'   events,
#'   grouping_var = "PentaClass", resolution = "monthly", directed = FALSE,
#'   start = "2018-01-01", end = "2019-12-31"
#' )
#' @export
aggregate_events <- function(
  data,
  grouping_var,
  resolution = c("yearly", "monthly", "weekly", "daily"),
  directed = TRUE,
  start = NULL,
  end = NULL,
  class_order = NULL,
  actor1 = "Actor1CountryCode",
  actor2 = "Actor2CountryCode",
  date = "SQLDATE"
) {
  resolution <- match.arg(resolution)
  fn <- "aggregate_events()"

  validate_bilatr_events(data, grouping_var, actor1, actor2, date)
  slim <- .slim_event_table(data, grouping_var, actor1, actor2, date)

  if (nrow(slim) == 0) {
    stop(fn, ": `data` is empty; nothing to aggregate.", call. = FALSE)
  }
  data_range <- range(slim$.date)
  start <- .as_window_date(start, "start")
  end <- .as_window_date(end, "end")
  # Only when both are given: a one-sided window beyond the data is an
  # empty window (reported below with the data's range), not a bad pair.
  if (!is.null(start) && !is.null(end) && start > end) {
    stop(sprintf("%s: `start` (%s) is after `end` (%s).", fn, start, end), call. = FALSE)
  }
  start <- start %||% data_range[1]
  end <- end %||% data_range[2]

  # Window first, as assemble_stan_data() does: the class set is built
  # from in-window rows only, but the dropped-class message compares
  # against every class in `data`.
  all_classes <- unique(slim$.class)
  in_window <- slim$.date >= start & slim$.date <= end
  n_dropped <- sum(!in_window)
  if (n_dropped > 0) {
    message(sprintf("%s: dropping %d row(s) outside `start`/`end`.", fn, n_dropped))
  }
  slim <- slim[in_window, , drop = FALSE]
  if (nrow(slim) == 0) {
    stop(sprintf(
      "%s: no events fall within `start`/`end` (%s to %s); the data's dates span %s to %s.",
      fn, start, end, data_range[1], data_range[2]
    ), call. = FALSE)
  }
  zero_window_classes <- setdiff(all_classes, unique(slim$.class))
  if (length(zero_window_classes) > 0) {
    message(
      fn, ": the following class(es) of `", grouping_var,
      "` have zero in-window events and get no column: ",
      paste(sort(zero_window_classes), collapse = ", ")
    )
  }

  classes <- .resolve_class_order(class_order, slim$.class)
  time_col <- c(yearly = "year", monthly = "month_id", weekly = "week_id", daily = "day_id")[[resolution]]
  clashing <- intersect(classes, c("dyad", "actor1", "actor2", time_col, "w_send"))
  if (any(classes == "") || length(clashing) > 0) {
    stop(
      fn, ": event classes become column names, so a class may not be empty or named ",
      "dyad/actor1/actor2/", time_col, "/w_send. Offending class(es): ",
      paste0("\"", c(classes[classes == ""], clashing), "\"", collapse = ", "),
      call. = FALSE
    )
  }

  origin <- .time_origin(resolution, min(start, slim$.date))
  slim$.time <- .time_id(slim$.date, resolution, origin, slim$.year, slim$.month)
  time_range <- seq.int(.time_id(start, resolution, origin), .time_id(end, resolution, origin))

  if (directed) {
    slim$actor1 <- slim$.actor1
    slim$actor2 <- slim$.actor2
  } else {
    pair <- .order_pair_c_locale(slim$.actor1, slim$.actor2)
    slim$actor1 <- pair$actor_a
    slim$actor2 <- pair$actor_b
  }
  slim$dyad <- paste(slim$actor1, slim$actor2, sep = "_")

  # Per-dyad actors and sender share, computed exactly as
  # .dyad_period_from_slim() does (1 by construction when directed).
  dyad_info <- slim %>%
    dplyr::mutate(is_side_a_sender = .data$.actor1 == .data$actor1) %>%
    dplyr::group_by(.data$dyad) %>%
    dplyr::summarise(
      actor1 = dplyr::first(.data$actor1),
      actor2 = dplyr::first(.data$actor2),
      w_send = mean(.data$is_side_a_sender),
      .groups = "drop"
    )

  counts <- slim %>%
    dplyr::count(.data$dyad, .data$.time, .data$.class) %>%
    tidyr::pivot_wider(names_from = ".class", values_from = "n", values_fill = 0L)

  result <- tidyr::expand_grid(dyad_info, .time = time_range) %>%
    dplyr::left_join(counts, by = c("dyad", ".time")) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(classes), ~ tidyr::replace_na(.x, 0L))) %>%
    dplyr::arrange(.data$dyad, .data$.time) %>%
    dplyr::select("dyad", "actor1", "actor2", !!time_col := ".time", "w_send", dplyr::all_of(classes))

  attr(result, "time_labels") <- .time_labels(time_range, resolution, origin)
  attr(result, "time_col") <- time_col
  attr(result, "time_origin") <- origin
  attr(result, "resolution") <- resolution
  attr(result, "grouping_var") <- grouping_var
  attr(result, "directed") <- directed
  attr(result, "action_cols") <- classes
  result
}

#' Parse a `start`/`end` window bound
#'
#' @param x `NULL`, a `Date`, or a string `as.Date()` parses.
#' @param arg Argument name, for the error message.
#' @return `NULL` or a length-one `Date`.
#' @keywords internal
.as_window_date <- function(x, arg) {
  if (is.null(x)) {
    return(NULL)
  }
  out <- if (length(x) != 1) {
    NA
  } else if (inherits(x, "Date")) {
    x
  } else if (is.character(x)) {
    tryCatch(as.Date(x), error = function(e) as.Date(NA))
  } else {
    NA
  }
  if (is.na(out)) {
    stop(
      "aggregate_events(): `", arg, "` must be NULL, a single Date, or a single ",
      "date string such as \"1991-01-01\".",
      call. = FALSE
    )
  }
  out
}

#' Check or default the order of event-class columns
#'
#' @param class_order `NULL`, or a permutation of the in-window classes.
#' @param classes The in-window class values (character).
#' @return The class order, as character.
#' @keywords internal
.resolve_class_order <- function(class_order, classes) {
  observed <- unique(classes)
  if (is.null(class_order)) {
    return(order_event_classes(observed, NULL))
  }
  class_order <- as.character(class_order)
  missing <- setdiff(observed, class_order)
  extra <- setdiff(class_order, observed)
  duplicated_classes <- unique(class_order[duplicated(class_order)])
  if (length(missing) > 0 || length(extra) > 0 || length(duplicated_classes) > 0) {
    show <- function(x) if (length(x) == 0) "none" else paste(x, collapse = ", ")
    stop(
      "aggregate_events(): `class_order` must be a permutation of the in-window classes.\n",
      "- missing: ", show(order_event_classes(missing)), "\n",
      "- not in the (in-window) data: ", show(extra), "\n",
      "- duplicated: ", show(duplicated_classes),
      call. = FALSE
    )
  }
  class_order
}

#' The origin date for weekly/daily time IDs
#'
#' 1946-01-01, unless `earliest` is before that, in which case 1 January
#' of `earliest`'s year, with a warning. `NA` for yearly/monthly.
#'
#' @param resolution The matched resolution.
#' @param earliest The earliest date in the window (a `Date`).
#' @return A `Date` (possibly `NA`).
#' @keywords internal
.time_origin <- function(resolution, earliest) {
  if (resolution %in% c("yearly", "monthly")) {
    return(as.Date(NA))
  }
  default <- as.Date("1946-01-01")
  if (earliest >= default) {
    return(default)
  }
  origin <- as.Date(sprintf("%04d-01-01", as.integer(format(earliest, "%Y"))))
  warning(
    "aggregate_events(): the window starts before 1946-01-01, so ", resolution,
    " time IDs count from ", origin, " instead; they are not comparable with ",
    "tables built on the default 1946-01-01 origin.",
    call. = FALSE
  )
  origin
}

#' Convert dates to absolute period IDs
#'
#' @param date A `Date` vector.
#' @param resolution The matched resolution.
#' @param origin The origin `Date` (weekly/daily only).
#' @param year,month Optional precomputed integer year/month of `date`,
#'   to avoid formatting a long date vector again.
#' @return An integer vector of period IDs; see [aggregate_events()].
#' @keywords internal
.time_id <- function(date, resolution, origin, year = NULL, month = NULL) {
  switch(resolution,
    yearly = year %||% as.integer(format(date, "%Y")),
    monthly = {
      year <- year %||% as.integer(format(date, "%Y"))
      month <- month %||% as.integer(format(date, "%m"))
      year * 12L + month - 1L
    },
    weekly = as.integer(unclass(date) - unclass(origin)) %/% 7L,
    daily = as.integer(unclass(date) - unclass(origin))
  )
}

#' Calendar labels for a range of period IDs
#'
#' @param time_range Integer vector of period IDs.
#' @param resolution The matched resolution.
#' @param origin The origin `Date` (weekly/daily only).
#' @return A tibble keyed by the time column; see [aggregate_events()].
#' @keywords internal
.time_labels <- function(time_range, resolution, origin) {
  switch(resolution,
    yearly = tibble::tibble(year = time_range),
    monthly = tibble::tibble(
      month_id = time_range,
      year = time_range %/% 12L,
      month = time_range %% 12L + 1L
    ),
    weekly = tibble::tibble(week_id = time_range, week_start = origin + time_range * 7L),
    daily = tibble::tibble(day_id = time_range, date = origin + time_range)
  )
}
