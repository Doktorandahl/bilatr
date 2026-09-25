#' The bilatr event-data contract
#'
#' Every function in this package that consumes event-level data (
#' [grouped_events_to_dyad_period()], [assemble_stan_data()]) expects a
#' data frame in the format documented here, and validates it with
#' [validate_bilatr_events()] before doing anything else. This topic is the
#' single place the contract is spelled out; every other function's
#' `@param data` links here instead of repeating it.
#'
#' | Role | Default column | Requirement |
#' | --- | --- | --- |
#' | Actor 1 (sender, for directed dyads) | `Actor1CountryCode` | character, factor or integer (used as character); no `NA`; no empty string; must not contain `"_"` (the dyad-key separator) |
#' | Actor 2 (target) | `Actor2CountryCode` | same; and `actor1 != actor2` in every row (bilatr models *bilateral* relations) |
#' | Event date | `SQLDATE` | a `Date`; a `POSIXt` (converted in its own time zone); a whole-number `YYYYMMDD` (integer, numeric, or 8-digit character/factor, e.g. `20200101`); or an ISO `YYYY-MM-DD` character/factor (e.g. `"2020-01-01"`); no `NA`; anything else (short/long digit runs, trailing text, an invalid calendar date) fails to parse and is reported as `NA` |
#' | Event class | `grouping_var` (no default) | any atomic type, used as character; no `NA` |
#'
#' Actor codes may have any length and any characters apart from `"_"`:
#' ISO3 (`"USA"`), ISO2 (`"US"`), COW numeric codes (`"2"`, `"365"`), and
#' free-text actor labels all work. Every column not named above is
#' ignored by these functions and never modified -- in particular, a
#' user's own `dyad`/`date`/`year`/`month`/`event_type` columns (if
#' present) are left exactly as they are; the internal aggregate these
#' functions build uses its own column names.
#'
#' [validate_bilatr_events()] runs against **every** row of `data`, including
#' rows outside the analysis `years` window (see [assemble_stan_data()]) --
#' so a data problem in an out-of-window row (e.g. an unparseable date) is
#' still caught. Windowing happens afterwards: rows outside `years` are
#' dropped before the set of observed event classes is determined, and a
#' `message()` reports how many rows were dropped.
#'
#' Two helpers produce data in this format: [read_gdelt()]/[download_gdelt()]
#' (GDELT) and [recode_cameo()] (attaches a CAMEO event-class column such
#' as `QuadClass`/`PentaClass`/`ModifiedRootCode` to an existing event
#' table).
#'
#' @examples
#' # A 10-row event table with COW-style numeric actor codes and a Date
#' # column -- neither ISO3 codes nor a GDELT-style integer SQLDATE are
#' # required.
#' events <- tibble::tibble(
#'   sender = c("2", "2", "365", "365", "710", "710", "2", "365", "710", "2"),
#'   target = c("365", "710", "2", "710", "2", "365", "710", "2", "365", "365"),
#'   when = as.Date("2020-01-01") + 0:9,
#'   action = c(1, 1, 2, 2, 3, 3, 4, 4, 1, 2)
#' )
#'
#' stan_data <- assemble_stan_data(
#'   events,
#'   years = 2020,
#'   resolution = "monthly",
#'   grouping_var = "action",
#'   actor1 = "sender",
#'   actor2 = "target",
#'   date = "when"
#' )
#'
#' @name bilatr_event_data
NULL

#' Parse an event date column into a `Date` vector
#'
#' A `Date` column is returned as-is. A `POSIXt` column is converted with
#' `as.Date(format(x, "\%Y-\%m-\%d"))`, taking the calendar date in the
#' object's own time zone rather than converting to UTC. Numeric input
#' must be a whole number; it is formatted as an integer
#' (`formatC(x, format = "d", big.mark = "")`, so no value goes through
#' scientific notation) and then parsed as character. Character input
#' (including numeric-as-character and factor levels) is accepted only
#' if it matches exactly `^[0-9]{8}$` (parsed as `\%Y\%m\%d`) or exactly
#' `^[0-9]{4}-[0-9]{2}-[0-9]{2}$` (parsed as `\%Y-\%m-\%d`); anything else,
#' including a non-whole numeric value, an 8-digit string with trailing
#' text, or an invalid calendar date (e.g. month 13), becomes `NA`, which
#' is what [validate_bilatr_events()]'s date rule flags. This is
#' deliberately stricter than plain `as.Date(x, format = "\%Y\%m\%d")`
#' (`strptime` underneath), which allows one- or two-digit months/days
#' and silently ignores trailing text.
#'
#' @param x A `Date`, `POSIXt`, numeric, character, or factor vector.
#' @return A `Date` vector, the same length as `x`.
#' @keywords internal
.parse_event_date <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }
  if (inherits(x, "POSIXt")) {
    return(as.Date(format(x, "%Y-%m-%d")))
  }
  if (is.numeric(x)) {
    whole <- !is.na(x) & x == trunc(x)
    chr <- rep(NA_character_, length(x))
    chr[whole] <- formatC(x[whole], format = "d", big.mark = "")
    return(.parse_event_date_strict(chr))
  }
  .parse_event_date_strict(as.character(x))
}

#' Parse a character vector as a strict `YYYYMMDD` or `YYYY-MM-DD` date
#'
#' The shared strict-parsing core of [.parse_event_date()]: values not
#' matching one of the two accepted forms exactly become `NA` before
#' `as.Date()` is even called, so `strptime`'s leniency (short digit
#' runs, trailing text) never applies.
#'
#' @param x A character vector.
#' @return A `Date` vector, the same length as `x`.
#' @keywords internal
.parse_event_date_strict <- function(x) {
  out <- as.Date(rep(NA_character_, length(x)))
  ymd <- !is.na(x) & grepl("^[0-9]{8}$", x)
  iso <- !is.na(x) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", x)
  out[ymd] <- as.Date(x[ymd], format = "%Y%m%d")
  out[iso] <- as.Date(x[iso], format = "%Y-%m-%d")
  out
}

#' Build the slim four-column working table from a validated event table
#'
#' Shared by [grouped_events_to_dyad_period()] and [assemble_stan_data()]
#' so that, between the two of them, a `data` table is validated exactly
#' once and its date column parsed exactly once per call (0.9.1, audit
#' 0f) -- previously `assemble_stan_data()` validated `data` itself and
#' then again inside [grouped_events_to_dyad_period()], and parsed dates
#' up to four times. Callers must validate `data` (via
#' [validate_bilatr_events()]) before calling this; it does no validation
#' of its own.
#'
#' @param data A data frame; see `?bilatr_event_data`.
#' @param grouping_var,actor1,actor2,date Column names; see
#'   [validate_bilatr_events()].
#' @return A tibble with `.actor1`, `.actor2`, `.date` (parsed), `.class`
#'   (as character), `.year`, and `.month`.
#' @keywords internal
.slim_event_table <- function(data, grouping_var, actor1, actor2, date) {
  slim <- tibble::tibble(
    .actor1 = as.character(data[[actor1]]),
    .actor2 = as.character(data[[actor2]]),
    .date = .parse_event_date(data[[date]]),
    .class = as.character(data[[grouping_var]])
  )
  slim$.year <- as.integer(format(slim$.date, "%Y"))
  slim$.month <- as.integer(format(slim$.date, "%m"))
  slim
}

#' Deterministically order an unordered pair of actor codes
#'
#' Ranks `a`/`b` against a radix-sorted level set rather than comparing
#' them directly with `pmin()`/`pmax()`/`<`, which are locale-dependent
#' for character input (the same concern [order_event_classes()]'s
#' `locale = "C"` and [assemble_stan_data()]'s country index address for
#' the action-class and country indices). This is the single place an
#' unordered actor pair is built, used both for undirected dyads'
#' side A/B assignment and for the undirected dyad key `dyad2`.
#'
#' @param a,b Character vectors of actor codes, same length.
#' @return A list with `actor_a` (whichever of `a`/`b` sorts first in
#'   C-locale byte order, element-wise) and `actor_b` (the other).
#' @keywords internal
.order_pair_c_locale <- function(a, b) {
  lev <- sort(unique(c(a, b)), method = "radix")
  a_first <- match(a, lev) < match(b, lev)
  list(
    actor_a = dplyr::if_else(a_first, a, b),
    actor_b = dplyr::if_else(a_first, b, a)
  )
}

#' Validate a data frame against the bilatr event-data contract
#'
#' Checks `data` against every rule in `?bilatr_event_data`, collecting
#' every broken rule rather than stopping at the first, and `stop()`s once
#' with a bulleted message naming each broken rule, the number of
#' offending rows, and up to five example row numbers. Called first thing
#' by [grouped_events_to_dyad_period()] and [assemble_stan_data()].
#'
#' @param data A data frame to validate.
#' @param grouping_var Name of the event-class column, or `NULL` to skip
#'   that rule (e.g. when validating a table before choosing/deriving a
#'   class column).
#' @param actor1,actor2 Names of the sender/target actor columns.
#' @param date Name of the event-date column.
#' @return `data`, invisibly and unmodified, if every rule passes.
#' @export
validate_bilatr_events <- function(
  data,
  grouping_var = NULL,
  actor1 = "Actor1CountryCode",
  actor2 = "Actor2CountryCode",
  date = "SQLDATE"
) {
  required_cols <- c(actor1 = actor1, actor2 = actor2, date = date)
  if (!is.null(grouping_var)) {
    required_cols <- c(required_cols, grouping_var = grouping_var)
  }
  missing_cols <- required_cols[!required_cols %in% names(data)]
  if (length(missing_cols) > 0) {
    stop(
      "validate_bilatr_events(): missing column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  n <- nrow(data)
  a1 <- as.character(data[[actor1]])
  a2 <- as.character(data[[actor2]])
  event_date <- .parse_event_date(data[[date]])

  example_rows <- function(idx) paste(utils::head(which(idx), 5), collapse = ", ")

  bullets <- character(0)
  add_rule <- function(label, idx) {
    if (any(idx)) {
      bullets <<- c(bullets, sprintf(
        "- %s (%d row(s); e.g. row(s) %s)", label, sum(idx), example_rows(idx)
      ))
    }
  }

  add_rule(sprintf("`%s` has missing (NA) values", actor1), is.na(a1))
  add_rule(sprintf("`%s` has empty-string values", actor1), !is.na(a1) & a1 == "")
  add_rule(sprintf("`%s` contains \"_\", the dyad-key separator", actor1), !is.na(a1) & grepl("_", a1, fixed = TRUE))
  add_rule(
    sprintf("`%s` has leading/trailing whitespace (e.g. \" USA\" and \"USA\" would otherwise become two actors)", actor1),
    !is.na(a1) & grepl("^\\s|\\s$", a1)
  )
  add_rule(sprintf("`%s` has missing (NA) values", actor2), is.na(a2))
  add_rule(sprintf("`%s` has empty-string values", actor2), !is.na(a2) & a2 == "")
  add_rule(sprintf("`%s` contains \"_\", the dyad-key separator", actor2), !is.na(a2) & grepl("_", a2, fixed = TRUE))
  add_rule(
    sprintf("`%s` has leading/trailing whitespace (e.g. \" USA\" and \"USA\" would otherwise become two actors)", actor2),
    !is.na(a2) & grepl("^\\s|\\s$", a2)
  )
  add_rule(
    sprintf("`%s` equals `%s` (a self-dyad; bilatr models bilateral relations)", actor1, actor2),
    !is.na(a1) & !is.na(a2) & a1 == a2
  )
  add_rule(sprintf("`%s` could not be parsed as a date", date), is.na(event_date))

  if (!is.null(grouping_var)) {
    class_na <- is.na(data[[grouping_var]])
    if (any(class_na)) {
      suggestion <- sprintf(
        "`%s` has missing (NA) values; drop them first with dplyr::filter(!is.na(%s))",
        grouping_var, grouping_var
      )
      if (identical(grouping_var, "ModifiedRootCode")) {
        suggestion <- paste0(
          suggestion,
          ". Bare two-digit CAMEO root codes (e.g. \"01\", \"04\") have no ",
          "ModifiedRootCode (see ?recode_cameo)."
        )
      }
      add_rule(suggestion, class_na)
    }
  }

  if (length(bullets) > 0) {
    stop(
      "validate_bilatr_events(): the following rule(s) were broken:\n",
      paste(bullets, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(data)
}
