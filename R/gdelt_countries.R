# GDELT country codes: names, entity types, state-system membership and
# ISO3 harmonization (1.0.0f). The reviewed table is `gdelt_countries`
# (data/, built by data-raw/gdelt_country_codes.R); the year-level periods
# and day-level spells behind gdelt_is_state()/gdelt_country_iso3()/
# gdelt_state_spells() are internal (R/sysdata.rda, same script).

#' Look up GDELT country codes in [gdelt_countries]
#'
#' @param code A character vector of codes (`NA` and `""` allowed).
#' @param caller The calling function, for the warning.
#' @return Integer row positions in [gdelt_countries] (`NA` for `NA`, `""`
#'   and unknown codes). Warns once, naming the unknown codes.
#' @noRd
.gdelt_country_index <- function(code, caller) {
  code <- .as_country_code(code, caller)
  idx <- match(code, bilatr::gdelt_countries$gdelt_code)
  unknown <- unique(code[is.na(idx) & !is.na(code) & code != ""])
  if (length(unknown) > 0) {
    warning(
      caller, ": ", length(unknown), " code(s) not in `gdelt_countries`, returned as NA: ",
      .show_codes(unknown), ".",
      call. = FALSE
    )
  }
  idx
}

#' Check and coerce a country-code argument
#' @noRd
.as_country_code <- function(code, caller) {
  if (is.factor(code)) {
    code <- as.character(code)
  }
  if (!is.character(code) && !(is.logical(code) && all(is.na(code)))) {
    stop(caller, ": `code` must be a character vector.", call. = FALSE)
  }
  as.character(code)
}

#' Format a code vector for a message: at most 10, then "and N more"
#' @noRd
.show_codes <- function(x, n = 10) {
  out <- paste(utils::head(x, n), collapse = ", ")
  if (length(x) > n) {
    out <- paste0(out, sprintf(" and %d more", length(x) - n))
  }
  out
}

#' Check `year` and recycle `code` and `year` to a common length
#'
#' Either may have length 1; otherwise the lengths must match.
#'
#' @return A list: `code` (character) and `year` (integer).
#' @noRd
.recycle_code_year <- function(code, year, caller) {
  if (!is.numeric(year) && !(is.logical(year) && all(is.na(year)))) {
    stop(caller, ": `year` must be NULL or a vector of whole-number years.", call. = FALSE)
  }
  if (!all(is.na(year) | .is_whole_number(as.numeric(year)))) {
    stop(caller, ": `year` must be NULL or a vector of whole-number years.", call. = FALSE)
  }
  n <- max(length(code), length(year))
  if (!length(code) %in% c(1L, n) || !length(year) %in% c(1L, n) ||
      length(code) == 0 || length(year) == 0) {
    if (length(code) == 0 || length(year) == 0) {
      n <- 0L
    } else {
      stop(
        caller, ": `code` (length ", length(code), ") and `year` (length ", length(year),
        ") must have the same length, or one of them length 1.",
        call. = FALSE
      )
    }
  }
  list(code = rep_len(code, n), year = rep_len(as.integer(year), n))
}

#' Resolve (code, year) pairs to the harmonized ISO3 code of a state
#'
#' Works on the unique pairs, so it is fast on long event tables.
#'
#' @param code Character vector of GDELT codes.
#' @param year Integer vector of years, the same length.
#' @return Character vector: the ISO3 code where `code` is a state in
#'   `year`, else `NA`.
#' @noRd
.gdelt_state_iso3 <- function(code, year) {
  key <- paste(code, year, sep = "\r")
  first <- !duplicated(key)
  u_code <- code[first]
  u_year <- year[first]
  iso <- rep(NA_character_, length(u_code))
  periods <- .gdelt_state_periods
  periods <- periods[periods$gdelt_code %in% u_code, , drop = FALSE]
  for (i in seq_len(nrow(periods))) {
    to <- periods$to_year[i]
    hit <- !is.na(u_code) & !is.na(u_year) & u_code == periods$gdelt_code[i] &
      u_year >= periods$from_year[i] & (is.na(to) | u_year <= to)
    iso[hit] <- periods$iso3c[i]
  }
  iso[match(key, key[first])]
}

#' Names of GDELT country codes
#'
#' Returns a readable name for each `Actor1CountryCode`/`Actor2CountryCode`
#' value, from [gdelt_countries]: GDELT's own label, with a few corrections,
#' or an added name for a code that is not in GDELT's lookup.
#'
#' @param code A character vector of GDELT country codes. `NA` and `""`
#'   (GDELT's empty country code, read as `NA` by [read_gdelt()]) give `NA`.
#' @return A character vector of names, the same length as `code`. Codes
#'   not in [gdelt_countries] give `NA`, with one warning naming them.
#' @seealso [gdelt_country_type()], [gdelt_is_state()],
#'   [gdelt_country_iso3()], [filter_gdelt_states()].
#' @examples
#' gdelt_country_name(c("USA", "ROM", "EUR", "WSB", NA))
#' @export
gdelt_country_name <- function(code) {
  idx <- .gdelt_country_index(code, "gdelt_country_name()")
  bilatr::gdelt_countries$name[idx]
}

#' Entity types of GDELT country codes
#'
#' Classifies each GDELT country code as a state, a former state, a
#' non-state territory, a region or other grouping, or something else; see
#' the `entity_type` column of [gdelt_countries] for the definitions. The
#' type is per code, not per year: use [gdelt_is_state()] with `year` for
#' state-system membership at the time of an event.
#'
#' @inheritParams gdelt_country_name
#' @return A character vector, the same length as `code`, of `"state"`,
#'   `"former_state"`, `"territory"`, `"region"` or `"other"`. `NA` for
#'   `NA`, `""` and unknown codes, with one warning naming the unknown
#'   codes.
#' @seealso [gdelt_countries], [gdelt_is_state()].
#' @examples
#' gdelt_country_type(c("USA", "EUR", "WAF", "WSB", "SJM", "GME"))
#' @export
gdelt_country_type <- function(code) {
  idx <- .gdelt_country_index(code, "gdelt_country_type()")
  bilatr::gdelt_countries$entity_type[idx]
}

#' Is a GDELT country code a state?
#'
#' With `year`, tests whether each code is a member of the Correlates of
#' War state system in that year: a member at any point during the year
#' (see [gdelt_countries], "States"). This is the rule
#' [filter_gdelt_states()] applies to events.
#'
#' Without `year`, tests whether the code is a state at any time from 1946
#' on (an `entity_type` of `"state"` or `"former_state"`). That is looser:
#' `gdelt_is_state("SVK")` is `TRUE`, but Slovakia is a state only from
#' 1993, so `gdelt_is_state("SVK", 1992)` is `FALSE`.
#'
#' @inheritParams gdelt_country_name
#' @param year `NULL` (default), or a vector of whole-number years (for
#'   example `SQLDATE %/% 10000`). `code` and `year` are recycled to a
#'   common length: they must have the same length, or one of them length
#'   1.
#' @return A logical vector, of the common length of `code` and `year`. `FALSE` for known
#'   non-state codes and for `NA`/`""`; `NA` for unknown codes (with one
#'   warning naming them) and for `NA` years.
#' @seealso [gdelt_country_iso3()], [filter_gdelt_states()],
#'   [gdelt_state_spells()].
#' @examples
#' gdelt_is_state(c("USA", "EUR", "SVK"))
#' gdelt_is_state("SVK", c(1992, 1993))
#'
#' # Serbia is one state (COW 345) in every year
#' gdelt_is_state("SRB", c(1995, 2004, 2010))
#' @export
gdelt_is_state <- function(code, year = NULL) {
  caller <- "gdelt_is_state()"
  code <- .as_country_code(code, caller)
  if (!is.null(year)) {
    cy <- .recycle_code_year(code, year, caller)
    code <- cy$code
    year <- cy$year
  }
  idx <- .gdelt_country_index(code, caller)
  known <- !is.na(idx)
  missing_code <- is.na(code) | code == ""
  if (is.null(year)) {
    out <- bilatr::gdelt_countries$entity_type[idx] %in% c("state", "former_state")
  } else {
    out <- !is.na(.gdelt_state_iso3(code, year))
    out[known & is.na(year)] <- NA
  }
  out[!known & !missing_code] <- NA
  out[missing_code] <- FALSE
  out
}

#' Harmonized ISO 3166-1 alpha-3 codes for GDELT country codes
#'
#' Maps each GDELT code that is a state to its ISO 3166-1 alpha-3 code,
#' for example `ROM` (CAMEO's code for Romania) to `ROU` and `TMP` to
#' `TLS`. With `year`, the code must be a state in that year, and
#' successor cases are resolved by year: `CZE` is `CSK` (Czechoslovakia)
#' before 1993. See [gdelt_countries], "Harmonization".
#'
#' @inheritParams gdelt_is_state
#' @return A character vector, of the common length of `code` and `year`
#'   (see [gdelt_is_state()]): the ISO3 code,
#'   or `NA` for a code that is not a state (in `year`, if given), for
#'   `NA`/`""`, and for unknown codes (with one warning naming them).
#'   Without `year`, a state code gives the `iso3c` of its row in
#'   [gdelt_countries], the code's latest ISO3 code.
#' @seealso [gdelt_is_state()], [filter_gdelt_states()].
#' @examples
#' gdelt_country_iso3(c("ROM", "TMP", "USA", "EUR"))
#' gdelt_country_iso3("CZE", c(1992, 1993))
#' gdelt_country_iso3("TMP", c(2001, 2002))
#' @export
gdelt_country_iso3 <- function(code, year = NULL) {
  caller <- "gdelt_country_iso3()"
  code <- .as_country_code(code, caller)
  if (!is.null(year)) {
    cy <- .recycle_code_year(code, year, caller)
    code <- cy$code
    year <- cy$year
  }
  idx <- .gdelt_country_index(code, caller)
  if (is.null(year)) {
    is_state <- bilatr::gdelt_countries$entity_type[idx] %in% c("state", "former_state")
    out <- bilatr::gdelt_countries$iso3c[idx]
    out[!is_state] <- NA_character_
    return(out)
  }
  .gdelt_state_iso3(code, year)
}

#' State-system membership spells
#'
#' Returns the Correlates of War state-system membership spells of every
#' state in [gdelt_countries], keyed by the harmonized ISO3 code that
#' [filter_gdelt_states()] writes, one row per spell. This is the
#' `actor_periods` input of [aggregate_events()], which marks the
#' dyad-periods in which both countries are members (`active`).
#'
#' Spells come from COW v2024, which ends on 2024-12-31: a spell open at
#' that date is a current membership and has `end = NA`, taken as
#' continuing through 2025 and after. Spells ending before 1946 are left
#' out; a spell that runs across 1946 keeps its true start date. Two
#' spells of one COW state number are merged when they are contiguous at
#' year level (the second starts no later than the calendar year after the
#' first ends), the same rule [gdelt_is_state()] uses, so the two always
#' agree on which years a country is a state in. Some ISO3 codes have two
#' spells under different COW numbers: `DEU` is West Germany (COW 260)
#' until 1990-10-02 and Germany (255) from 1990-10-03, and `YEM` is North
#' Yemen (678) until 1990-05-21 and Yemen (679) from 1990-05-22.
#'
#' [filter_gdelt_states()] works at year level: it keeps an event when
#' both countries are members at any point in the event's year. So with
#' [aggregate_events()] at a resolution finer than yearly, an event in a
#' country's entry or exit year can fall in a period before the spell
#' starts or after it ends (South Sudan in February 2011, say, before its
#' 9 July independence), which [aggregate_events()] treats as an error.
#' Use `whole_years = TRUE` to widen every spell to whole calendar years
#' and match [filter_gdelt_states()] exactly.
#'
#' Kosovo (COW 347, a member from 2008-02-20) is a state under this rule
#' but has no GDELT code, so it has no spell here.
#'
#' @param whole_years `FALSE` (default) gives COW's exact dates. `TRUE`
#'   moves every `start` to 1 January and every non-missing `end` to 31
#'   December of its year (see Details).
#' @return A tibble with one row per spell, sorted by `iso3c` then
#'   `start`: `iso3c` (character), `cown` (integer, the COW state number),
#'   `start` and `end` (`Date`; `end` is `NA` for a current member).
#' @source Correlates of War Project (2024). State System Membership
#'   (v2024). <https://correlatesofwar.org/data-sets/state-system-membership/>
#' @seealso [aggregate_events()]'s `actor_periods`; [gdelt_countries].
#' @examples
#' spells <- gdelt_state_spells()
#' spells[spells$iso3c %in% c("CSK", "CZE", "SVK", "DEU"), ]
#' @export
gdelt_state_spells <- function(whole_years = FALSE) {
  if (!rlang::is_bool(whole_years)) {
    stop("gdelt_state_spells(): `whole_years` must be TRUE or FALSE.", call. = FALSE)
  }
  out <- tibble::as_tibble(.gdelt_state_spells)
  if (whole_years) {
    out$start <- as.Date(paste0(format(out$start, "%Y"), "-01-01"))
    has_end <- !is.na(out$end)
    out$end[has_end] <- as.Date(paste0(format(out$end[has_end], "%Y"), "-12-31"))
  }
  out
}

#' Keep state-to-state events, with harmonized country codes
#'
#' Keeps the rows of an event table in which both actors are members of
#' the Correlates of War state system in the year of the event (see
#' [gdelt_is_state()] and [gdelt_countries]), and by default rewrites both
#' actor columns to harmonized ISO 3166-1 alpha-3 codes (see
#' [gdelt_country_iso3()]).
#'
#' In order:
#'
#' 1. rows whose date does not parse are dropped (the date column is read
#'    as in `?bilatr_event_data`, e.g. GDELT's `SQLDATE` as `YYYYMMDD`);
#' 2. rows in which either actor's code is not a state in the event's
#'    year are dropped: regions, territories, other codes, `NA`/`""`, codes
#'    not in [gdelt_countries], and states outside their membership years
#'    (e.g. `SVK` in 1992, `TMP` before 2002);
#' 3. with `harmonize = TRUE`, both actor columns are overwritten with the
#'    ISO3 codes, and rows that this turns into self-dyads (two codes for
#'    one state, e.g. `ROM` and `ROU`) are dropped.
#'
#' One `message()` reports the rows kept and dropped, the most frequent
#' dropped codes (counted per actor side; the ten most frequent, then how
#' many others), and the codes harmonization changed.
#'
#' [read_gdelt()]'s `states_only = TRUE` applies the same filter to each
#' chunk as the files are read, so that rows that are dropped never take
#' up memory.
#'
#' @param data A data frame of event-level records, e.g. from
#'   [read_gdelt()].
#' @param actor1,actor2 Names of the two actor country-code columns.
#' @param date Name of the event-date column.
#' @param harmonize `TRUE` (default) rewrites the actor columns with ISO3
#'   codes and drops harmonized self-dyads; `FALSE` keeps GDELT's codes.
#' @return `data`'s kept rows, with the same columns in the same order
#'   (and the same class).
#' @seealso [gdelt_countries]; [gdelt_state_spells()] and
#'   [aggregate_events()]'s `actor_periods` for the dyad-periods outside
#'   membership.
#' @examples
#' events <- tibble::tibble(
#'   Actor1CountryCode = c("USA", "ROM", "EUR", "SVK", "SVK", "ROM"),
#'   Actor2CountryCode = c("CHN", "HUN", "RUS", "CZE", "CZE", "ROU"),
#'   SQLDATE = c(20200101, 19950601, 20200101, 19920301, 19930301, 20000101)
#' )
#' filter_gdelt_states(events)
#' @export
filter_gdelt_states <- function(
  data,
  actor1 = "Actor1CountryCode",
  actor2 = "Actor2CountryCode",
  date = "SQLDATE",
  harmonize = TRUE
) {
  caller <- "filter_gdelt_states()"
  res <- .gdelt_filter_states(data, actor1, actor2, date, harmonize, caller)
  message(.gdelt_states_message(res$stats, caller))
  res$data
}

#' Worker for [filter_gdelt_states()]: filter, and return counts
#'
#' @inheritParams filter_gdelt_states
#' @param caller The calling function, for errors.
#' @return A list: `data` (the kept rows) and `stats` (see
#'   [.gdelt_states_message()]).
#' @noRd
.gdelt_filter_states <- function(data, actor1, actor2, date, harmonize, caller) {
  if (!is.data.frame(data)) {
    stop(caller, ": `data` must be a data frame.", call. = FALSE)
  }
  if (!rlang::is_bool(harmonize)) {
    stop(caller, ": `harmonize` must be TRUE or FALSE.", call. = FALSE)
  }
  for (arg in c("actor1", "actor2", "date")) {
    if (!rlang::is_string(get(arg))) {
      stop(caller, ": `", arg, "` must be a single column name.", call. = FALSE)
    }
  }
  missing_cols <- setdiff(c(actor1, actor2, date), names(data))
  if (length(missing_cols) > 0) {
    stop(
      caller, ": column(s) not found in `data`: ", paste0("`", missing_cols, "`", collapse = ", "),
      ". Set `actor1`, `actor2` and `date` to the actor country-code and event-date columns.",
      call. = FALSE
    )
  }

  a1 <- as.character(data[[actor1]])
  a2 <- as.character(data[[actor2]])
  year <- as.integer(format(.parse_event_date(data[[date]]), "%Y"))
  iso1 <- .gdelt_state_iso3(a1, year)
  iso2 <- .gdelt_state_iso3(a2, year)

  bad_date <- is.na(year)
  non_state <- !bad_date & (is.na(iso1) | is.na(iso2))
  self <- if (harmonize) !bad_date & !non_state & iso1 == iso2 else rep(FALSE, nrow(data))
  keep <- !bad_date & !non_state & !self

  # Dropped codes, per actor side, among rows dropped for a non-state actor.
  dropped_codes <- c(a1[non_state & is.na(iso1)], a2[non_state & is.na(iso2)])
  dropped_codes[is.na(dropped_codes) | dropped_codes == ""] <- "<empty>"

  harmonized <- character(0)
  if (harmonize) {
    changed1 <- keep & a1 != iso1
    changed2 <- keep & a2 != iso2
    harmonized <- c(paste(a1[changed1], iso1[changed1], sep = " -> "),
                    paste(a2[changed2], iso2[changed2], sep = " -> "))
  }

  out <- data[keep, , drop = FALSE]
  if (harmonize) {
    out[[actor1]] <- iso1[keep]
    out[[actor2]] <- iso2[keep]
  }

  list(
    data = out,
    stats = list(
      n_in = nrow(data),
      n_kept = sum(keep),
      n_bad_date = sum(bad_date),
      n_non_state = sum(non_state),
      n_self = sum(self),
      dropped_codes = .count_values(dropped_codes),
      harmonized = .count_values(harmonized)
    )
  )
}

#' Named integer counts of a character vector, most frequent first
#' @noRd
.count_values <- function(x) {
  if (length(x) == 0) {
    return(stats::setNames(integer(0), character(0)))
  }
  tab <- table(x)
  out <- stats::setNames(as.integer(tab), names(tab))
  out[order(-out, names(out), method = "radix")]
}

#' Add two named count vectors
#' @noRd
.add_counts <- function(a, b) {
  keys <- union(names(a), names(b))
  out <- stats::setNames(integer(length(keys)), keys)
  out[names(a)] <- out[names(a)] + a
  out[names(b)] <- out[names(b)] + b
  out[order(-out, names(out), method = "radix")]
}

#' Combine the counts of two state-filter passes (e.g. two chunks)
#' @noRd
.gdelt_merge_state_stats <- function(a, b) {
  if (is.null(a)) {
    return(b)
  }
  list(
    n_in = a$n_in + b$n_in,
    n_kept = a$n_kept + b$n_kept,
    n_bad_date = a$n_bad_date + b$n_bad_date,
    n_non_state = a$n_non_state + b$n_non_state,
    n_self = a$n_self + b$n_self,
    dropped_codes = .add_counts(a$dropped_codes, b$dropped_codes),
    harmonized = .add_counts(a$harmonized, b$harmonized)
  )
}

#' The one message of [filter_gdelt_states()] (and of `states_only`)
#'
#' @param stats Counts from [.gdelt_filter_states()].
#' @param caller The calling function.
#' @return A single string.
#' @noRd
.gdelt_states_message <- function(stats, caller) {
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  n_dropped <- stats$n_in - stats$n_kept
  reasons <- c(
    if (stats$n_non_state > 0) sprintf("%s with an actor that is not a state in the event's year", fmt(stats$n_non_state)),
    if (stats$n_self > 0) sprintf("%s that harmonization made self-dyads", fmt(stats$n_self)),
    if (stats$n_bad_date > 0) sprintf("%s with an unparseable date", fmt(stats$n_bad_date))
  )
  msg <- sprintf(
    "%s: kept %s of %s rows; dropped %s%s.",
    caller, fmt(stats$n_kept), fmt(stats$n_in), fmt(n_dropped),
    if (length(reasons) > 0) paste0(" (", paste(reasons, collapse = "; "), ")") else ""
  )
  codes <- stats$dropped_codes
  if (length(codes) > 0) {
    top <- utils::head(codes, 10)
    rest <- codes[-seq_along(top)]
    msg <- paste0(
      msg, "\nDropped codes (actor sides): ",
      paste(names(top), fmt(top), collapse = ", "),
      if (length(rest) > 0) sprintf(", and %d other code(s) (%s)", length(rest), fmt(sum(rest))) else "",
      "."
    )
  }
  if (length(stats$harmonized) > 0) {
    msg <- paste0(
      msg, "\nHarmonized codes (actor sides): ",
      paste0(names(stats$harmonized), " (", fmt(stats$harmonized), ")", collapse = ", "),
      "."
    )
  }
  msg
}

#' Fill empty actor country codes from the actor-code prefix
#'
#' GDELT derives `Actor1CountryCode`/`Actor2CountryCode` from the leading
#' three characters of `Actor1Code`/`Actor2Code`, and leaves the country
#' code empty when its country table lacks that segment (e.g. Romania's
#' `ROU...`). Where the country code is `NA` or `""` and the prefix is a
#' code in [gdelt_countries], the prefix is written in.
#'
#' @param chunk A data frame with the four actor code columns.
#' @return `chunk`, with the country codes filled and two extra columns,
#'   `.filled1`/`.filled2`, holding the filled code (or `NA`), for the
#'   caller to count and drop.
#' @noRd
.gdelt_recover_country_codes <- function(chunk) {
  codes <- bilatr::gdelt_countries$gdelt_code
  fill <- function(country, actor) {
    prefix <- substr(actor, 1, 3)
    (is.na(country) | country == "") & !is.na(prefix) & prefix %in% codes
  }
  f1 <- fill(chunk$Actor1CountryCode, chunk$Actor1Code)
  f2 <- fill(chunk$Actor2CountryCode, chunk$Actor2Code)
  chunk$.filled1 <- ifelse(f1, substr(chunk$Actor1Code, 1, 3), NA_character_)
  chunk$.filled2 <- ifelse(f2, substr(chunk$Actor2Code, 1, 3), NA_character_)
  chunk$Actor1CountryCode[f1] <- chunk$.filled1[f1]
  chunk$Actor2CountryCode[f2] <- chunk$.filled2[f2]
  chunk
}

#' The message for `recover_country_codes = TRUE`
#'
#' @param filled Named integer counts of filled actor sides, by code, in
#'   the rows returned.
#' @param caller The calling function.
#' @return A single string.
#' @noRd
.gdelt_recover_message <- function(filled, caller) {
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  if (length(filled) == 0) {
    return(paste0(
      caller, ": `recover_country_codes` filled no actor country codes in the rows kept."
    ))
  }
  sprintf(
    "%s: `recover_country_codes` filled %s empty actor country code(s) from the actor-code prefix, in the rows kept (actor sides): %s.",
    caller, fmt(sum(filled)), paste(names(filled), fmt(filled), collapse = ", ")
  )
}
