#' GDELT 1.0 event file layout
#'
#' GDELT 1.0's event archive (`https://data.gdeltproject.org/events/`)
#' ships three file types, chosen by date, verified 2026-09-25 against the
#' live index at that URL and the GDELT 1.0 Data Format Codebook
#' (`https://data.gdeltproject.org/documentation/GDELT-Data_Format_Codebook.pdf`):
#'
#' - yearly `YYYY.zip`, 1979-01-01 through 2005-12-31 (confirmed: `1979.zip`
#'   through `2005.zip` are the only 4-digit-named zips on the index);
#' - monthly `YYYYMM.zip`, 2006-01-01 through 2013-03-31 (confirmed:
#'   `200601.zip` through `201303.zip` are the only 6-digit-named zips);
#' - daily `YYYYMMDD.export.CSV.zip`, from 2013-04-01 (confirmed: the
#'   earliest 8-digit-named file is `20130401.export.CSV.zip`, and the
#'   codebook states the Historical Backfile -- 57 fields -- "runs January
#'   1, 1979 through March 31, 2013" while the Daily Updates collection --
#'   58 fields -- "begins April 1, 2013").
#'
#' Yearly/monthly backfiles have 57 tab-separated columns; daily files have
#' 58, the extra one being `SOURCEURL` (confirmed against the codebook,
#' which states `SOURCEURL` "is only present in the daily event stream
#' files beginning April 1, 2013"). All files are tab-delimited with no
#' header row, despite the `.csv`/`.CSV` extension.
#'
#' Backfile rows are grouped by event date (`SQLDATE`); daily files hold
#' the events *added* that day (`DATEADDED`), whose `SQLDATE` can be
#' earlier (the codebook: `DATEADDED` "gives the date the event was added
#' to the database... news coverage published today could add events from
#' the distant past, which would result in the SQLDATE and other event
#' date fields containing the date the event actually took place, while
#' the DATEADDED field... will carry today's date"). This matters for
#' `date_range` filtering in [read_gdelt()], which filters on `SQLDATE`.
#'
#' The index also links two verification files, confirmed by fetching them
#' directly: `md5sums` (`https://data.gdeltproject.org/events/md5sums`,
#' one `<md5>␣␣<filename>` pair per line, standard `md5sum` format) and
#' `filesizes` (`https://data.gdeltproject.org/events/filesizes`, one
#' `<bytes>␣<filename>` pair per line).
#'
#' `https://data.gdeltproject.org/events/` serves the files directly
#' (confirmed: `curl -sS https://data.gdeltproject.org/events/index.html`
#' returns 200); a plain `http://` request to the same host redirects to
#' `https://` (confirmed via `curl -L`, one redirect). The base URL is
#' overridable with `getOption("bilatr.gdelt_base_url", <default>)`, which
#' is what makes [download_gdelt()] testable offline against `file://`
#' fixtures (see `tests/testthat/helper_gdelt.R`).
#'
#' @keywords internal
#' @name bilatr_gdelt_layout
NULL

.GDELT_YEARLY_START <- as.Date("1979-01-01")
.GDELT_YEARLY_END <- as.Date("2005-12-31")
.GDELT_MONTHLY_START <- as.Date("2006-01-01")
.GDELT_MONTHLY_END <- as.Date("2013-03-31")
.GDELT_DAILY_START <- as.Date("2013-04-01")

.GDELT_MD5SUMS_FILE <- "md5sums"
.GDELT_FILESIZES_FILE <- "filesizes"

#' The base URL GDELT downloads/reads are served from
#'
#' Overridable with `getOption("bilatr.gdelt_base_url", ...)`, so tests
#' can point it at a `file://` fixture directory instead of the live
#' server. See `?bilatr_gdelt_layout`.
#'
#' @return A single string, always ending in `"/"`.
#' @keywords internal
.gdelt_base_url <- function() {
  base <- getOption("bilatr.gdelt_base_url", "https://data.gdeltproject.org/events/")
  if (!grepl("/$", base)) base <- paste0(base, "/")
  base
}

# ---------------------------------------------------------------------------
# Part 2: schema
# ---------------------------------------------------------------------------

#' The full 58-column GDELT 1.0 schema, in file order
#'
#' Column types are read off the GDELT 1.0 Data Format Codebook (see
#' `?bilatr_gdelt_layout`), made explicit rather than guessed per file.
#' `GLOBALEVENTID` and `DATEADDED` are stored as `"double"` rather than the
#' codebook's `(integer)`, since a 32-bit R integer can overflow for large
#' identifiers/timestamps and `double` is exact to 2^53. All `*Code`
#' fields (actor, event, type, known-group, ethnic, religion,
#' `*Geo_CountryCode`, `*Geo_ADM1Code`), `*Name`/`*FullName` fields,
#' `*Geo_FeatureID` fields, and `SOURCEURL` are `"character"`.
#'
#' @return A tibble with `name` and `type` (`"character"`, `"integer"`, or
#'   `"double"`), 58 rows, in file order (the daily-file layout; the
#'   backfile layout is the first 57).
#' @keywords internal
.gdelt_schema <- function() {
  actor_fields <- function(prefix) {
    c(
      paste0(prefix, "Code"), paste0(prefix, "Name"), paste0(prefix, "CountryCode"),
      paste0(prefix, "KnownGroupCode"), paste0(prefix, "EthnicCode"),
      paste0(prefix, "Religion1Code"), paste0(prefix, "Religion2Code"),
      paste0(prefix, "Type1Code"), paste0(prefix, "Type2Code"), paste0(prefix, "Type3Code")
    )
  }
  geo_fields <- function(prefix) {
    c(
      paste0(prefix, "Geo_Type"), paste0(prefix, "Geo_FullName"),
      paste0(prefix, "Geo_CountryCode"), paste0(prefix, "Geo_ADM1Code"),
      paste0(prefix, "Geo_Lat"), paste0(prefix, "Geo_Long"), paste0(prefix, "Geo_FeatureID")
    )
  }

  name <- c(
    "GLOBALEVENTID", "SQLDATE", "MonthYear", "Year", "FractionDate",
    actor_fields("Actor1"),
    actor_fields("Actor2"),
    "IsRootEvent", "EventCode", "EventBaseCode", "EventRootCode", "QuadClass",
    "GoldsteinScale", "NumMentions", "NumSources", "NumArticles", "AvgTone",
    geo_fields("Actor1"),
    geo_fields("Actor2"),
    geo_fields("Action"),
    "DATEADDED", "SOURCEURL"
  )

  character_names <- c(
    actor_fields("Actor1"), actor_fields("Actor2"),
    "EventCode", "EventBaseCode", "EventRootCode",
    "Actor1Geo_FullName", "Actor1Geo_CountryCode", "Actor1Geo_ADM1Code", "Actor1Geo_FeatureID",
    "Actor2Geo_FullName", "Actor2Geo_CountryCode", "Actor2Geo_ADM1Code", "Actor2Geo_FeatureID",
    "ActionGeo_FullName", "ActionGeo_CountryCode", "ActionGeo_ADM1Code", "ActionGeo_FeatureID",
    "SOURCEURL"
  )
  integer_names <- c(
    "SQLDATE", "MonthYear", "Year", "IsRootEvent", "QuadClass",
    "NumMentions", "NumSources", "NumArticles",
    "Actor1Geo_Type", "Actor2Geo_Type", "ActionGeo_Type"
  )
  double_names <- c(
    "GLOBALEVENTID", "DATEADDED", "FractionDate", "GoldsteinScale", "AvgTone",
    "Actor1Geo_Lat", "Actor1Geo_Long", "Actor2Geo_Lat", "Actor2Geo_Long",
    "ActionGeo_Lat", "ActionGeo_Long"
  )

  type <- dplyr::case_when(
    name %in% character_names ~ "character",
    name %in% integer_names ~ "integer",
    name %in% double_names ~ "double",
    TRUE ~ NA_character_
  )

  tibble::tibble(name = name, type = type)
}

#' Named presets of GDELT columns
#'
#' Name-vector presets over the full 58-column schema (`?bilatr_gdelt_layout`),
#' for use as [read_gdelt()]'s `columns` argument.
#'
#' GDELT's own `QuadClass` is deliberately left out of `"core"`:
#' [recode_cameo()] attaches its own `QuadClass` derived from `EventCode`,
#' and if both are present it keeps the existing column and warns --
#' correct, but easy to trip over if you expect `"core"` to include
#' GDELT's `QuadClass`. Request it explicitly (or via `"all"`) if you want
#' GDELT's own value instead of `recode_cameo()`'s.
#'
#' @param set One of `"core"` (the essentials: identifiers, both actors'
#'   `Code`/`CountryCode`/`Type1-3Code`, the event-code hierarchy,
#'   `IsRootEvent`, and the numeric/impact fields), `"actors"` (`"core"`
#'   plus every other `Actor1*`/`Actor2*` field: `*Name`,
#'   `*KnownGroupCode`, `*EthnicCode`, `*Religion1Code`,
#'   `*Religion2Code`), `"geo"` (`"core"` plus every `*Geo_*` field), or
#'   `"all"` (every column, including `SOURCEURL`).
#' @return A character vector of column names, in file order.
#' @examples
#' gdelt_columns("core")
#' @export
gdelt_columns <- function(set = c("core", "actors", "geo", "all")) {
  set <- match.arg(set)
  schema <- .gdelt_schema()
  all_names <- schema$name

  core <- c(
    "GLOBALEVENTID", "SQLDATE",
    "Actor1Code", "Actor1CountryCode", "Actor1Type1Code", "Actor1Type2Code", "Actor1Type3Code",
    "Actor2Code", "Actor2CountryCode", "Actor2Type1Code", "Actor2Type2Code", "Actor2Type3Code",
    "IsRootEvent", "EventCode", "EventBaseCode", "EventRootCode",
    "GoldsteinScale", "NumMentions", "NumSources", "NumArticles", "AvgTone"
  )

  if (set == "all") {
    return(all_names)
  }
  if (set == "core") {
    return(all_names[all_names %in% core])
  }
  if (set == "actors") {
    actor_all <- grep("^Actor[12](?!Geo_)", all_names, value = TRUE, perl = TRUE)
    return(all_names[all_names %in% union(core, actor_all)])
  }
  # geo
  geo_all <- grep("Geo_", all_names, value = TRUE, fixed = TRUE)
  all_names[all_names %in% union(core, geo_all)]
}

#' CAMEO actor/role type codes
#'
#' The vocabulary [read_gdelt()]'s `actor_types` argument is checked
#' against. Taken from the CAMEO 1.1b3 Manual's Table 3.1 ("Generic
#' Domestic Role Codes": primary, secondary, and tertiary role codes) and
#' Table 3.2 ("International/Transnational Generic Codes") -- the source
#' the GDELT codebook itself cites for "the complete available taxonomy"
#' of the `*Type1-3Code` fields, since the codebook lists examples but not
#' the full vocabulary itself. Verified 2026-09-25 against
#' `http://gdeltproject.org/data/documentation/CAMEO.Manual.1.1b3.pdf`.
#'
#' @return A character vector of 3-character codes.
#' @keywords internal
.gdelt_actor_type_codes <- function() {
  c(
    # Table 3.1, primary role codes
    "COP", "GOV", "INS", "JUD", "MIL", "OPP", "REB", "SEP", "SPY", "UAF",
    # Table 3.1, secondary role codes
    "AGR", "BUS", "CRM", "CVL", "DEV", "EDU", "ELI", "ENV", "HLH", "HRI",
    "LAB", "LEG", "MED", "REF",
    # Table 3.1, tertiary role codes
    "MOD", "RAD",
    # Table 3.2, international/transnational generic codes
    "IGO", "IMG", "INT", "MNC", "NGM", "NGO", "UIS"
  )
}

# ---------------------------------------------------------------------------
# Part 3a: gdelt_files()
# ---------------------------------------------------------------------------

#' List the GDELT files covering a date range
#'
#' A pure function (no network): rolls `[start, end]` up to the distinct
#' yearly/monthly/daily files that cover it, per the boundaries at
#' `?bilatr_gdelt_layout`. A range crossing a boundary returns a mix of
#' types.
#'
#' @param start,end Anything [.parse_event_date()] accepts (`Date`,
#'   `POSIXt`, `YYYYMMDD`, or `YYYY-MM-DD`). `end` defaults to `start`
#'   (a single period). Both must fall within 1979-01-01 and today.
#' @return A tibble with one row per file: `period` (the file's date
#'   stamp, `"YYYY"`/`"YYYYMM"`/`"YYYYMMDD"`), `type`
#'   (`"yearly"`/`"monthly"`/`"daily"`), `file` (the basename), `url`
#'   (under [.gdelt_base_url()]), and `n_cols` (57 or 58).
#' @examples
#' gdelt_files("2020-01-01", "2020-01-03")
#' gdelt_files("1990-06-15") # one yearly file
#' @export
gdelt_files <- function(start, end = start) {
  start_date <- .parse_event_date(start)
  end_date <- .parse_event_date(end)
  if (is.na(start_date)) {
    stop("gdelt_files(): `start` could not be parsed as a date.", call. = FALSE)
  }
  if (is.na(end_date)) {
    stop("gdelt_files(): `end` could not be parsed as a date.", call. = FALSE)
  }
  if (end_date < start_date) {
    stop("gdelt_files(): `end` must not be before `start`.", call. = FALSE)
  }

  today <- Sys.Date()
  if (start_date < .GDELT_YEARLY_START || end_date > today) {
    stop(sprintf(
      "gdelt_files(): dates must be between %s and today (%s); got start = %s, end = %s.",
      .GDELT_YEARLY_START, today, start_date, end_date
    ), call. = FALSE)
  }

  dates <- seq(start_date, end_date, by = "day")
  type <- dplyr::case_when(
    dates <= .GDELT_YEARLY_END ~ "yearly",
    dates <= .GDELT_MONTHLY_END ~ "monthly",
    TRUE ~ "daily"
  )
  period <- dplyr::case_when(
    type == "yearly" ~ format(dates, "%Y"),
    type == "monthly" ~ format(dates, "%Y%m"),
    TRUE ~ format(dates, "%Y%m%d")
  )

  out <- dplyr::distinct(tibble::tibble(period = period, type = type))
  out$file <- ifelse(out$type == "daily", paste0(out$period, ".export.CSV.zip"), paste0(out$period, ".zip"))
  out$url <- paste0(.gdelt_base_url(), out$file)
  out$n_cols <- ifelse(out$type == "daily", 58L, 57L)
  out
}

# ---------------------------------------------------------------------------
# Shared internals for read_gdelt() / count_gdelt_actor_types()
# ---------------------------------------------------------------------------

#' Infer a GDELT file's column count from its name
#' @keywords internal
.gdelt_filename_ncols <- function(path) {
  if (grepl("^[0-9]{8}\\.export\\.CSV\\.zip$", basename(path))) 58L else 57L
}

#' Read the first line of a (possibly zipped) GDELT file
#' @keywords internal
.gdelt_peek_first_line <- function(path) {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    inner <- utils::unzip(path, list = TRUE)$Name[1]
    con <- unz(path, inner, open = "rt")
  } else {
    con <- file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  readLines(con, n = 1, warn = FALSE)
}

#' Resolve and validate a `date_range` argument
#' @keywords internal
.gdelt_resolve_date_range <- function(date_range) {
  if (is.null(date_range)) {
    return(NULL)
  }
  if (length(date_range) != 2) {
    stop("`date_range` must have length 2 (start, end).", call. = FALSE)
  }
  parsed <- .parse_event_date(date_range)
  if (anyNA(parsed)) {
    stop("`date_range` could not be parsed as dates.", call. = FALSE)
  }
  sort(parsed)
}

#' Build an empty tibble with the right column types
#' @keywords internal
.gdelt_empty_tibble <- function(schema, cols) {
  make_col <- function(type) {
    switch(type, character = character(0), integer = integer(0), double = double(0))
  }
  types <- schema$type[match(cols, schema$name)]
  tibble::as_tibble(stats::setNames(purrr::map(types, make_col), cols))
}

#' Apply `read_gdelt()`'s row filters to one chunk
#' @keywords internal
.gdelt_apply_filters <- function(chunk, actor_types, actor_type_match,
                                  countries, country_match,
                                  cross_border, root_events_only,
                                  event_codes, date_bounds) {
  keep <- rep(TRUE, nrow(chunk))

  if (isTRUE(cross_border)) {
    keep <- keep & !is.na(chunk$Actor1CountryCode) & !is.na(chunk$Actor2CountryCode) &
      chunk$Actor1CountryCode != chunk$Actor2CountryCode
  }

  if (!is.null(actor_types)) {
    a1 <- (chunk$Actor1Type1Code %in% actor_types) |
      (chunk$Actor1Type2Code %in% actor_types) |
      (chunk$Actor1Type3Code %in% actor_types)
    a2 <- (chunk$Actor2Type1Code %in% actor_types) |
      (chunk$Actor2Type2Code %in% actor_types) |
      (chunk$Actor2Type3Code %in% actor_types)
    keep <- keep & (if (identical(actor_type_match, "both")) a1 & a2 else a1 | a2)
  }

  if (!is.null(countries)) {
    a1 <- chunk$Actor1CountryCode %in% countries
    a2 <- chunk$Actor2CountryCode %in% countries
    keep <- keep & (if (identical(country_match, "both")) a1 & a2 else a1 | a2)
  }

  if (isTRUE(root_events_only)) {
    keep <- keep & !is.na(chunk$IsRootEvent) & chunk$IsRootEvent == 1L
  }

  if (!is.null(event_codes)) {
    matches_any <- Reduce(
      `|`,
      lapply(event_codes, function(p) startsWith(chunk$EventCode, p)),
      init = rep(FALSE, nrow(chunk))
    )
    keep <- keep & !is.na(chunk$EventCode) & matches_any
  }

  if (!is.null(date_bounds)) {
    event_date <- .parse_event_date(chunk$SQLDATE)
    keep <- keep & !is.na(event_date) & event_date >= date_bounds[1] & event_date <= date_bounds[2]
  }

  chunk[keep, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Part 3b: read_gdelt()
# ---------------------------------------------------------------------------

#' Read and filter one or more raw GDELT export zips
#'
#' Reads GDELT 1.0 event files (as listed by [gdelt_files()] or downloaded
#' by [download_gdelt()]), applying an explicit column spec (`?bilatr_gdelt_layout`)
#' and row filters. Filters run on the full row **before** `columns` is
#' applied, so a filter column does not need to be requested in `columns`.
#' Memory is kept bounded with `readr::read_tsv_chunked()` (filtered chunk
#' by chunk, so only surviving rows are ever held); this reads each `.zip`
#' directly (no `unzip()` to a temp directory needed -- confirmed
#' empirically for this readr version).
#'
#' For each file, the column count is inferred from the filename (a
#' daily-pattern name -- `YYYYMMDD.export.CSV.zip` -- means 58 columns,
#' anything else 57) and then confirmed by counting the tab-separated
#' fields on the file's first line; a disagreement is trusted to the file
#' (not the name) and produces a warning. Reads use `quote = ""` and
#' `na = ""`: GDELT's free-text fields (`Actor1Name`, `*Geo_FullName`) can
#' contain literal `"` characters, and readr's default quote handling
#' mis-splits such rows.
#'
#' GDELT 2.0 (15-minute-interval files, the Mentions/GKG tables) is not
#' supported by this function or by [gdelt_files()]/[download_gdelt()];
#' this package reads GDELT 1.0's yearly/monthly/daily event files only.
#'
#' @param files Character vector of paths to GDELT export zips (or, for a
#'   non-zipped file, a plain tab-separated path). Must be non-empty.
#' @param columns Character vector of column names to keep in the result
#'   (see [gdelt_columns()]), or `"all"`. Unknown names error, listing the
#'   valid ones. If a 57-column (backfile) file is read with `"SOURCEURL"`
#'   requested, that file's rows get `NA` for it, so output from mixed
#'   eras binds cleanly.
#' @param actor_types Character vector of CAMEO actor/role type codes
#'   (`?bilatr_gdelt_layout`'s `.gdelt_actor_type_codes()`); `NULL` turns
#'   the filter off. Unknown codes error, listing the valid ones.
#' @param actor_type_match `"both"` (default; matches
#'   `extract_all_relevant_gdelt()`'s pre-0.9.1 semantics): both actors
#'   must carry one of `actor_types` in one of their `Type1-3Code`
#'   fields. `"either"`: one actor is enough (covers state/non-state
#'   interactions without a second code path).
#' @param countries Character vector matched against
#'   `Actor1CountryCode`/`Actor2CountryCode`, or `NULL` (default) for no
#'   country filter.
#' @param country_match `"either"` (default): keep events involving any
#'   listed country. `"both"`: keep only events between two listed
#'   countries (the closed set).
#' @param cross_border `TRUE` (default): keep only rows where
#'   `Actor1CountryCode`/`Actor2CountryCode` are both non-`NA` and
#'   different. `FALSE`: no country-presence filter at all.
#' @param root_events_only `FALSE` (default), or `TRUE` to keep only
#'   `IsRootEvent == 1` rows.
#' @param event_codes Character vector of CAMEO event-code prefixes, or
#'   `NULL` (default). A row is kept if its `EventCode` **starts with**
#'   any given code -- e.g. `"19"` keeps every code under root `19`, while
#'   `"0211"` keeps only that exact code.
#' @param date_range Length-2 vector (anything [.parse_event_date()]
#'   accepts), or `NULL` (default). Filters on `SQLDATE`. For daily files
#'   (organised by `DATEADDED`, not `SQLDATE`; see `?bilatr_gdelt_layout`),
#'   this keeps events *dated* in the range among the files actually read
#'   -- it does not recover events dated in the range but added (and so
#'   filed) outside those files.
#' @return A tibble of event-level rows, `columns` plus `source_file` (the
#'   file's basename), row-bound across `files`.
#' @examples
#' \dontrun{
#' events <- read_gdelt("data/gdelt_raw/20200101.export.CSV.zip")
#' }
#' @export
read_gdelt <- function(
  files,
  columns = gdelt_columns("core"),
  actor_types = c("GOV", "MIL", "SPY"),
  actor_type_match = c("both", "either"),
  countries = NULL,
  country_match = c("either", "both"),
  cross_border = TRUE,
  root_events_only = FALSE,
  event_codes = NULL,
  date_range = NULL
) {
  actor_type_match <- match.arg(actor_type_match)
  country_match <- match.arg(country_match)

  if (length(files) == 0) {
    stop("read_gdelt(): `files` must be a non-empty character vector of paths.", call. = FALSE)
  }

  schema <- .gdelt_schema()
  if (identical(columns, "all")) {
    columns <- schema$name
  }
  unknown_cols <- setdiff(columns, schema$name)
  if (length(unknown_cols) > 0) {
    stop(
      "read_gdelt(): unknown column(s) in `columns`: ", paste(unknown_cols, collapse = ", "),
      ". Valid columns: ", paste(schema$name, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.null(actor_types)) {
    unknown_types <- setdiff(actor_types, .gdelt_actor_type_codes())
    if (length(unknown_types) > 0) {
      stop(
        "read_gdelt(): unknown `actor_types`: ", paste(unknown_types, collapse = ", "),
        ". Valid CAMEO actor/role codes: ", paste(.gdelt_actor_type_codes(), collapse = ", "),
        call. = FALSE
      )
    }
  }

  date_bounds <- .gdelt_resolve_date_range(date_range)

  filter_cols <- c(
    "Actor1CountryCode", "Actor2CountryCode",
    "Actor1Type1Code", "Actor1Type2Code", "Actor1Type3Code",
    "Actor2Type1Code", "Actor2Type2Code", "Actor2Type3Code",
    "IsRootEvent", "EventCode", "SQLDATE"
  )

  results <- purrr::map(files, function(f) {
    .read_gdelt_one(
      f,
      schema = schema, columns = columns, filter_cols = filter_cols,
      actor_types = actor_types, actor_type_match = actor_type_match,
      countries = countries, country_match = country_match,
      cross_border = cross_border, root_events_only = root_events_only,
      event_codes = event_codes, date_bounds = date_bounds
    )
  })

  dplyr::bind_rows(results)
}

#' Read and filter one GDELT file (the per-file worker for [read_gdelt()])
#' @keywords internal
.read_gdelt_one <- function(file, schema, columns, filter_cols,
                             actor_types, actor_type_match,
                             countries, country_match,
                             cross_border, root_events_only,
                             event_codes, date_bounds) {
  expected_n <- .gdelt_filename_ncols(file)
  first_line <- .gdelt_peek_first_line(file)
  actual_n <- length(strsplit(first_line, "\t", fixed = TRUE)[[1]])
  n_cols <- expected_n
  if (actual_n != expected_n && actual_n > 0) {
    warning(sprintf(
      "read_gdelt(): '%s' looks like a %d-column file by name, but its first line has %d field(s); using %d (the file's own field count).",
      basename(file), expected_n, actual_n, actual_n
    ), call. = FALSE)
    n_cols <- actual_n
  }
  n_cols <- min(n_cols, nrow(schema))

  file_names <- schema$name[seq_len(n_cols)]
  type_letter <- c(character = "c", integer = "i", double = "d")
  col_types_str <- paste(type_letter[schema$type[seq_len(n_cols)]], collapse = "")

  read_cols <- intersect(union(columns, filter_cols), file_names)

  out <- NULL
  readr::read_tsv_chunked(
    file,
    callback = readr::DataFrameCallback$new(function(chunk, pos) {
      filtered <- .gdelt_apply_filters(
        chunk, actor_types, actor_type_match, countries, country_match,
        cross_border, root_events_only, event_codes, date_bounds
      )
      if (nrow(filtered) > 0) {
        out <<- dplyr::bind_rows(out, filtered[read_cols])
      }
    }),
    col_names = file_names,
    col_types = col_types_str,
    quote = "",
    na = "",
    progress = FALSE
  )

  if (is.null(out)) {
    out <- .gdelt_empty_tibble(schema, read_cols)
  }

  # Every requested column present, even one absent from this particular
  # file (e.g. SOURCEURL for a 57-column backfile) -- filled NA, so output
  # from mixed eras binds cleanly (see @param columns).
  missing_requested <- setdiff(columns, names(out))
  for (mc in missing_requested) {
    mc_type <- schema$type[schema$name == mc]
    out[[mc]] <- switch(mc_type, character = NA_character_, integer = NA_integer_, double = NA_real_)
  }

  out <- out[columns]
  out$source_file <- basename(file)
  out
}

# ---------------------------------------------------------------------------
# Part 3c: download_gdelt()
# ---------------------------------------------------------------------------

#' Fetch and parse GDELT's `md5sums` index
#' @return A named character vector (name = filename, value = md5), or
#'   `NULL` if the index could not be fetched.
#' @keywords internal
.gdelt_fetch_md5sums <- function() {
  url <- paste0(.gdelt_base_url(), .GDELT_MD5SUMS_FILE)
  tryCatch(
    {
      lines <- readLines(url, warn = FALSE)
      lines <- lines[nzchar(trimws(lines))]
      parts <- strsplit(trimws(lines), "\\s+")
      md5 <- vapply(parts, `[[`, character(1), 1)
      fname <- vapply(parts, function(p) paste(p[-1], collapse = " "), character(1))
      stats::setNames(md5, fname)
    },
    error = function(e) NULL,
    warning = function(w) NULL
  )
}

#' Download one GDELT file to `destfile`, verifying and staging via `.part`
#'
#' Exactly one request: `utils::download.file(mode = "wb")` to
#' `<destfile>.part`, renamed to `destfile` only on success, so a failed
#' or interrupted download never leaves a file at the final path. Never
#' `try()`-and-return-the-path (audit G2): a failed/missing/corrupt result
#' becomes a status, with no warning raised here (the caller raises one
#' summary warning listing every failure).
#'
#' @return A list with `status` (`"downloaded"`, `"missing"`, `"failed"`,
#'   or `"corrupt"`) and `path` (`destfile` if downloaded, else `NA`).
#' @keywords internal
.gdelt_download_one <- function(url, destfile, md5_lookup, verify) {
  part <- paste0(destfile, ".part")
  on.exit(if (file.exists(part)) unlink(part), add = TRUE)

  dl_status <- tryCatch(
    {
      utils::download.file(url, destfile = part, mode = "wb", quiet = TRUE)
      "ok"
    },
    warning = function(w) {
      if (grepl("404", conditionMessage(w), fixed = TRUE)) "missing" else "failed"
    },
    error = function(e) {
      if (grepl("404", conditionMessage(e), fixed = TRUE) ||
        grepl("No such file", conditionMessage(e), fixed = TRUE)) {
        "missing"
      } else {
        "failed"
      }
    }
  )

  if (!identical(dl_status, "ok")) {
    return(list(status = dl_status, path = NA_character_))
  }
  if (!file.exists(part)) {
    return(list(status = "failed", path = NA_character_))
  }

  if (isTRUE(verify) && !is.null(md5_lookup)) {
    expected <- unname(md5_lookup[basename(destfile)])
    if (!is.na(expected)) {
      actual <- unname(tools::md5sum(part))
      if (!identical(actual, expected)) {
        return(list(status = "corrupt", path = NA_character_))
      }
    }
  }

  file.rename(part, destfile)
  list(status = "downloaded", path = destfile)
}

#' Raise one summary warning listing every failed/missing/corrupt file
#' @keywords internal
.gdelt_warn_failures <- function(plan, statuses) {
  bad <- statuses %in% c("failed", "missing", "corrupt")
  if (any(bad)) {
    warning(
      "download_gdelt(): ", sum(bad), " of ", length(statuses), " file(s) not available: ",
      paste(sprintf("%s (%s)", plan$period[bad], statuses[bad]), collapse = ", "),
      call. = FALSE
    )
  }
}

#' Human-readable byte count for progress messages
#' @keywords internal
.gdelt_format_bytes <- function(bytes) {
  if (is.na(bytes)) {
    return("")
  }
  if (bytes >= 1e6) {
    return(sprintf("%.1fMB", bytes / 1e6))
  }
  if (bytes >= 1e3) {
    return(sprintf("%.0fKB", bytes / 1e3))
  }
  sprintf("%dB", as.integer(bytes))
}

#' `download_gdelt()`, temp mode: download, read, discard, per file
#' @keywords internal
.download_gdelt_temp <- function(plan, md5_lookup, verify, ...) {
  statuses <- character(nrow(plan))
  all_events <- vector("list", nrow(plan))

  for (i in seq_len(nrow(plan))) {
    tmp <- file.path(tempdir(), plan$file[i])
    res <- .gdelt_download_one(plan$url[i], tmp, md5_lookup = md5_lookup, verify = verify)
    statuses[i] <- res$status
    size_txt <- if (identical(res$status, "downloaded")) paste0(" (", .gdelt_format_bytes(file.size(res$path)), ")") else ""
    message(sprintf("download_gdelt(): %s (%s)%s -- %s", plan$period[i], plan$type[i], size_txt, res$status))

    if (identical(res$status, "downloaded")) {
      all_events[[i]] <- tryCatch(
        read_gdelt(res$path, ...),
        finally = unlink(res$path)
      )
    }
  }

  .gdelt_warn_failures(plan, statuses)

  result <- dplyr::bind_rows(all_events)
  attr(result, "files") <- dplyr::mutate(plan, status = statuses)
  result
}

#' `download_gdelt()`, cache mode: download to `dest_dir`, keep the zips
#' @keywords internal
.download_gdelt_cache <- function(plan, dest_dir, overwrite, md5_lookup, verify) {
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }

  paths <- file.path(dest_dir, plan$file)
  statuses <- character(nrow(plan))

  for (i in seq_len(nrow(plan))) {
    if (file.exists(paths[i]) && !overwrite) {
      if (isTRUE(verify) && !is.null(md5_lookup)) {
        expected <- unname(md5_lookup[plan$file[i]])
        if (!is.na(expected) && !identical(unname(tools::md5sum(paths[i])), expected)) {
          statuses[i] <- "corrupt"
          unlink(paths[i])
          message(sprintf("download_gdelt(): %s (%s) -- corrupt", plan$period[i], plan$type[i]))
          next
        }
      }
      statuses[i] <- "cached"
      size_txt <- paste0(" (", .gdelt_format_bytes(file.size(paths[i])), ")")
      message(sprintf("download_gdelt(): %s (%s)%s -- cached", plan$period[i], plan$type[i], size_txt))
      next
    }

    res <- .gdelt_download_one(plan$url[i], paths[i], md5_lookup = md5_lookup, verify = verify)
    statuses[i] <- res$status
    size_txt <- if (identical(res$status, "downloaded")) paste0(" (", .gdelt_format_bytes(file.size(res$path)), ")") else ""
    message(sprintf("download_gdelt(): %s (%s)%s -- %s", plan$period[i], plan$type[i], size_txt, res$status))
  }

  .gdelt_warn_failures(plan, statuses)

  out <- plan
  out$path <- ifelse(statuses %in% c("downloaded", "cached"), paths, NA_character_)
  out$status <- statuses
  invisible(out)
}

#' Download raw GDELT export zips
#'
#' Lists the files covering `[start, end]` with [gdelt_files()], then
#' downloads each (this fixes audit G2/G3: one request per file, under a
#' locally-raised, generous `timeout`, never a swallowed failure).
#'
#' **`dest_dir = NULL`** (default): each file is downloaded to a
#' temporary path, read with `read_gdelt(<tmp>, ...)` (`...` passed
#' through), and deleted -- so this mode needs no folder and leaves no
#' zips behind. Returns the row-bound, filtered events; the per-file
#' status table (`gdelt_files()`'s columns plus `status`) is attached as
#' `attr(result, "files")`.
#'
#' **`dest_dir = "some/dir"`**: creates the folder if needed and keeps the
#' zips as a cache. A file already present is skipped (verified against
#' `md5sums` too, if `verify = TRUE`) unless `overwrite = TRUE`. Data is
#' **not** read in this mode -- call `read_gdelt(status$path[status$status
#' %in% c("downloaded", "cached")], ...)` yourself. Returns (invisibly)
#' `gdelt_files()`'s columns plus `path` (`NA` for anything not on disk)
#' and `status` (one of `"downloaded"`, `"cached"`, `"missing"` -- HTTP
#' 404, no such file on the server -- `"failed"`, or `"corrupt"`).
#'
#' A failed/missing/corrupt file never aborts the whole call: it becomes
#' a status, and one summary `warning()` at the end lists every such
#' file. Downloads run sequentially (GDELT is a shared public server), with
#' a one-line progress `message()` per file.
#'
#' @param start,end See [gdelt_files()].
#' @param dest_dir `NULL` (default, temp mode) or a directory path (cache
#'   mode).
#' @param overwrite In cache mode, re-download files already present.
#'   Ignored in temp mode.
#' @param timeout Seconds before a single download times out, raised via
#'   a locally-scoped `options(timeout = timeout)` (restored on exit) --
#'   R's default 60s truncates large backfiles.
#' @param verify Verify each downloaded (and, in cache mode, each
#'   already-cached) file's MD5 against GDELT's `md5sums` index. If the
#'   index itself can't be fetched, warns once and proceeds unverified.
#' @param ... In temp mode, passed through to `read_gdelt()`. Ignored in
#'   cache mode.
#' @return See "Details" above.
#' @examples
#' \dontrun{
#' # temp mode: a quick look, nothing kept on disk
#' events <- download_gdelt("2020-01-01")
#'
#' # cache mode: for real work -- keep the zips, read them yourself
#' status <- download_gdelt("2020-01-01", "2020-01-31", dest_dir = "data/gdelt_raw")
#' events <- read_gdelt(status$path[status$status %in% c("downloaded", "cached")])
#' }
#' @export
download_gdelt <- function(
  start,
  end = start,
  dest_dir = NULL,
  overwrite = FALSE,
  timeout = 3600,
  verify = TRUE,
  ...
) {
  plan <- gdelt_files(start, end)

  old_timeout <- getOption("timeout")
  options(timeout = timeout)
  on.exit(options(timeout = old_timeout), add = TRUE)

  md5_lookup <- NULL
  if (isTRUE(verify)) {
    md5_lookup <- .gdelt_fetch_md5sums()
    if (is.null(md5_lookup)) {
      warning("download_gdelt(): could not fetch the md5sums index; proceeding unverified.", call. = FALSE)
    }
  }

  if (is.null(dest_dir)) {
    .download_gdelt_temp(plan, md5_lookup = md5_lookup, verify = verify, ...)
  } else {
    .download_gdelt_cache(plan, dest_dir, overwrite = overwrite, md5_lookup = md5_lookup, verify = verify)
  }
}

# ---------------------------------------------------------------------------
# Part 3d: count_gdelt_actor_types()
# ---------------------------------------------------------------------------

#' Count one chunk's actor-type-code combinations, over all three slots per side
#' @keywords internal
.gdelt_count_actor_types_chunk <- function(chunk) {
  chunk$.row_id <- seq_len(nrow(chunk))

  side <- function(prefix, value_name) {
    tidyr::pivot_longer(
      chunk[c(".row_id", paste0(prefix, "Type1Code"), paste0(prefix, "Type2Code"), paste0(prefix, "Type3Code"))],
      cols = -".row_id",
      values_to = value_name
    ) %>%
      dplyr::filter(!is.na(.data[[value_name]])) %>%
      dplyr::distinct(.data$.row_id, .data[[value_name]])
  }

  a1 <- side("Actor1", "actor1_type")
  a2 <- side("Actor2", "actor2_type")

  dplyr::inner_join(a1, a2, by = ".row_id", relationship = "many-to-many") %>%
    dplyr::count(.data$actor1_type, .data$actor2_type, name = "n")
}

#' Read one file's contribution to `count_gdelt_actor_types()`
#' @keywords internal
.count_gdelt_actor_types_one <- function(file, schema,
                                          countries, country_match,
                                          cross_border, root_events_only,
                                          date_bounds) {
  expected_n <- .gdelt_filename_ncols(file)
  first_line <- .gdelt_peek_first_line(file)
  actual_n <- length(strsplit(first_line, "\t", fixed = TRUE)[[1]])
  n_cols <- if (actual_n != expected_n && actual_n > 0) actual_n else expected_n
  n_cols <- min(n_cols, nrow(schema))

  file_names <- schema$name[seq_len(n_cols)]
  type_letter <- c(character = "c", integer = "i", double = "d")
  col_types_str <- paste(type_letter[schema$type[seq_len(n_cols)]], collapse = "")

  out <- NULL
  readr::read_tsv_chunked(
    file,
    callback = readr::DataFrameCallback$new(function(chunk, pos) {
      filtered <- .gdelt_apply_filters(
        chunk,
        actor_types = NULL, actor_type_match = "both",
        countries = countries, country_match = country_match,
        cross_border = cross_border, root_events_only = root_events_only,
        event_codes = NULL, date_bounds = date_bounds
      )
      if (nrow(filtered) > 0) {
        out <<- dplyr::bind_rows(out, .gdelt_count_actor_types_chunk(filtered))
      }
    }),
    col_names = file_names,
    col_types = col_types_str,
    quote = "",
    na = "",
    progress = FALSE
  )

  if (is.null(out)) {
    out <- tibble::tibble(actor1_type = character(0), actor2_type = character(0), n = integer(0))
  }
  out
}

#' Count GDELT actor-type-code combinations
#'
#' Diagnostic helper for tuning [read_gdelt()]'s `actor_types` argument:
#' counts how often each `actor1_type` x `actor2_type` combination occurs,
#' over **all three** type slots per side (a row contributes each
#' distinct type code it carries on each side, so a row with `GOV`/`MIL`
#' on Actor1 and `SPY` on Actor2 contributes to both `(GOV, SPY)` and
#' `(MIL, SPY)`). Applies the same `countries`/`cross_border`/
#' `root_events_only`/`date_range` filters as [read_gdelt()], but not
#' `actor_types`, since its purpose is to tune that argument.
#'
#' @inheritParams read_gdelt
#' @return A tibble of `actor1_type`, `actor2_type`, `n`, sorted by `n`
#'   descending.
#' @examples
#' \dontrun{
#' count_gdelt_actor_types("data/gdelt_raw/20200101.export.CSV.zip")
#' }
#' @export
count_gdelt_actor_types <- function(
  files,
  countries = NULL,
  country_match = c("either", "both"),
  cross_border = TRUE,
  root_events_only = FALSE,
  date_range = NULL
) {
  country_match <- match.arg(country_match)
  if (length(files) == 0) {
    stop("count_gdelt_actor_types(): `files` must be a non-empty character vector of paths.", call. = FALSE)
  }

  schema <- .gdelt_schema()
  date_bounds <- .gdelt_resolve_date_range(date_range)

  per_file <- purrr::map(files, function(f) {
    .count_gdelt_actor_types_one(
      f,
      schema = schema,
      countries = countries, country_match = country_match,
      cross_border = cross_border, root_events_only = root_events_only,
      date_bounds = date_bounds
    )
  })

  dplyr::bind_rows(per_file) %>%
    dplyr::group_by(.data$actor1_type, .data$actor2_type) %>%
    dplyr::summarise(n = sum(.data$n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$n))
}
