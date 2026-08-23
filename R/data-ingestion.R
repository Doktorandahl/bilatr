#' GDELT raw event file column names
#'
#' The fixed column layout of GDELT's raw daily/monthly/yearly export
#' files, which ship without a header row.
#'
#' @return A character vector of column names, in file order.
#' @keywords internal
get_gdelt_column_names <- function() {
  c(
    "GLOBALEVENTID", "SQLDATE", "MonthYear", "Year", "FractionDate",
    "Actor1Code", "Actor1Name", "Actor1CountryCode", "Actor1KnownGroupCode",
    "Actor1EthnicCode", "Actor1Religion1Code", "Actor1Religion2Code",
    "Actor1Type1Code", "Actor1Type2Code", "Actor1Type3Code",
    "Actor2Code", "Actor2Name", "Actor2CountryCode", "Actor2KnownGroupCode",
    "Actor2EthnicCode", "Actor2Religion1Code", "Actor2Religion2Code",
    "Actor2Type1Code", "Actor2Type2Code", "Actor2Type3Code",
    "IsRootEvent", "EventCode", "EventBaseCode", "EventRootCode", "QuadClass",
    "GoldsteinScale", "NumMentions", "NumSources", "NumArticles", "AvgTone",
    "Actor1Geo_Type", "Actor1Geo_FullName", "Actor1Geo_CountryCode",
    "Actor1Geo_ADM1Code", "Actor1Geo_Lat", "Actor1Geo_Long", "Actor1Geo_FeatureID",
    "Actor2Geo_Type", "Actor2Geo_FullName", "Actor2Geo_CountryCode",
    "Actor2Geo_ADM1Code", "Actor2Geo_Lat", "Actor2Geo_Long", "Actor2Geo_FeatureID",
    "ActionGeo_Type", "ActionGeo_FullName", "ActionGeo_CountryCode",
    "ActionGeo_ADM1Code", "ActionGeo_Lat", "ActionGeo_Long", "ActionGeo_FeatureID",
    "DATEADDED"
  )
}

#' Actor-type codes treated as "relevant" state/security actors in GDELT
#'
#' @return A character vector of GDELT actor-type codes.
#' @keywords internal
relevant_actors <- function() {
  c("GOV", "MIL", "SPY")
}

#' Actor-sector name fragments treated as "relevant" state/security
#' actors in ICEWS
#'
#' @return A character vector of ICEWS sector-name substrings.
#' @keywords internal
relevant_actors_icews <- function() {
  c(
    "Air Force", "Army", "Cabinet", "Coast Guard", "Ministry",
    "Security Mi", "Execut", "Gov", "Legisl", "Lower House", "Milit",
    "Supreme Court", "Navy", "Police", "Unicameral", "Upper House"
  )
}

#' Download a raw GDELT export zip
#'
#' Downloads a single daily, monthly, or yearly GDELT events export to
#' `local_folder`, without unzipping or processing it. Intended to be
#' called once per date/period; use [furrr::future_walk()] (after setting
#' a parallel [future::plan()]) to download many periods concurrently,
#' since each download is independent I/O-bound work.
#'
#' @param date A date (or date-like string, `"YYYY-MM-DD"`) identifying
#'   the period to download.
#' @param local_folder Directory to save the zip file into. Must already
#'   exist.
#' @param type One of `"daily"`, `"monthly"`, or `"yearly"`.
#' @return Invisibly, the path to the downloaded zip file, or `NULL`
#'   (with a warning) if GDELT has no export for `date`.
#' @examples
#' \dontrun{
#' dir.create("data/gdelt_raw", recursive = TRUE, showWarnings = FALSE)
#' download_gdelt_raw_zip("2020-01-01", "data/gdelt_raw", type = "daily")
#'
#' # many dates in parallel
#' future::plan(future::multisession, workers = 8)
#' furrr::future_walk(
#'   as.character(seq(as.Date("2020-01-01"), as.Date("2020-01-31"), by = "day")),
#'   download_gdelt_raw_zip, local_folder = "data/gdelt_raw"
#' )
#' }
#' @export
download_gdelt_raw_zip <- function(date, local_folder = "data/gdelt_raw", type = c("daily", "monthly", "yearly")) {
  type <- match.arg(type)
  urls <- gdelt_export_url(date, local_folder, type)

  url_check <- httr::HEAD(urls$remote)
  if (httr::status_code(url_check) != 200L) {
    warning("GDELT data not available for ", date, call. = FALSE)
    return(invisible(NULL))
  }
  try(download.file(urls$remote, destfile = urls$local, mode = "wb", quiet = TRUE))
  invisible(urls$local)
}

#' Build the remote/local file paths for a GDELT export
#' @keywords internal
gdelt_export_url <- function(date, local_folder, type) {
  stamp <- switch(
    type,
    daily = stringr::str_remove_all(date, "-"),
    monthly = stringr::str_sub(stringr::str_remove_all(date, "-"), 1, 6),
    yearly = stringr::str_sub(stringr::str_remove_all(date, "-"), 1, 4)
  )
  suffix <- if (type == "daily") ".export.CSV.zip" else ".zip"
  list(
    remote = paste0("http://data.gdeltproject.org/events/", stamp, suffix),
    local = paste0(local_folder, "/", stamp, ".zip")
  )
}

#' Read and filter a raw GDELT zip to relevant dyadic events
#'
#' Reads a raw GDELT export zip (as downloaded by
#' [download_gdelt_raw_zip()]), and filters to cross-country dyadic events
#' where both actors are state/security actors (per [relevant_actors()]).
#' Retains event-level records (one row per event) rather than
#' aggregating, so the result can be recoded via [recode_cameo()] and fed
#' to [grouped_events_to_dyad_period()] for any CAMEO-derived grouping
#' variable.
#'
#' @param file Path to a raw GDELT export zip file.
#' @return A data frame of event-level records, filtered to relevant
#'   cross-country dyadic events.
#' @examples
#' \dontrun{
#' events <- extract_all_relevant_gdelt("data/gdelt_raw/20200101.zip")
#'
#' # many files at once
#' future::plan(future::multisession, workers = 8)
#' events <- furrr::future_map_dfr(
#'   list.files("data/gdelt_raw", full.names = TRUE),
#'   extract_all_relevant_gdelt
#' )
#' }
#' @export
extract_all_relevant_gdelt <- function(file) {
  rel_act <- relevant_actors()

  raw <- readr::read_delim(
    unz(file, unzip(file, list = TRUE)$Name[1]),
    num_threads = 1,
    col_names = get_gdelt_column_names(),
    progress = FALSE,
    show_col_types = FALSE
  )

  raw %>%
    dplyr::filter(
      !is.na(Actor1CountryCode) & !is.na(Actor2CountryCode),
      Actor1CountryCode != Actor2CountryCode,
      (Actor1Type1Code %in% rel_act | Actor1Type2Code %in% rel_act | Actor1Type3Code %in% rel_act) &
        (Actor2Type1Code %in% rel_act | Actor2Type2Code %in% rel_act | Actor2Type3Code %in% rel_act)
    )
}

#' Count actor-type-code combinations in a raw GDELT file
#'
#' Diagnostic helper for tuning [relevant_actors()] on new data: counts
#' how often each `Actor1Type1Code` x `Actor2Type1Code` combination occurs
#' among cross-country dyadic events in a raw GDELT export.
#'
#' @param file Path to a raw GDELT export zip file.
#' @return A data frame with `Actor1Type1Code`, `Actor2Type1Code`, `n`,
#'   and `date` (the file's date stamp).
#' @examples
#' \dontrun{
#' get_actor_combos("data/gdelt_raw/20200101.zip")
#' }
#' @export
get_actor_combos <- function(file) {
  raw <- readr::read_delim(
    unz(file, unzip(file, list = TRUE)$Name[1]),
    num_threads = 1,
    col_names = get_gdelt_column_names(),
    progress = FALSE,
    show_col_types = FALSE
  )

  raw %>%
    dplyr::filter(!is.na(Actor1CountryCode) & !is.na(Actor2CountryCode)) %>%
    dplyr::count(Actor1Type1Code, Actor2Type1Code) %>%
    dplyr::mutate(date = stringr::str_remove_all(basename(file), ".zip"))
}

#' Normalize a possibly leading-zero-stripped CAMEO code
#'
#' ICEWS event exports frequently store `CAMEO Code` as a numeric type,
#' which silently drops leading zeros (e.g. `"044"` becomes `44`). Pads
#' codes shorter than 3 characters back out with a leading zero.
#'
#' @param code Character or numeric vector of CAMEO codes.
#' @return Character vector, left-padded to at least 3 characters.
#' @keywords internal
normalize_cameo_code <- function(code) {
  code <- as.character(code)
  dplyr::if_else(stringr::str_length(code) < 3, paste0("0", code), code)
}

#' Read and filter ICEWS event export zips to relevant dyadic events
#'
#' Reads one or more ICEWS event export zips, filters to cross-country
#' dyadic events where both actors are state/security sectors (per
#' [relevant_actors_icews()]), and recodes CAMEO codes to QuadClass/
#' PentaClass via the package's own [cameo_lookup] (matching a code that
#' initially fails to match after a first round of leading-zero
#' normalization is retried with a second round, since ICEWS's numeric
#' CAMEO code field can drop more than one leading zero for some codes).
#' Output columns are renamed to match [extract_all_relevant_gdelt()]'s
#' schema (`Actor1CountryCode`, `Actor2CountryCode`, `SQLDATE`), so both
#' data sources can be aggregated with the same
#' [grouped_events_to_dyad_period()] / [assemble_stan_data()] pipeline.
#'
#' @param files Character vector of paths to ICEWS event export zip
#'   files.
#' @param relevant_sectors Character vector of sector-name substrings
#'   identifying relevant actors. Defaults to [relevant_actors_icews()].
#' @return A data frame of event-level records with `Actor1CountryCode`,
#'   `Actor2CountryCode`, `SQLDATE`, `QuadClass`, `PentaClass`,
#'   `PentaClass_modified`, and the original ICEWS columns.
#' @examples
#' \dontrun{
#' future::plan(future::multisession, workers = 8)
#' events <- ingest_icews(list.files("data/icews", pattern = "zip$", full.names = TRUE))
#' }
#' @export
ingest_icews <- function(files, relevant_sectors = relevant_actors_icews()) {
  raw <- furrr::future_map_dfr(files, function(f) {
    readr::read_delim(
      unz(f, unzip(f, list = TRUE)$Name[1]),
      num_threads = 1,
      progress = FALSE,
      show_col_types = FALSE
    )
  })

  sector_pattern <- paste(relevant_sectors, collapse = "|")
  filtered <- raw %>%
    dplyr::filter(`Source Country` != `Target Country`) %>%
    dplyr::filter(
      stringr::str_detect(`Source Sectors`, sector_pattern) &
        stringr::str_detect(`Target Sectors`, sector_pattern)
    ) %>%
    dplyr::mutate(cameo = normalize_cameo_code(`CAMEO Code`))

  matched <- recode_cameo(filtered, code_col = "cameo")
  still_unmatched <- matched %>% dplyr::filter(is.na(QuadClass))
  rematched <- still_unmatched %>%
    dplyr::select(-QuadClass, -PentaClass, -PentaClass_modified, -CAMEOLabel, -GoldsteinScore) %>%
    dplyr::mutate(cameo = normalize_cameo_code(cameo)) %>%
    recode_cameo(code_col = "cameo")

  dplyr::bind_rows(dplyr::filter(matched, !is.na(QuadClass)), rematched) %>%
    dplyr::rename(Actor1CountryCode = `Source Country`, Actor2CountryCode = `Target Country`) %>%
    dplyr::mutate(SQLDATE = as.integer(format(`Event Date`, "%Y%m%d")))
}
