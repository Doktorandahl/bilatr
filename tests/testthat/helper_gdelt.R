# Fixtures for tests/testthat/test_gdelt.R: two tiny hand-built GDELT
# zips (one daily-layout, 58 columns; one monthly-backfile-layout, 57
# columns) plus a matching md5sums index, written to a temp directory so
# download_gdelt()/read_gdelt() can be exercised offline via
# options(bilatr.gdelt_base_url = "file://<dir>/").

#' Twenty hand-built GDELT event rows covering every case Part 4 asks for
#'
#' Base row: a cross-border, root, GOV/GOV USA->CHN event. Rows 2-10
#' override specific fields to hit each required case (see inline
#' comments); rows 11-20 are filler, distinct only in `GLOBALEVENTID`/
#' `SOURCEURL`. Columns are in `.gdelt_schema()` order (58, including
#' `SOURCEURL`); callers drop the last column for a 57-column backfile.
#'
#' @param sqldate An 8-digit integer `SQLDATE` shared by every row.
#' @param id_offset Added to `GLOBALEVENTID`, so the two fixture files
#'   don't share event IDs.
#' @keywords internal
.make_gdelt_fixture_rows <- function(sqldate, id_offset) {
  n <- 20L
  id <- seq_len(n)

  base <- tibble::tibble(
    GLOBALEVENTID = id_offset + id,
    SQLDATE = sqldate,
    MonthYear = as.integer(sqldate %/% 100L),
    Year = as.integer(sqldate %/% 10000L),
    FractionDate = 2020.0027,
    Actor1Code = "USAGOV",
    Actor1Name = "UNITED STATES",
    Actor1CountryCode = "USA",
    Actor1KnownGroupCode = NA_character_,
    Actor1EthnicCode = NA_character_,
    Actor1Religion1Code = NA_character_,
    Actor1Religion2Code = NA_character_,
    Actor1Type1Code = "GOV",
    Actor1Type2Code = NA_character_,
    Actor1Type3Code = NA_character_,
    Actor2Code = "CHNGOV",
    Actor2Name = "CHINA",
    Actor2CountryCode = "CHN",
    Actor2KnownGroupCode = NA_character_,
    Actor2EthnicCode = NA_character_,
    Actor2Religion1Code = NA_character_,
    Actor2Religion2Code = NA_character_,
    Actor2Type1Code = "GOV",
    Actor2Type2Code = NA_character_,
    Actor2Type3Code = NA_character_,
    IsRootEvent = 1L,
    EventCode = "043",
    EventBaseCode = "043",
    EventRootCode = "04",
    QuadClass = 1L,
    GoldsteinScale = 1.9,
    NumMentions = 3L,
    NumSources = 2L,
    NumArticles = 2L,
    AvgTone = 0.5,
    Actor1Geo_Type = 1L,
    Actor1Geo_FullName = "United States",
    Actor1Geo_CountryCode = "US",
    Actor1Geo_ADM1Code = "US",
    Actor1Geo_Lat = 39.0,
    Actor1Geo_Long = -95.0,
    Actor1Geo_FeatureID = NA_character_,
    Actor2Geo_Type = 1L,
    Actor2Geo_FullName = "China",
    Actor2Geo_CountryCode = "CH",
    Actor2Geo_ADM1Code = "CH",
    Actor2Geo_Lat = 35.0,
    Actor2Geo_Long = 105.0,
    Actor2Geo_FeatureID = NA_character_,
    ActionGeo_Type = 1L,
    ActionGeo_FullName = "United States",
    ActionGeo_CountryCode = "US",
    ActionGeo_ADM1Code = "US",
    ActionGeo_Lat = 39.0,
    ActionGeo_Long = -95.0,
    ActionGeo_FeatureID = NA_character_,
    DATEADDED = as.numeric(sqldate) * 1e6,
    SOURCEURL = "http://example.com/1"
  )
  base <- base[rep(1L, n), ]
  base$GLOBALEVENTID <- id_offset + id
  base$SOURCEURL <- paste0("http://example.com/", id)

  # Row 2: cross-border RUS -> USA, MIL/SPY (Type1), non-root.
  base$Actor1CountryCode[2] <- "RUS"
  base$Actor2CountryCode[2] <- "USA"
  base$Actor1Type1Code[2] <- "MIL"
  base$Actor2Type1Code[2] <- "SPY"
  base$IsRootEvent[2] <- 0L
  base$EventCode[2] <- "190"
  base$EventBaseCode[2] <- "190"
  base$EventRootCode[2] <- "19"
  base$QuadClass[2] <- 4L

  # Row 3: domestic (Actor1CountryCode == Actor2CountryCode); leading-zero
  # EventCode "010".
  base$Actor1CountryCode[3] <- "USA"
  base$Actor2CountryCode[3] <- "USA"
  base$EventCode[3] <- "010"
  base$EventBaseCode[3] <- "010"
  base$EventRootCode[3] <- "01"

  # Row 4: NA Actor1CountryCode.
  base$Actor1CountryCode[4] <- NA_character_
  base$Actor1Code[4] <- "REB"
  base$Actor1Type1Code[4] <- "REB"

  # Row 5: NA Actor2CountryCode.
  base$Actor2CountryCode[5] <- NA_character_
  base$Actor2Code[5] <- "UIS"
  base$Actor2Type1Code[5] <- NA_character_

  # Row 6: a non-state Type1Code (BUS) on Actor1; second leading-zero
  # EventCode example, "0211".
  base$Actor1CountryCode[6] <- "USA"
  base$Actor2CountryCode[6] <- "RUS"
  base$Actor1Type1Code[6] <- "BUS"
  base$Actor2Type1Code[6] <- "GOV"
  base$EventCode[6] <- "0211"
  base$EventBaseCode[6] <- "021"
  base$EventRootCode[6] <- "02"
  base$IsRootEvent[6] <- 0L

  # Row 7: GOV in Type2Code (Actor1); MIL in Type3Code (Actor2).
  base$Actor1CountryCode[7] <- "CHN"
  base$Actor2CountryCode[7] <- "RUS"
  base$Actor1Type1Code[7] <- "ELI"
  base$Actor1Type2Code[7] <- "GOV"
  base$Actor2Type1Code[7] <- "BUS"
  base$Actor2Type2Code[7] <- "ELI"
  base$Actor2Type3Code[7] <- "MIL"
  base$EventCode[7] <- "051"
  base$EventBaseCode[7] <- "051"
  base$EventRootCode[7] <- "05"

  # Row 8: SPY in Type3Code (Actor1); non-root.
  base$Actor1CountryCode[8] <- "USA"
  base$Actor2CountryCode[8] <- "CHN"
  base$Actor1Type1Code[8] <- "ELI"
  base$Actor1Type2Code[8] <- "BUS"
  base$Actor1Type3Code[8] <- "SPY"
  base$Actor2Type1Code[8] <- "MIL"
  base$EventCode[8] <- "130"
  base$EventBaseCode[8] <- "130"
  base$EventRootCode[8] <- "13"
  base$IsRootEvent[8] <- 0L
  base$QuadClass[8] <- 3L

  # Row 9: explicit non-root, GOV/GOV.
  base$Actor1CountryCode[9] <- "RUS"
  base$Actor2CountryCode[9] <- "CHN"
  base$IsRootEvent[9] <- 0L
  base$EventCode[9] <- "036"
  base$EventBaseCode[9] <- "036"
  base$EventRootCode[9] <- "03"

  # Row 10: Actor1Name with an embedded double-quote (the readr
  # quote-mis-splitting fixture row).
  base$Actor1Name[10] <- 'John "Big Jim" Smith'
  base$Actor1CountryCode[10] <- "USA"
  base$Actor2CountryCode[10] <- "RUS"
  base$EventCode[10] <- "057"
  base$EventBaseCode[10] <- "057"
  base$EventRootCode[10] <- "05"

  base
}

#' Write GDELT fixture rows to a tab-separated file inside a zip
#'
#' @param rows Output of [.make_gdelt_fixture_rows()].
#' @param zip_path Destination `.zip` path.
#' @param n_cols 57 or 58 (drops `SOURCEURL` for 57).
#' @keywords internal
.write_gdelt_fixture_zip <- function(rows, zip_path, n_cols) {
  schema <- bilatr:::.gdelt_schema()
  cols <- schema$name[seq_len(n_cols)]
  types <- schema$type[seq_len(n_cols)]
  df <- rows[cols]

  formatted <- purrr::map2(df, types, function(x, type) {
    out <- switch(
      type,
      character = x,
      integer = sprintf("%.0f", x),
      double = sprintf("%.4f", x)
    )
    out[is.na(x)] <- ""
    out
  })
  lines <- do.call(paste, c(formatted, sep = "\t"))

  is_daily <- n_cols == 58L
  inner_name <- if (is_daily) {
    sub("\\.zip$", ".export.CSV", basename(zip_path))
  } else {
    sub("\\.zip$", ".csv", basename(zip_path))
  }

  tmp_dir <- tempfile("gdelt_fixture_")
  dir.create(tmp_dir)
  tmp_file <- file.path(tmp_dir, inner_name)
  writeLines(lines, tmp_file, useBytes = TRUE)

  old_wd <- setwd(tmp_dir)
  on.exit(setwd(old_wd), add = TRUE)
  status <- utils::zip(zip_path, inner_name, flags = "-q")
  if (!identical(status, 0L)) {
    stop("zip() failed building the GDELT test fixture; is a zip binary on PATH?")
  }
  invisible(zip_path)
}

#' Build the daily + monthly GDELT fixture zips and an md5sums index
#'
#' Skips the calling test (with a clear reason) if no `zip` binary is on
#' `PATH`. Writes into a fresh temp directory and returns its path plus
#' the two files' basenames and their row tibbles (pre-filter, full
#' 58-column shape, for hand-computing expected results in tests).
#'
#' @return A list: `dir`, `daily_file`, `monthly_file`, `daily_rows`,
#'   `monthly_rows`.
#' @keywords internal
.build_gdelt_fixtures <- function() {
  if (!nzchar(Sys.which("zip"))) {
    testthat::skip("no `zip` binary on PATH; cannot build GDELT test fixtures")
  }

  dir <- tempfile("gdelt_fixtures_")
  dir.create(dir)

  daily_rows <- .make_gdelt_fixture_rows(sqldate = 20200101L, id_offset = 900000000)
  monthly_rows <- .make_gdelt_fixture_rows(sqldate = 20100115L, id_offset = 800000000)

  daily_file <- "20200101.export.CSV.zip"
  monthly_file <- "201001.zip"

  .write_gdelt_fixture_zip(daily_rows, file.path(dir, daily_file), n_cols = 58L)
  .write_gdelt_fixture_zip(monthly_rows, file.path(dir, monthly_file), n_cols = 57L)

  md5 <- tools::md5sum(c(file.path(dir, daily_file), file.path(dir, monthly_file)))
  writeLines(
    sprintf("%s  %s", unname(md5), basename(names(md5))),
    file.path(dir, "md5sums")
  )

  list(
    dir = dir,
    daily_file = daily_file,
    monthly_file = monthly_file,
    daily_rows = daily_rows,
    monthly_rows = monthly_rows
  )
}
