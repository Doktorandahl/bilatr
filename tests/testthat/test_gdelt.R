# Tests for the 0.9.1 GDELT rewrite (dev/claude_code_prompt_0.9.1_gdelt.md,
# Parts 1-4): gdelt_files(), read_gdelt(), download_gdelt(),
# count_gdelt_actor_types(), gdelt_columns(). Offline only: no network, no
# CmdStan. Fixtures are built once per test via .build_gdelt_fixtures()
# (tests/testthat/helper_gdelt.R), which skip()s if no `zip` binary is on
# PATH.

# ---------------------------------------------------------------------------
# 1. gdelt_files()
# ---------------------------------------------------------------------------

test_that("gdelt_files() gets the yearly/monthly boundary right", {
  files <- gdelt_files("2005-12-31", "2006-01-01")
  expect_equal(files$period, c("2005", "200601"))
  expect_equal(files$type, c("yearly", "monthly"))
  expect_equal(files$file, c("2005.zip", "200601.zip"))
  expect_equal(files$n_cols, c(57L, 57L))
})

test_that("gdelt_files() gets the monthly/daily boundary right", {
  files <- gdelt_files("2013-03-31", "2013-04-01")
  expect_equal(files$period, c("201303", "20130401"))
  expect_equal(files$type, c("monthly", "daily"))
  expect_equal(files$file, c("201303.zip", "20130401.export.CSV.zip"))
  expect_equal(files$n_cols, c(57L, 58L))
})

test_that("gdelt_files() maps a range inside 1990 to a single yearly file", {
  files <- gdelt_files("1990-01-01", "1990-12-31")
  expect_equal(nrow(files), 1L)
  expect_equal(files$period, "1990")
  expect_equal(files$type, "yearly")
})

test_that("gdelt_files() maps a range inside 2008 to the covering monthly files", {
  files <- gdelt_files("2008-01-15", "2008-03-15")
  expect_equal(files$period, c("200801", "200802", "200803"))
  expect_true(all(files$type == "monthly"))
})

test_that("gdelt_files() handles a range spanning all three file types", {
  files <- gdelt_files("2005-12-01", "2013-04-02")
  expect_setequal(files$type, c("yearly", "monthly", "daily"))
  expect_equal(files$type[files$period == "2005"], "yearly")
  expect_equal(files$type[files$period == "200601"], "monthly")
  expect_equal(files$type[files$period == "20130401"], "daily")
})

test_that("gdelt_files() errors outside 1979-01-01..today", {
  expect_error(gdelt_files("1978-12-31"), "between")
  expect_error(gdelt_files(Sys.Date() + 1), "between")
  expect_error(gdelt_files("2020-06-01", "2019-01-01"), "before")
})

# ---------------------------------------------------------------------------
# 2. read_gdelt(): both layouts, types, quoting, binding
# ---------------------------------------------------------------------------

test_that("read_gdelt() reads the 58-column daily layout with correct types", {
  fx <- .build_gdelt_fixtures()
  events <- read_gdelt(
    file.path(fx$dir, fx$daily_file),
    columns = "all", actor_types = NULL, cross_border = FALSE
  )

  expect_equal(nrow(events), 20L)
  expect_type(events$EventCode, "character")
  expect_equal(events$EventCode[3], "010") # leading zero preserved
  expect_equal(events$EventCode[6], "0211")
  expect_type(events$GLOBALEVENTID, "double")
  expect_type(events$SQLDATE, "integer")
  expect_type(events$GoldsteinScale, "double")
  expect_false(any(is.na(events$SOURCEURL))) # daily file: SOURCEURL present
})

test_that("read_gdelt() reads the 57-column monthly layout, SOURCEURL all NA", {
  fx <- .build_gdelt_fixtures()
  events <- read_gdelt(
    file.path(fx$dir, fx$monthly_file),
    columns = "all", actor_types = NULL, cross_border = FALSE
  )

  expect_equal(nrow(events), 20L)
  expect_true(all(is.na(events$SOURCEURL)))
  expect_equal(events$EventCode[3], "010")
})

test_that("read_gdelt() survives an embedded double-quote in a free-text field", {
  fx <- .build_gdelt_fixtures()
  events <- read_gdelt(
    file.path(fx$dir, fx$daily_file),
    columns = "all", actor_types = NULL, cross_border = FALSE
  )
  expect_equal(events$Actor1Name[10], 'John "Big Jim" Smith')
})

test_that("read_gdelt() binds the 57- and 58-column layouts cleanly", {
  fx <- .build_gdelt_fixtures()
  both <- read_gdelt(
    c(file.path(fx$dir, fx$daily_file), file.path(fx$dir, fx$monthly_file)),
    columns = "all", actor_types = NULL, cross_border = FALSE
  )
  expect_equal(nrow(both), 40L)
  expect_setequal(both$source_file, c(fx$daily_file, fx$monthly_file))
  expect_equal(sum(is.na(both$SOURCEURL)), 20L)
})

# ---------------------------------------------------------------------------
# 3. Filters and *_match modes
# ---------------------------------------------------------------------------

test_that("cross_border TRUE drops domestic and NA-country rows; FALSE keeps everything", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  cb_true <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = TRUE)
  cb_false <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = FALSE)

  expect_equal(nrow(cb_false), 20L)
  # rows 3 (domestic), 4 (NA actor1), 5 (NA actor2) dropped
  expect_false(900000003 %in% cb_true$GLOBALEVENTID)
  expect_false(900000004 %in% cb_true$GLOBALEVENTID)
  expect_false(900000005 %in% cb_true$GLOBALEVENTID)
  expect_equal(nrow(cb_true), 17L)
})

test_that("actor_type_match 'both' requires both sides to match; 'either' needs only one", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  # Row 6: Actor1Type1Code = BUS (non-state), Actor2Type1Code = GOV.
  both <- read_gdelt(f, columns = "all", cross_border = TRUE, actor_type_match = "both")
  either <- read_gdelt(f, columns = "all", cross_border = TRUE, actor_type_match = "either")

  expect_false(900000006 %in% both$GLOBALEVENTID)
  expect_true(900000006 %in% either$GLOBALEVENTID)
})

test_that("actor_types checks all three type slots per side, not just Type1Code", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  # Row 7: GOV in Actor1Type2Code, MIL in Actor2Type3Code (both slots > 1).
  # Row 8: SPY in Actor1Type3Code, MIL in Actor2Type1Code.
  out <- read_gdelt(f, columns = "all", cross_border = TRUE, actor_type_match = "both")
  expect_true(900000007 %in% out$GLOBALEVENTID)
  expect_true(900000008 %in% out$GLOBALEVENTID)
})

test_that("actor_types = NULL turns the actor-type filter off", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)
  out <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = TRUE)
  expect_equal(nrow(out), 17L) # only the cross-border filter applies
})

test_that("countries + country_match 'either'/'both' behave as documented", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  either <- read_gdelt(
    f,
    columns = "all", actor_types = NULL, cross_border = FALSE,
    countries = "RUS", country_match = "either"
  )
  both <- read_gdelt(
    f,
    columns = "all", actor_types = NULL, cross_border = FALSE,
    countries = c("RUS", "USA"), country_match = "both"
  )

  # "either": every row involving RUS on either side.
  expect_true(all(either$Actor1CountryCode == "RUS" | either$Actor2CountryCode == "RUS"))
  expect_true(nrow(either) > 0)
  # "both": the closed set -- both actors from {RUS, USA} (which, since
  # cross_border = FALSE here, also admits the USA<->USA domestic row).
  expect_true(all(both$Actor1CountryCode %in% c("RUS", "USA") & both$Actor2CountryCode %in% c("RUS", "USA")))
  expect_true(any(both$Actor1CountryCode != both$Actor2CountryCode)) # at least one real RUS<->USA pair
})

test_that("root_events_only keeps IsRootEvent == 1 rows only", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)
  out <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = FALSE, root_events_only = TRUE)
  expect_true(all(out$IsRootEvent == 1L))
  expect_false(900000002 %in% out$GLOBALEVENTID) # row 2 is non-root
})

test_that("event_codes keeps rows whose EventCode starts with any given prefix", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  root19 <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = FALSE, event_codes = "19")
  expect_true(all(startsWith(root19$EventCode, "19")))
  expect_true(900000002 %in% root19$GLOBALEVENTID) # EventCode "190"

  exact <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = FALSE, event_codes = "0211")
  expect_equal(exact$EventCode, "0211")
})

test_that("date_range filters on SQLDATE", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$monthly_file) # all rows SQLDATE = 20100115
  in_range <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = FALSE, date_range = c("2010-01-01", "2010-01-31"))
  out_of_range <- read_gdelt(f, columns = "all", actor_types = NULL, cross_border = FALSE, date_range = c("2011-01-01", "2011-01-31"))
  expect_equal(nrow(in_range), 20L)
  expect_equal(nrow(out_of_range), 0L)
})

test_that("filtering on a column not requested in `columns` still works", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)
  # Actor1CountryCode/Actor2CountryCode aren't in this `columns` set, but
  # cross_border filtering must still apply.
  out <- read_gdelt(f, columns = c("GLOBALEVENTID", "EventCode"), actor_types = NULL, cross_border = TRUE)
  expect_equal(names(out), c("GLOBALEVENTID", "EventCode", "source_file"))
  expect_equal(nrow(out), 17L)
})

test_that("unknown `columns`/`actor_types`/empty `files` all error clearly", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)
  expect_error(read_gdelt(f, columns = c("NotAColumn")), "unknown column")
  expect_error(read_gdelt(f, actor_types = c("ZZZ")), "unknown `actor_types`")
  expect_error(read_gdelt(character(0)), "non-empty")
})

# ---------------------------------------------------------------------------
# 4. Parity with 0.9.0's extract_all_relevant_gdelt()
# ---------------------------------------------------------------------------

test_that("the default read_gdelt() call matches extract_all_relevant_gdelt()'s (0.9.0) filter exactly", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  # extract_all_relevant_gdelt() (0.9.0, R/data_ingestion.R, now removed):
  # cross-border (both actor country codes non-NA and different) AND both
  # actors carry GOV/MIL/SPY in one of Type1/2/3Code -- computed by hand
  # here from the fixture's own rows, since the old function is gone.
  relevant_actors <- c("GOV", "MIL", "SPY")
  expected <- fx$daily_rows %>%
    dplyr::filter(
      !is.na(Actor1CountryCode) & !is.na(Actor2CountryCode),
      Actor1CountryCode != Actor2CountryCode,
      (Actor1Type1Code %in% relevant_actors | Actor1Type2Code %in% relevant_actors | Actor1Type3Code %in% relevant_actors) &
        (Actor2Type1Code %in% relevant_actors | Actor2Type2Code %in% relevant_actors | Actor2Type3Code %in% relevant_actors)
    ) %>%
    dplyr::pull(GLOBALEVENTID)

  actual <- read_gdelt(f, columns = "all")$GLOBALEVENTID
  expect_setequal(actual, expected)
})

# ---------------------------------------------------------------------------
# 5. download_gdelt(), via a file:// fixture base URL
# ---------------------------------------------------------------------------

test_that("download_gdelt() temp mode returns events and leaves tempdir() untouched", {
  fx <- .build_gdelt_fixtures()
  withr_base_url <- getOption("bilatr.gdelt_base_url")
  options(bilatr.gdelt_base_url = paste0("file://", fx$dir, "/"))
  on.exit(options(bilatr.gdelt_base_url = withr_base_url), add = TRUE)

  events <- download_gdelt("2020-01-01", actor_types = NULL, cross_border = FALSE)

  expect_equal(nrow(events), 20L)
  expect_equal(attr(events, "files")$status, "downloaded")
  # download_gdelt() names its temp file after the real GDELT filename
  # (so read_gdelt()'s daily/backfile column-count sniffing works); that
  # exact path must not survive the call.
  expect_false(file.exists(file.path(tempdir(), "20200101.export.CSV.zip")))
  expect_false(file.exists(file.path(tempdir(), "20200101.export.CSV.zip.part")))
})

test_that("download_gdelt() cache mode writes zips, and a second call reports 'cached'", {
  fx <- .build_gdelt_fixtures()
  withr_base_url <- getOption("bilatr.gdelt_base_url")
  options(bilatr.gdelt_base_url = paste0("file://", fx$dir, "/"))
  on.exit(options(bilatr.gdelt_base_url = withr_base_url), add = TRUE)

  cache_dir <- tempfile("gdelt_cache_")
  status1 <- download_gdelt("2010-01-15", dest_dir = cache_dir)
  expect_equal(status1$status, "downloaded")
  expect_true(file.exists(status1$path))

  status2 <- download_gdelt("2010-01-15", dest_dir = cache_dir)
  expect_equal(status2$status, "cached")
})

test_that("download_gdelt() reports 'missing' with one warning for a period with no fixture", {
  fx <- .build_gdelt_fixtures()
  withr_base_url <- getOption("bilatr.gdelt_base_url")
  options(bilatr.gdelt_base_url = paste0("file://", fx$dir, "/"))
  on.exit(options(bilatr.gdelt_base_url = withr_base_url), add = TRUE)

  cache_dir <- tempfile("gdelt_cache_")
  warnings <- character(0)
  status <- withCallingHandlers(
    download_gdelt("1985-06-15", dest_dir = cache_dir),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(status$status, "missing")
  expect_true(is.na(status$path))
  expect_length(warnings, 1L)
  expect_match(warnings, "not available")
})

test_that("download_gdelt() reports 'corrupt' for a fixture with a wrong md5", {
  fx <- .build_gdelt_fixtures()
  bad_dir <- tempfile("gdelt_bad_")
  dir.create(bad_dir)
  file.copy(file.path(fx$dir, fx$monthly_file), file.path(bad_dir, fx$monthly_file))
  writeLines("not the right content", file.path(bad_dir, fx$monthly_file))
  writeLines(readLines(file.path(fx$dir, "md5sums")), file.path(bad_dir, "md5sums"))

  withr_base_url <- getOption("bilatr.gdelt_base_url")
  options(bilatr.gdelt_base_url = paste0("file://", bad_dir, "/"))
  on.exit(options(bilatr.gdelt_base_url = withr_base_url), add = TRUE)

  cache_dir <- tempfile("gdelt_cache_")
  status <- suppressWarnings(download_gdelt("2010-01-15", dest_dir = cache_dir))
  expect_equal(status$status, "corrupt")
  expect_true(is.na(status$path))
})

test_that("download_gdelt() never leaves a .part file behind, success or failure", {
  fx <- .build_gdelt_fixtures()
  withr_base_url <- getOption("bilatr.gdelt_base_url")
  options(bilatr.gdelt_base_url = paste0("file://", fx$dir, "/"))
  on.exit(options(bilatr.gdelt_base_url = withr_base_url), add = TRUE)

  cache_dir <- tempfile("gdelt_cache_")
  suppressWarnings(download_gdelt("2010-01-15", dest_dir = cache_dir)) # succeeds
  suppressWarnings(download_gdelt("1985-06-15", dest_dir = cache_dir)) # missing

  part_files <- list.files(cache_dir, pattern = "\\.part$", full.names = TRUE)
  expect_length(part_files, 0L)
})

# ---------------------------------------------------------------------------
# 6. count_gdelt_actor_types()
# ---------------------------------------------------------------------------

test_that("count_gdelt_actor_types() counts combinations over all three type slots per side", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  counts <- count_gdelt_actor_types(f, cross_border = FALSE)
  expect_true(all(c("actor1_type", "actor2_type", "n") %in% names(counts)))
  expect_true(is.unsorted(-counts$n) == FALSE) # sorted descending by n

  # Row 7 (Actor1 in {ELI, GOV}, Actor2 in {BUS, ELI, MIL}) and Row 8
  # (Actor1 in {ELI, BUS, SPY}, Actor2 in {MIL}) both contribute an
  # (ELI, MIL) pair -- only visible if all three slots per side are used.
  eli_mil <- counts$n[counts$actor1_type == "ELI" & counts$actor2_type == "MIL"]
  expect_equal(eli_mil, 2L)
})

# ---------------------------------------------------------------------------
# 7. End-to-end: read_gdelt() -> recode_cameo() -> validate/assemble
# ---------------------------------------------------------------------------

test_that("read_gdelt() output passes validate_bilatr_events() after recode_cameo(), and assembles", {
  fx <- .build_gdelt_fixtures()
  f <- file.path(fx$dir, fx$daily_file)

  events <- read_gdelt(f, columns = gdelt_columns("core"), actor_types = NULL, cross_border = TRUE)
  events <- recode_cameo(events, code_col = "EventCode")
  expect_true(all(!is.na(events$PentaClass)))

  expect_silent(validate_bilatr_events(events, grouping_var = "PentaClass"))

  sd <- assemble_stan_data(
    events,
    years = 2020, resolution = "yearly", grouping_var = "PentaClass", min_n_events = 1
  )
  expect_gt(sd$D, 0L)
})
