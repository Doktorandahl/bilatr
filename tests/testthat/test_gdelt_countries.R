# GDELT country codes (1.0.0f): gdelt_countries, the lookup helpers,
# filter_gdelt_states(), gdelt_state_spells(), and read_gdelt()'s
# recover_country_codes/states_only/harmonize_codes. No network.

# --- names ------------------------------------------------------------------

test_that("gdelt_country_name() names states, regions and territories", {
  expect_identical(
    gdelt_country_name(c("USA", "SWE", "ROM", "COL", "EUR", "WAF", "WSB", "HKG", "ROU", "SSD")),
    c("United States", "Sweden", "Romania", "Colombia", "Europe", "West Africa",
      "West Bank", "Hong Kong", "Romania", "South Sudan")
  )
})

test_that("unknown codes give NA with one warning naming them; NA and \"\" give NA silently", {
  expect_warning(
    out <- gdelt_country_name(c("USA", "XXX", "ZZZ", "XXX")),
    "2 code\\(s\\) not in `gdelt_countries`.*XXX, ZZZ"
  )
  expect_identical(out, c("United States", NA, NA, NA))
  expect_no_warning(out <- gdelt_country_name(c(NA, "", "USA")))
  expect_identical(out, c(NA, NA, "United States"))
  expect_identical(gdelt_country_name(character(0)), character(0))
  expect_identical(gdelt_country_name(NA), NA_character_)
  expect_identical(gdelt_country_name(factor("SWE")), "Sweden")
  expect_error(gdelt_country_name(1), "must be a character vector")
})

# --- types ------------------------------------------------------------------

test_that("gdelt_country_type() classifies regions, territories, states and others", {
  expect_identical(
    gdelt_country_type(c("EUR", "WAF", "WSB", "SJM", "USA", "GME", "VAT", NA)),
    c("region", "region", "territory", "territory", "state", "former_state", "other", NA)
  )
  expect_warning(expect_identical(gdelt_country_type("QQQ"), NA_character_), "QQQ")
})

# --- time-varying states ------------------------------------------------------

expect_state_years <- function(code, years, state, iso3) {
  expect_identical(gdelt_is_state(code, years), state, label = paste(code, "is_state"))
  expect_identical(gdelt_country_iso3(code, years), iso3, label = paste(code, "iso3"))
}

test_that("Serbia is SRB, and a state, in every year (COW 345)", {
  expect_state_years("SRB", c(1995, 2004, 2010), c(TRUE, TRUE, TRUE), c("SRB", "SRB", "SRB"))
  expect_state_years("SCG", c(1991, 1992, 2006, 2007), c(FALSE, TRUE, TRUE, FALSE), c(NA, "SRB", "SRB", NA))
})

test_that("successor sequences give the right state and ISO3 either side of each transition", {
  # Czechoslovakia -> Czech Republic / Slovakia
  expect_state_years("CZE", c(1992, 1993), c(TRUE, TRUE), c("CSK", "CZE"))
  expect_state_years("SVK", c(1992, 1993), c(FALSE, TRUE), c(NA, "SVK"))
  # Germany
  expect_state_years("DEU", c(1989, 1990, 1991), c(TRUE, TRUE, TRUE), c("DEU", "DEU", "DEU"))
  expect_state_years("GMW", c(1990, 1991), c(TRUE, FALSE), c("DEU", NA))
  expect_state_years("GME", c(1990, 1991), c(TRUE, FALSE), c("DDR", NA))
  # Yemen
  expect_state_years("YEM", c(1989, 1990), c(TRUE, TRUE), c("YEM", "YEM"))
  # USSR -> Russia (one COW state, 365)
  expect_state_years("RUS", c(1991, 1992), c(TRUE, TRUE), c("RUS", "RUS"))
  # Yugoslav successors
  expect_state_years("HRV", c(1991, 1992), c(FALSE, TRUE), c(NA, "HRV"))
  expect_state_years("MKD", c(1992, 1993), c(FALSE, TRUE), c(NA, "MKD"))
  # Sudan / South Sudan; Eritrea; Timor-Leste; Namibia; Palau
  expect_state_years("SDN", c(2010, 2012), c(TRUE, TRUE), c("SDN", "SDN"))
  expect_state_years("SSD", c(2010, 2011), c(FALSE, TRUE), c(NA, "SSD"))
  expect_state_years("ERI", c(1992, 1993), c(FALSE, TRUE), c(NA, "ERI"))
  expect_state_years("TMP", c(2001, 2002), c(FALSE, TRUE), c(NA, "TLS"))
  expect_state_years("NAM", c(1989, 1990), c(FALSE, TRUE), c(NA, "NAM"))
  expect_state_years("PLW", c(1993, 1994), c(FALSE, TRUE), c(NA, "PLW"))
  # Taiwan is a COW member (D-f3)
  expect_state_years("TWN", c(1991, 2025), c(TRUE, TRUE), c("TWN", "TWN"))
})

test_that("the six originally missing states are states in their membership years, with ISO3 codes", {
  years <- 1991:2025
  expect_true(all(gdelt_is_state("SRB", years)))
  expect_true(all(gdelt_country_iso3("SRB", years) == "SRB"))
  expect_true(all(gdelt_is_state("ROM", years)))
  expect_true(all(gdelt_country_iso3(c("ROM", "ROU"), 2000) == "ROU"))
  expect_identical(gdelt_is_state("BIH", years), years >= 1992)
  expect_identical(gdelt_is_state("SVN", years), years >= 1992)
  expect_identical(gdelt_is_state("MTN", years), years >= 2006)
  expect_identical(gdelt_country_iso3("MTN", 2010), "MNE")
  expect_identical(gdelt_is_state("TMP", years), years >= 2002)
  expect_identical(gdelt_country_iso3("TMP", 2010), "TLS")
})

test_that("gdelt_is_state()/gdelt_country_iso3() without year, and argument handling", {
  expect_identical(gdelt_is_state(c("USA", "EUR", "SVK", "GME", NA, "")), c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE))
  expect_identical(gdelt_country_iso3(c("ROM", "EUR", "CZE")), c("ROU", NA, "CZE"))
  expect_warning(expect_identical(gdelt_is_state(c("USA", "QQQ"), 2000), c(TRUE, NA)), "QQQ")
  expect_identical(gdelt_is_state("USA", c(2000, NA)), c(TRUE, NA))
  expect_identical(gdelt_is_state(c("USA", "SVK"), 1992), c(TRUE, FALSE))
  expect_error(gdelt_is_state(c("USA", "SVK", "CZE"), c(1992, 1993)), "same length")
  expect_error(gdelt_is_state("USA", 1992.5), "whole-number years")
  expect_error(gdelt_is_state("USA", "1992"), "whole-number years")
})

# --- filter_gdelt_states() ----------------------------------------------------

state_events <- function() {
  tibble::tibble(
    GLOBALEVENTID = 1:8,
    Actor1CountryCode = c("USA", "ROM", "EUR", "SVK", "SVK", "ROM", NA, "USA"),
    Actor2CountryCode = c("CHN", "HUN", "RUS", "CZE", "CZE", "ROU", "USA", "WSB"),
    SQLDATE = c(20200101L, 19950601L, 20200101L, 19920301L, 19930301L, 20000101L, 20200101L, 20200101L),
    EventCode = "043"
  )
}

test_that("filter_gdelt_states() keeps state pairs in the event year, harmonizes, and drops self-dyads", {
  ev <- state_events()
  expect_message(out <- filter_gdelt_states(ev), "kept 3 of 8 rows; dropped 5")
  expect_identical(out$GLOBALEVENTID, c(1L, 2L, 5L))
  expect_identical(out$Actor1CountryCode, c("USA", "ROU", "SVK"))
  expect_identical(names(out), names(ev))
  expect_s3_class(out, "tbl_df")
  expect_identical(lapply(out, class), lapply(ev, class))
})

test_that("filter_gdelt_states()'s message reports reasons, dropped codes and harmonized codes", {
  msg <- testthat::capture_messages(filter_gdelt_states(state_events()))
  expect_length(msg, 1)
  expect_match(msg, "4 with an actor that is not a state in the event's year")
  expect_match(msg, "1 that harmonization made self-dyads")
  expect_match(msg, "Dropped codes \\(actor sides\\): <empty> 1, EUR 1, SVK 1, WSB 1\\.")
  expect_match(msg, "Harmonized codes \\(actor sides\\): ROM -> ROU \\(1\\)\\.")
})

test_that("the dropped-code list shows the ten most frequent codes, then how many others", {
  regions <- c("AFR", "ASA", "BLK", "CAS", "CAU", "CEU", "CFR", "CRB", "EAF", "EEU", "EIN", "EUR")
  ev <- tibble::tibble(
    Actor1CountryCode = c(rep("AFR", 3), regions[-1]),
    Actor2CountryCode = "USA",
    SQLDATE = 20200101L
  )
  msg <- testthat::capture_messages(filter_gdelt_states(ev))
  expect_match(msg, "AFR 3, ASA 1")
  expect_match(msg, "and 2 other code\\(s\\) \\(2\\)")
})

test_that("harmonize = FALSE keeps GDELT's codes and keeps two-code dyads of one state", {
  out <- suppressMessages(filter_gdelt_states(state_events(), harmonize = FALSE))
  expect_identical(out$GLOBALEVENTID, c(1L, 2L, 5L, 6L))
  expect_identical(out$Actor1CountryCode, c("USA", "ROM", "SVK", "ROM"))
})

test_that("filter_gdelt_states() works on a data.frame, custom columns and Date dates; drops bad dates", {
  df <- data.frame(a = c("ROM", "USA", "USA"), b = c("USA", "CHN", "CHN"),
                   d = as.Date(c("2000-01-01", "2000-01-01", NA)))
  expect_message(out <- filter_gdelt_states(df, actor1 = "a", actor2 = "b", date = "d"), "1 with an unparseable date")
  expect_s3_class(out, "data.frame")
  expect_identical(out$a, c("ROU", "USA"))
})

test_that("filter_gdelt_states() errors on missing columns and bad arguments", {
  expect_error(filter_gdelt_states(state_events(), actor1 = "nope"), "column\\(s\\) not found in `data`: `nope`")
  expect_error(filter_gdelt_states(dplyr::select(state_events(), -"SQLDATE")), "`SQLDATE`")
  expect_error(filter_gdelt_states(list()), "must be a data frame")
  expect_error(filter_gdelt_states(state_events(), harmonize = NA), "TRUE or FALSE")
})

# --- read_gdelt() integration ---------------------------------------------------

test_that("read_gdelt(states_only = TRUE) equals filter_gdelt_states() on read_gdelt(states_only = FALSE)", {
  fx <- .build_gdelt_country_fixtures()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  for (recover in c(FALSE, TRUE)) {
    all_rows <- suppressMessages(read_gdelt(fx$files, recover_country_codes = recover))
    expected <- suppressMessages(filter_gdelt_states(all_rows))
    msg <- testthat::capture_messages(
      got <- read_gdelt(fx$files, recover_country_codes = recover, states_only = TRUE)
    )
    expect_length(msg, 1L + recover)
    expect_match(msg, "read_gdelt\\(\\): kept", all = FALSE)
    expect_identical(got, expected)
  }
  expect_identical(
    suppressMessages(read_gdelt(fx$files, states_only = TRUE, harmonize_codes = FALSE)),
    suppressMessages(filter_gdelt_states(read_gdelt(fx$files), harmonize = FALSE))
  )
})

test_that("a cross-border ROM/ROU event passes cross_border on raw codes and is dropped after harmonization", {
  fx <- .build_gdelt_country_fixtures()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  raw <- read_gdelt(fx$files[1])
  expect_true(any(raw$Actor1CountryCode == "ROM" & raw$Actor2CountryCode %in% "ROU"))
  msg <- testthat::capture_messages(states <- read_gdelt(fx$files[1], states_only = TRUE))
  expect_match(msg, "1 that harmonization made self-dyads", all = FALSE)
  expect_false(any(states$Actor1CountryCode == states$Actor2CountryCode))
  expect_false(any(states$GLOBALEVENTID == 600000024))
  expect_true(600000023 %in% states$GLOBALEVENTID) # ROM -> HUN, harmonized
  expect_identical(states$Actor1CountryCode[states$GLOBALEVENTID == 600000023], "ROU")
})

test_that("states_only drops regions, territories, non-members and out-of-membership years", {
  fx <- .build_gdelt_country_fixtures()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  daily <- suppressMessages(read_gdelt(fx$files[1], states_only = TRUE))
  monthly <- suppressMessages(read_gdelt(fx$files[2], states_only = TRUE))
  # EUR, WSB, VAT out in both; TMP (TLS from 2002) in both; SSD (2011) only in 2020
  for (id in c(21, 22, 28)) {
    expect_false((600000000 + id) %in% daily$GLOBALEVENTID)
    expect_false((500000000 + id) %in% monthly$GLOBALEVENTID)
  }
  expect_true(600000026 %in% daily$GLOBALEVENTID)
  expect_true(500000026 %in% monthly$GLOBALEVENTID)
  expect_identical(daily$Actor1CountryCode[daily$GLOBALEVENTID == 600000026], "TLS")
  expect_true(600000027 %in% daily$GLOBALEVENTID)
  expect_false(500000027 %in% monthly$GLOBALEVENTID)
})

test_that("recover_country_codes fills empty codes from a known actor-code prefix, and reports them", {
  fx <- .build_gdelt_country_fixtures()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  plain <- read_gdelt(fx$files[1])
  expect_false(any(plain$GLOBALEVENTID %in% c(600000025, 600000029)))
  msg <- testthat::capture_messages(rec <- read_gdelt(fx$files[1], recover_country_codes = TRUE))
  expect_length(msg, 1)
  expect_match(msg, "filled 2 empty actor country code\\(s\\).*ROU 1, TWN 1")
  expect_identical(rec$Actor1CountryCode[rec$GLOBALEVENTID == 600000025], "ROU")
  expect_identical(rec$Actor2CountryCode[rec$GLOBALEVENTID == 600000029], "TWN")
  # a role-only actor code is not a country prefix
  expect_false(600000030 %in% rec$GLOBALEVENTID)
  expect_false(any(c(".filled1", ".filled2") %in% names(rec)))
  expect_identical(names(rec), names(plain))
  # recovery runs before cross_border, so recovered rows count as cross-border
  expect_identical(nrow(rec), nrow(plain) + 2L)
})

test_that("count_gdelt_actor_types() takes the same country-code steps", {
  fx <- .build_gdelt_country_fixtures()
  on.exit(unlink(fx$dir, recursive = TRUE), add = TRUE)
  msg <- testthat::capture_messages(
    state_counts <- count_gdelt_actor_types(fx$files[1], recover_country_codes = TRUE, states_only = TRUE)
  )
  expect_match(msg, "count_gdelt_actor_types\\(\\): kept", all = FALSE)
  expect_match(msg, "filled 2 empty", all = FALSE)
  rows <- suppressMessages(read_gdelt(
    fx$files[1], columns = "all", actor_types = NULL, recover_country_codes = TRUE, states_only = TRUE
  ))
  expected <- bilatr:::.gdelt_count_actor_types_chunk(rows) %>%
    dplyr::group_by(.data$actor1_type, .data$actor2_type) %>%
    dplyr::summarise(n = sum(.data$n), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$n))
  expect_identical(
    dplyr::arrange(state_counts, actor1_type, actor2_type),
    dplyr::arrange(expected, actor1_type, actor2_type)
  )
  expect_error(read_gdelt(fx$files[1], states_only = "yes"), "TRUE or FALSE")
})

# --- gdelt_state_spells() -----------------------------------------------------

test_that("gdelt_state_spells() returns one row per spell keyed by ISO3", {
  sp <- gdelt_state_spells()
  expect_named(sp, c("iso3c", "cown", "start", "end"))
  expect_s3_class(sp$start, "Date")
  expect_s3_class(sp$end, "Date")
  expect_type(sp$cown, "integer")
  expect_false(anyNA(sp$start))
  deu <- sp[sp$iso3c == "DEU", ]
  expect_identical(deu$cown, c(260L, 255L))
  expect_identical(deu$end[1], as.Date("1990-10-02"))
  expect_identical(deu$start[2], as.Date("1990-10-03"))
  expect_identical(sp$end[sp$iso3c == "CSK"], as.Date("1992-12-31"))
  expect_identical(sp$start[sp$iso3c == "SVK"], as.Date("1993-01-01"))
  expect_identical(sp$start[sp$iso3c == "HRV"], as.Date("1992-01-15"))
  expect_true(is.na(sp$end[sp$iso3c == "USA"]))
  # COW 345 is one spell (YUG and SRB rows merged)
  expect_identical(nrow(sp[sp$iso3c == "SRB", ]), 1L)
  expect_false("XKX" %in% sp$iso3c)
  wy <- gdelt_state_spells(whole_years = TRUE)
  expect_identical(wy$start[wy$iso3c == "HRV"], as.Date("1992-01-01"))
  expect_identical(wy$end[wy$iso3c == "DDR"], as.Date("1990-12-31"))
  expect_error(gdelt_state_spells(NA), "TRUE or FALSE")
})

test_that("the spells and the year-level rule agree on which years each state code is a state in", {
  sp <- gdelt_state_spells()
  codes <- gdelt_countries$gdelt_code[gdelt_countries$entity_type %in% c("state", "former_state")]
  years <- 1946:2025
  grid <- expand.grid(code = codes, year = years, stringsAsFactors = FALSE)
  grid$iso3 <- gdelt_country_iso3(grid$code, grid$year)
  pairs <- unique(grid[!is.na(grid$iso3), c("iso3", "year")])
  spell_years <- do.call(rbind, lapply(seq_len(nrow(sp)), function(i) {
    from <- max(as.integer(format(sp$start[i], "%Y")), 1946L)
    to <- if (is.na(sp$end[i])) 2025L else as.integer(format(sp$end[i], "%Y"))
    data.frame(iso3 = sp$iso3c[i], year = from:to)
  }))
  spell_years <- unique(spell_years[spell_years$year %in% years, ])
  key <- function(d) sort(paste(d$iso3, d$year))
  expect_identical(key(pairs), key(spell_years))
})

# --- table integrity ------------------------------------------------------------

test_that("gdelt_countries is internally consistent", {
  gc <- gdelt_countries
  expect_named(gc, c("gdelt_code", "name", "name_source", "entity_type", "iso3c", "cown",
                     "cow_start", "cow_end", "note"))
  expect_false(anyDuplicated(gc$gdelt_code) > 0)
  expect_true(all(grepl("^[A-Z]{3}$", gc$gdelt_code)))
  expect_setequal(unique(gc$entity_type), c("state", "former_state", "territory", "region", "other"))
  is_state <- gc$entity_type %in% c("state", "former_state")
  expect_false(anyNA(gc$iso3c[is_state]))
  expect_false(anyNA(gc$cown[is_state]))
  expect_false(anyNA(gc$cow_start[is_state]))
  expect_true(all(is.na(gc$cow_end[gc$entity_type == "state"])))
  expect_false(anyNA(gc$cow_end[gc$entity_type == "former_state"]))
  expect_true(all(is.na(gc$cown[!is_state])))
  expect_true(all(is.na(gc$iso3c[!is_state])))
  expect_false(anyNA(gc$name))
  expect_identical(sum(gc$entity_type == "state"), 195L)
  expect_true(all(c("SRB", "ROU", "TLS") %in% gc$iso3c))
  expect_identical(gc$entity_type[gc$gdelt_code == "VAT"], "other")
  # every state code's years are reachable through the year-level periods
  expect_true(all(gdelt_is_state(gc$gdelt_code[is_state], gc$cow_start[is_state])))
})
