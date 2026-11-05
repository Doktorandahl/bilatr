# Classify every GDELT 1.0 country code (1.0.0f, Part 1b).
#
# Builds one row per GDELT country code -- the union of the official CAMEO
# country lookup, every Actor1/Actor2CountryCode value observed in the full
# GDELT 1.0 archive 1991-2025, and the country-like actor-code prefixes
# that GDELT leaves without a country code -- plus a long table of state
# periods used to resolve successor cases by year.
#
# Sources (read at build time):
#   * CAMEO country lookup, GDELT Project:
#     https://www.gdeltproject.org/data/lookups/CAMEO.country.txt
#     (CODE, LABEL; 261 codes; CRLF, tab-separated).
#   * Correlates of War Project, State System Membership, v2024
#     (statelist2024.csv inside States2024.zip):
#     https://correlatesofwar.org/wp-content/uploads/States2024.zip
#     Membership runs through 2024-12-31 for current members; a state whose
#     last spell ends on that date is treated as a current member (cow_end
#     NA). v2024 splits code 345 into Yugoslavia (1944-10-20 .. 2006-06-03)
#     and Serbia (from 2006-06-03); states2016 did not.
#   * countrycode::codelist (countrycode >= 1.6.1): iso3c, cown,
#     country.name.en. Used only here; not a bilatr dependency.
#   * data-raw/gdelt_country_codes/{code,prefix}_counts.csv, written by
#     data-raw/gdelt_country_codes_inventory.R from the full archive.
#
# Rules (decisions D-f1..D-f3, 2026-11-05):
#   * D-f1: a code is a state in year t when it maps to a COW system member
#     at any point during t ("in the year of the event", the year-level
#     convention of COW's own system-year file).
#   * D-f2: state codes are harmonized to ISO 3166-1 alpha-3. One ISO3
#     code per COW state number and period; where COW keeps one number
#     through a renaming (345: Yugoslavia -> Serbia and Montenegro ->
#     Serbia; 365: USSR -> Russia) the current code is used throughout.
#   * D-f3: contested cases follow D-f1 (Taiwan and Kosovo are COW members;
#     Palestine, Western Sahara and the Vatican are not). Kosovo (COW 347,
#     from 2008-02-20) has no GDELT code, so it has no row here.
#   * In a transition year the period holding on 31 December wins (e.g.
#     Germany 1990 is COW 255, not 260).
#   * cow_end NA = a current member. COW v2024 ends on 2024-12-31; its
#     current members are treated as members through 2025 (and after).
#
# Review decisions (dev/reply_1.0.0f_review.md, 2026-11-05): R1 (actor-
# prefix recovery in read_gdelt()), E1 (00b drops the class-0 codes), Y1
# (SRB for COW 345 in every year), C1 (CZE before 1993 -> CSK, COW 315),
# K1 (Kosovo is a COW state; docs only), S4 (cow_end NA = ongoing). All
# were accepted as proposed in dev/review_1.0.0f_country_codes.md, so the
# rules below are unchanged by the review.
#
# Writes:
#   * data-raw/gdelt_country_codes/gdelt_countries.csv (with row counts),
#     gdelt_state_periods.csv and gdelt_state_spells.csv;
#   * data/gdelt_countries.rda (the shipped table, ?gdelt_countries);
#   * R/sysdata.rda: .gdelt_state_periods (year-level code -> ISO3/COW
#     periods, for gdelt_is_state()/gdelt_country_iso3()) and
#     .gdelt_state_spells (day-level spells by ISO3, for
#     gdelt_state_spells()).
# Run from the package root.

library(dplyr)
library(tibble)
library(readr)
library(tidyr)

out_dir <- file.path("data-raw", "gdelt_country_codes")
first_year <- 1946L

# --- sources ----------------------------------------------------------------

lookup_url <- "https://www.gdeltproject.org/data/lookups/CAMEO.country.txt"
lookup <- utils::read.delim(lookup_url, colClasses = "character", strip.white = TRUE) |>
  as_tibble() |>
  rename(gdelt_code = CODE, label = LABEL)
stopifnot(nrow(lookup) == 261, !anyDuplicated(lookup$gdelt_code))

cow_url <- "https://correlatesofwar.org/wp-content/uploads/States2024.zip"
cow_zip <- tempfile(fileext = ".zip")
utils::download.file(cow_url, cow_zip, mode = "wb", quiet = TRUE)
cow_csv <- grep("^[^_].*/statelist2024\\.csv$", utils::unzip(cow_zip, list = TRUE)$Name, value = TRUE)
stopifnot(length(cow_csv) == 1)
cow <- utils::read.csv(unz(cow_zip, cow_csv)) |> as_tibble()
cow_data_end <- max(cow$endyear)
stopifnot(cow_data_end == 2024)

cl <- countrycode::codelist |>
  select(iso3c, cown, country.name.en) |>
  filter(!is.na(iso3c))

code_counts <- read_csv(file.path(out_dir, "code_counts.csv"),
                        na = "", show_col_types = FALSE)
prefix_counts <- read_csv(file.path(out_dir, "prefix_counts.csv"),
                          na = "", show_col_types = FALSE)

# --- COW spells from 1946 ---------------------------------------------------

# Day-level COW spells, for gdelt_state_spells(): spells ending before
# 1946 are dropped; a spell running to the end of the data is open (NA).
cow_day <- cow |>
  transmute(
    cown = as.integer(ccode),
    start = as.Date(sprintf("%04d-%02d-%02d", styear, stmonth, stday)),
    end = as.Date(sprintf("%04d-%02d-%02d", endyear, endmonth, endday))
  ) |>
  filter(end >= as.Date(sprintf("%d-01-01", first_year)))
stopifnot(!anyNA(cow_day$start), !anyNA(cow_day$end))

spells <- cow |>
  filter(endyear >= first_year) |>
  transmute(
    cown = as.integer(ccode), cow_abb = stateabb, cow_name = statenme,
    start = pmax(as.integer(styear), first_year),
    end = if_else(endyear == cow_data_end & endmonth == 12 & endday == 31,
                  NA_integer_, as.integer(endyear))
  )

# Post-1946 spells of the same state number, merged when the gap between
# them is under a calendar year (none in v2024 after 1946), since membership is
# counted per year.
spells <- spells |>
  arrange(cown, start) |>
  group_by(cown) |>
  mutate(new_spell = c(TRUE, start[-1] > lag(coalesce(end, 9999L))[-1] + 1L),
         spell_id = cumsum(new_spell)) |>
  group_by(cown, spell_id) |>
  summarise(cow_abb = last(cow_abb), cow_name = last(cow_name),
            start = min(start),
            end = if (anyNA(end)) NA_integer_ else max(end),
            .groups = "drop") |>
  select(-spell_id)

# --- observed codes ---------------------------------------------------------

observed <- code_counts |>
  filter(!is.na(code)) |>
  group_by(gdelt_code = code) |>
  summarise(n_country_code = sum(n), .groups = "drop")

# Country-like actor-code prefixes on rows with an empty country code
# (inventory: ROU, TWN, SSD, GGY, GME, GMW, SCG); keep those that are not
# role/religion codes.
country_like <- union(cl$iso3c, c(lookup$gdelt_code, "GME", "GMW", "SCG"))
prefix_only <- prefix_counts |>
  filter(prefix %in% country_like, year >= 1979) |>
  group_by(gdelt_code = prefix) |>
  summarise(n_actor_prefix = sum(n), .groups = "drop") |>
  filter(n_actor_prefix >= 50)

codes <- tibble(gdelt_code = sort(unique(c(lookup$gdelt_code, observed$gdelt_code,
                                           prefix_only$gdelt_code)))) |>
  left_join(lookup, by = "gdelt_code") |>
  left_join(observed, by = "gdelt_code") |>
  left_join(prefix_only, by = "gdelt_code") |>
  mutate(across(c(n_country_code, n_actor_prefix), \(x) coalesce(x, 0L)))

# --- hand rules -------------------------------------------------------------

regions <- c(
  AFR = "Africa", ASA = "Asia", BLK = "Balkans", CAS = "Central Asia",
  CAU = "Caucasus", CEU = "Central Europe", CFR = "Central Africa",
  CRB = "Caribbean", EAF = "Eastern Africa", EEU = "Eastern Europe",
  EIN = "East Indies", EUR = "Europe", LAM = "Latin America",
  MDT = "Mediterranean", MEA = "Middle East", NAF = "North Africa",
  NMR = "North America", PGS = "Persian Gulf", SAF = "Southern Africa",
  SAM = "South America", SAS = "South Asia", SCN = "Scandinavia",
  SEA = "Southeast Asia", WAF = "West Africa", WST = "The West"
)

# Names: the lookup label, except for these corrections.
name_fixes <- c(
  COL = "Colombia",                 # lookup: "Columbia"
  REU = "Reunion",                  # lookup: "Runion" (lost accent)
  COG = "Republic of the Congo"     # lookup: "People's Republic of the Congo"
)
# Names for codes not in the lookup.
name_added <- c(
  SSD = "South Sudan", ROU = "Romania", TWN = "Taiwan",
  GGY = "Guernsey", GME = "East Germany", GMW = "West Germany",
  SCG = "Serbia and Montenegro"
)

# Codes whose COW mapping is set by hand rather than via countrycode:
# gdelt_code -> periods (from, to, cown, iso3c). `to = NA` means ongoing.
hand_periods <- tribble(
  ~gdelt_code, ~from, ~to, ~cown, ~iso3c,
  # CAMEO-only codes for current states
  "ROM", 1946L, NA, 360L, "ROU",
  "ROU", 1946L, NA, 360L, "ROU",
  "TMP", 2002L, NA, 860L, "TLS",
  "MTN", 2006L, NA, 341L, "MNE",
  # COW 345 (Yugoslavia 1944-2006, Serbia since): one state, one code
  "SRB", 1946L, NA, 345L, "SRB",
  "SCG", 1992L, 2006L, 345L, "SRB",
  # Germany: GDELT writes DEU for West Germany before 1990 (COW 260)
  "DEU", 1955L, 1989L, 260L, "DEU",
  "DEU", 1990L, NA, 255L, "DEU",
  "GMW", 1955L, 1990L, 260L, "DEU",
  "GME", 1954L, 1990L, 265L, "DDR",
  # Czechoslovakia: GDELT writes CZE (Prague) before 1993 too (COW 315)
  "CZE", 1946L, 1992L, 315L, "CSK",
  "CZE", 1993L, NA, 316L, "CZE",
  # Yemen: YEM before unification is North Yemen (COW 678)
  "YEM", 1946L, 1989L, 678L, "YEM",
  "YEM", 1990L, NA, 679L, "YEM"
)
not_states <- c(VAT = "other", PSE = "territory", ESH = "territory",
                WSB = "territory", GZS = "territory", BAG = "other")

notes <- c(
  ROM = "CAMEO code for Romania; harmonized to ISO3 ROU. GDELT writes ROM only rarely (about 11,000 rows 1991-2025): Romanian actors are coded ROU... in Actor1/2Code, a segment GDELT's country table lacks, so their country code is usually empty.",
  ROU = "Not a GDELT country code: the ISO3 code for Romania, which GDELT writes as the actor-code prefix (ROU...) but leaves out of Actor1/2CountryCode. Listed so harmonized data and actor-prefix recovery resolve.",
  TMP = "CAMEO code for East Timor; harmonized to ISO3 TLS. A COW state from 2002 (independence 2002-05-20; COW member from 2002-09-27); before that an occupied/UN-administered territory, so not a state.",
  MTN = "CAMEO code for Montenegro; harmonized to ISO3 MNE. In the lookup but never written by GDELT 1.0 (no Montenegrin actor segment in its dictionaries). A COW state from 2006 (COW 341).",
  SRB = "Serbia. COW treats Yugoslavia (1944-2006), Serbia and Montenegro, and Serbia as one state, 345; GDELT writes SRB (actor codes SCGSRB...) throughout. Harmonized to SRB in every year, so the series is continuous. Before 1992, SRB (the Serbian republic) is mapped to COW 345 (SFR Yugoslavia).",
  SCG = "Not a GDELT country code: the actor-code prefix for Serbia and Montenegro (SCG...) when no SRB segment follows, which leaves the country code empty. COW 345 (FRY 1992-2003, Serbia and Montenegro 2003-2006); harmonized to SRB; not a state outside 1992-2006.",
  DEU = "Germany. Before 1990 GDELT's DEU refers to West Germany (COW 260, member 1955-1990); from 1990 COW 255. ISO3 DEU throughout.",
  GMW = "Not a GDELT country code: the actor-code prefix for West Germany, with an empty country code. COW 260 (1955-1990); harmonized to DEU.",
  GME = "Not a GDELT country code: the actor-code prefix for East Germany, with an empty country code. COW 265 (1954-1990); ISO3 DDR (historical).",
  CZE = "Czech Republic from 1993 (COW 316). GDELT has no Czechoslovakia code and writes CZE for Czechoslovak actors (e.g. 'Prague'), so CZE before 1993 is Czechoslovakia (COW 315, ISO3 CSK).",
  SVK = "Slovakia, a COW state from 1993 (COW 317); not a state before.",
  YEM = "Yemen. Before unification (1990-05-22) GDELT's YEM is taken as North Yemen (COW 678), the predecessor in place; South Yemen (COW 680) has no GDELT code.",
  RUS = "Russia. COW 365 covers the USSR and Russia as one state; GDELT has no USSR code (Soviet actors are coded RUS). ISO3 RUS throughout.",
  VNM = "Vietnam, COW 816 (North Vietnam from 1954); GDELT has no code for the Republic of Vietnam (COW 817, 1954-1975).",
  SDN = "Sudan, COW 625, continuous through South Sudan's secession.",
  SSD = "Not in the CAMEO lookup, but written by GDELT 1.0's daily files (and as an actor-code prefix, with an empty country code, in the backfiles). A COW state from 2011 (COW 626).",
  ERI = "Eritrea, a COW state from 1993 (COW 531); part of Ethiopia before.",
  TWN = "Taiwan, a COW state (713) under D-f3. In the backfiles (to March 2013) Taiwanese actors usually have an empty country code and only the actor-code prefix TWN.",
  BIH = "Bosnia and Herzegovina, COW 346 from 1992. In the lookup but never written by GDELT 1.0: its dictionaries have no Bosnian actor segment (Bosnian actors are coded by role only, e.g. GOV, MIL).",
  SVN = "Slovenia, COW 349 from 1992. In the lookup but never written by GDELT 1.0 (no Slovenian actor segment).",
  HRV = "Croatia, a COW state from 1992 (COW 344, member from 1992-01-15); not a state in 1991.",
  MKD = "North Macedonia (lookup: 'Macedonia'), a COW state from 1993 (COW 343).",
  COD = "Democratic Republic of the Congo; GDELT also uses COD for Zaire (1971-1997). COW 490.",
  MMR = "Myanmar; GDELT also uses MMR for Burma. COW 775.",
  COG = "Republic of the Congo; the lookup's label 'People's Republic of the Congo' is the 1970-1991 name. COW 484.",
  VAT = "Vatican City: sovereign, but not a COW system member (D-f3). countrycode's cown 327 is the Papal States (1816-1860) and is not used.",
  PSE = "Occupied Palestinian Territory: not a COW member (D-f3).",
  ESH = "Western Sahara: not a COW member (D-f3). Never written by GDELT 1.0.",
  WSB = "West Bank. CAMEO code; never written by GDELT 1.0.",
  GZS = "Gaza Strip. CAMEO code; never written by GDELT 1.0.",
  BAG = "Baghdad: a city code in the CAMEO lookup; never written by GDELT 1.0.",
  HKG = "Hong Kong: a special administrative region of China, not a COW member.",
  MAC = "Macao: a special administrative region of China, not a COW member.",
  IMY = "CAMEO code for the Isle of Man (ISO3 IMN); a Crown dependency. Never written by GDELT 1.0.",
  ANT = "Netherlands Antilles (dissolved 2010; ISO3 ANT withdrawn). Never written by GDELT 1.0.",
  GGY = "Not a GDELT country code: the actor-code prefix GGY (Guernsey, ISO3) with an empty country code.",
  SJM = "Svalbard and Jan Mayen: Norwegian territory.",
  PRI = "Puerto Rico: US territory.",
  GRL = "Greenland: part of the Kingdom of Denmark. Never written by GDELT 1.0.",
  SWZ = "Eswatini (lookup: 'Swaziland'), COW 572.",
  PLW = "Palau, a COW state from 1994.",
  NAM = "Namibia, a COW state from 1990."
)

# --- assemble ----------------------------------------------------------------

auto <- codes |>
  filter(!gdelt_code %in% c(hand_periods$gdelt_code, names(regions), names(not_states))) |>
  left_join(cl, by = c("gdelt_code" = "iso3c")) |>
  mutate(cown = as.integer(cown))

auto_periods <- auto |>
  filter(!is.na(cown)) |>
  inner_join(spells, by = "cown", relationship = "many-to-many") |>
  transmute(gdelt_code, from = start, to = end, cown, iso3c = gdelt_code)

# countrycode cowns that are not COW members after 1946 would be a silent
# error: there must be none apart from VAT (handled above).
stopifnot(all(auto$cown[!is.na(auto$cown)] %in% spells$cown))

periods <- bind_rows(hand_periods, auto_periods) |>
  arrange(gdelt_code, from)

table_state <- periods |>
  group_by(gdelt_code) |>
  summarise(
    iso3c = last(iso3c), cown = last(cown),
    cow_start = min(from),
    cow_end = if (anyNA(to)) NA_integer_ else max(to),
    .groups = "drop"
  ) |>
  mutate(entity_type = if_else(is.na(cow_end), "state", "former_state"))

gdelt_countries <- codes |>
  left_join(table_state, by = "gdelt_code") |>
  left_join(select(cl, iso3c_cc = iso3c, country.name.en), by = c("gdelt_code" = "iso3c_cc")) |>
  mutate(
    entity_type = case_when(
      !is.na(entity_type) ~ entity_type,
      gdelt_code %in% names(regions) ~ "region",
      gdelt_code %in% names(not_states) ~ not_states[gdelt_code],
      gdelt_code %in% c("IMY", "ANT", "GGY") ~ "territory",
      gdelt_code %in% cl$iso3c ~ "territory",
      TRUE ~ "other"
    ),
    name = case_when(
      gdelt_code %in% names(name_fixes) ~ name_fixes[gdelt_code],
      !is.na(label) ~ label,
      gdelt_code %in% names(name_added) ~ name_added[gdelt_code],
      TRUE ~ country.name.en
    ),
    name_source = case_when(
      gdelt_code %in% names(name_fixes) ~ "cameo_lookup_corrected",
      !is.na(label) ~ "cameo_lookup",
      TRUE ~ "added"
    ),
    note = unname(notes[gdelt_code]),
    note = case_when(
      !is.na(note) ~ note,
      entity_type == "region" ~ "Region or supranational grouping (CAMEO).",
      entity_type == "territory" ~ "Non-state territory: an ISO 3166-1 code with no COW state.",
      TRUE ~ NA_character_
    )
  ) |>
  select(gdelt_code, name, name_source, entity_type, iso3c, cown,
         cow_start, cow_end, note, n_country_code, n_actor_prefix) |>
  arrange(gdelt_code)

gdelt_state_periods <- periods |>
  rename(from_year = from, to_year = to) |>
  arrange(gdelt_code, from_year)

# Day-level spells by harmonized ISO3 code: every (cown, iso3c) pair the
# periods use, with that state number's COW spells. Spells of one state
# number are merged when they are contiguous at year level (the next
# starts no later than the calendar year after the previous ends), the
# same rule as the year-level merge above, so the year-level and
# day-level views agree; the YUG/SRB rows of COW 345 (which share
# 2006-06-03) become one spell.
data_end <- as.Date(sprintf("%d-12-31", cow_data_end))
gdelt_state_spells <- periods |>
  distinct(iso3c, cown) |>
  inner_join(cow_day, by = "cown", relationship = "many-to-many") |>
  arrange(iso3c, cown, start) |>
  group_by(iso3c, cown) |>
  mutate(
    prev_end_year = lag(as.integer(format(end, "%Y"))),
    new_spell = is.na(prev_end_year) | as.integer(format(start, "%Y")) > prev_end_year + 1L,
    spell_id = cumsum(new_spell)
  ) |>
  group_by(iso3c, cown, spell_id) |>
  summarise(start = min(start), end = max(end), .groups = "drop") |>
  mutate(end = if_else(end == data_end, as.Date(NA), end)) |>
  select(iso3c, cown, start, end) |>
  arrange(iso3c, start)

# --- checks ------------------------------------------------------------------

stopifnot(
  !anyDuplicated(gdelt_countries$gdelt_code),
  all(gdelt_countries$entity_type %in% c("state", "former_state", "territory", "region", "other")),
  all(!is.na(gdelt_countries$iso3c[gdelt_countries$entity_type %in% c("state", "former_state")])),
  all(!is.na(gdelt_countries$cown[gdelt_countries$entity_type %in% c("state", "former_state")])),
  all(is.na(gdelt_countries$cown[!gdelt_countries$entity_type %in% c("state", "former_state")])),
  all(grepl("^[A-Z]{3}$", gdelt_countries$gdelt_code)),
  # periods do not overlap within a code
  gdelt_state_periods |>
    group_by(gdelt_code) |>
    summarise(ok = all(head(coalesce(to_year, 9999L), -1) < tail(from_year, -1))) |>
    pull(ok) |> all(),
  # every state ISO3 code has a spell, and every spell an ISO3 code in use
  setequal(gdelt_state_spells$iso3c, gdelt_state_periods$iso3c),
  !anyNA(gdelt_state_spells$start)
)

# The year-level periods and the day-level spells must agree on which
# years each ISO3 code is a state in (S2: a member at any point in the
# year), over 1946..2025.
years <- first_year:2025L
year_view <- gdelt_state_periods |>
  rowwise() |>
  reframe(iso3c, year = intersect(years, from_year:coalesce(to_year, 2025L))) |>
  distinct()
spell_view <- gdelt_state_spells |>
  mutate(from = pmax(as.integer(format(start, "%Y")), first_year),
         to = coalesce(as.integer(format(end, "%Y")), 2025L)) |>
  rowwise() |>
  reframe(iso3c, year = intersect(years, from:to)) |>
  distinct()
year_mismatch <- full_join(mutate(year_view, periods = TRUE),
                           mutate(spell_view, spells = TRUE),
                           by = c("iso3c", "year")) |>
  filter(is.na(periods) | is.na(spells))
if (nrow(year_mismatch) > 0) {
  print(year_mismatch, n = Inf)
  stop("year-level periods and day-level spells disagree (above)")
}

write_csv(gdelt_countries, file.path(out_dir, "gdelt_countries.csv"), na = "")
write_csv(gdelt_state_periods, file.path(out_dir, "gdelt_state_periods.csv"), na = "")
write_csv(gdelt_state_spells, file.path(out_dir, "gdelt_state_spells.csv"), na = "")

# --- ship ---------------------------------------------------------------------

# R/ must stay ASCII: the shipped strings are checked here.
ship <- gdelt_countries |>
  select(gdelt_code, name, name_source, entity_type, iso3c, cown,
         cow_start, cow_end, note) |>
  mutate(across(c(cown, cow_start, cow_end), as.integer)) |>
  as.data.frame() |>
  tibble::as_tibble()
stopifnot(all(!grepl("[^ -~]", unlist(ship[c("name", "note")]), perl = TRUE) |
                is.na(unlist(ship[c("name", "note")]))))
gdelt_countries <- ship
save(gdelt_countries, file = file.path("data", "gdelt_countries.rda"),
     compress = "xz", version = 2)

.gdelt_state_periods <- gdelt_state_periods |>
  mutate(across(c(from_year, to_year, cown), as.integer)) |>
  as.data.frame()
.gdelt_state_spells <- gdelt_state_spells |>
  mutate(cown = as.integer(cown)) |>
  as.data.frame()
.gdelt_cow_version <- "Correlates of War State System Membership v2024 (data end 2024-12-31)"
save(.gdelt_state_periods, .gdelt_state_spells, .gdelt_cow_version,
     file = file.path("R", "sysdata.rda"), compress = "xz", version = 2)

message(sprintf("%d codes: %s", nrow(gdelt_countries),
                paste(names(table(gdelt_countries$entity_type)),
                      table(gdelt_countries$entity_type), sep = " = ", collapse = ", ")))
