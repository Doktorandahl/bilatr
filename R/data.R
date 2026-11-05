#' CAMEO event code lookup table
#'
#' A reference table mapping every CAMEO event code to its Goldstein
#' conflict-cooperation score, a human-readable label, and the
#' QuadClass/PentaClass/ModifiedRootCode recoding used throughout the
#' package. Built by `data-raw/build_cameo_lookup.R`; see
#' [recode_cameo()] for how it is applied to event data.
#'
#' @format A data frame with one row per CAMEO event code and columns:
#' \describe{
#'   \item{CAMEOEVENTCODE}{Character. The CAMEO event code, e.g. `"0211"`.}
#'   \item{CAMEOLabel}{Character. Human-readable description of the event
#'     code.}
#'   \item{GoldsteinScore}{Numeric. Goldstein conflict-cooperation score,
#'     ranging from -10 (most conflictual) to 10 (most cooperative).}
#'   \item{QuadClass}{Integer, 1-4. The same mapping as GDELT's own
#'     `QuadClass` column: 1 = verbal cooperation (roots 01-05),
#'     2 = material cooperation (06-09), 3 = verbal conflict (10-14),
#'     4 = material conflict (15-20).}
#'   \item{QuadClassName}{Character. Human-readable label for
#'     `QuadClass`.}
#'   \item{QuadClassEventCodes}{Character. Comma-separated CAMEO root
#'     codes underlying each `QuadClass` class.}
#'   \item{PentaClass}{Integer, 0-4. As QuadClass, except that roots
#'     01-02 (make public statement, appeal) form class 0, root 14
#'     (protest) moves to class 4 (material conflict), and root 16
#'     (reduce relations) moves to class 3 (verbal conflict): 0 = make
#'     statement (roots 01-02), 1 = verbal cooperation (03-05),
#'     2 = material cooperation (06-09), 3 = verbal conflict (10-13, 16),
#'     4 = material conflict (14, 15, 17-20).}
#'   \item{PentaClassName}{Character. Human-readable label for
#'     `PentaClass`.}
#'   \item{PentaClassEventCodes}{Character. Comma-separated CAMEO root
#'     codes underlying each `PentaClass` class.}
#'   \item{ModifiedRootCode}{Integer, 1-18. A regrouping of the CAMEO
#'     root codes that moves `016` (deny responsibility) into Reject
#'     (class 13) and `018`/`019` (make empathetic comment, express
#'     accord) into diplomatic cooperation (class 7); splits root `04`
#'     (consult) into unspecified (`040`), meet, discuss, or visit
#'     (`041`-`044`) and negotiate or mediate (`045`-`046`); merges
#'     Investigate (`09`) and Demand (`10`); merges Threaten (`13`) and
#'     Exhibit force posture (`15`); and merges `18`-`20` into "Assault,
#'     fight, or mass violence".
#'     Blank (`NA`) on the top-level (two-digit) row for a root whose own
#'     event codes don't all belong to the same `ModifiedRootCode` class
#'     (roots `01` and `04`); populated on that root's individual
#'     sub-codes as normal.}
#'   \item{ModifiedRootCodeName}{Character. Human-readable label for
#'     `ModifiedRootCode`. `NA`
#'     wherever `ModifiedRootCode` is blank.}
#'   \item{ModifiedRootCodeEventCodes}{Character. Comma-separated CAMEO
#'     root/event codes underlying each `ModifiedRootCode` class. `NA`
#'     wherever `ModifiedRootCode`
#'     is blank.}
#' }
#'
#' [event_class_labels()] gives a class's human-readable label from a
#' fitted model's `stan_data`, for each scheme with a `*Name` column
#' here (`QuadClass`, `PentaClass`, `ModifiedRootCode`).
#' @source Built by `data-raw/build_cameo_lookup.R` in the source
#'   repository. CAMEO event codes and labels from Schrodt (2012), *CAMEO
#'   Conflict and Mediation Event Observations Event and Actor Codebook*,
#'   version 1.1b3. Goldstein scale from Goldstein (1992), "A
#'   Conflict-Cooperation Scale for WEIS Events Data".
#' @examples
#' # the root-level rows of the verbal-conflict PentaClass
#' subset(cameo_lookup, PentaClass == 3 & nchar(CAMEOEVENTCODE) == 2,
#'   select = c(CAMEOEVENTCODE, CAMEOLabel, GoldsteinScore))
#'
#' # how the 310 codes spread over QuadClass and PentaClass
#' table(cameo_lookup$QuadClassName, cameo_lookup$PentaClass)
"cameo_lookup"

#' GDELT country codes: names, entity types and state-system membership
#'
#' One row per country code that GDELT 1.0 can write in
#' `Actor1CountryCode`/`Actor2CountryCode`: every code in GDELT's CAMEO
#' country lookup, every code observed in the full GDELT 1.0 event archive
#' for 1991-2025, and seven actor-code prefixes that GDELT leaves without a
#' country code (see "Codes GDELT does not write"). The helpers
#' [gdelt_country_name()], [gdelt_country_type()], [gdelt_is_state()],
#' [gdelt_country_iso3()], [filter_gdelt_states()] and
#' [gdelt_state_spells()] are built on it.
#'
#' @format A tibble with 267 rows and 9 columns:
#' \describe{
#'   \item{gdelt_code}{Character. The code as GDELT writes it (three
#'     upper-case letters); unique.}
#'   \item{name}{Character. GDELT's label, cleaned; see `name_source`.}
#'   \item{name_source}{Character. `"cameo_lookup"` (GDELT's label as
#'     is), `"cameo_lookup_corrected"` (three corrections: `COL`
#'     "Columbia" to "Colombia", `REU` "Runion" to "Reunion", `COG`
#'     "People's Republic of the Congo" to "Republic of the Congo") or
#'     `"added"` (a code not in the lookup).}
#'   \item{entity_type}{Character; see "Entity types".}
#'   \item{iso3c}{Character. The ISO 3166-1 alpha-3 code a state code is
#'     harmonized to (its latest one, for a code with successor periods),
#'     or `NA` for a non-state code.}
#'   \item{cown}{Integer. The Correlates of War (COW) state number for the
#'     code's latest period, or `NA`.}
#'   \item{cow_start, cow_end}{Integer. The first and last year, from 1946
#'     on, in which the code is a state system member; `cow_end` is `NA`
#'     for a current member. A code can map to different states in
#'     different years (see "Successor rules"); [gdelt_is_state()] and
#'     [gdelt_country_iso3()] resolve them by year.}
#'   \item{note}{Character. The reason for every decision that is not
#'     automatic; `NA` for plain cases.}
#' }
#'
#' @section Entity types:
#' * `"state"` (195 codes): a current member of the COW state system. The
#'   195 codes cover 194 of COW's 195 current members; Kosovo has no GDELT
#'   code.
#' * `"former_state"` (3): a code that is a state only in earlier years:
#'   `GME` (East Germany, 1954-1990), `GMW` (West Germany, 1955-1990) and
#'   `SCG` (Serbia and Montenegro, 1992-2006).
#' * `"territory"` (42): an ISO 3166-1 code with no COW state, e.g. `HKG`,
#'   `MAC`, `PRI`, `SJM`, `GRL`, and the Palestinian codes `PSE`, `WSB`
#'   (West Bank) and `GZS` (Gaza Strip), and Western Sahara (`ESH`).
#' * `"region"` (25): a continent, sub-region or other grouping, e.g. `EUR`
#'   (Europe), `AFR`, `WAF` (West Africa), `MEA` (Middle East), `WST` (the
#'   West).
#' * `"other"` (2): `VAT` (Vatican City, sovereign but not a COW member)
#'   and `BAG` (Baghdad, a city code in the lookup).
#'
#' @section States:
#' A code is a state in year \eqn{t} when it maps to a COW state system
#' member at any point during \eqn{t}, so the rule is time-varying:
#' Slovakia (`SVK`) is a state from 1993, Timor-Leste (`TMP`) from 2002,
#' South Sudan (`SSD`) from 2011. Two spells of one COW state number are
#' merged when they are contiguous at year level. In a year with two
#' candidate periods for one code, the period holding on 31 December wins
#' (Germany in 1990 is COW 255, not 260). The COW data (v2024) end on
#' 2024-12-31; a state that is a member on that date is a current member
#' (`cow_end = NA`), taken as a member through 2025 and after.
#'
#' @section Harmonization:
#' State codes are harmonized to ISO 3166-1 alpha-3, so that dyad keys
#' match the [theta_trajectories] tables: `ROM` (CAMEO's code for Romania)
#' becomes `ROU`, `TMP` becomes `TLS` and `MTN` becomes `MNE`; every other
#' current-state code already is its ISO3 code. One code is used per COW
#' state number and period. Where COW keeps one number through a change
#' of name, the current code is used throughout: `SRB` for COW 345
#' (Yugoslavia, Serbia and Montenegro, Serbia) in every year, and `RUS`
#' for COW 365 (the Soviet Union and Russia).
#'
#' @section Successor rules:
#' * **Serbia** (`SRB`, `SCG`): COW 345 is one state from 1944 on; both
#'   codes harmonize to `SRB`, so Serbia's series is continuous. `SCG` is a
#'   state only in 1992-2006.
#' * **Czechoslovakia** (`CZE`): GDELT has no Czechoslovakia code and
#'   writes `CZE` for Czechoslovak actors too, so `CZE` before 1993 is
#'   Czechoslovakia (COW 315, ISO3 `CSK`) and from 1993 the Czech Republic
#'   (316). Slovakia (`SVK`, 317) is a state from 1993.
#' * **Germany** (`DEU`, `GMW`, `GME`): before 1990 GDELT's `DEU` is West
#'   Germany (COW 260); from 1990, Germany (255). Both harmonize to `DEU`.
#'   `GMW` (West Germany, 260) is `DEU` and `GME` (East Germany, 265) is
#'   `DDR`, through 1990.
#' * **Yemen** (`YEM`): before 1990, North Yemen (COW 678); from 1990,
#'   Yemen (679). South Yemen (680) has no GDELT code.
#' * **Russia** (`RUS`): COW 365 throughout; GDELT has no Soviet Union
#'   code.
#' * **Others:** Croatia (`HRV`) and Bosnia and Herzegovina (`BIH`) are
#'   states from 1992, North Macedonia (`MKD`) and Eritrea (`ERI`) from
#'   1993, Palau (`PLW`) from 1994, Montenegro (`MTN`) from 2006; Sudan
#'   (`SDN`) is continuous through South Sudan's secession; `COD` covers
#'   Zaire and `MMR` covers Burma.
#'
#' @section Contested cases:
#' Contested and borderline cases follow the COW rule. **Taiwan** (`TWN`,
#' COW 713) is a state. **Kosovo** (COW 347, a member from 2008-02-20) is a
#' state, but GDELT 1.0 has no code for it. **Palestine** (`PSE`, `WSB`,
#' `GZS`) and **Western Sahara** (`ESH`) are territories. **Vatican City**
#' (`VAT`) is not a COW member; the COW number 327 that some crosswalks
#' give it is the Papal States (1816-1860).
#'
#' @section Codes GDELT does not write:
#' GDELT 1.0 derives the country code from the leading country segment of
#' the actor code (`Actor1Code`), and leaves the country code empty when
#' its own country table lacks that segment. Seven such prefixes are
#' listed here, with `name_source = "added"` where GDELT has no label:
#' `ROU` (Romania: actors are coded `ROU...`, but GDELT's table has
#' `ROM`), `TWN` (Taiwan, in the backfiles to March 2013), `SSD` (South
#' Sudan, in the backfiles), `SCG`, `GME`, `GMW` and `GGY` (Guernsey).
#' [read_gdelt()]'s `recover_country_codes = TRUE` fills an empty country
#' code from such a prefix.
#'
#' Four states have no code GDELT 1.0 ever writes, as a country code or
#' as an actor-code prefix, in any year: **Bosnia and Herzegovina**,
#' **Slovenia** and **Montenegro** (`BIH`, `SVN` and `MTN` are in the
#' lookup but never occur) and **Kosovo** (no code at all). Their actors
#' are coded by role only (e.g. `GOV`, `MIL`), so their events cannot be
#' assigned to them from the country codes, and state-to-state event
#' tables built from GDELT 1.0 do not include them.
#'
#' @source Built by `data-raw/gdelt_country_codes.R` in the source
#'   package, from:
#'   * the GDELT Project's CAMEO country lookup,
#'     <https://gdeltproject.org/data/lookups/CAMEO.country.txt>
#'     (261 codes, accessed 2026-11-05);
#'   * Correlates of War Project (2024), State System Membership, v2024,
#'     <https://correlatesofwar.org/data-sets/state-system-membership/>;
#'   * Arel-Bundock, V., Enevoldsen, N. and Yetman, C. (2018).
#'     countrycode: An R package to convert country names and country
#'     codes. *Journal of Open Source Software* 3(28), 848
#'     (`countrycode::codelist`, version 1.6.1, for ISO3 codes, COW
#'     numbers and names);
#'   * the country codes in the full GDELT 1.0 event archive, 1991-2025
#'     (`data-raw/gdelt_country_codes_inventory.R`).
#' @seealso [filter_gdelt_states()], [gdelt_state_spells()].
#' @examples
#' table(gdelt_countries$entity_type)
#' subset(gdelt_countries, entity_type == "former_state", gdelt_code:cow_end)
#' @keywords datasets
"gdelt_countries"

#' Estimated bilateral-relations trajectories, 1991-2025
#'
#' Yearly latent trajectories (\eqn{\theta}) of the level of bilateral
#' relations for every country pair with at least one recorded event,
#' estimated with bilatr from GDELT 1.0 event data. The data are not
#' shipped with the package: [download_bilatr_theta()] downloads them from a
#' versioned data release (see "Getting the data"). A small extract,
#' [theta_sample], ships for examples.
#'
#' There are two tables. `download_bilatr_theta("directed")` holds the directed
#' model (one series per sender-receiver pair); `download_bilatr_theta("undirected")`
#' holds the undirected model (one series per unordered pair, events in
#' both directions pooled). Each row gives both the smoothed estimate --
#' the posterior of \eqn{\theta_{dt}} given all years, 1991-2025 -- and the
#' forward-filtered estimate, which uses only events up to and including
#' year \eqn{t} and so is the quantity to use when no look-ahead is allowed
#' (e.g. as a predictor in forecasting).
#'
#' **Higher values of \eqn{\theta} mean more cooperative relations.** The
#' reference category is `ModifiedRootCode` 5, "Consult: meet, discuss, or
#' visit", a cooperative class, and \eqn{\theta} is oriented so that higher
#' values make the reference category more likely relative to the other
#' action types (see the `reference_category` argument of
#' [prepare_stan_data()]). The sign is fixed by the model's orientation
#' constraint (\eqn{\alpha_1 \ge 0} for the reference class), the location
#' by the \eqn{\theta_0} prior (mean pinned at 0 across dyads) and the scale
#' by normalizing the discrimination parameters \eqn{\alpha} to
#' root-mean-square 1. Because the `stable_gamma` model includes
#' country-level category offsets, \eqn{\theta} is measured net of each
#' country's baseline mix of event types.
#'
#' Dyad-years before a dyad's first recorded event, and dyad-years in which
#' either country was not a member of the state system, are dropped: the
#' likelihood is only evaluated in years with events, so \eqn{\theta}
#' before the first event is the random-walk prior smoothed backwards
#' from that event and carries no information from the dyad's own data,
#' and a year outside membership has no events by construction (see
#' [aggregate_events()], "Actor periods"). The membership rule applies from
#' the next data update (`v2027.01`); `v2026.11` used a membership test
#' without years (see "Known gaps"), so it still has, for example,
#' Slovakia's dyad-years for 1991-1992.
#' After the first event every year through 2025 is kept, including
#' years without events (where \eqn{\theta} is interpolated or, at the
#' end of a series, extrapolated by the random walk), so each dyad's
#' series has no gaps and ends in 2025. All numeric columns are rounded
#' to 3 significant digits.
#'
#' @format Tibbles with one row per dyad-year, sorted by `dyad` then
#'   `year`. In data version `v2026.11`, the directed table has 355,821
#'   rows (16,750 directed dyads) and the undirected table 195,106 rows
#'   (9,139 undirected dyads). Both cover the same 189 countries (see
#'   "Known gaps"). Columns:
#' \describe{
#'   \item{dyad}{Character. Two three-letter country codes joined by `_`.
#'     The codes are GDELT's `Actor1CountryCode`/`Actor2CountryCode`
#'     values, which are CAMEO country codes; every code in these tables
#'     is also the country's ISO 3166-1 alpha-3 code. Directed:
#'     `"SENDER_RECEIVER"`, Actor1 then Actor2, so `"SWE_NOR"` is Sweden
#'     toward Norway. Undirected: the two codes in alphabetical (C-locale)
#'     order, so `"NOR_SWE"` only.}
#'   \item{year}{Integer, 1991-2025.}
#'   \item{mean}{Numeric. Posterior mean of the smoothed \eqn{\theta}.}
#'   \item{\code{5\%}, \code{95\%}}{Numeric. 5th and 95th percentiles of
#'     the smoothed \eqn{\theta} posterior (a 90% credible interval).}
#'   \item{filtered_mean}{Numeric. Posterior mean of the forward-filtered
#'     \eqn{\theta}: for each posterior draw, a West-Harrison/
#'     Fisher-scoring linear-Bayes filter is run forward from 1991
#'     conditional on that draw's hyperparameters, and the draw-wise
#'     filtered state means are averaged.}
#'   \item{\code{filtered_5\%}, \code{filtered_95\%}}{Numeric. 5th and
#'     95th percentiles, across posterior draws, of the draw-wise filtered
#'     state mean. These reflect hyperparameter uncertainty only, not the
#'     filter's own state uncertainty (`theta_filtered_sd`), which is why
#'     they are much narrower than the smoothed interval and stay constant
#'     through years without events. They are not a credible interval for
#'     the filtered \eqn{\theta}.}
#' }
#'
#' The interval columns have non-syntactic names, so they must be quoted
#' with backticks in R code.
#'
#' @details
#' **Events.** GDELT 1.0 events (the historical backfiles and the daily
#' updates) dated 1991-01-01 through 2025-12-31, accessed 2026-01-31. An
#' event is kept when both actors carry `GOV`, `MIL` or `SPY` in any of
#' their `Type1Code`/`Type2Code`/`Type3Code` fields and the two actors'
#' country codes differ. This is
#' `read_gdelt(actor_types = c("GOV", "MIL", "SPY"), actor_type_match =
#' "both", cross_border = TRUE)`, the defaults. All events are used, root
#' events or not, except that events with CAMEO codes 010-015, 017, 024x,
#' 025x, 040 and 041 (`BilatrClass` 0 in the original coding: neutral and
#' low-intensity statements) are excluded, which is why two
#' `ModifiedRootCode` classes do not occur. Events were further restricted
#' to those where both actors' country codes map to a Correlates of War
#' state system member (see "Known gaps"); [filter_gdelt_states()] applies
#' this restriction, by the year of each event, and the next data update
#' (`v2027.01`) is built with it (`read_gdelt(states_only = TRUE)`).
#'
#' **Model.** The `stable_gamma` model, at yearly resolution, with events
#' grouped by `ModifiedRootCode` (16 of its 18 classes occur in the data,
#' so \eqn{A = 16}), reference category 5 ("Consult: meet, discuss, or
#' visit"), `min_n_events = 1`, and 4 chains of 1,000 warmup and 1,000
#' sampling iterations each. The Stan program and the data preparation
#' used are the ones in this version of bilatr.
#'
#' @section Getting the data:
#' `download_bilatr_theta("directed")` and
#' `download_bilatr_theta("undirected")` download a table once, verify its
#' checksum, and cache it in `tools::R_user_dir("bilatr", "data")`; later
#' calls read the cache without using the internet.
#' The data are published as versioned releases of
#' <https://github.com/doktorandahl/bilatr-data> under the CC BY 4.0
#' license, archived on Zenodo. Each bilatr release is pinned to one data
#' version (`v2026.11` for this one), so the same code returns the same
#' data; `download_bilatr_theta(version = "latest")` asks for the newest
#' release, and [bilatr_theta_versions()] lists them. [bilatr_clear_data_cache()]
#' deletes the cache.
#'
#' Without bilatr, the same files can be read directly, e.g.
#' `readRDS(url(".../releases/download/v2026.11/theta_directed.rds"))`, or
#' `read.csv()` on `theta_directed.csv.gz`; see the data repository's
#' README.
#'
#' Please cite the data release by the DOI of the version you use, for
#' these data (`v2026.11`) \doi{10.5281/zenodo.23951819}; the DOI
#' \doi{10.5281/zenodo.23951818} refers to all versions of the data
#' release. Please also cite the GDELT references below. Cite the bilatr
#' package (`citation("bilatr")`) for the method and the software; the
#' entry for the data release is there too.
#'
#' @section Known gaps:
#' * **Correlates of War restriction.** Events were further restricted to
#'   those where both actors' country codes map to a Correlates of War
#'   state system member, so GDELT's regional and other non-state codes do
#'   not occur. The restriction was applied when the event table was
#'   built; [read_gdelt()] cannot reproduce it on its own.
#' * **Six missing states.** Serbia, Romania, Bosnia and Herzegovina,
#'   Montenegro, Slovenia and Timor-Leste (`SRB`, `ROU`, `BIH`, `MNE`,
#'   `SVN`, `TLS`) are absent from both tables, although they are
#'   Correlates of War system members. Their CAMEO country codes were not
#'   matched in the country-code crosswalk used to build the event table.
#'   So their dyads are missing, and their neighbors' other dyads include
#'   no events with them.
#'
#' This will be corrected in the next data update (`v2027.01`), expected in
#' January 2027.
#'
#' @source GDELT 1.0 event data, The GDELT Project,
#'   <https://gdeltproject.org/>, accessed 2026-01-31.
#' @references
#' Leetaru, K. and Schrodt, P. A. (2013). "GDELT: Global Data on Events,
#' Location and Tone, 1979--2012." Paper presented at the International
#' Studies Association annual convention, San Francisco.
#'
#' The GDELT Project. Global Database of Events, Language, and Tone
#' (GDELT). Accessed 31 January 2026, <https://gdeltproject.org/>.
#' @seealso [download_bilatr_theta()] to download the tables; [theta_sample]
#'   for the extract that ships with the package; [cameo_lookup] for the
#'   `ModifiedRootCode` classes; [diagnose_and_extract_bilatr()] for how
#'   the theta summaries are extracted from a fit.
#' @examples
#' # Norway-Sweden, both directions, from the shipped extract
#' subset(theta_sample, type == "directed" & dyad %in% c("NOR_SWE", "SWE_NOR"))
#'
#' # the same pair in the undirected model
#' subset(theta_sample, type == "undirected" & dyad == "NOR_SWE")
#' @examplesIf interactive()
#' # the full tables (downloaded once, then cached)
#' theta_directed <- download_bilatr_theta("directed")
#' theta_undirected <- download_bilatr_theta("undirected")
#' @keywords datasets
#' @name theta_trajectories
NULL

#' A sample of the estimated theta trajectories
#'
#' An extract of the trajectories described in [theta_trajectories], for
#' examples and vignettes: ten country pairs, both directions of each from
#' the directed table and the pair from the undirected table, all years.
#' The full tables are downloaded with [download_bilatr_theta()].
#'
#' The pairs are United States-China, United States-Russia,
#' Russia-Ukraine, China-Japan, South Korea-North Korea, India-Pakistan,
#' Israel-Iran, United Kingdom-France, Germany-France and Norway-Sweden.
#'
#' @format A tibble with 1,038 rows (690 directed and 348 undirected
#'   dyad-years; 20 directed and 10 undirected dyads), sorted by `type`,
#'   `dyad` and `year`. The first column, `type`, is `"directed"` or
#'   `"undirected"`; the other columns are those of [theta_trajectories],
#'   with the same values (data version `v2026.11`).
#' @source Data version `v2026.11` of the bilatr theta trajectories; see
#'   [theta_trajectories].
#' @seealso [theta_trajectories], [download_bilatr_theta()].
#' @examples
#' # United States toward China and China toward the United States
#' us_cn <- subset(theta_sample, type == "directed" & dyad %in% c("USA_CHN", "CHN_USA"))
#' tail(us_cn[c("dyad", "year", "mean", "filtered_mean")])
#'
#' # the dyads in the sample
#' table(theta_sample$type)
#' unique(theta_sample$dyad[theta_sample$type == "undirected"])
"theta_sample"
