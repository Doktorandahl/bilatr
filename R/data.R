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

#' Estimated bilateral-relations trajectories, 1991-2025
#'
#' Yearly latent trajectories (\eqn{\theta}) of the level of bilateral
#' relations for every country pair with at least one recorded event,
#' estimated with bilatr from GDELT 1.0 event data. The data are not
#' shipped with the package: [bilatr_theta()] downloads them from a
#' versioned data release (see "Getting the data"). A small extract,
#' [theta_sample], ships for examples.
#'
#' There are two tables. `bilatr_theta("directed")` holds the directed
#' model (one series per sender-receiver pair); `bilatr_theta("undirected")`
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
#' Dyad-years before a dyad's first recorded event are dropped: the
#' likelihood is only evaluated in years with events, so \eqn{\theta}
#' before the first event is the random-walk prior smoothed backwards
#' from that event and carries no information from the dyad's own data.
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
#' events or not, with no filter on event codes. Events were further
#' restricted to those where both actors' country codes map to a
#' Correlates of War state system member (see "Known gaps").
#'
#' **Model.** The `stable_gamma` model, at yearly resolution, with events
#' grouped by `ModifiedRootCode` (16 of its 18 classes occur in the data,
#' so \eqn{A = 16}), reference category 5 ("Consult: meet, discuss, or
#' visit"), `min_n_events = 1`, and 4 chains of 1,000 warmup and 1,000
#' sampling iterations each. The Stan program and the data preparation
#' used are the ones in this version of bilatr.
#'
#' @section Getting the data:
#' `bilatr_theta("directed")` and `bilatr_theta("undirected")` download a
#' table once, verify its checksum, and cache it in
#' `tools::R_user_dir("bilatr", "data")`; later calls read the cache.
#' The data are published as versioned releases of
#' <https://github.com/doktorandahl/bilatr-data> under the CC BY 4.0
#' license, archived on Zenodo. Each bilatr release is pinned to one data
#' version (`v2026.11` for this one), so the same code returns the same
#' data; `bilatr_theta(version = "latest")` asks for the newest release,
#' and [bilatr_theta_versions()] lists them. [bilatr_clear_data_cache()]
#' deletes the cache.
#'
#' Without bilatr, the same files can be read directly, e.g.
#' `readRDS(url(".../releases/download/v2026.11/theta_directed.rds"))`, or
#' `read.csv()` on `theta_directed.csv.gz`; see the data repository's
#' README.
#'
#' Please cite the data release by its Zenodo DOI (Zenodo DOI: to be
#' added once the record is published), the bilatr package
#' (`citation("bilatr")`), and the GDELT references below.
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
#' @seealso [bilatr_theta()] to download the tables; [theta_sample] for
#'   the extract that ships with the package; [cameo_lookup] for the
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
#' theta_directed <- bilatr_theta("directed")
#' theta_undirected <- bilatr_theta("undirected")
#' @keywords datasets
#' @name theta_trajectories
NULL

#' A sample of the estimated theta trajectories
#'
#' An extract of the trajectories described in [theta_trajectories], for
#' examples and vignettes: ten country pairs, both directions of each from
#' the directed table and the pair from the undirected table, all years.
#' The full tables are downloaded with [bilatr_theta()].
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
#' @seealso [theta_trajectories], [bilatr_theta()].
#' @examples
#' # United States toward China and China toward the United States
#' us_cn <- subset(theta_sample, type == "directed" & dyad %in% c("USA_CHN", "CHN_USA"))
#' tail(us_cn[c("dyad", "year", "mean", "filtered_mean")])
#'
#' # the dyads in the sample
#' table(theta_sample$type)
#' unique(theta_sample$dyad[theta_sample$type == "undirected"])
"theta_sample"
