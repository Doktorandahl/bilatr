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
#'     `QuadClass`; see [quadclass_name()].}
#'   \item{QuadClassEventCodes}{Character. Comma-separated CAMEO root
#'     codes underlying each `QuadClass` class; see
#'     [quadclass_eventcodes()].}
#'   \item{PentaClass}{Integer, 0-4. As QuadClass, except that roots
#'     01-02 (make public statement, appeal) form class 0, root 14
#'     (protest) moves to class 4 (material conflict), and root 16
#'     (reduce relations) moves to class 3 (verbal conflict): 0 = make
#'     statement (roots 01-02), 1 = verbal cooperation (03-05),
#'     2 = material cooperation (06-09), 3 = verbal conflict (10-13, 16),
#'     4 = material conflict (14, 15, 17-20).}
#'   \item{PentaClassName}{Character. Human-readable label for
#'     `PentaClass`; see [pentaclass_name()].}
#'   \item{PentaClassEventCodes}{Character. Comma-separated CAMEO root
#'     codes underlying each `PentaClass` class; see
#'     [pentaclass_eventcodes()].}
#'   \item{PentaClass_modified}{Integer, 0-4. As PentaClass, except that
#'     every code with `QuadClass` 1 and `GoldsteinScore` at most 1 is
#'     class 0. Class 0 is therefore roots 01-02 plus the codes of roots
#'     03-05 scoring at most 1 (in this table `04`, `040` and `041`).}
#'   \item{ModifiedRootCode}{Integer, 1-18. A regrouping of the CAMEO
#'     root codes that moves `016` (deny responsibility) into Reject
#'     (class 13) and `018`/`019` (make empathetic comment, express
#'     accord) into diplomatic cooperation (class 7); splits root `04`
#'     (consult) into unspecified (`040`), meet, discuss, or visit
#'     (`041`-`044`) and negotiate or mediate (`045`-`046`); merges
#'     Investigate (`09`) and Demand (`10`); merges Threaten (`13`) and
#'     Exhibit force posture (`15`); and merges `18`-`20` into "Assault,
#'     fight, or mass violence". See [assign_modified_root_code()].
#'     Blank (`NA`) on the top-level (two-digit) row for a root whose own
#'     event codes don't all belong to the same `ModifiedRootCode` class
#'     (roots `01` and `04`); populated on that root's individual
#'     sub-codes as normal.}
#'   \item{ModifiedRootCodeName}{Character. Human-readable label for
#'     `ModifiedRootCode`; see [modified_root_code_name()]. `NA`
#'     wherever `ModifiedRootCode` is blank.}
#'   \item{ModifiedRootCodeEventCodes}{Character. Comma-separated CAMEO
#'     root/event codes underlying each `ModifiedRootCode` class; see
#'     [modified_root_code_eventcodes()]. `NA` wherever `ModifiedRootCode`
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
"cameo_lookup"

#' Estimated bilateral-relations trajectories, 1991-2025
#'
#' Yearly latent trajectories (\eqn{\theta}) of the level of bilateral
#' relations for every country pair with at least one recorded event,
#' estimated with bilatr from GDELT 1.0 event data. `theta_directed` holds
#' the directed model (one series per sender-receiver pair);
#' `theta_undirected` holds the undirected model (one series per unordered
#' pair, events in both directions pooled). Each row gives both the
#' smoothed estimate -- the posterior of \eqn{\theta_{dt}} given all
#' years, 1991-2025 -- and the forward-filtered estimate, which uses only
#' events up to and including year \eqn{t} and so is the quantity to use
#' when no look-ahead is allowed (e.g. as a predictor in forecasting).
#'
#' **Higher values of \eqn{\theta} mean more cooperative relations.** The
#' reference category is `ModifiedRootCode` 5, "Consult: meet, discuss, or
#' visit", a cooperative class, and \eqn{\theta} is oriented so that higher
#' values make the reference category more likely relative to the other
#' action types (see the `reference_category` argument of
#' [prepare_stan_data()]). The sign is fixed by the model's orientation
#' fold (\eqn{\alpha_1 \ge 0}; see the `stable_gamma` Stan program's
#' header), the location by the \eqn{\theta_0} prior (mean exactly 0
#' across dyads) and the scale by normalising the discrimination
#' parameters \eqn{\alpha} to root-mean-square 1. Because the
#' `stable_gamma` model includes country-level category offsets,
#' \eqn{\theta} is measured net of each country's baseline mix of event
#' types.
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
#'   `year`. `theta_directed` has 355,821 rows (16,750 directed dyads);
#'   `theta_undirected` has 195,106 rows (9,139 undirected dyads). Both
#'   cover the same 189 countries. Columns:
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
#' events or not, with no filter on event codes. In addition, both actors'
#' country codes map to a Correlates of War state, so GDELT's regional and
#' other non-state codes (such as `"AFR"`, `"EUR"`, `"PSE"` or `"HKG"`) do
#' not occur.
#'
#' **Model.** The `stable_gamma` Stan program, at yearly resolution, with
#' events grouped by `ModifiedRootCode` (16 of its 18 classes occur in the
#' data, so \eqn{A = 16}), reference category 5 ("Consult: meet, discuss,
#' or visit"), `min_n_events = 1`, and 4 chains of 1,000 warmup and 1,000
#' sampling iterations each. Estimated with an earlier development version
#' of bilatr; the Stan program and the data preparation are unchanged in
#' this version.
#'
#' Please cite the bilatr methods article (forthcoming); see
#' `citation("bilatr")`.
#'
#' @source GDELT 1.0 event data, The GDELT Project,
#'   <https://www.gdeltproject.org/>, accessed 2026-01-31.
#' @references
#' Leetaru, K. and Schrodt, P. A. (2013). "GDELT: Global Data on Events,
#' Location and Tone, 1979–2012." Paper presented at the International
#' Studies Association annual convention, San Francisco.
#'
#' The GDELT Project. Global Database of Events, Language, and Tone
#' (GDELT). Accessed 31 January 2026, <https://gdeltproject.org/>.
#' @seealso [cameo_lookup] for the `ModifiedRootCode` classes;
#'   [diagnose_and_extract_bilatr()] for how the theta summaries are
#'   extracted from a fit.
#' @examples
#' # Norway-Sweden, both directions
#' subset(theta_directed, dyad %in% c("NOR_SWE", "SWE_NOR"))
#'
#' # same pair in the undirected model
#' subset(theta_undirected, dyad == "NOR_SWE")
#' @rdname theta_trajectories
"theta_directed"

#' @rdname theta_trajectories
#' @format NULL
"theta_undirected"
