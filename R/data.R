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
#'   \item{QuadClass}{Integer, 1-4. 1 = verbal cooperation, 2 = material
#'     cooperation, 3 = verbal conflict, 4 = material conflict.}
#'   \item{QuadClassName}{Character. Human-readable label for
#'     `QuadClass`; see [quadclass_name()].}
#'   \item{QuadClassEventCodes}{Character. Comma-separated CAMEO root
#'     codes underlying each `QuadClass` class; see
#'     [quadclass_eventcodes()].}
#'   \item{PentaClass}{Integer, 0-4. As QuadClass, but with verbal
#'     cooperation (0) split out from material cooperation (1), and
#'     protest/reduce-relations events reassigned to their own classes.}
#'   \item{PentaClassName}{Character. Human-readable label for
#'     `PentaClass`; see [pentaclass_name()].}
#'   \item{PentaClassEventCodes}{Character. Comma-separated CAMEO root
#'     codes underlying each `PentaClass` class; see
#'     [pentaclass_eventcodes()].}
#'   \item{PentaClass_modified}{Integer, 0-4. As PentaClass, but low-
#'     intensity verbal cooperation (Goldstein score <= 1) is folded into
#'     class 0, useful as a near-neutral reference category.}
#'   \item{ModifiedRootCode}{Integer, 1-18. A regrouping of the CAMEO
#'     root codes that relocates `016` into the Reject class and
#'     `018`/`019` into the diplomatic-cooperation class, splits root
#'     `04` ("Consult") by kind of consultation, splits root `13`
#'     ("Threaten") off from force posture only where noted, and merges
#'     Investigate (root `09`) and Demand (root `10`) into a single
#'     "Investigate or demand" class; see
#'     [assign_modified_root_code()]. Blank (`NA`) on the top-level
#'     (two-digit) row for a root whose own event codes don't all belong
#'     to the same `ModifiedRootCode` class (currently roots `01` and
#'     `04`); populated on that root's individual sub-codes as normal.
#'     This scheme was previously named `ERC16NZ`.}
#'   \item{ModifiedRootCodeName}{Character. Human-readable label for
#'     `ModifiedRootCode`; see [modified_root_code_name()]. `NA`
#'     wherever `ModifiedRootCode` is blank.}
#'   \item{ModifiedRootCodeEventCodes}{Character. Comma-separated CAMEO
#'     root/event codes underlying each `ModifiedRootCode` class; see
#'     [modified_root_code_eventcodes()]. `NA` wherever `ModifiedRootCode`
#'     is blank.}
#' }
#'
#' A handful of other classification schemes previously shipped in this
#' table (`EventRootCode2`, `EventRootCode3`, `EventRootCode4`,
#' `BilatrClass`, `BilatrClass2`, each with their `*Name`/`*RootCodes`
#' columns) have been retired from it, and their `assign_*()`/`*_name()`
#' functions removed entirely (0.10.1; they are no longer attached by
#' [recode_cameo()], and reading pre-0.10.1 output that used them needs
#' an older bilatr install). [event_class_labels()] is the current,
#' data-driven way to get a class's human-readable label from a fitted
#' model's `stan_data`, for any scheme that still has a `*Name` column
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
#' estimated with bilatr from GDELT event data. `theta_directed` holds the
#' directed model (one series per sender-receiver pair); `theta_undirected`
#' holds the undirected model (one series per unordered pair, events in
#' both directions pooled). Each row gives both the smoothed estimate --
#' the posterior of \eqn{\theta_{dt}} given all years, 1991-2025 -- and
#' the forward-filtered estimate, which uses only events up to and
#' including year \eqn{t} and so is the quantity to use when no
#' look-ahead is allowed (e.g. as a predictor in forecasting).
#'
#' Higher values of \eqn{\theta} indicate ??? (more cooperative / more
#' conflictual) relations. The sign is fixed by the model's orientation
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
#'   \item{dyad}{Character. Three-letter country codes joined by `_`,
#'     e.g. `"SWE_NOR"`. Directed: `"SENDER_RECEIVER"` ???. Undirected:
#'     the two codes in alphabetical order, so `"NOR_SWE"` only.}
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
#'     through years without events. ???}
#' }
#'
#' The interval columns have non-syntactic names, so they must be quoted
#' with backticks in R code.
#'
#' @details Both models were fitted with bilatr 0.10.2 (commit `4a4ffd5`)
#'   via `runscripts/submit_bilatr_16cat.R`: the `stable_gamma` Stan
#'   program, yearly resolution, events grouped by `ModifiedRootCode`
#'   (16 of its 18 classes occur in the data, so \eqn{A = 16}) with
#'   reference category 5 ("Consult: meet, discuss, or visit"),
#'   `min_n_events = 1`, and 4 chains of 1,000 warmup and 1,000 sampling
#'   iterations each. The theta summaries were then processed with
#'   `runscripts/process_theta_16cat.R` and saved here by
#'   `data-raw/theta_trajectories.R`.
#'
#' @source Estimated from GDELT ??? (version, download date, actor and
#'   event filters) events between state actors of ??? countries, 1991-2025.
#'   Please cite ??? when using these data.
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
