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
#' columns) have been retired from it. Their `assign_*()`/`*_name()`
#' functions are still available internally in `R/cameo_recode.R` for
#' backward compatibility (see e.g. [assign_eventrootcode2()],
#' [assign_bilatr_class()]), but are no longer attached by
#' [recode_cameo()].
#' @source Built by `data-raw/build_cameo_lookup.R` in the source
#'   repository. CAMEO event codes and labels from Schrodt (2012), *CAMEO
#'   Conflict and Mediation Event Observations Event and Actor Codebook*,
#'   version 1.1b3. Goldstein scale from Goldstein (1992), "A
#'   Conflict-Cooperation Scale for WEIS Events Data".
"cameo_lookup"
