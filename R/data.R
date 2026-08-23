#' CAMEO event code lookup table
#'
#' A reference table mapping every CAMEO event code to its Goldstein
#' conflict-cooperation score, a human-readable label, and the
#' QuadClass/PentaClass recoding used throughout the package. Built by
#' `data-raw/build_cameo_lookup.R`; see [recode_cameo()] for how it is
#' applied to event data.
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
#'   \item{PentaClass}{Integer, 0-4. As QuadClass, but with verbal
#'     cooperation (0) split out from material cooperation (1), and
#'     protest/reduce-relations events reassigned to their own classes.}
#'   \item{PentaClass_modified}{Integer, 0-4. As PentaClass, but low-
#'     intensity verbal cooperation (Goldstein score <= 1) is folded into
#'     class 0, useful as a near-neutral reference category.}
#' }
#' @source GDELT/CAMEO event taxonomy; Goldstein scale from Goldstein
#'   (1992), "A Conflict-Cooperation Scale for WEIS Events Data".
"cameo_lookup"
