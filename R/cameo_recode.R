#' Extract the two-digit CAMEO root code
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Character vector of two-digit root codes.
#' @keywords internal
get_root <- function(code) {
  stringr::str_sub(stringr::str_pad(as.character(code), width = 2, pad = "0"), 1, 2)
}

#' Assign QuadClass from a CAMEO root code
#'
#' @param root Character vector of two-digit CAMEO root codes, as returned
#'   by [get_root()].
#' @return Integer vector with values 1 (verbal cooperation), 2 (material
#'   cooperation), 3 (verbal conflict), 4 (material conflict), or `NA`.
#' @keywords internal
assign_quad <- function(root) {
  root_num <- as.integer(root)
  dplyr::case_when(
    root_num >= 1 & root_num <= 5 ~ 1L,
    root_num >= 6 & root_num <= 9 ~ 2L,
    root_num >= 10 & root_num <= 14 ~ 3L,
    root_num >= 15 & root_num <= 20 ~ 4L,
    TRUE ~ NA_integer_
  )
}

#' Assign PentaClass from a CAMEO root code and QuadClass
#'
#' Refines QuadClass into a five-level scheme by splitting out verbal
#' cooperation (root 01-02) as its own class and carving protest (root 14)
#' and reduce-relations (root 16) out of the conflict classes.
#'
#' @param root Character vector of two-digit CAMEO root codes.
#' @param quad Integer vector of QuadClass values, as returned by
#'   [assign_quad()].
#' @return Integer vector with values 0-4.
#' @keywords internal
assign_penta <- function(root, quad) {
  root_num <- as.integer(root)
  dplyr::case_when(
    root_num %in% c(1, 2) ~ 0L, # verbal cooperation
    root_num == 14 ~ 4L, # protest (override quad = 3)
    root_num == 16 ~ 3L, # reduce relations (override quad = 4)
    TRUE ~ quad
  )
}

#' Recode a data frame of CAMEO event codes to quad/penta classes
#'
#' Joins `data` against the package's built-in [cameo_lookup] table to
#' attach `QuadClass`, `PentaClass`, and `PentaClass_modified` columns.
#' `PentaClass_modified` folds low-intensity verbal cooperation (Goldstein
#' score <= 1) into its own class, which can be useful as a near-neutral
#' reference category.
#'
#' @param data A data frame containing a CAMEO event code column.
#' @param code_col Name of the column in `data` holding CAMEO event codes
#'   (as a string). Defaults to `"EventCode"`, matching the raw GDELT
#'   column name.
#' @return `data` with `QuadClass`, `PentaClass`, `PentaClass_modified`,
#'   `GoldsteinScore`, and `CAMEOLabel` columns attached.
#' @examples
#' \dontrun{
#' events <- data.frame(EventCode = c("01", "190", "0862"))
#' recode_cameo(events)
#' }
#' @export
recode_cameo <- function(data, code_col = "EventCode") {
  data %>%
    dplyr::left_join(
      bilatr::cameo_lookup,
      by = rlang::set_names("CAMEOEVENTCODE", code_col)
    )
}
