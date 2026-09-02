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

#' Regroup a CAMEO event code into the bilatr "EventRootCode2" scheme
#'
#' `EventRootCode2` is a coarser regrouping of the 20 CAMEO root codes
#' used by this project, motivated by how the underlying event types
#' behave in dyadic conflict data:
#'
#' \itemize{
#'   \item Root 04 ("Consult") is *split* by the kind of consultation:
#'     bare `04`/`040` stay `"040"`; meeting/visiting/phoning
#'     (`041`-`044`) become `"044"`; mediating/negotiating (`045`-`046`)
#'     become `"046"`.
#'   \item Root 09 ("Investigate") is folded into 10 ("Demand") -> `"10"`.
#'   \item Roots 14 ("Protest") and 15 ("Exhibit Force Posture") are
#'     folded into 13 ("Threaten") -> `"13"`.
#'   \item Roots 18 ("Assault") and 20 ("Use Unconventional Mass
#'     Violence") are folded into 19 ("Fight") -> `"19"`.
#'   \item Every other root code maps to itself.
#' }
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Character vector of `EventRootCode2` values.
#' @keywords internal
assign_eventrootcode2 <- function(code) {
  code <- stringr::str_pad(as.character(code), width = 2, pad = "0")
  root <- get_root(code)
  sub3 <- stringr::str_sub(code, 1, 3)
  dplyr::case_when(
    root == "04" & sub3 %in% c("041", "042", "043", "044") ~ "044",
    root == "04" & sub3 %in% c("045", "046") ~ "046",
    root == "04" ~ "040",
    root == "09" ~ "10",
    root %in% c("14", "15") ~ "13",
    root %in% c("18", "20") ~ "19",
    TRUE ~ root
  )
}

#' Human-readable label for a bilatr event class
#'
#' @param class Integer vector of `BilatrClass` values (0-10).
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 0-10.
#' @keywords internal
bilatr_class_name <- function(class) {
  names <- c(
    "0"  = "Neutral / low-intensity statement",
    "1"  = "Express intent to cooperate",
    "2"  = "Consult: meet, discuss, or visit",
    "3"  = "Consult: negotiate or mediate",
    "4"  = "Engage in diplomatic cooperation",
    "5"  = "Engage in material cooperation or provide aid",
    "6"  = "Yield",
    "7"  = "Investigate, demand, reject, or reduce relations",
    "8"  = "Disapprove",
    "9"  = "Threaten or coerce",
    "10" = "Assault, fight, or mass violence"
  )
  unname(names[as.character(class)])
}

#' Human-readable label for a coarsened bilatr event class
#'
#' @param class Integer vector of `BilatrClass2` values (0-8).
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 0-8.
#' @keywords internal
bilatr_class2_name <- function(class) {
  names <- c(
    "0" = "Neutral / low-intensity statement",
    "1" = "Express intent to cooperate",
    "2" = "Consult: meet, discuss, or visit",
    "3" = "Consult: negotiate or mediate",
    "4" = "Engage in diplomatic cooperation",
    "5" = "Engage in material cooperation or provide aid",
    "6" = "Yield",
    "7" = "Disapprove, demand, reject, or reduce relations",
    "8" = "Threaten, coerce, or use force"
  )
  unname(names[as.character(class)])
}

#' Coarsen `BilatrClass` into the 9-level `BilatrClass2` scheme
#'
#' `BilatrClass2` merges the two adjacent pairs of hostile `BilatrClass`
#' levels (see [assign_bilatr_class()]) that behave similarly in dyadic
#' conflict data, giving a 9-level (0-8) scheme:
#'
#' \itemize{
#'   \item `BilatrClass` 0-6 are unchanged (`BilatrClass2` 0-6).
#'   \item `BilatrClass` 7 ("Investigate, demand, reject, or reduce
#'     relations") and 8 ("Disapprove") merge into `BilatrClass2` 7
#'     ("Disapprove, demand, reject, or reduce relations").
#'   \item `BilatrClass` 9 ("Threaten or coerce") and 10 ("Assault,
#'     fight, or mass violence") merge into `BilatrClass2` 8 ("Threaten,
#'     coerce, or use force").
#' }
#'
#' See [bilatr_class2_name()] for the level labels.
#'
#' @param bilatr_class Integer vector of `BilatrClass` values (0-10), as
#'   returned by [assign_bilatr_class()].
#' @return Integer vector of `BilatrClass2` values (0-8), `NA` where
#'   `bilatr_class` is `NA` or outside 0-10.
#' @keywords internal
assign_bilatr_class2 <- function(bilatr_class) {
  bilatr_class <- as.integer(bilatr_class)
  dplyr::case_when(
    bilatr_class >= 0L & bilatr_class <= 6L ~ bilatr_class,
    bilatr_class %in% c(7L, 8L) ~ 7L,
    bilatr_class %in% c(9L, 10L) ~ 8L,
    TRUE ~ NA_integer_
  )
}

#' Assign the bilatr event class from a CAMEO event code
#'
#' `BilatrClass` is an 11-level (0-10) collapse of the CAMEO taxonomy
#' used as the default action-class scheme for the bilatr model. It is
#' mostly a function of [assign_eventrootcode2()] (see that function for
#' the root regrouping and [bilatr_class_name()] for the level labels),
#' with a handful of per-code refinements where a specific event type
#' fits a different class than its regrouped root:
#'
#' \itemize{
#'   \item `016` "Deny responsibility" -> 8 (Disapprove).
#'   \item `018` "Make empathetic comment" -> 4 (diplomatic cooperation).
#'   \item `019` "Express accord" -> 1 (express intent to cooperate).
#'   \item Appeals for political reform / to yield (`024*`, `025*`) ->
#'     0 (neutral statement).
#'   \item Appeals to others to meet / settle / mediate (`026`, `027`,
#'     `028`) -> 3 (consult: negotiate or mediate).
#'   \item `041` "Discuss by telephone" -> 0 (neutral statement).
#' }
#'
#' @param code Character or numeric vector of CAMEO event codes.
#' @param eventrootcode2 Character vector of `EventRootCode2` values, as
#'   returned by [assign_eventrootcode2()]. Defaults to computing it from
#'   `code`.
#' @return Integer vector of `BilatrClass` values (0-10), `NA` for
#'   unrecognized codes.
#' @keywords internal
assign_bilatr_class <- function(code, eventrootcode2 = assign_eventrootcode2(code)) {
  code <- stringr::str_pad(as.character(code), width = 2, pad = "0")
  dplyr::case_when(
    code == "016" ~ 8L,
    code == "018" ~ 4L,
    code == "019" ~ 1L,
    stringr::str_starts(code, "024") ~ 0L,
    stringr::str_starts(code, "025") ~ 0L,
    code %in% c("026", "027", "028") ~ 3L,
    code == "041" ~ 0L,
    eventrootcode2 == "01" ~ 0L,
    eventrootcode2 == "02" ~ 1L,
    eventrootcode2 == "03" ~ 1L,
    eventrootcode2 == "040" ~ 0L,
    eventrootcode2 == "044" ~ 2L,
    eventrootcode2 == "046" ~ 3L,
    eventrootcode2 == "05" ~ 4L,
    eventrootcode2 == "06" ~ 5L,
    eventrootcode2 == "07" ~ 5L,
    eventrootcode2 == "08" ~ 6L,
    eventrootcode2 == "10" ~ 7L,
    eventrootcode2 == "11" ~ 8L,
    eventrootcode2 == "12" ~ 7L,
    eventrootcode2 == "13" ~ 9L,
    eventrootcode2 == "16" ~ 7L,
    eventrootcode2 == "17" ~ 9L,
    eventrootcode2 == "19" ~ 10L,
    TRUE ~ NA_integer_
  )
}

#' Recode a data frame of CAMEO event codes to quad/penta/bilatr classes
#'
#' Joins `data` against the package's built-in [cameo_lookup] table to
#' attach `CAMEOLabel`, `GoldsteinScore`, `QuadClass`, `PentaClass`,
#' `PentaClass_modified`, `EventRootCode2`, `BilatrClass`,
#' `BilatrClassName`, `BilatrClass2`, and `BilatrClass2Name` columns.
#' `PentaClass_modified` folds low-intensity verbal cooperation (Goldstein
#' score <= 1) into its own class, which can be useful as a near-neutral
#' reference category. `EventRootCode2` (see [assign_eventrootcode2()]) is
#' a coarser regrouping of the CAMEO root codes; `BilatrClass` /
#' `BilatrClassName` (see [assign_bilatr_class()]) is the 11-level action
#' scheme used as the model's default; `BilatrClass2` / `BilatrClass2Name`
#' (see [assign_bilatr_class2()]) is a 9-level coarsening that merges the
#' two adjacent pairs of hostile levels.
#'
#' Any of these columns that already exist in `data` are left as they are
#' (not overwritten, no `.x`/`.y` suffixing), with a warning naming them.
#'
#' @param data A data frame containing a CAMEO event code column.
#' @param code_col Name of the column in `data` holding CAMEO event codes
#'   (as a string). Defaults to `"EventCode"`, matching the raw GDELT
#'   column name.
#' @return `data` with the recode columns above attached (minus any that
#'   were already present).
#' @examples
#' \dontrun{
#' events <- data.frame(EventCode = c("01", "190", "0862"))
#' recode_cameo(events)
#' }
#' @export
recode_cameo <- function(data, code_col = "EventCode") {
  lookup <- bilatr::cameo_lookup

  already_present <- intersect(
    setdiff(names(lookup), "CAMEOEVENTCODE"),
    names(data)
  )
  if (length(already_present) > 0) {
    warning(
      "recode_cameo(): ", length(already_present),
      " recode column(s) already present in `data`; leaving them untouched: ",
      paste(already_present, collapse = ", "), ".",
      call. = FALSE
    )
    lookup <- dplyr::select(lookup, -dplyr::all_of(already_present))
  }

  dplyr::left_join(
    data,
    lookup,
    by = rlang::set_names("CAMEOEVENTCODE", code_col)
  )
}
