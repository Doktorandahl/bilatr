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

#' Human-readable label for a `QuadClass` value
#'
#' @param class Integer vector of `QuadClass` values, as returned by
#'   [assign_quad()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 1-4.
#' @keywords internal
quadclass_name <- function(class) {
  names <- c(
    "1" = "Verbal cooperation",
    "2" = "Material cooperation",
    "3" = "Verbal conflict",
    "4" = "Material conflict"
  )
  unname(names[as.character(class)])
}

#' Underlying CAMEO root codes for a `QuadClass` value
#'
#' @param class Integer vector of `QuadClass` values, as returned by
#'   [assign_quad()].
#' @return Character vector of comma-separated CAMEO root codes, `NA` for
#'   values outside 1-4.
#' @keywords internal
quadclass_eventcodes <- function(class) {
  codes <- c(
    "1" = "01, 02, 03, 04, 05",
    "2" = "06, 07, 08, 09",
    "3" = "10, 11, 12, 13, 14",
    "4" = "15, 16, 17, 18, 19, 20"
  )
  unname(codes[as.character(class)])
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

#' Human-readable label for a `PentaClass` value
#'
#' @param class Integer vector of `PentaClass` values, as returned by
#'   [assign_penta()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 0-4.
#' @keywords internal
pentaclass_name <- function(class) {
  names <- c(
    "0" = "Make statement",
    "1" = "Verbal cooperation",
    "2" = "Material cooperation",
    "3" = "Verbal conflict",
    "4" = "Material conflict"
  )
  unname(names[as.character(class)])
}

#' Underlying CAMEO root codes for a `PentaClass` value
#'
#' @param class Integer vector of `PentaClass` values, as returned by
#'   [assign_penta()].
#' @return Character vector of comma-separated CAMEO root codes, `NA` for
#'   values outside 0-4.
#' @keywords internal
pentaclass_eventcodes <- function(class) {
  codes <- c(
    "0" = "01, 02",
    "1" = "03, 04, 05",
    "2" = "06, 07, 08, 09",
    "3" = "10, 11, 12, 13, 16",
    "4" = "14, 15, 17, 18, 19, 20"
  )
  unname(codes[as.character(class)])
}

#' Regroup a CAMEO event code into the bilatr "EventRootCode2" scheme
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup]. Superseded by `ModifiedRootCode`; see
#' [assign_modified_root_code()].
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

#' Human-readable label for an `EventRootCode2` value
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
#'
#' Preliminary labels, one per [assign_eventrootcode2()] output value.
#' Where an `EventRootCode2` value corresponds exactly to a `BilatrClass`
#' level (the two root-04 splits `"044"`/`"046"`), the label text matches
#' [bilatr_class_name()]'s for that level, since they denote the same
#' category; every other label is new here, since `EventRootCode2` is
#' finer-grained than `BilatrClass` everywhere else (e.g. `"06"`/`"07"`
#' get distinct labels despite both folding into `BilatrClass`'s single
#' "material cooperation or provide aid").
#'
#' @param class Character vector of `EventRootCode2` values, as returned
#'   by [assign_eventrootcode2()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values [assign_eventrootcode2()] does not produce.
#' @keywords internal
eventrootcode2_name <- function(class) {
  names <- c(
    "01"  = "Make a public statement",
    "02"  = "Appeal for action",
    "03"  = "Express intent to cooperate",
    "040" = "Consult, unspecified",
    "044" = "Consult: meet, discuss, or visit",
    "046" = "Consult: negotiate or mediate",
    "05"  = "Engage in diplomatic cooperation",
    "06"  = "Engage in material cooperation",
    "07"  = "Provide aid",
    "08"  = "Yield",
    "10"  = "Investigate or demand",
    "11"  = "Disapprove",
    "12"  = "Reject",
    "13"  = "Threaten, protest, or exhibit force posture",
    "16"  = "Reduce relations",
    "17"  = "Coerce",
    "19"  = "Assault, fight, or mass violence"
  )
  unname(names[as.character(class)])
}

#' Regroup a CAMEO event code into the bilatr "EventRootCode3" scheme
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup]. Superseded by `ModifiedRootCode`; see
#' [assign_modified_root_code()].
#'
#' `EventRootCode3` is a further refinement of [assign_eventrootcode2()],
#' motivated by wanting a few of its coarser groupings split back out
#' while relocating individual codes that behave more like a different
#' category than their root suggests:
#'
#' \itemize{
#'   \item `016` ("Deny responsibility") moves out of root 01 into the
#'     Reject class (with root 12).
#'   \item `018` ("Make empathetic comment") and `019` ("Express accord")
#'     move out of root 01 into the diplomatic-cooperation class (with
#'     root 05).
#'   \item `EventRootCode2`'s `"10"` (roots 09 + 10, "Investigate or
#'     demand") is split back into Investigate (root 09) and Demand
#'     (root 10).
#'   \item `EventRootCode2`'s `"13"` (roots 13 + 14 + 15, "Threaten,
#'     protest, or exhibit force posture") is split so that root 14
#'     ("Protest") becomes its own class, while roots 13 and 15
#'     ("Threaten" and "Exhibit force posture") stay merged.
#'   \item Every other root maps as in [assign_eventrootcode2()].
#' }
#'
#' Classes are numbered 1-19 in ascending root-code order; see
#' [eventrootcode3_name()] for labels and
#' [eventrootcode3_rootcodes()] for the underlying CAMEO root/event
#' codes each class draws from.
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Integer vector of `EventRootCode3` values (1-19).
#' @keywords internal
assign_eventrootcode3 <- function(code) {
  code <- stringr::str_pad(as.character(code), width = 2, pad = "0")
  root <- get_root(code)
  sub3 <- stringr::str_sub(code, 1, 3)
  dplyr::case_when(
    code == "016" ~ 14L,
    code == "018" ~ 7L,
    code == "019" ~ 7L,
    root == "01" ~ 1L,
    root == "02" ~ 2L,
    root == "03" ~ 3L,
    root == "04" & sub3 %in% c("041", "042", "043", "044") ~ 5L,
    root == "04" & sub3 %in% c("045", "046") ~ 6L,
    root == "04" ~ 4L,
    root == "05" ~ 7L,
    root == "06" ~ 8L,
    root == "07" ~ 9L,
    root == "08" ~ 10L,
    root == "09" ~ 11L,
    root == "10" ~ 12L,
    root == "11" ~ 13L,
    root == "12" ~ 14L,
    root %in% c("13", "15") ~ 15L,
    root == "14" ~ 16L,
    root == "16" ~ 17L,
    root == "17" ~ 18L,
    root %in% c("18", "19", "20") ~ 19L,
    TRUE ~ NA_integer_
  )
}

#' Human-readable label for an `EventRootCode3` value
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
#'
#' Preliminary labels, one per [assign_eventrootcode3()] output value.
#'
#' @param class Integer vector of `EventRootCode3` values, as returned
#'   by [assign_eventrootcode3()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 1-19.
#' @keywords internal
eventrootcode3_name <- function(class) {
  names <- c(
    "1"  = "Make a public statement",
    "2"  = "Appeal for action",
    "3"  = "Express intent to cooperate",
    "4"  = "Consult, unspecified",
    "5"  = "Consult: meet, discuss, or visit",
    "6"  = "Consult: negotiate or mediate",
    "7"  = "Engage in diplomatic cooperation",
    "8"  = "Engage in material cooperation",
    "9"  = "Provide aid",
    "10" = "Yield",
    "11" = "Investigate",
    "12" = "Demand",
    "13" = "Disapprove",
    "14" = "Reject",
    "15" = "Threaten or exhibit force posture",
    "16" = "Protest",
    "17" = "Reduce relations",
    "18" = "Coerce",
    "19" = "Assault, fight, or mass violence"
  )
  unname(names[as.character(class)])
}

#' Underlying CAMEO root/event codes for an `EventRootCode3` value
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
#'
#' For most `EventRootCode3` classes this is just the single two-digit
#' CAMEO root code the class was built from. Where a class merges
#' several roots (e.g. "Threaten or exhibit force posture"), splits a
#' root into finer sub-codes (e.g. "Consult: meet, discuss, or visit"),
#' or relocates an individual code out of its root (e.g. "Reject", which
#' picks up `016` from root 01), this lists every code that feeds into
#' it, comma-separated in the order given in [assign_eventrootcode3()]'s
#' docs.
#'
#' @param class Integer vector of `EventRootCode3` values, as returned
#'   by [assign_eventrootcode3()].
#' @return Character vector of comma-separated CAMEO root/event codes,
#'   `NA` for values outside 1-19.
#' @keywords internal
eventrootcode3_rootcodes <- function(class) {
  codes <- c(
    "1"  = "01",
    "2"  = "02",
    "3"  = "03",
    "4"  = "04, 040",
    "5"  = "041, 042, 043, 044",
    "6"  = "045, 046",
    "7"  = "05, 018, 019",
    "8"  = "06",
    "9"  = "07",
    "10" = "08",
    "11" = "09",
    "12" = "10",
    "13" = "11",
    "14" = "12, 016",
    "15" = "13, 15",
    "16" = "14",
    "17" = "16",
    "18" = "17",
    "19" = "18, 19, 20"
  )
  unname(codes[as.character(class)])
}

#' Regroup a CAMEO event code into the bilatr "EventRootCode4" scheme
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup]. Superseded by `ModifiedRootCode`; see
#' [assign_modified_root_code()].
#'
#' `EventRootCode4` is a further refinement of [assign_eventrootcode3()],
#' splitting three of its remaining merged/aggregated groupings back out
#' into their individual root/sub-codes, while keeping its relocations of
#' `016`/`018`/`019`:
#'
#' \itemize{
#'   \item `EventRootCode3`'s "Consult: meet, discuss, or visit" class
#'     (`041`-`044`) is split so that `041` ("Discuss by telephone")
#'     becomes its own class, separate from `042`-`044`.
#'   \item `EventRootCode3`'s "Threaten or exhibit force posture" class
#'     (roots 13 + 15) is split back into Threaten (root 13) and Exhibit
#'     force posture (root 15); Protest (root 14) was already its own
#'     class in [assign_eventrootcode3()] and is unaffected.
#'   \item `EventRootCode3`'s "Assault, fight, or mass violence" class
#'     (roots 18 + 19 + 20) is split back into Assault (root 18), Fight
#'     (root 19), and Use unconventional mass violence (root 20).
#'   \item `016` ("Deny responsibility") still moves out of root 01 into
#'     the Reject class (with root 12); `018`/`019` ("Make empathetic
#'     comment"/"Express accord") still move out of root 01 into the
#'     diplomatic-cooperation class (with root 05) -- same relocations
#'     as [assign_eventrootcode3()].
#'   \item Every other root maps as in [assign_eventrootcode3()].
#' }
#'
#' Classes are numbered 1-23 in ascending root-code order; see
#' [eventrootcode4_name()] for labels and [eventrootcode4_rootcodes()]
#' for the underlying CAMEO root/event codes each class draws from.
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Integer vector of `EventRootCode4` values (1-23).
#' @keywords internal
assign_eventrootcode4 <- function(code) {
  code <- stringr::str_pad(as.character(code), width = 2, pad = "0")
  root <- get_root(code)
  sub3 <- stringr::str_sub(code, 1, 3)
  dplyr::case_when(
    code == "016" ~ 15L,
    code == "018" ~ 8L,
    code == "019" ~ 8L,
    root == "01" ~ 1L,
    root == "02" ~ 2L,
    root == "03" ~ 3L,
    root == "04" & sub3 == "041" ~ 5L,
    root == "04" & sub3 %in% c("042", "043", "044") ~ 6L,
    root == "04" & sub3 %in% c("045", "046") ~ 7L,
    root == "04" ~ 4L,
    root == "05" ~ 8L,
    root == "06" ~ 9L,
    root == "07" ~ 10L,
    root == "08" ~ 11L,
    root == "09" ~ 12L,
    root == "10" ~ 13L,
    root == "11" ~ 14L,
    root == "12" ~ 15L,
    root == "13" ~ 16L,
    root == "14" ~ 17L,
    root == "15" ~ 18L,
    root == "16" ~ 19L,
    root == "17" ~ 20L,
    root == "18" ~ 21L,
    root == "19" ~ 22L,
    root == "20" ~ 23L,
    TRUE ~ NA_integer_
  )
}

#' Human-readable label for an `EventRootCode4` value
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
#'
#' Preliminary labels, one per [assign_eventrootcode4()] output value.
#'
#' @param class Integer vector of `EventRootCode4` values, as returned
#'   by [assign_eventrootcode4()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 1-23.
#' @keywords internal
eventrootcode4_name <- function(class) {
  names <- c(
    "1"  = "Make a public statement",
    "2"  = "Appeal for action",
    "3"  = "Express intent to cooperate",
    "4"  = "Consult, unspecified",
    "5"  = "Discuss by telephone",
    "6"  = "Consult: meet, discuss, or visit",
    "7"  = "Consult: negotiate or mediate",
    "8"  = "Engage in diplomatic cooperation",
    "9"  = "Engage in material cooperation",
    "10" = "Provide aid",
    "11" = "Yield",
    "12" = "Investigate",
    "13" = "Demand",
    "14" = "Disapprove",
    "15" = "Reject",
    "16" = "Threaten",
    "17" = "Protest",
    "18" = "Exhibit force posture",
    "19" = "Reduce relations",
    "20" = "Coerce",
    "21" = "Assault",
    "22" = "Fight",
    "23" = "Use unconventional mass violence"
  )
  unname(names[as.character(class)])
}

#' Underlying CAMEO root/event codes for an `EventRootCode4` value
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
#'
#' For most `EventRootCode4` classes this is just the single two-digit
#' CAMEO root code the class was built from. Where a class merges
#' several roots (e.g. "Engage in diplomatic cooperation", which picks
#' up `018`/`019` from root 01), splits a root into finer sub-codes
#' (e.g. "Consult: meet, discuss, or visit"), or relocates an individual
#' code out of its root (e.g. "Reject", which picks up `016` from root
#' 01), this lists every code that feeds into it, comma-separated in the
#' order given in [assign_eventrootcode4()]'s docs.
#'
#' @param class Integer vector of `EventRootCode4` values, as returned
#'   by [assign_eventrootcode4()].
#' @return Character vector of comma-separated CAMEO root/event codes,
#'   `NA` for values outside 1-23.
#' @keywords internal
eventrootcode4_rootcodes <- function(class) {
  codes <- c(
    "1"  = "01",
    "2"  = "02",
    "3"  = "03",
    "4"  = "04, 040",
    "5"  = "041",
    "6"  = "042, 043, 044",
    "7"  = "045, 046",
    "8"  = "05, 018, 019",
    "9"  = "06",
    "10" = "07",
    "11" = "08",
    "12" = "09",
    "13" = "10",
    "14" = "11",
    "15" = "12, 016",
    "16" = "13",
    "17" = "14",
    "18" = "15",
    "19" = "16",
    "20" = "17",
    "21" = "18",
    "22" = "19",
    "23" = "20"
  )
  unname(codes[as.character(class)])
}

#' Regroup a CAMEO event code into the bilatr "ModifiedRootCode" scheme
#'
#' `ModifiedRootCode` is a variant of [assign_eventrootcode3()] that merges
#' its "Investigate" (root 09) and "Demand" (root 10) classes back into a
#' single "Investigate or demand" class. Every other class, including the
#' `016`/`018`/`019` relocations, is unchanged from `EventRootCode3`.
#'
#' Classes are numbered 1-18 in ascending root-code order (so
#' `EventRootCode3` classes 1-11 keep their values, with 11 now covering
#' both roots 09 and 10, and `EventRootCode3` classes 13-19 shift down by
#' one to 12-18); see [modified_root_code_name()] for labels. This scheme
#' was previously named `ERC16NZ`.
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Integer vector of `ModifiedRootCode` values (1-18).
#' @keywords internal
assign_modified_root_code <- function(code) {
  erc3 <- assign_eventrootcode3(code)
  dplyr::if_else(erc3 >= 12L, erc3 - 1L, erc3)
}

#' Human-readable label for a `ModifiedRootCode` value
#'
#' Preliminary labels, one per [assign_modified_root_code()] output value.
#'
#' @param class Integer vector of `ModifiedRootCode` values, as returned by
#'   [assign_modified_root_code()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 1-18.
#' @keywords internal
modified_root_code_name <- function(class) {
  names <- c(
    "1"  = "Make a public statement",
    "2"  = "Appeal for action",
    "3"  = "Express intent to cooperate",
    "4"  = "Consult, unspecified",
    "5"  = "Consult: meet, discuss, or visit",
    "6"  = "Consult: negotiate or mediate",
    "7"  = "Engage in diplomatic cooperation",
    "8"  = "Engage in material cooperation",
    "9"  = "Provide aid",
    "10" = "Yield",
    "11" = "Investigate or demand",
    "12" = "Disapprove",
    "13" = "Reject",
    "14" = "Threaten or exhibit force posture",
    "15" = "Protest",
    "16" = "Reduce relations",
    "17" = "Coerce",
    "18" = "Assault, fight, or mass violence"
  )
  unname(names[as.character(class)])
}

#' Underlying CAMEO event codes for a `ModifiedRootCode` value
#'
#' For most `ModifiedRootCode` classes this is just the single two-digit
#' CAMEO root code the class was built from. Where a class merges several
#' roots (e.g. "Investigate or demand"), splits a root into finer
#' sub-codes (e.g. "Consult: meet, discuss, or visit"), or relocates an
#' individual code out of its root (e.g. "Reject", which picks up `016`
#' from root 01), this lists every code that feeds into it, comma-
#' separated. Root `01`'s own entry notes the codes carved out of it
#' (`016`/`018`/`019`) as an exception, since the rest of root `01` stays
#' in this class. The bare two-digit root `04` is omitted from the
#' "Consult, unspecified" entry (`"040"` only, not `"04, 040"`), since
#' root `04` itself splits across several `ModifiedRootCode` classes and
#' so is blanked at the top level; see [cameo_lookup]'s docs.
#'
#' @param class Integer vector of `ModifiedRootCode` values, as returned
#'   by [assign_modified_root_code()].
#' @return Character vector of comma-separated CAMEO root/event codes,
#'   `NA` for values outside 1-18.
#' @keywords internal
modified_root_code_eventcodes <- function(class) {
  codes <- c(
    "1"  = "01 except 016, 018, 019",
    "2"  = "02",
    "3"  = "03",
    "4"  = "040",
    "5"  = "041, 042, 043, 044",
    "6"  = "045, 046",
    "7"  = "05, 018, 019",
    "8"  = "06",
    "9"  = "07",
    "10" = "08",
    "11" = "09, 10",
    "12" = "11",
    "13" = "12, 016",
    "14" = "13, 15",
    "15" = "14",
    "16" = "16",
    "17" = "17",
    "18" = "18, 19, 20"
  )
  unname(codes[as.character(class)])
}

#' Human-readable label for a bilatr event class
#'
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
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
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
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
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
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
#' Deprecated: retained internally for backward compatibility; no longer
#' part of [cameo_lookup].
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

#' Recode a data frame of CAMEO event codes to quad/penta/modified-root
#' classes
#'
#' Joins `data` against the package's built-in [cameo_lookup] table to
#' attach `CAMEOLabel`, `GoldsteinScore`, `QuadClass`, `QuadClassName`,
#' `QuadClassEventCodes`, `PentaClass`, `PentaClassName`,
#' `PentaClassEventCodes`, `PentaClass_modified`, `ModifiedRootCode`,
#' `ModifiedRootCodeName`, and `ModifiedRootCodeEventCodes` columns.
#' `PentaClass_modified` folds low-intensity verbal cooperation
#' (Goldstein score <= 1) into its own class, which can be useful as a
#' near-neutral reference category. `ModifiedRootCode` /
#' `ModifiedRootCodeName` / `ModifiedRootCodeEventCodes` (see
#' [assign_modified_root_code()], [modified_root_code_name()], and
#' [modified_root_code_eventcodes()]) is an 18-level regrouping of the
#' CAMEO root codes (formerly named `ERC16NZ`).
#'
#' A handful of other classification schemes
#' (`EventRootCode2`/`EventRootCode3`/`EventRootCode4`, `BilatrClass`/
#' `BilatrClass2`) previously shipped in [cameo_lookup] have been retired
#' from it; their `assign_*()`/`*_name()` functions are still available
#' internally (see [assign_eventrootcode2()], [assign_bilatr_class()],
#' etc.) for backward compatibility, but are no longer attached by this
#' function.
#'
#' Any of these columns that already exist in `data` are left as they are
#' (not overwritten, no `.x`/`.y` suffixing), with a warning naming them.
#'
#' `code_col` must be character or factor. A numeric `EventCode` has
#' already lost any leading zero (`"010"` -> `10`) by the time it reaches
#' this function, and that loss is not recoverable here -- `stop()`s with
#' a message to that effect rather than failing inside `left_join()` with
#' an opaque type-mismatch error. After the join, `message()`s the number
#' of rows whose code did not match [cameo_lookup] at all, and separately
#' the number that matched but got `NA` `ModifiedRootCode` because they
#' are a bare top-level root (e.g. `"01"`, `"04"`; see [cameo_lookup]'s
#' docs) -- one line each, nothing when both are zero.
#'
#' @param data A data frame containing a CAMEO event code column.
#' @param code_col Name of the column in `data` holding CAMEO event codes
#'   (as a string). Defaults to `"EventCode"`, matching the raw GDELT
#'   column name. Must be character or factor (see Details).
#' @return `data` with the recode columns above attached (minus any that
#'   were already present).
#' @examples
#' \dontrun{
#' events <- data.frame(EventCode = c("01", "190", "0862"))
#' recode_cameo(events)
#' }
#' @export
recode_cameo <- function(data, code_col = "EventCode") {
  code <- data[[code_col]]
  if (is.factor(code)) {
    data[[code_col]] <- as.character(code)
  } else if (!is.character(code)) {
    stop(
      "recode_cameo(): `", code_col, "` must be character (or factor), ",
      "not ", class(code)[1], ". A numeric CAMEO code has already lost ",
      "any leading zero (e.g. \"010\" -> 10) by the time it reaches this ",
      "function, and that is not recoverable here -- convert the raw code ",
      "to a zero-padded character string first.",
      call. = FALSE
    )
  }

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

  out <- dplyr::left_join(
    data,
    lookup,
    by = rlang::set_names("CAMEOEVENTCODE", code_col)
  )

  if ("QuadClass" %in% names(out)) {
    n_unmatched <- sum(is.na(out$QuadClass))
    if (n_unmatched > 0) {
      message(
        "recode_cameo(): ", n_unmatched,
        " row(s) did not match any code in cameo_lookup."
      )
    }
  }
  if ("ModifiedRootCode" %in% names(out) && "QuadClass" %in% names(out)) {
    n_bare_root <- sum(!is.na(out$QuadClass) & is.na(out$ModifiedRootCode))
    if (n_bare_root > 0) {
      message(
        "recode_cameo(): ", n_bare_root,
        " row(s) matched a bare top-level CAMEO root code with no ",
        "ModifiedRootCode (its sub-codes span several classes)."
      )
    }
  }

  out
}
