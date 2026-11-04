#' Extract the two-digit CAMEO root code
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Character vector of two-digit root codes.
#' @noRd
get_root <- function(code) {
  stringr::str_sub(stringr::str_pad(as.character(code), width = 2, pad = "0"), 1, 2)
}

#' Assign QuadClass from a CAMEO root code
#'
#' The same mapping as GDELT's own `QuadClass` column: roots 01-05 are
#' 1 (verbal cooperation), 06-09 are 2 (material cooperation), 10-14 are
#' 3 (verbal conflict) and 15-20 are 4 (material conflict). Checked
#' against raw GDELT 1.0 daily, monthly and yearly event files.
#'
#' @param root Character vector of two-digit CAMEO root codes, as returned
#'   by [get_root()].
#' @return Integer vector with values 1 (verbal cooperation), 2 (material
#'   cooperation), 3 (verbal conflict), 4 (material conflict), or `NA`.
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' As QuadClass, except that roots 01-02 (make public statement, appeal)
#' form class 0 "Make statement", root 14 (protest) moves to class 4
#' (material conflict), and root 16 (reduce relations) moves to class 3
#' (verbal conflict). The classes are 0 "Make statement", 1 "Verbal
#' cooperation" (roots 03-05), 2 "Material cooperation" (06-09),
#' 3 "Verbal conflict" (10-13, 16) and 4 "Material conflict" (14, 15,
#' 17-20).
#'
#' @param root Character vector of two-digit CAMEO root codes.
#' @param quad Integer vector of QuadClass values, as returned by
#'   [assign_quad()].
#' @return Integer vector with values 0-4.
#' @noRd
assign_penta <- function(root, quad) {
  root_num <- as.integer(root)
  dplyr::case_when(
    # Make statement: make public statement (01) and appeal (02).
    root_num %in% c(1, 2) ~ 0L,
    # Protest (14) is material conflict, not QuadClass's verbal conflict.
    root_num == 14 ~ 4L,
    # Reduce relations (16) is verbal conflict, not QuadClass's material conflict.
    root_num == 16 ~ 3L,
    TRUE ~ quad
  )
}

#' Human-readable label for a `PentaClass` value
#'
#' @param class Integer vector of `PentaClass` values, as returned by
#'   [assign_penta()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 0-4.
#' @noRd
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
#' @noRd
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

#' Regroup a CAMEO event code into the bilatr "ModifiedRootCode" scheme
#'
#' An 18-class regrouping of the 20 CAMEO root codes. Relative to the
#' roots, it:
#'
#' * moves `016` (deny responsibility) into Reject (class 13), and `018`
#'   (make empathetic comment) and `019` (express accord) into diplomatic
#'   cooperation (class 7);
#' * splits root `04` (consult) into unspecified (`040`, class 4),
#'   meet, discuss, or visit (`041`-`044`, class 5) and negotiate or
#'   mediate (`045`-`046`, class 6);
#' * merges Investigate (`09`) and Demand (`10`) into "Investigate or
#'   demand" (class 11);
#' * merges Threaten (`13`) and Exhibit force posture (`15`) into
#'   "Threaten or exhibit force posture" (class 14);
#' * merges `18`, `19` and `20` into "Assault, fight, or mass violence"
#'   (class 18).
#'
#' Classes are numbered 1-18 in ascending root-code order; see
#' [modified_root_code_name()] for labels.
#'
#' @param code Character or numeric vector of CAMEO event codes (e.g.
#'   `"0211"`, `"19"`).
#' @return Integer vector of `ModifiedRootCode` values (1-18).
#' @noRd
assign_modified_root_code <- function(code) {
  code <- stringr::str_pad(as.character(code), width = 2, pad = "0")
  root <- get_root(code)
  sub3 <- stringr::str_sub(code, 1, 3)
  dplyr::case_when(
    code == "016" ~ 13L,
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
    root == "10" ~ 11L,
    root == "11" ~ 12L,
    root == "12" ~ 13L,
    root %in% c("13", "15") ~ 14L,
    root == "14" ~ 15L,
    root == "16" ~ 16L,
    root == "17" ~ 17L,
    root %in% c("18", "19", "20") ~ 18L,
    TRUE ~ NA_integer_
  )
}

#' Human-readable label for a `ModifiedRootCode` value
#'
#' Preliminary labels, one per [assign_modified_root_code()] output value.
#'
#' @param class Integer vector of `ModifiedRootCode` values, as returned by
#'   [assign_modified_root_code()].
#' @return Character vector of the corresponding class names, `NA` for
#'   values outside 1-18.
#' @noRd
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
#' @noRd
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

#' Recode a data frame of CAMEO event codes to quad/penta/modified-root
#' classes
#'
#' Joins `data` against the package's built-in [cameo_lookup] table to
#' attach 11 columns: `CAMEOLabel`, `GoldsteinScore`, `QuadClass`,
#' `QuadClassName`, `QuadClassEventCodes`, `PentaClass`, `PentaClassName`,
#' `PentaClassEventCodes`, `ModifiedRootCode`, `ModifiedRootCodeName`,
#' and `ModifiedRootCodeEventCodes`. See
#' [cameo_lookup] for the schemes:
#'
#' * `QuadClass` is the same mapping as GDELT's own `QuadClass` (roots
#'   01-05 verbal cooperation, 06-09 material cooperation, 10-14 verbal
#'   conflict, 15-20 material conflict);
#' * `PentaClass` is QuadClass, except that roots 01-02 form class 0
#'   "Make statement", root 14 (protest) moves to class 4 (material
#'   conflict), and root 16 (reduce relations) moves to class 3 (verbal
#'   conflict);
#' * `ModifiedRootCode` is an 18-class regrouping of the root codes that
#'   splits Consult and merges some conflict roots (see [cameo_lookup]).
#'
#' [event_class_labels()] labels the classes of a fitted model's
#' `stan_data` from these columns.
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
#' events <- data.frame(EventCode = c("01", "190", "0862"))
#' recode_cameo(events)
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

#' Human-readable labels for a Stan data object's event classes
#'
#' Reads the labels from [cameo_lookup]'s `*Name` column for `scheme`,
#' keyed by event class value, never by position in `event_classes`, so
#' the result is correct wherever `reference_category` put a class in that
#' order. The other post-fit functions use it for their `class_labels`
#' default.
#'
#' `scheme` names a [cameo_lookup] column with a matching `*Name` column:
#' currently `"QuadClass"` (`QuadClassName`), `"PentaClass"`
#' (`PentaClassName`), and `"ModifiedRootCode"` (`ModifiedRootCodeName`).
#' For any other `scheme` -- `NULL`, an unrecognized name, or a
#' non-CAMEO grouping variable from the caller's own data (e.g. a
#' `"sector"` column you coded yourself) -- the event
#' classes' own codes are used as their labels, since there is no lookup
#' to consult; if `scheme` was non-`NULL` but unrecognized, a `message()`
#' names the schemes that do have labels.
#'
#' If the event classes come from a renamed CAMEO column (for example,
#' `ModifiedRootCode` copied into a column of another name before
#' aggregating), pass `scheme` explicitly, since the default
#' `attr(stan_data, "grouping_var")` is then not a name [cameo_lookup]
#' recognizes.
#'
#' @param stan_data A Stan data list from [prepare_stan_data()] (or
#'   `attr(fit, "stan_data")`; must carry its `"event_classes"`
#'   attribute).
#' @param scheme Name of a [cameo_lookup] column with a matching `*Name`
#'   column (see Details), or anything else to fall back to the raw event
#'   classes. Defaults to `attr(stan_data, "grouping_var")`.
#' @return A named character vector: names are `attr(stan_data,
#'   "event_classes")`, in order; values are the corresponding labels. An
#'   event class missing from the lookup (possible for a real CAMEO
#'   `scheme` if the classes came from a non-`cameo_lookup` source) keeps
#'   its code as its own label, with a `warning()` naming it.
#' @examples
#' # the Stan data of the example fit, aggregated by PentaClass
#' csv <- bilatr_example_fit()
#' stan_data <- readRDS(file.path(dirname(csv[1]), "example_stan_data.rds"))
#' event_class_labels(stan_data)
#'
#' # a grouping variable with no labels in cameo_lookup: the codes are used
#' event_class_labels(stan_data, scheme = "sector")
#' @export
event_class_labels <- function(stan_data, scheme = attr(stan_data, "grouping_var")) {
  event_classes <- attr(stan_data, "event_classes")
  if (is.null(event_classes)) {
    stop(
      "event_class_labels(): `stan_data` has no \"event_classes\" ",
      "attribute -- pass a Stan data list from prepare_stan_data(), or ",
      "attr(fit, \"stan_data\") for a fit from fit_panel()/fit_dyad_ts().",
      call. = FALSE
    )
  }
  key <- as.character(event_classes)

  name_col <- if (!is.null(scheme)) {
    switch(scheme,
      QuadClass = "QuadClassName",
      PentaClass = "PentaClassName",
      ModifiedRootCode = "ModifiedRootCodeName",
      NULL
    )
  } else {
    NULL
  }

  if (is.null(name_col)) {
    if (!is.null(scheme)) {
      message(
        "event_class_labels(): scheme '", scheme, "' has no labels in ",
        "cameo_lookup; using the raw event classes as labels. Schemes ",
        "with labels: QuadClass, PentaClass, ModifiedRootCode."
      )
    }
    return(rlang::set_names(key, key))
  }

  lookup <- bilatr::cameo_lookup
  pairs <- unique(stats::na.omit(lookup[c(scheme, name_col)]))
  lookup_vec <- rlang::set_names(as.character(pairs[[name_col]]), as.character(pairs[[scheme]]))

  labels <- unname(lookup_vec[key])
  missing <- key[is.na(labels)]
  if (length(missing) > 0) {
    warning(
      "event_class_labels(): event class(es) not found in cameo_lookup's ",
      name_col, " column, keeping their code as the label: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
    labels[is.na(labels)] <- missing
  }

  rlang::set_names(labels, key)
}

#' Resolve a `class_labels` argument to one label per event class, in order
#'
#' Shared by [diagnose_category_merges()], [extract_gamma()],
#' [icc_curves()], and [check_compositional_residuals()]. Always looks
#' labels up by event class VALUE, never by position, so a
#' `class_labels` vector (or the default from [event_class_labels()])
#' need not share `event_classes`' order.
#'
#' @param class_labels `NULL`, a character vector named by event class
#'   (must cover every value in `event_classes`; errors listing what's
#'   missing otherwise), or an unnamed character vector the same length
#'   as `event_classes`, in `event_classes` order (trusted as-is, since
#'   there is no class value to key an unnamed vector by).
#' @param event_classes Character vector of event classes, in the order
#'   the caller's output rows use.
#' @param stan_data Optional. When `class_labels` is `NULL` and this is
#'   supplied, defaults to `event_class_labels(stan_data)` (see
#'   [event_class_labels()]). Falls back to `event_classes` themselves
#'   when `stan_data` is `NULL` (possible for [icc_curves()]), or when it
#'   lacks a usable `"event_classes"` attribute (e.g. a hand-built
#'   `stan_data` list) --
#'   `class_labels`/`event_classes` are the caller-supplied fallback for
#'   exactly that case. Only the missing-attribute case is treated this
#'   way; a genuine error raised by [event_class_labels()] itself (once
#'   the attribute is present) is not caught here and propagates.
#' @return Character vector, `length(event_classes)`, in `event_classes`
#'   order.
#' @noRd
.resolve_class_labels <- function(class_labels, event_classes, stan_data = NULL) {
  if (is.null(class_labels)) {
    class_labels <- if (!is.null(stan_data) && !is.null(attr(stan_data, "event_classes"))) {
      event_class_labels(stan_data)
    } else {
      NULL
    }
    if (is.null(class_labels)) {
      return(as.character(event_classes))
    }
  }

  if (is.null(names(class_labels))) {
    if (length(class_labels) != length(event_classes)) {
      stop(
        "`class_labels` has length ", length(class_labels), " but there ",
        "are ", length(event_classes), " event classes; supply one label ",
        "per class (in event_classes order) or a vector named by event ",
        "class.",
        call. = FALSE
      )
    }
    return(as.character(class_labels))
  }

  missing <- setdiff(as.character(event_classes), names(class_labels))
  if (length(missing) > 0) {
    stop(
      "`class_labels` is missing a label for event class(es): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  unname(as.character(class_labels[as.character(event_classes)]))
}
