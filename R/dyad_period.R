#' Validate a reference-category argument against the observed classes
#'
#' @param value The candidate reference value, or `NULL`.
#' @param classes Vector of observed event-class values.
#' @param arg_name Name of the argument, used in the error message.
#' @param null_hint What `NULL` anchors on instead, for the error message.
#'   The default describes [order_event_classes()]'s C-locale numeric
#'   sort; [prepare_stan_data()] passes "the first action column".
#' @return `value` (coerced to character) if it is present in `classes`.
#'   `NULL` is returned as-is: the only way to say "no preference", which
#'   in [prepare_stan_data()] anchors on the first action column.
#' @keywords internal
validate_reference_class <- function(value, classes, arg_name,
                                     null_hint = "the first class instead (C-locale numeric sort)") {
  if (is.null(value)) {
    return(NULL)
  }
  value <- as.character(value)
  observed <- unique(as.character(classes))
  if (!(value %in% observed)) {
    stop(
      arg_name, " = '", value, "' was requested, but is not present in the ",
      "(in-window) data. Classes present: ",
      paste(sort(observed), collapse = ", "),
      ". `", arg_name, " = NULL` anchors on ", null_hint, ".",
      call. = FALSE
    )
  }
  value
}

#' Order event classes with the reference category first
#'
#' Puts `reference_category` first (if supplied and present), with all other
#' classes after it in C-locale, numeric-aware order (`"2"` before `"10"`,
#' digits before letters, upper case before lower case). This ordering is what implements
#' the model's identification constraint on the R side: `stable`/`ou` fold
#' `alpha[1]`'s (the first column's) sign into the reported `alpha`/`theta`
#' (see each `.stan` file's header, "IDENTIFICATION: ORIENTATION FOLD"),
#' so the reference class anchors the sign of the latent scale (see the
#' `reference_category` argument of [prepare_stan_data()]). All
#' remaining `alpha[2:A]` are freely estimated.
#'
#' @param classes Vector of observed event-class values.
#' @param reference_category Value to place first, or `NULL`.
#' @return Character vector of unique classes in anchor order.
#' @keywords internal
order_event_classes <- function(classes, reference_category = NULL) {
  # locale = "C" (0.7.1): str_sort()'s default locale is the system's
  # ("en" here), which is locale-COLLATION-dependent -- and this ordering
  # *defines* what alpha/mu_intercept/gamma's action-index columns mean
  # (the same class of concern as prepare_stan_data()'s country index;
  # see R/stan_data.R). "C" gives byte/ASCII-order collation, so a
  # stan_data.rds re-derived on a machine with a different locale can't
  # silently relabel action classes; `numeric = TRUE` still gives the
  # intended natural/numeric-aware ordering ("2" before "10") under "C".
  classes <- stringr::str_sort(unique(as.character(classes)), numeric = TRUE, locale = "C")
  middle <- setdiff(classes, reference_category)
  c(reference_category, middle)
}

#' Build the `dyad_ids` attribute for [prepare_stan_data()]
#'
#' Recovers the `dyad_id` <-> `dyad` (and undirected `dyad2`) mapping, so
#' posterior draws indexed by `dyad_id` can be joined back to readable
#' identifiers (see [extract_theta()]). Builds the table directly from the
#' (already sorted, already `min_n_events`-filtered) dyads and the full
#' time range, so its rows are in the same order as `Y`'s: dyad
#' (`dyad_id`), then time. Until 0.12.0 this was `.make_dyad_ids_generic()`;
#' the `make_dyad_ids()` of bilatr 0.11.0 and earlier built the old
#' `actor_a`/`actor_b` table.
#'
#' @param dyads Character vector of retained dyads, in `dyad_id` order.
#' @param actor1,actor2 Character vectors of each dyad's actors (aligned
#'   with `dyads`), or `NULL` when the data has no actor columns.
#' @param time_values Integer vector of the time-column values, `T` long.
#' @param time_col Name of the time column.
#' @param time_labels Optional data frame keyed by `time_col`; its other
#'   columns are left-joined on.
#' @return A tibble with `dyad_id`, `time_index`, `dyad`, `dyad2`
#'   (the undirected dyad key), `actor1`, `actor2` (`dyad2`/`actor1`/
#'   `actor2` are `NA` when `actor1` is `NULL`), the time column, and any
#'   `time_labels` columns.
#' @keywords internal
make_dyad_ids <- function(dyads, actor1, actor2, time_values, time_col, time_labels = NULL) {
  D <- length(dyads)
  Tn <- length(time_values)
  if (is.null(actor1)) {
    actor1 <- actor2 <- dyad2 <- rep(NA_character_, D)
  } else {
    pair <- .order_pair_c_locale(actor1, actor2)
    dyad2 <- paste(pair$actor_a, pair$actor_b, sep = "_")
  }

  ids <- tibble::tibble(
    dyad_id = rep(seq_len(D), each = Tn),
    time_index = rep(seq_len(Tn), times = D),
    dyad = rep(dyads, each = Tn),
    dyad2 = rep(dyad2, each = Tn),
    actor1 = rep(actor1, each = Tn),
    actor2 = rep(actor2, each = Tn),
    !!time_col := rep(time_values, times = D)
  )

  if (is.null(time_labels)) {
    return(ids)
  }
  if (!is.data.frame(time_labels) || !time_col %in% names(time_labels)) {
    stop(
      "prepare_stan_data(): `time_labels` must be a data frame with a `", time_col,
      "` column (or NULL).",
      call. = FALSE
    )
  }
  clashing <- setdiff(intersect(names(time_labels), names(ids)), time_col)
  if (length(clashing) > 0) {
    stop(
      "prepare_stan_data(): `time_labels` columns clash with `dyad_ids` columns: ",
      paste(clashing, collapse = ", "),
      call. = FALSE
    )
  }
  dplyr::left_join(ids, time_labels, by = time_col, relationship = "many-to-one")
}
