#' The bilatr dyad-period data contract
#'
#' [prepare_stan_data()] takes an already aggregated dyad x time-period
#' table of action counts, from any source, in the format documented
#' here, and validates it with [validate_bilatr_dyad_period()] before
#' doing anything else. [aggregate_events()] produces a table in this
#' format from event-level data (see `?bilatr_event_data`); any other
#' source works too, so long as it meets the rules below. This topic is
#' the single place the contract is spelled out.
#'
#' | Role | Default column | Requirement |
#' | --- | --- | --- |
#' | Dyad | `dyad` | no `NA`; no empty string; no leading/trailing whitespace. Optional if both actor columns are given: it is then built as `paste(actor1, actor2, sep = "_")` |
#' | Actor 1 / Actor 2 | `actor1`, `actor2` | optional (both or neither); no `NA`; no empty string; no leading/trailing whitespace; must not contain `"_"`; `actor1 != actor2`; constant within a dyad (each dyad maps to exactly one actor pair). Needed for the country-level offset of the `stable_gamma` model |
#' | Time | chosen by the user (`time`) | integer-valued (integer, or a whole-number double); one unit = one period; no `NA`. Gaps are allowed: every dyad is filled over the whole `min(time):max(time)` range, missing periods becoming unobserved |
#' | Sender share | `w_send` | optional; numeric in `[0, 1]`; constant within a dyad. The share of the dyad's events sent by `actor1`; only used by `stable_gamma` for undirected dyads |
#' | Action counts | `action_cols` | at least two columns; non-negative whole numbers; no `NA` |
#'
#' There must be at most one row per dyad and time period: duplicates are
#' an error, never summed.
#'
#' `action_cols` takes a character vector or a tidyselect expression
#' (e.g. `c(coop, conflict)`, `starts_with("n_")`). If it is not given
#' and `data` carries an `action_cols` attribute (as [aggregate_events()]
#' output does) naming only columns still present, those columns are
#' used silently. Otherwise every column except the dyad, actor, time
#' and `w_send` columns is used,
#' with a warning naming the columns chosen; a column picked up this way
#' that is not a non-negative whole-number count is an error. When
#' `action_cols` is given explicitly, every other column is ignored.
#'
#' The time column holds absolute period IDs, so tables built separately
#' can be compared. [aggregate_events()] writes `year`, `month_id`
#' (`year * 12 + month - 1`), `week_id` or `day_id` (whole weeks or days
#' since an origin date, by default 1946-01-01), but any integer period
#' index works.
#'
#' @examples
#' # Three dyads with arbitrary labels and no actor columns, observed in
#' # periods 1-6 with a gap, and three action types.
#' dp <- tibble::tibble(
#'   dyad = c("a-b", "a-b", "a-b", "a-c", "a-c", "b-c"),
#'   period = c(1L, 2L, 6L, 1L, 3L, 4L),
#'   talks = c(3, 0, 2, 1, 4, 2),
#'   threats = c(0, 1, 0, 2, 0, 1),
#'   aid = c(1, 1, 0, 0, 2, 0)
#' )
#' validate_bilatr_dyad_period(dp, time = "period", action_cols = c(talks, threats, aid))
#'
#' @name bilatr_dyad_period_data
NULL

#' Validate a data frame against the bilatr dyad-period data contract
#'
#' Checks `data` against every rule in `?bilatr_dyad_period_data`,
#' collecting every broken rule rather than stopping at the first, and
#' `stop()`s once with a bulleted message naming each one. Missing
#' columns (other than the optional ones) are reported up front, before
#' any other rule is checked. The actor and `w_send` columns are
#' validated only if present. Called first thing by [prepare_stan_data()].
#'
#' @param data A data frame to validate.
#' @param time Name of the time-period column.
#' @param action_cols The action-count columns: a character vector of
#'   column names or a tidyselect expression. `NULL` (default) means every
#'   column except `dyad`, `actor1`, `actor2`, `time` and `w_send`, with a
#'   warning listing the columns chosen, unless `data` carries a usable
#'   `action_cols` attribute; see `?bilatr_dyad_period_data`.
#' @param dyad Name of the dyad column. It may be absent if both actor
#'   columns are present (it is then built from them).
#' @param actor1,actor2 Names of the optional actor columns.
#' @param w_send Name of the optional sender-share column.
#' @return Invisibly, if every rule passes, a list of the resolved column
#'   names: `dyad`, `time`, `action_cols` (character, in `data`'s
#'   selection order), and `actor1`, `actor2`, `w_send` (each `NULL` if
#'   the column is absent).
#' @export
validate_bilatr_dyad_period <- function(
  data,
  time,
  action_cols = NULL,
  dyad = "dyad",
  actor1 = "actor1",
  actor2 = "actor2",
  w_send = "w_send"
) {
  .validate_dyad_period(data, time, rlang::enquo(action_cols), dyad, actor1, actor2, w_send)
}

#' Worker for [validate_bilatr_dyad_period()]
#'
#' Takes `action_cols` as an already captured quosure, so that
#' [prepare_stan_data()] can forward its own `action_cols` argument
#' unevaluated.
#'
#' @inheritParams validate_bilatr_dyad_period
#' @param action_quo A quosure of the `action_cols` expression.
#' @return See [validate_bilatr_dyad_period()].
#' @keywords internal
.validate_dyad_period <- function(data, time, action_quo, dyad, actor1, actor2, w_send) {
  fn <- "validate_bilatr_dyad_period()"
  if (!is.data.frame(data)) {
    stop(fn, ": `data` must be a data frame.", call. = FALSE)
  }
  if (!rlang::is_string(time)) {
    stop(fn, ": `time` must be a single column name.", call. = FALSE)
  }

  has_actor1 <- actor1 %in% names(data)
  has_actor2 <- actor2 %in% names(data)
  has_actors <- has_actor1 && has_actor2
  has_dyad <- dyad %in% names(data)
  has_w_send <- w_send %in% names(data)

  structural <- character(0)
  if (!time %in% names(data)) {
    structural <- c(structural, sprintf("- missing time column `%s`", time))
  }
  if (has_actor1 != has_actor2) {
    structural <- c(structural, sprintf(
      "- only one actor column is present (`%s`); give both `%s` and `%s`, or neither",
      if (has_actor1) actor1 else actor2, actor1, actor2
    ))
  }
  if (!has_dyad && !has_actors) {
    structural <- c(structural, sprintf(
      "- missing dyad column `%s` (it can be omitted only if both actor columns `%s` and `%s` are present)",
      dyad, actor1, actor2
    ))
  }
  if (length(structural) > 0) {
    stop(fn, ": ", paste(structural, collapse = "\n"), call. = FALSE)
  }

  role_cols <- c(dyad, actor1, actor2, time, w_send)
  action <- .resolve_action_cols(action_quo, data, role_cols)
  action_cols <- action$cols

  bullets <- character(0)
  add_rule <- function(label, idx) {
    if (any(idx)) {
      bullets <<- c(bullets, sprintf(
        "- %s (%d row(s); e.g. row(s) %s)",
        label, sum(idx), paste(utils::head(which(idx), 5), collapse = ", ")
      ))
    }
  }
  add_dyad_rule <- function(label, bad_dyads) {
    if (length(bad_dyads) > 0) {
      bullets <<- c(bullets, sprintf(
        "- %s (%d dyad(s); e.g. %s)",
        label, length(bad_dyads), paste(utils::head(bad_dyads, 5), collapse = ", ")
      ))
    }
  }
  # Dyads whose rows carry more than one distinct value of `x`.
  nonconstant_dyads <- function(dyad_vec, x) {
    tibble::tibble(dyad = dyad_vec, x = x) %>%
      dplyr::distinct() %>%
      dplyr::count(.data$dyad) %>%
      dplyr::filter(.data$n > 1) %>%
      dplyr::pull(.data$dyad)
  }

  if (has_actors) {
    a1 <- as.character(data[[actor1]])
    a2 <- as.character(data[[actor2]])
    for (role in list(list(actor1, a1), list(actor2, a2))) {
      add_rule(sprintf("`%s` has missing (NA) values", role[[1]]), is.na(role[[2]]))
      add_rule(sprintf("`%s` has empty-string values", role[[1]]), !is.na(role[[2]]) & role[[2]] == "")
      add_rule(
        sprintf("`%s` contains \"_\", the dyad-key separator", role[[1]]),
        !is.na(role[[2]]) & grepl("_", role[[2]], fixed = TRUE)
      )
      add_rule(
        sprintf("`%s` has leading/trailing whitespace (e.g. \" USA\" and \"USA\" would otherwise become two actors)", role[[1]]),
        !is.na(role[[2]]) & grepl("^\\s|\\s$", role[[2]])
      )
    }
    add_rule(
      sprintf("`%s` equals `%s` (a self-dyad; bilatr models bilateral relations)", actor1, actor2),
      !is.na(a1) & !is.na(a2) & a1 == a2
    )
  }

  dyad_vec <- if (has_dyad) as.character(data[[dyad]]) else paste(a1, a2, sep = "_")
  if (has_dyad) {
    add_rule(sprintf("`%s` has missing (NA) values", dyad), is.na(dyad_vec))
    add_rule(sprintf("`%s` has empty-string values", dyad), !is.na(dyad_vec) & dyad_vec == "")
    add_rule(
      sprintf("`%s` has leading/trailing whitespace (e.g. \" USA_CHN\" and \"USA_CHN\" would otherwise become two dyads)", dyad),
      !is.na(dyad_vec) & grepl("^\\s|\\s$", dyad_vec)
    )
    if (has_actors) {
      add_dyad_rule(
        sprintf("`%s`/`%s` are not constant within a dyad (each dyad must map to exactly one actor pair)", actor1, actor2),
        nonconstant_dyads(dyad_vec, paste(a1, a2, sep = "\r"))
      )
    }
  }

  time_vec <- data[[time]]
  if (!is.numeric(time_vec)) {
    bullets <- c(bullets, sprintf(
      "- `%s` must be integer-valued (an integer or whole-number double column), not %s",
      time, class(time_vec)[1]
    ))
  } else {
    add_rule(sprintf("`%s` has missing (NA) values", time), is.na(time_vec))
    add_rule(
      sprintf("`%s` is not integer-valued", time),
      !is.na(time_vec) & !.is_whole_number(time_vec)
    )
  }

  if (has_w_send) {
    w <- data[[w_send]]
    if (!is.numeric(w)) {
      bullets <- c(bullets, sprintf("- `%s` must be numeric, not %s", w_send, class(w)[1]))
    } else {
      add_rule(sprintf("`%s` has missing (NA) values", w_send), is.na(w))
      add_rule(sprintf("`%s` is outside [0, 1]", w_send), !is.na(w) & (w < 0 | w > 1))
      add_dyad_rule(sprintf("`%s` is not constant within a dyad", w_send), nonconstant_dyads(dyad_vec, w))
    }
  }

  if (length(action_cols) < 2) {
    bullets <- c(bullets, sprintf(
      "- at least two action-count columns are needed; `action_cols` selected %d%s",
      length(action_cols),
      if (length(action_cols) > 0) paste0(" (", paste(action_cols, collapse = ", "), ")") else ""
    ))
  }
  clashing <- intersect(action_cols, role_cols)
  if (length(clashing) > 0) {
    bullets <- c(bullets, sprintf(
      "- `action_cols` includes the dyad/actor/time/w_send column(s) %s",
      paste0("`", clashing, "`", collapse = ", ")
    ))
  }
  for (col in setdiff(action_cols, clashing)) {
    x <- data[[col]]
    bad_type <- !is.numeric(x)
    bad_na <- if (bad_type) logical(0) else is.na(x)
    bad_value <- if (bad_type) logical(0) else !is.na(x) & (!.is_whole_number(x) | x < 0)
    if (action$default) {
      if (bad_type || any(bad_na) || any(bad_value)) {
        bullets <- c(bullets, sprintf(
          paste0(
            "- `%s` was picked up by the default `action_cols` (every column except ",
            "the dyad/actor/time/w_send columns), but is not a column of non-negative ",
            "whole-number counts without NA; set `action_cols` explicitly"
          ),
          col
        ))
      }
    } else if (bad_type) {
      bullets <- c(bullets, sprintf("- action column `%s` must be numeric, not %s", col, class(x)[1]))
    } else {
      add_rule(sprintf("action column `%s` has missing (NA) values", col), bad_na)
      add_rule(sprintf("action column `%s` has values that are not non-negative whole numbers", col), bad_value)
    }
  }

  dup <- tibble::tibble(dyad = dyad_vec, time = time_vec) %>%
    dplyr::count(.data$dyad, .data$time) %>%
    dplyr::filter(.data$n > 1)
  add_dyad_rule(
    sprintf("duplicate dyad x `%s` rows (one row per dyad and period is required; counts are never summed)", time),
    unique(dup$dyad)
  )

  if (length(bullets) > 0) {
    stop(
      fn, ": the following rule(s) were broken:\n",
      paste(bullets, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(list(
    dyad = dyad,
    time = time,
    action_cols = action_cols,
    actor1 = if (has_actors) actor1,
    actor2 = if (has_actors) actor2,
    w_send = if (has_w_send) w_send
  ))
}

#' Resolve `action_cols` to column names
#'
#' `NULL` first looks for an `action_cols` attribute on `data` (set by
#' [aggregate_events()]): if it names only columns still present in
#' `data`, it is used silently. Otherwise `NULL` selects every column of
#' `data` not in `role_cols`, with a `warning()` listing them. Anything else goes through
#' [tidyselect::eval_select()]; a bare symbol that is not a column of
#' `data` but names a character vector in the calling environment is
#' wrapped in `all_of()` first, so passing a variable holding column
#' names does not trigger tidyselect's "external vector" deprecation.
#'
#' @param action_quo A quosure of the `action_cols` expression.
#' @param data The data frame.
#' @param role_cols The dyad/actor/time/w_send column names.
#' @return A list: `cols` (character) and `default` (`TRUE` if selected by
#'   the `NULL` default).
#' @keywords internal
.resolve_action_cols <- function(action_quo, data, role_cols) {
  if (rlang::quo_is_null(action_quo)) {
    from_attr <- attr(data, "action_cols")
    if (is.character(from_attr) && length(from_attr) > 0 && all(from_attr %in% names(data))) {
      return(list(cols = from_attr, default = FALSE))
    }
    cols <- setdiff(names(data), role_cols)
    warning(
      "`action_cols` not set; using every column except the dyad/actor/time/w_send ",
      "columns as action counts: ", paste(cols, collapse = ", "),
      ". Set `action_cols` explicitly to silence this warning.",
      call. = FALSE
    )
    return(list(cols = cols, default = TRUE))
  }
  if (rlang::quo_is_symbol(action_quo) &&
      !rlang::as_string(rlang::quo_get_expr(action_quo)) %in% names(data)) {
    value <- tryCatch(rlang::eval_tidy(action_quo), error = function(e) NULL)
    if (is.character(value)) {
      action_quo <- rlang::quo(dplyr::all_of(!!value))
    }
  }
  pos <- tidyselect::eval_select(action_quo, data, allow_rename = FALSE)
  list(cols = names(data)[pos], default = FALSE)
}

#' Whole-number test for a numeric vector
#'
#' `TRUE` where `x` is finite, a whole number, and within integer range
#' (so `as.integer()` is lossless). `NA` gives `FALSE` for a double and
#' `TRUE` for an integer; callers test for `NA` separately.
#'
#' @param x A numeric vector.
#' @return A logical vector.
#' @keywords internal
.is_whole_number <- function(x) {
  if (is.integer(x)) {
    return(rep(TRUE, length(x)))
  }
  is.finite(x) & x == trunc(x) & abs(x) <= .Machine$integer.max
}
