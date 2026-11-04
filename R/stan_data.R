#' Prepare Stan-ready data from a dyad-period table
#'
#' Packages a dyad x time-period table of action counts (see
#' `?bilatr_dyad_period_data`), from [aggregate_events()] or any other
#' source, as the data list expected by the package's Stan models (see
#' [compile_bilatr_model()], "Models").
#'
#' In order:
#'
#' 1. validates `data` with [validate_bilatr_dyad_period()], building the
#'    dyad column from the actor columns if it is absent;
#' 2. drops action columns that are zero in every row, with a `message()`
#'    naming them;
#' 3. orders the action columns: `reference_category` first, then the
#'    rest in `data`'s column order;
#' 4. fills every dyad over the whole `time_range` (by default
#'    `min(time):max(time)` of the table), missing periods becoming
#'    all-zero and unobserved (`is_obs = 0`; a period is observed when its
#'    counts sum to more than zero);
#' 5. drops dyads with fewer than `min_n_events` events in total;
#' 6. if both actor columns are present, builds the country index
#'    (`n_countries`, `ctry_a`, `ctry_b`) and `w_send`, with a `message()`
#'    reporting `n_countries` and the dyads-per-country distribution;
#' 7. resolves `filter_dyads` and builds the list and its attributes.
#'
#' @section Size:
#' `Y` is a dense `D x T x A` integer array, about `D * T * A * 4` bytes
#' (the filled grid has `D * T` rows). `T` counts every period in
#' `time_range`, observed or not, so finer time units grow it fast:
#' 16,750 dyads over 35 years with 16 action types is about 37 MB yearly,
#' 450 MB monthly and 13.7 GB daily. Daily (or similarly fine) resolution
#' is only practical for few dyads or short windows; use `time_range` to
#' pin `T`.
#'
#' @param data A dyad-period data frame; see `?bilatr_dyad_period_data`.
#' @param time Name of the time-period column. `NULL` (default) takes
#'   `attr(data, "time_col")` (set by [aggregate_events()]); failing that,
#'   whichever single one of `year`, `month_id`, `week_id` or `day_id` is
#'   present, with a `message()`; otherwise errors.
#' @param action_cols The action-count columns: a character vector or a
#'   tidyselect expression. `NULL` (default) means every column except
#'   the dyad/actor/time/`w_send` columns, with a warning; see
#'   `?bilatr_dyad_period_data`.
#' @param dyad,actor1,actor2,w_send Names of the dyad, actor and
#'   sender-share columns; see `?bilatr_dyad_period_data`. The actor and
#'   `w_send` columns are optional.
#' @param directed Whether the dyads are directed. Defaults to
#'   `attr(data, "directed")` if set (by [aggregate_events()]), else
#'   `TRUE`. Only matters for `w_send`: with `directed = TRUE` a `w_send`
#'   column must be all `1`, and a missing one is taken as `1`; with
#'   `directed = FALSE` a missing `w_send` is taken as `0.5` and the
#'   `w_send_imputed` attribute is set to `TRUE`.
#' @param reference_category Name of the action column that anchors the
#'   sign of the latent scale. \eqn{\theta} is oriented so that higher
#'   values make `reference_category` more likely relative to the other
#'   action types. With a cooperative reference category, higher
#'   \eqn{\theta} means more cooperative relations.
#'
#'   The reference category is placed first, and every bilatr model folds
#'   the sign of its discrimination `alpha[1]` into the reported `alpha`
#'   and `theta`, so the reported `alpha[1]` is never negative (see the
#'   header of each
#'   model's Stan program, "Sign: the ORIENTATION FOLD"). Since
#'   `alpha` sums to zero, `alpha[1]` is at least the average
#'   discrimination of the other classes, and the reference class gains
#'   probability relative to them as \eqn{\theta} rises.
#'
#'   `NULL` (default) anchors on the first action column, with a
#'   `message()` naming it: the sign of \eqn{\theta} then depends on
#'   whichever class that is, so set a cooperative class explicitly (for
#'   example `"1"`, verbal cooperation, for `PentaClass`). Errors if it is
#'   not one of the (non-all-zero) action columns.
#' @param min_n_events Minimum total events for a dyad to be retained: a
#'   single non-negative whole number.
#' @param time_range `NULL` (default), or two whole numbers `c(min, max)`
#'   on the time column's scale. Rows outside it are dropped (with a
#'   `message()` giving the count), and every dyad is filled over
#'   `min:max`, so `T = max - min + 1` whatever periods the data covers.
#'   `NULL` means `min:max` of the data. Errors if it is not two whole
#'   numbers, if `min > max`, or if no rows fall inside it.
#' @param time_labels A data frame keyed by the time column, whose other
#'   columns (e.g. `year`, `month`) are joined into the `dyad_ids`
#'   attribute. Defaults to `attr(data, "time_labels")` (set by
#'   [aggregate_events()]); `NULL` adds none.
#' @param grouping_var Stored as the `grouping_var` attribute, for
#'   [event_class_labels()]. Defaults to `attr(data, "grouping_var")`.
#' @param chunk_size `reduce_sum` grainsize used to chunk the likelihood
#'   across dyads: a single whole number, at least 1. The default, 1, is
#'   conservative and fine for single-threaded fits; with several threads
#'   per chain a larger value is usually faster for large panels. See
#'   [fit_panel()]'s `threads_per_chain` and `vignette("tuning", package = "bilatr")`.
#' @param rho_prior_a,rho_prior_b Shape parameters (single positive
#'   numbers) of the `beta(rho_prior_a, rho_prior_b)` prior on the
#'   OU/AR(1) persistence parameter `rho`, used only by the experimental
#'   `"ou"` model (see [compile_bilatr_model()]); ignored by `"stable"`,
#'   which keeps a random-walk `theta` with no persistence parameter.
#'   Default `8, 2` (weighted toward strong persistence).
#' @param compute_log_lik `0` (default) or `1` (`FALSE`/`TRUE` also
#'   work; stored as an integer). Data flag consumed by
#'   every registered Stan model, gating a per-dyad-period `log_lik` in
#'   `generated quantities`; left off by default since it is `D x T x`
#'   draws and file size already scales with dyad count.
#' @param prior_only `0` (default) or `1` (`FALSE`/`TRUE` also work).
#'   Data flag consumed by every
#'   registered Stan model, gating the `reduce_sum` likelihood call in
#'   `model` -- `1` fits the prior alone, useful for prior-predictive
#'   checks and for validating [alpha_prior_moments()] against an actual
#'   fit.
#' @param compute_theta_filtered `0` (default) or `1` (`FALSE`/`TRUE`
#'   also work). Data flag consumed
#'   by every registered Stan model, gating
#'   `theta_filtered`/`theta_filtered_sd` in `generated quantities`: the
#'   forward-filtered state (conditional on each draw's hyperparameters,
#'   a West-Harrison/Fisher-scoring linear-Bayes update; see each `.stan`
#'   file's `generated quantities` block) as opposed to the smoothed
#'   `theta` every fit already returns. Left off by default: computing it
#'   for every dyad is `D * T` iterations with `A` `digamma()` calls each,
#'   single-threaded (no `reduce_sum` in `generated quantities`), on the
#'   order of a minute or two per chain for a panel of about 17,000 dyads
#'   over 35 years, and adds two more `D x T`-scaled Tier 3 blocks to the
#'   output.
#' @param filter_dyads Optional character vector of dyad names (matched
#'   against the `dyad_ids` attribute's `dyad` column) naming which dyads
#'   to compute `theta_filtered` for. Ignored, with a warning, if
#'   `compute_theta_filtered` is `0`. `NULL` (default) means all dyads, once
#'   `compute_theta_filtered = 1` is set -- pass a small subset (e.g.
#'   fifty dyads) to compute the filter cheaply for just those. Errors,
#'   listing the offending names, if any requested dyad was dropped by
#'   `min_n_events` or never existed.
#' @return A named list suitable as the `data` argument to
#'   `cmdstanr::CmdStanModel$sample()`: `D`, `T`, `A`, `C` (integers),
#'   `is_obs`,
#'   `Y`, `rho_prior_a`, `rho_prior_b`, `compute_log_lik`, `prior_only`,
#'   `compute_theta_filtered`, `n_filter_dyads`, `filter_dyads`, and, only
#'   when both actor columns are present, `n_countries`, `ctry_a`,
#'   `ctry_b` and `w_send` (consumed only by the experimental
#'   `stable_gamma` variant). `Y` is an integer `D x T x A` array, in
#'   `dyad_id` order (dyads sorted in C-locale byte order) and time order.
#'   Attributes:
#'
#'   * `dyad_ids`: one row per dyad-period, with `dyad_id`, `time_index`,
#'     `dyad`, `dyad2` (the undirected dyad key), `actor1`, `actor2`
#'     (`dyad2`/`actor1`/`actor2` are `NA` without actor columns), the time
#'     column, and any `time_labels` columns; see [extract_theta()];
#'   * `event_classes`: the action-column names, in `Y`'s order;
#'   * `country_codes`: the countries in `ctry_a`/`ctry_b` index order
#'     (`NULL` without actor columns); see [extract_gamma()];
#'   * `grouping_var`, `time_col`, `time_origin` (from
#'     `attr(data, "time_origin")`, else `NA`), `directed`;
#'   * `w_send_imputed`: `TRUE` if `w_send` was absent and taken as `0.5`
#'     for undirected dyads, else `FALSE`.
#' @examples
#' # A non-event source: arbitrary dyad labels, no actor columns, a weekly
#' # period index with gaps, and three action types.
#' dp <- tibble::tibble(
#'   dyad = c("a-b", "a-b", "a-b", "a-c", "a-c", "b-c"),
#'   week = c(1L, 2L, 6L, 1L, 3L, 4L),
#'   talks = c(3, 0, 2, 1, 4, 2),
#'   threats = c(0, 1, 0, 2, 0, 1),
#'   aid = c(1, 1, 0, 0, 2, 0)
#' )
#' stan_data <- prepare_stan_data(
#'   dp,
#'   time = "week", action_cols = c(talks, threats, aid),
#'   reference_category = "aid"
#' )
#' dim(stan_data$Y)
#'
#' # From event data
#' events <- tibble::tibble(
#'   Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN"),
#'   Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA"),
#'   SQLDATE = c(20180105L, 20180610L, 20190302L, 20190401L, 20200101L),
#'   PentaClass = c(0, 1, 0, 1, 3)
#' )
#' agg <- aggregate_events(events, grouping_var = "PentaClass")
#' stan_data <- prepare_stan_data(
#'   agg,
#'   action_cols = c("0", "1", "3"),
#'   reference_category = "1" # verbal cooperation
#' )
#' @export
prepare_stan_data <- function(
  data,
  time = NULL,
  action_cols = NULL,
  dyad = "dyad",
  actor1 = "actor1",
  actor2 = "actor2",
  w_send = "w_send",
  directed = attr(data, "directed") %||% TRUE,
  reference_category = NULL,
  min_n_events = 1,
  time_range = NULL,
  chunk_size = 1,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  time_labels = attr(data, "time_labels"),
  grouping_var = attr(data, "grouping_var")
) {
  fn <- "prepare_stan_data()"
  if (!rlang::is_bool(directed)) {
    stop(fn, ": `directed` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!.is_count(min_n_events, min = 0)) {
    stop(fn, ": `min_n_events` must be a single non-negative whole number.", call. = FALSE)
  }
  if (!.is_count(chunk_size, min = 1)) {
    stop(fn, ": `chunk_size` must be a single whole number, at least 1.", call. = FALSE)
  }
  for (arg in c("rho_prior_a", "rho_prior_b")) {
    val <- get(arg)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val) || val <= 0) {
      stop(fn, ": `", arg, "` must be a single positive number.", call. = FALSE)
    }
  }
  compute_log_lik <- .as_flag(compute_log_lik, "compute_log_lik", fn)
  prior_only <- .as_flag(prior_only, "prior_only", fn)
  compute_theta_filtered <- .as_flag(compute_theta_filtered, "compute_theta_filtered", fn)
  if (!is.null(filter_dyads) && compute_theta_filtered == 0L) {
    warning(
      fn, ": `filter_dyads` is ignored because `compute_theta_filtered = 0`. ",
      "Set `compute_theta_filtered = 1` to compute theta_filtered for these dyads.",
      call. = FALSE
    )
  }
  time <- .resolve_time_col(data, time)

  # 1. validate; resolves action_cols exactly once
  cols <- .validate_dyad_period(data, time, rlang::enquo(action_cols), dyad, actor1, actor2, w_send)
  has_actors <- !is.null(cols$actor1)
  dyad_vec <- if (dyad %in% names(data)) {
    as.character(data[[dyad]])
  } else {
    paste(data[[actor1]], data[[actor2]], sep = "_")
  }
  time_vec <- as.integer(data[[time]])
  w_send_vec <- if (!is.null(cols$w_send)) as.double(data[[w_send]])
  if (directed && !is.null(w_send_vec) && any(w_send_vec != 1)) {
    stop(
      fn, ": `directed = TRUE`, but `", w_send, "` is not 1 in every row (a directed ",
      "dyad's sender is always actor1). Set `directed = FALSE` for undirected data.",
      call. = FALSE
    )
  }
  counts <- purrr::map(data[cols$action_cols], as.integer)
  actor1_vec <- if (has_actors) as.character(data[[actor1]])
  actor2_vec <- if (has_actors) as.character(data[[actor2]])

  # 1b. time_range: drop rows outside it; the grid is filled over it below
  if (is.null(time_range)) {
    time_range <- range(time_vec)
  } else {
    if (!is.numeric(time_range) || length(time_range) != 2 || anyNA(time_range) ||
        !all(.is_whole_number(time_range))) {
      stop(
        fn, ": `time_range` must be two whole numbers, c(min, max), on the scale of `",
        time, "`.",
        call. = FALSE
      )
    }
    if (time_range[1] > time_range[2]) {
      stop(
        fn, ": `time_range` min (", time_range[1], ") is greater than its max (",
        time_range[2], ").",
        call. = FALSE
      )
    }
    in_range <- time_vec >= time_range[1] & time_vec <= time_range[2]
    if (!any(in_range)) {
      stop(
        fn, ": no rows fall within `time_range` (", time_range[1], " to ", time_range[2],
        "); `", time, "` spans ", min(time_vec), " to ", max(time_vec), ".",
        call. = FALSE
      )
    }
    if (!all(in_range)) {
      message(sprintf("%s: dropping %d row(s) outside `time_range`.", fn, sum(!in_range)))
      dyad_vec <- dyad_vec[in_range]
      time_vec <- time_vec[in_range]
      w_send_vec <- w_send_vec[in_range]
      counts <- purrr::map(counts, ~ .x[in_range])
      actor1_vec <- actor1_vec[in_range]
      actor2_vec <- actor2_vec[in_range]
    }
  }
  time_values <- seq.int(as.integer(time_range[1]), as.integer(time_range[2]))

  # 2. drop all-zero action columns, before the min_n_events filter
  all_zero <- purrr::map_lgl(counts, ~ all(.x == 0L))
  if (any(all_zero)) {
    message(
      fn, ": dropping action column(s) with zero events in every row: ",
      paste(names(counts)[all_zero], collapse = ", ")
    )
    counts <- counts[!all_zero]
  }
  if (length(counts) < 2) {
    stop(
      fn, ": fewer than two action columns have any events (",
      paste(names(counts), collapse = ", "), "); the model needs at least two.",
      call. = FALSE
    )
  }

  # 3. reference first, then the user's column order
  reference_category <- validate_reference_class(
    reference_category, names(counts), "reference_category",
    null_hint = "the first action column"
  )
  if (is.null(reference_category)) {
    anchor <- names(counts)[1]
    message(
      fn, ": reference_category not set; anchoring on the first action column, \"",
      anchor, "\". Higher theta means more \"", anchor, "\" relative to other ",
      "action types. Set reference_category to a cooperative class to read ",
      "higher theta as more cooperative relations."
    )
  }
  event_classes <- c(reference_category, setdiff(names(counts), reference_category))
  counts <- counts[event_classes]
  A <- length(event_classes)

  # 4. fill the dyad x time grid. Counts get positional internal names
  # (.k1, .k2, ...) so no action column can collide with the key columns.
  count_names <- paste0(".k", seq_len(A))
  long <- tibble::tibble(.dyad = dyad_vec, .time = time_vec, !!!rlang::set_names(counts, count_names))
  # C-locale (radix) sort: the same order dplyr::group_by() gives, so
  # dyad_id does not depend on the locale.
  dyads <- sort(unique(dyad_vec), method = "radix")
  grid <- tidyr::expand_grid(.dyad = dyads, .time = time_values) %>%
    dplyr::left_join(long, by = c(".dyad", ".time")) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(count_names), ~ tidyr::replace_na(.x, 0L)))
  total <- rowSums(as.matrix(grid[count_names]))

  # 5. min_n_events filter
  dyad_totals <- tapply(total, factor(grid$.dyad, levels = dyads), sum)
  keep <- dyad_totals >= min_n_events
  dyads <- dyads[keep]
  D <- length(dyads)
  if (D == 0) {
    stop(
      fn, ": no dyads have at least `min_n_events` (", min_n_events, ") total events. ",
      "Lower `min_n_events` or check the input data.",
      call. = FALSE
    )
  }
  keep_rows <- grid$.dyad %in% dyads
  grid <- grid[keep_rows, , drop = FALSE]
  total <- total[keep_rows]
  Tn <- length(time_values)

  is_obs <- matrix(as.integer(total > 0), nrow = D, byrow = TRUE)
  Y <- array(
    unlist(purrr::map(count_names, ~ matrix(grid[[.x]], nrow = D, byrow = TRUE)), use.names = FALSE),
    dim = c(D, Tn, A)
  )

  # 6. country index, only with actor columns
  dyad_info <- tibble::tibble(
    dyad = dyad_vec,
    actor1 = actor1_vec %||% NA_character_,
    actor2 = actor2_vec %||% NA_character_,
    w_send = w_send_vec %||% NA_real_
  ) %>%
    dplyr::distinct(.data$dyad, .keep_all = TRUE)
  dyad_info <- dplyr::left_join(tibble::tibble(dyad = dyads), dyad_info, by = "dyad")

  # 7. w_send (1 when directed; 0.5, flagged, when undirected and absent)
  w_send_imputed <- is.null(w_send_vec) && !directed
  if (is.null(w_send_vec)) {
    dyad_info$w_send <- rep(if (directed) 1 else 0.5, D)
  }

  country_fields <- NULL
  countries_present <- NULL
  if (has_actors) {
    # radix = C-locale byte order: this order defines what gamma's columns
    # mean, so it must not depend on locale.
    countries_present <- sort(unique(c(dyad_info$actor1, dyad_info$actor2)), method = "radix")
    n_countries <- length(countries_present)
    country_index <- stats::setNames(seq_len(n_countries), countries_present)
    ctry_a <- unname(country_index[dyad_info$actor1])
    ctry_b <- unname(country_index[dyad_info$actor2])
    dyads_per_country <- table(c(ctry_a, ctry_b))
    message(sprintf(
      "%s: n_countries = %d (dyads per country -- min %d, median %s, max %d)",
      fn, n_countries, min(dyads_per_country),
      format(stats::median(as.numeric(dyads_per_country)), nsmall = 1),
      max(dyads_per_country)
    ))
    country_fields <- list(
      n_countries = n_countries,
      ctry_a = ctry_a,
      ctry_b = ctry_b,
      w_send = dyad_info$w_send
    )
  }

  dyad_ids <- make_dyad_ids(
    dyads,
    actor1 = if (has_actors) dyad_info$actor1,
    actor2 = if (has_actors) dyad_info$actor2,
    time_values = time_values,
    time_col = time,
    time_labels = time_labels
  )

  # 8. filter_dyads, the list, and attributes
  if (compute_theta_filtered == 1L) {
    if (is.null(filter_dyads)) {
      filter_dyad_ids <- seq_len(D)
    } else {
      filter_dyad_ids <- match(filter_dyads, dyads)
      missing_dyads <- filter_dyads[is.na(filter_dyad_ids)]
      if (length(missing_dyads) > 0) {
        stop(
          fn, ": `filter_dyads` names not found (dropped by min_n_events, or never existed): ",
          paste(missing_dyads, collapse = ", "),
          call. = FALSE
        )
      }
    }
  } else {
    filter_dyad_ids <- integer(0)
  }

  stan_data <- c(
    list(
      D = D,
      T = Tn,
      A = A,
      C = as.integer(chunk_size),
      is_obs = is_obs,
      Y = Y,
      rho_prior_a = rho_prior_a,
      rho_prior_b = rho_prior_b,
      compute_log_lik = compute_log_lik,
      prior_only = prior_only,
      compute_theta_filtered = compute_theta_filtered,
      n_filter_dyads = length(filter_dyad_ids),
      filter_dyads = filter_dyad_ids
    ),
    # Only stable_gamma declares these; CmdStan ignores undeclared fields.
    country_fields
  )

  attr(stan_data, "dyad_ids") <- dyad_ids
  attr(stan_data, "event_classes") <- event_classes
  attr(stan_data, "country_codes") <- countries_present
  attr(stan_data, "grouping_var") <- grouping_var
  attr(stan_data, "time_col") <- time
  attr(stan_data, "time_origin") <- attr(data, "time_origin") %||% as.Date(NA)
  attr(stan_data, "directed") <- directed
  attr(stan_data, "w_send_imputed") <- w_send_imputed

  stan_data
}

#' Resolve [prepare_stan_data()]'s `time = NULL`
#'
#' @param data The dyad-period data frame.
#' @param time The `time` argument.
#' @return A single column name.
#' @noRd
.resolve_time_col <- function(data, time) {
  if (!is.null(time)) {
    return(time)
  }
  from_attr <- attr(data, "time_col")
  if (!is.null(from_attr)) {
    return(from_attr)
  }
  known <- intersect(c("year", "month_id", "week_id", "day_id"), names(data))
  if (length(known) == 1) {
    message("prepare_stan_data(): using `", known, "` as the time column.")
    return(known)
  }
  stop(
    "prepare_stan_data(): cannot tell which column is the time period",
    if (length(known) > 1) paste0(" (found ", paste(known, collapse = ", "), ")"),
    "; set `time`.",
    call. = FALSE
  )
}

#' Is `x` a single whole number of at least `min`?
#'
#' @param x Anything.
#' @param min The smallest allowed value.
#' @return `TRUE`/`FALSE`.
#' @noRd
.is_count <- function(x, min = 0) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) && .is_whole_number(x) && x >= min
}

#' Validate a 0/1 data flag and store it as an integer
#'
#' Accepts `0`/`1` and `FALSE`/`TRUE`.
#'
#' @param x The flag.
#' @param arg Its argument name, for the error.
#' @param fn The calling function, for the error.
#' @return `0L` or `1L`.
#' @noRd
.as_flag <- function(x, arg, fn) {
  if (length(x) != 1 || is.na(x) || !(is.logical(x) || (is.numeric(x) && x %in% c(0, 1)))) {
    stop(fn, ": `", arg, "` must be 0 or 1 (or FALSE or TRUE).", call. = FALSE)
  }
  as.integer(x)
}
