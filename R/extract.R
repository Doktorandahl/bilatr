#' Extract theta trajectories with dyad identifiers reattached
#'
#' Pulls posterior summaries of the latent bilateral-relations trajectory
#' `theta`
#' out of a fitted model, indexed by dyad and time period, and reattaches
#' the human-readable dyad identifiers (`dyad`, the undirected `dyad2`,
#' the actors, the time column and any time labels such as `year`/`month`)
#' via the `dyad_ids` attribute that [prepare_stan_data()] attaches to its
#' output. Works the same way for
#' single-dyad ([fit_dyad_ts()]) and panel ([fit_panel()]) fits.
#'
#' Every model folds `alpha[1]`'s sign into the reported `alpha` and
#' `theta`, so `theta` is already oriented, with no post-hoc relabeling
#' needed: higher values make `reference_category` more likely relative to the
#' other action types, and so mean more cooperative relations when the
#' reference category is cooperative. See the `reference_category`
#' argument of [prepare_stan_data()].
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths
#' (one per chain), as [diagnose_convergence()] does: `theta` is then read
#' and summarized from the CSV files in memory-bounded chunks, so the full
#' set of draws is never in memory at once. **This trades memory for
#' time**: each chunk costs about one pass over the largest chain file,
#' however few columns it keeps, so the total time grows with the number
#' of chunks. Use an in-memory `fit` when you have one, and otherwise set
#' `max_memory_mb` as large as you can afford (see "Large fits from CSV
#' files" in [diagnose_convergence()]). This mode supports only the
#' default `probs` (the 5%, 50% and 95% quantiles); pass an in-memory
#' `fit` for other quantiles.
#'
#' @section Inactive periods:
#' A dyad-period in which one of the actors does not exist (a country
#' before it enters the state system or after it leaves; see
#' [aggregate_events()], "Actor periods") has no events, so its
#' \eqn{\theta} is a random-walk extrapolation with no data, not an
#' estimate. With `active_only = TRUE`, the default, such periods (those
#' with `active = FALSE` in the `dyad_ids` attribute) are left out. Stan
#' data from before `active` existed has no such column; every period
#' then counts as active.
#'
#' @param fit A `CmdStanMCMC` fit object from [fit_dyad_ts()] or
#'   [fit_panel()], or a character vector of raw CmdStan CSV file paths
#'   (see Details).
#' @param stan_data The Stan data list used to produce `fit`, from
#'   [prepare_stan_data()]. `NULL` (default) uses the one [fit_panel()]/
#'   [fit_dyad_ts()] attached to `fit` (`attr(fit, "stan_data")`), or, when
#'   `fit` is CSV paths, the single `*_stan_data.rds` file next to them;
#'   errors if neither is found.
#'   It must carry its `dyad_ids` attribute.
#' @param probs Posterior quantiles to report alongside the mean. Only
#'   the default is supported when `fit` is CSV file paths (see Details).
#' @param active_only `TRUE` (default) drops inactive dyad-periods; `FALSE`
#'   keeps every period. See "Inactive periods".
#' @param max_memory_mb,chunk_size,n_cores Only used when
#'   `fit` is CSV file paths; the same as [diagnose_convergence()]'s
#'   arguments of the same name (including the message when
#'   `max_memory_mb` is left at its default), applied to `theta` alone.
#' @return A tibble with one row per dyad-period: `dyad_id`,
#'   `time_index`, the other `dyad_ids` columns (`dyad`, `dyad2`, `actor1`,
#'   `actor2`, the time column, `active` and the time labels), the posterior `mean`
#'   of theta, and one column per requested quantile.
#' @examples
#' # CmdStan CSVs from a completed run, never read into this R session (here,
#' # the example fit shipped with bilatr); the *_stan_data.rds that
#' # fit_panel() wrote next to the CSVs is found
#' csv_files <- bilatr_example_fit()
#' theta <- extract_theta(csv_files, max_memory_mb = 4096)
#' theta
#'
#' # one dyad's trajectory
#' theta[theta$dyad == "USA_CHN", c("year", "mean", "5%", "95%")]
#' @export
extract_theta <- function(
  fit, stan_data = NULL, probs = c(0.05, 0.5, 0.95), active_only = TRUE,
  max_memory_mb = 8192, chunk_size = NULL, n_cores = 1L
) {
  max_memory_mb_missing <- missing(max_memory_mb)
  if (!rlang::is_bool(active_only)) {
    stop("extract_theta(): `active_only` must be TRUE or FALSE.", call. = FALSE)
  }
  stan_data <- .stan_data_from_fit(fit, stan_data, "extract_theta()")

  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be a Stan data list from prepare_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }
  dyad_ids <- .dyad_ids_with_active(dyad_ids)

  if (is.character(fit)) {
    if (!identical(probs, c(0.05, 0.5, 0.95))) {
      stop(
        "extract_theta()'s CSV-file-path mode only supports the default ",
        "`probs` (0.05, 0.5, 0.95). Pass an in-memory `fit` object if you ",
        "need other quantiles.",
        call. = FALSE
      )
    }

    csv_files <- fit
    prepared <- .prepare_fast_csv_read(csv_files)

    var_tiers <- .classify_bilatr_tier(prepared$variables)
    theta_vars <- var_tiers$variable[var_tiers$tier == 3L & startsWith(var_tiers$variable, "theta[")]

    chunk_size_used <- .resolve_chunk_size_and_report(
      length(theta_vars), prepared, max_memory_mb, chunk_size, n_cores,
      max_memory_mb_missing
    )
    theta_summ <- .chunked_summarise_csv(prepared, theta_vars, chunk_size_used, n_cores) %>%
      dplyr::select(variable, mean, `5%` = q5, `50%` = median, `95%` = q95)
  } else {
    draws <- fit$draws(variables = "theta")
    theta_summ <- posterior::summarise_draws(
      draws,
      mean = mean,
      ~ stats::quantile(.x, probs = probs)
    )
  }

  theta_summ %>%
    dplyr::mutate(variable = stringr::str_remove_all(variable, "theta\\[|\\]")) %>%
    tidyr::separate(variable, into = c("dyad_id", "time_index"), sep = ",", convert = TRUE) %>%
    dplyr::left_join(dyad_ids, by = c("dyad_id", "time_index")) %>%
    .keep_active(active_only)
}

#' Add `active = TRUE` to a `dyad_ids` table that predates the column
#'
#' Stan data prepared before `active` existed (e.g. saved
#' `*_stan_data.rds` files) has none; every period then counts as active.
#'
#' @param dyad_ids The `dyad_ids` attribute.
#' @return `dyad_ids`, with a logical `active` column.
#' @noRd
.dyad_ids_with_active <- function(dyad_ids) {
  if (!"active" %in% names(dyad_ids)) {
    dyad_ids$active <- rep(TRUE, nrow(dyad_ids))
  }
  dyad_ids
}

#' Drop inactive dyad-periods from a table joined to `dyad_ids`
#'
#' @param x A table with an `active` column.
#' @param active_only Whether to drop rows with `active = FALSE`.
#' @return `x`, filtered when `active_only`.
#' @noRd
.keep_active <- function(x, active_only) {
  if (!active_only) {
    return(x)
  }
  x[x$active %in% TRUE, , drop = FALSE]
}

#' Read named variables from either an in-memory fit or raw CmdStan CSVs
#'
#' Shared one-line dispatch behind [extract_alpha()]/
#' [extract_mu_intercept()]'s `fit`-as-CSV-paths support: unlike
#' [extract_theta()]/[diagnose_convergence()]'s CSV-path branches, these
#' two already know exactly which (small) set of variable names they
#' want, so there is no tier enumeration or chunking to do -- just one
#' plain, unchunked read across all chain files, via
#' [.prepare_fast_csv_read]/[.fast_read_post_warmup_draws]. Not free in
#' absolute terms, though: `data.table::fread()` touches each file's
#' full byte range once per call regardless of how few columns are
#' requested (see [.fast_read_post_warmup_draws]'s docs), so this is the
#' minimum possible number of file touches (one), not an instant lookup.
#'
#' @param fit A `CmdStanMCMC`-like fit object, or a character vector of
#'   CmdStan CSV file paths.
#' @param variables Character vector of variable names to read.
#' @return A `posterior::draws_array`.
#' @noRd
.get_draws <- function(fit, variables) {
  if (is.character(fit)) {
    prepared <- .prepare_fast_csv_read(fit)
    .fast_read_post_warmup_draws(prepared, variables)
  } else {
    fit$draws(variables = variables)
  }
}

#' Extract discrimination parameters (alpha) with action-class labels
#'
#' Pulls posterior summaries of the action-type discrimination parameters
#' `alpha` out of a fitted model. `alpha` sums to exactly 0 and has
#' root-mean-square 1 in every model; `alpha[1]`, the discrimination of
#' the `reference_category` class, is always reported non-negative (the
#' models' orientation fold), so no post-hoc relabeling is needed. See
#' the identification notes in `vignette("dyad_time_series", package =
#' "bilatr")`.
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths
#' (one per chain), for when there is no fit object in memory, for
#' example after chains were run as separate jobs. `alpha` has only `A`
#' values, so the files are read in one pass, without chunking.
#'
#' @inheritParams extract_theta
#' @param event_classes Optional character vector of event-class labels,
#'   in the same order used to build the Stan data (i.e. its
#'   `"event_classes"` attribute, from [prepare_stan_data()]). If supplied,
#'   an `event_class` column is added alongside the raw action index.
#'   `NULL` (default) uses the `event_classes` of the Stan data attached to
#'   `fit` by [fit_panel()]/[fit_dyad_ts()], if any; otherwise the output
#'   is unlabeled.
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of alpha, and one column per requested quantile.
#' @examples
#' # cross-chain roll-up from saved CSVs, no in-memory fit (here, the
#' # example fit shipped with bilatr). A fit object from fit_panel() is
#' # labeled from attr(fit, "stan_data"); CSVs need `event_classes`.
#' csv_files <- bilatr_example_fit()
#' stan_data <- readRDS(file.path(dirname(csv_files[1]), "example_stan_data.rds"))
#' extract_alpha(csv_files, event_classes = attr(stan_data, "event_classes"))
#' @export
extract_alpha <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95)) {
  event_classes <- event_classes %||% .attached_event_classes(fit)
  draws <- .get_draws(fit, "alpha")

  out <- posterior::summarise_draws(
    draws,
    mean = mean,
    ~ stats::quantile(.x, probs = probs)
  ) %>%
    dplyr::mutate(
      action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])"))
    )

  if (!is.null(event_classes)) {
    out <- dplyr::mutate(out, event_class = event_classes[action_index])
  }
  out
}

#' Extract global action-type intercepts (mu_intercept) with action-class
#' labels
#'
#' Pulls posterior summaries of the global action-type intercepts
#' `mu_intercept` out of a fitted model. `mu_intercept` sums to exactly 0
#' (a `sum_to_zero_vector`, not a single fixed-to-0 element) for
#' `stable`/`ou`, so no residual location degree of freedom hides in the
#' softmax level-shift.
#'
#' `mu_intercept` never flips under the alpha/theta reflection symmetry
#' (`alpha .* theta` is invariant under the joint negation), and `stable`/
#' `ou`'s orientation fold leaves nothing else to correct here, so
#' `mu_intercept` is read alone.
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths,
#' for the same post-hoc/no-in-memory-fit case described in
#' [extract_alpha()]'s documentation (also Tier 1/2: small regardless of
#' dyad-set size, so no chunking here either).
#'
#' @inheritParams extract_alpha
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of mu_intercept, and one column per requested quantile.
#' @examples
#' extract_mu_intercept(bilatr_example_fit())
#' @export
extract_mu_intercept <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95)) {
  event_classes <- event_classes %||% .attached_event_classes(fit)
  draws <- .get_draws(fit, "mu_intercept")

  out <- posterior::summarise_draws(
    draws,
    mean = mean,
    ~ stats::quantile(.x, probs = probs)
  ) %>%
    dplyr::mutate(
      action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])"))
    )

  if (!is.null(event_classes)) {
    out <- dplyr::mutate(out, event_class = event_classes[action_index])
  }
  out
}

#' Extract country-level category offsets (`gamma`) with country/action
#' labels
#'
#' Pulls posterior summaries of the country-level category-offset matrix
#' `gamma` out of a `stable_gamma` fit, the only model with one (see
#' [compile_bilatr_model()], "Models"). `gamma` is identified by three
#' constraints built into the model: for each country it sums to zero
#' over action classes and is orthogonal to `alpha`, and for each action
#' class it sums to zero over countries. It is unaffected by the models'
#' orientation fold, so no relabeling applies to it.
#'
#' `fit` also accepts a character vector of raw CmdStan CSV file paths,
#' as [extract_alpha()] does. `gamma` has `A x n_countries` values, so
#' the files are read in one pass, without chunking.
#'
#' @param fit A `CmdStanMCMC` fit object from
#'   `fit_panel(stan_model = "stable_gamma")` or
#'   `fit_dyad_ts(stan_model = "stable_gamma")`, or a character vector of
#'   raw CmdStan CSV file paths.
#' @param stan_data The Stan data list used to produce `fit`, from
#'   [prepare_stan_data()]. `NULL` (default) uses the one [fit_panel()]/
#'   [fit_dyad_ts()] attached to `fit` (`attr(fit, "stan_data")`), or, when
#'   `fit` is CSV paths, the single `*_stan_data.rds` file next to them;
#'   errors if neither is found.
#'   It must carry its `country_codes` attribute (built only from actor
#'   columns).
#' @param probs Posterior quantiles to report alongside the mean.
#' @param event_classes Optional character vector of event-class labels,
#'   in the same order used to build `stan_data`. `NULL` (default) uses
#'   `stan_data`'s `"event_classes"` attribute; an explicit value wins.
#' @param class_labels Optional labels for the event classes: a character
#'   vector named by event class, or an unnamed one in `event_classes`
#'   order. Ignored when `event_classes` is `NULL` (nothing to key labels
#'   by). `NULL` (default) uses `event_class_labels(stan_data)` (see
#'   [event_class_labels()]).
#' @return A tibble (subclassed `bilatr_gamma`, with [autoplot()] and
#'   `plot()` methods) with one row per (country, action type):
#'   `country_code`, `country_index`, `action_index` (and
#'   `event_class`/`class_label` if supplied), the posterior `mean` of
#'   `gamma`, and one column per requested quantile.
#' @examples
#' gamma <- extract_gamma(bilatr_example_fit("stable_gamma"))
#' gamma
#' gamma_norm(gamma)
#' @export
extract_gamma <- function(
  fit, stan_data = NULL, probs = c(0.05, 0.5, 0.95),
  event_classes = NULL,
  class_labels = NULL
) {
  stan_data <- .stan_data_from_fit(fit, stan_data, "extract_gamma()")
  event_classes <- event_classes %||% attr(stan_data, "event_classes")
  country_codes <- attr(stan_data, "country_codes")
  if (is.null(country_codes)) {
    stop(
      "`stan_data` has no 'country_codes' attribute. prepare_stan_data() ",
      "builds it only from actor columns (`actor1`/`actor2`), which ",
      "stable_gamma needs.",
      call. = FALSE
    )
  }

  draws <- .get_draws(fit, "gamma")
  out <- posterior::summarise_draws(
    draws,
    mean = mean,
    ~ stats::quantile(.x, probs = probs)
  ) %>%
    dplyr::mutate(
      idx = stringr::str_match(variable, "\\[(\\d+),(\\d+)\\]"),
      action_index = as.integer(idx[, 2]),
      country_index = as.integer(idx[, 3]),
      country_code = country_codes[country_index]
    ) %>%
    dplyr::select(-idx, -variable)

  if (!is.null(event_classes)) {
    out <- dplyr::mutate(out, event_class = event_classes[action_index])
    labels <- .resolve_class_labels(class_labels, event_classes, stan_data = stan_data)
    out <- dplyr::mutate(out, class_label = labels[action_index])
  }

  meta_cols <- c("country_code", "country_index", "action_index", intersect(c("event_class", "class_label"), names(out)))
  out <- dplyr::relocate(out, dplyr::all_of(meta_cols))

  class(out) <- c("bilatr_gamma", class(out))
  out
}

#' Per-country RMS of `gamma`, the "how idiosyncratic is this country's
#' repertoire" summary
#'
#' `sqrt(mean(gamma_c^2))` across action categories, for each country --
#' the natural ranking key for which countries carry the largest
#' country-level offset (see [extract_gamma()]). Computed from
#' `x$mean` (posterior means), not per-draw: a summary of a summary,
#' adequate for ranking/ordering (e.g. [autoplot.bilatr_gamma()]'s
#' country axis) but not itself a posterior quantity -- no interval is
#' reported.
#'
#' @param x An [extract_gamma()] tibble.
#' @return A tibble, one row per country, `country_code`, `country_index`,
#'   `gamma_norm`, sorted descending by `gamma_norm`.
#' @examples
#' gamma <- extract_gamma(bilatr_example_fit("stable_gamma"))
#' gamma_norm(gamma)
#' @export
gamma_norm <- function(x) {
  stopifnot(inherits(x, "bilatr_gamma"))
  x %>%
    dplyr::group_by(.data$country_code, .data$country_index) %>%
    dplyr::summarise(gamma_norm = sqrt(mean(.data$mean^2)), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$gamma_norm))
}

.plot_gamma <- function(x) {
  fill_col <- if ("50%" %in% names(x)) "50%" else "mean"
  label_col <- if ("class_label" %in% names(x)) "class_label" else if ("event_class" %in% names(x)) "event_class" else "action_index"

  order <- gamma_norm(x)
  d <- x %>%
    dplyr::mutate(
      country_code = factor(.data$country_code, levels = order$country_code),
      .label = factor(.data[[label_col]], levels = unique(.data[[label_col]]))
    )

  ggplot2::ggplot(d, ggplot2::aes(x = .data$country_code, y = .data$.label, fill = .data[[fill_col]])) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(name = "gamma", low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0) +
    ggplot2::labs(
      x = "country (ordered by gamma_norm())", y = NULL,
      title = "Country-level category offsets (gamma)",
      subtitle = if (fill_col == "50%") "fill: posterior median" else "fill: posterior mean"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5))
}

#' Plot `extract_gamma()` output
#'
#' Countries on the x-axis (ordered by [gamma_norm()], largest first, so
#' the structure is visible rather than alphabetical), action categories
#' on the y-axis, posterior median as fill (falling back to the mean if
#' `x` was built with quantiles that don't include the median).
#'
#' @param object,x A `bilatr_gamma` tibble, as returned by
#'   [extract_gamma()].
#' @param ... Ignored; present for S3 consistency.
#' @return A `ggplot` object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   gamma <- extract_gamma(bilatr_example_fit("stable_gamma"))
#'   ggplot2::autoplot(gamma)
#'   plot(gamma)
#' }
#' @exportS3Method ggplot2::autoplot
autoplot.bilatr_gamma <- function(object, ...) {
  rlang::check_installed("ggplot2", "for autoplot.bilatr_gamma()")
  .plot_gamma(object)
}

#' @rdname autoplot.bilatr_gamma
#' @export
plot.bilatr_gamma <- function(x, ...) {
  rlang::check_installed("ggplot2", "for plot.bilatr_gamma()")
  print(ggplot2::autoplot(x))
  invisible(x)
}
