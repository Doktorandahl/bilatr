#' Extract theta trajectories with dyad identifiers reattached
#'
#' Pulls posterior summaries of the latent conflict trajectory `theta`
#' out of a fitted model, indexed by dyad and time period, and reattaches
#' the human-readable dyad identifiers (`dyad`, the undirected `dyad2`,
#' `year`, and `month` if applicable) via the `dyad_ids` attribute that
#' [assemble_stan_data()] attaches to its output. Works the same way for
#' single-dyad ([fit_dyad_ts()]) and panel ([fit_panel()]) fits.
#'
#' @param fit A `CmdStanMCMC` fit object from [fit_dyad_ts()] or
#'   [fit_panel()].
#' @param stan_data The Stan data list used to produce `fit`, as returned
#'   by [assemble_stan_data()] (must still carry its `dyad_ids`
#'   attribute).
#' @param probs Posterior quantiles to report alongside the mean.
#' @return A tibble with one row per dyad-period: `dyad_id`,
#'   `time_index`, `dyad`, `dyad2`, `year` (and `month`, if applicable),
#'   the posterior `mean` of theta, and one column per requested quantile.
#' @examples
#' \dontrun{
#' theta <- extract_theta(fit, stan_data)
#' }
#' @export
extract_theta <- function(fit, stan_data, probs = c(0.05, 0.5, 0.95)) {
  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be the output of assemble_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }

  fit$summary(
    variables = "theta",
    mean = mean,
    ~ stats::quantile(.x, probs = probs)
  ) %>%
    dplyr::mutate(variable = stringr::str_remove_all(variable, "theta\\[|\\]")) %>%
    tidyr::separate(variable, into = c("dyad_id", "time_index"), sep = ",", convert = TRUE) %>%
    dplyr::left_join(dyad_ids, by = c("dyad_id", "time_index"))
}

#' Extract discrimination parameters (alpha) with action-class labels
#'
#' Pulls posterior summaries of the action-type discrimination parameters
#' `alpha` out of a fitted model. `alpha[1]` is fixed at 1 (the model's
#' scale reference) and the last element is the negative, known-hostile
#' anchor; see the package's identification notes in
#' `vignette("dyad_time_series")`.
#'
#' @inheritParams extract_theta
#' @param event_classes Optional character vector of event-class labels,
#'   in the same order used to build `stan_data` (i.e. the
#'   `"event_classes"` attribute attached by [assemble_stan_data()]). If
#'   supplied, an `event_class` column is added alongside the raw action
#'   index.
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of alpha, and one column per requested quantile.
#' @examples
#' \dontrun{
#' alpha <- extract_alpha(fit, event_classes = attr(stan_data, "event_classes"))
#' }
#' @export
extract_alpha <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95)) {
  out <- fit$summary(
    variables = "alpha",
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
#' `mu_intercept` out of a fitted model. `mu_intercept[1]` is fixed at 0
#' (the model's softmax level-shift reference).
#'
#' @inheritParams extract_alpha
#' @return A tibble with one row per action type: `action_index`
#'   (and `event_class` if `event_classes` is supplied), the posterior
#'   `mean` of mu_intercept, and one column per requested quantile.
#' @examples
#' \dontrun{
#' mu_intercept <- extract_mu_intercept(fit, event_classes = attr(stan_data, "event_classes"))
#' }
#' @export
extract_mu_intercept <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95)) {
  out <- fit$summary(
    variables = "mu_intercept",
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
