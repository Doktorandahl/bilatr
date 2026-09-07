#' Extract theta trajectories with dyad identifiers reattached
#'
#' Pulls posterior summaries of the latent conflict trajectory `theta`
#' out of a fitted model, indexed by dyad and time period, and reattaches
#' the human-readable dyad identifiers (`dyad`, the undirected `dyad2`,
#' `year`, and `month` if applicable) via the `dyad_ids` attribute that
#' [assemble_stan_data()] attaches to its output. Works the same way for
#' single-dyad ([fit_dyad_ts()]) and panel ([fit_panel()]) fits.
#'
#' Draws are passed through [bilatr_orient()] before summarizing, so
#' `theta`'s sign always reflects the canonical orientation (higher theta
#' = better relations) regardless of which reflection-symmetry basin the
#' sampler actually landed in. This is a no-op for `stan_model`s other
#' than `"alphanorm"`/`"alphanorm_ou"`, so the default preserves existing
#' behavior for `stable`-model fits (its `alpha[1] = 1` is fixed, never
#' negative).
#'
#' @param fit A `CmdStanMCMC` fit object from [fit_dyad_ts()] or
#'   [fit_panel()].
#' @param stan_data The Stan data list used to produce `fit`, as returned
#'   by [assemble_stan_data()] (must still carry its `dyad_ids`
#'   attribute).
#' @param probs Posterior quantiles to report alongside the mean.
#' @param stan_model Name registered in `.bilatr_stan_models` identifying
#'   which model produced `fit`; see [bilatr_orient()]. Defaults to
#'   `.BILATR_DEFAULT_MODEL` (`"stable"`), matching what [fit_dyad_ts()]/
#'   [fit_panel()] always fit.
#' @return A tibble with one row per dyad-period: `dyad_id`,
#'   `time_index`, `dyad`, `dyad2`, `year` (and `month`, if applicable),
#'   the posterior `mean` of theta, and one column per requested quantile.
#' @examples
#' \dontrun{
#' theta <- extract_theta(fit, stan_data)
#' }
#' @export
extract_theta <- function(fit, stan_data, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL) {
  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be the output of assemble_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }

  draws <- fit$draws(variables = c("alpha[1]", "theta"))
  draws <- bilatr_orient(draws, stan_model = stan_model, variables = "theta")

  posterior::summarise_draws(
    draws,
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
#' scale/sign reference) for `stable`/`ou`; every other element is freely
#' estimated. See the package's identification notes in
#' `vignette("dyad_time_series")`.
#'
#' Draws are passed through [bilatr_orient()] before summarizing, so
#' `alpha`'s sign always reflects the canonical orientation regardless of
#' which reflection-symmetry basin the sampler actually landed in
#' (`stan_model = "alphanorm"`/`"alphanorm_ou"` only; a no-op otherwise).
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
extract_alpha <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL) {
  draws <- bilatr_orient(fit$draws(variables = "alpha"), stan_model = stan_model, variables = "alpha")

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
#' `mu_intercept` out of a fitted model. `mu_intercept[1]` is fixed at 0
#' (the model's softmax level-shift reference) for `stable`/`ou`.
#'
#' `mu_intercept` is unaffected by alphanorm/alphanorm_ou's reflection
#' symmetry (`alpha .* theta` is invariant under the joint negation, so
#' `mu_intercept` never needs to flip; see [bilatr_orient()]) -- `fit` is
#' passed through it anyway, for a single consistent code path across the
#' `extract_*()` functions, but it is always a no-op here.
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
extract_mu_intercept <- function(fit, event_classes = NULL, probs = c(0.05, 0.5, 0.95), stan_model = .BILATR_DEFAULT_MODEL) {
  draws <- fit$draws(variables = c("alpha[1]", "mu_intercept"))
  draws <- bilatr_orient(draws, stan_model = stan_model, variables = "mu_intercept")

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
