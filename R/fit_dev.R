# Internal development/debugging entry points for fitting bilatr models.
#
# These are NOT part of the public API (no @export) and their signature
# — in particular the `stan_model` argument — may change without notice.
# Use fit_dyad_ts()/fit_panel() for anything outside active model
# development; those always fit the stable model and are unaffected by
# whatever gets registered here.
#
# New Stan model variants are registered in R/model_registry.R
# (.bilatr_stan_models); once registered there, they're immediately
# fittable by name via fit_dyad_ts_dev()/fit_panel_dev().

#' Development entry point: fit a single dyad with a selectable Stan model
#'
#' Identical to [fit_dyad_ts()] except for the `stan_model` argument.
#' Not exported; for use during model development only.
#'
#' @inheritParams fit_dyad_ts
#' @param stan_model Name of a model registered in `.bilatr_stan_models`
#'   (see `R/model_registry.R`): `"stable"` (also the default used by the
#'   exported [fit_dyad_ts()]) or one of the experimental variants
#'   (`"alphanorm"`, `"ou"`, `"alphanorm_ou"`). Unrecognized names error
#'   immediately, before any compilation is attempted, listing the
#'   currently registered options.
#' @return A `CmdStanMCMC` fit object.
#' @keywords internal
fit_dyad_ts_dev <- function(
  stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 2,
  output_dir = NULL,
  stan_model = .BILATR_DEFAULT_MODEL,
  ...
) {
  .resolve_stan_model(stan_model)
  if (stan_data$D != 1) {
    stop(
      "fit_dyad_ts_dev() expects a single-dyad stan_data (D == 1), got D = ",
      stan_data$D,
      ". ",
      "Use fit_panel_dev() for multiple dyads.",
      call. = FALSE
    )
  }
  fit_bilatr(
    stan_data,
    chains,
    parallel_chains,
    threads_per_chain,
    iter_warmup,
    iter_sampling,
    seed,
    opt_level,
    output_dir,
    stan_model = stan_model,
    ...
  )
}

#' Development entry point: fit a panel with a selectable Stan model
#'
#' Identical to [fit_panel()] except for the `stan_model` argument. Not
#' exported; for use during model development only.
#'
#' @inheritParams fit_panel
#' @param stan_model Name of a model registered in `.bilatr_stan_models`
#'   (see `R/model_registry.R`): `"stable"` (also the default used by the
#'   exported [fit_panel()]) or one of the experimental variants
#'   (`"alphanorm"`, `"ou"`, `"alphanorm_ou"`). Unrecognized names error
#'   immediately, before any compilation is attempted, listing the
#'   currently registered options.
#' @return A `CmdStanMCMC` fit object.
#' @keywords internal
fit_panel_dev <- function(
  stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 16,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 3,
  output_dir = NULL,
  stan_model = .BILATR_DEFAULT_MODEL,
  ...
) {
  .resolve_stan_model(stan_model)
  if (stan_data$D < 2) {
    stop(
      "fit_panel_dev() expects multiple dyads (D >= 2), got D = ",
      stan_data$D,
      ". ",
      "Use fit_dyad_ts_dev() for a single dyad.",
      call. = FALSE
    )
  }
  fit_bilatr(
    stan_data,
    chains,
    parallel_chains,
    threads_per_chain,
    iter_warmup,
    iter_sampling,
    seed,
    opt_level,
    output_dir,
    stan_model = stan_model,
    ...
  )
}
