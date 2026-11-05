#' bilatr: Hierarchical Dynamic IRT Models for Latent Bilateral Relations
#'
#' bilatr estimates the latent state of bilateral relations between pairs
#' of actors (dyads) over time from event data. Each dyad's latent state
#' \eqn{\theta} follows a random walk, and the counts of its events in a
#' set of action classes (for example CAMEO-coded GDELT events grouped by
#' [recode_cameo()]) follow a Dirichlet-multinomial distribution whose
#' class probabilities depend on \eqn{\theta} through a discrimination
#' parameter per class, as in an item-response model. The model is
#' hierarchical across dyads and is fitted in 'Stan' through 'cmdstanr'.
#' The latent scale runs from conflict to cooperation: it is oriented by a
#' reference action class, so that higher \eqn{\theta} makes that class
#' more likely; with a cooperative reference class, higher values mean
#' more cooperative relations (see `reference_category` in
#' [prepare_stan_data()]).
#'
#' @section Workflow:
#' 1. Get dyadic event data: [download_gdelt()] and [read_gdelt()] for
#'    GDELT 1.0, or your own events in the format described in
#'    [bilatr_event_data].
#' 2. Group the events into action classes: [recode_cameo()] for CAMEO
#'    codes.
#' 3. Count events per dyad and period: [aggregate_events()] (the
#'    resulting table is described in [bilatr_dyad_period_data]).
#' 4. Fit the model: [fit_panel()] for many dyads, [fit_dyad_ts()] for one
#'    (both call [prepare_stan_data()] and [compile_bilatr_model()]).
#' 5. Extract the estimates: [extract_theta()], [extract_alpha()],
#'    [extract_mu_intercept()] and, for the `stable_gamma` model,
#'    [extract_gamma()].
#' 6. Check the fit: [diagnose_convergence()] for convergence by
#'    parameter tier, then [check_compositional_residuals()],
#'    [icc_curves()] and [diagnose_category_merges()].
#'
#' Fitting needs the 'cmdstanr' package and a CmdStan installation; check
#' with [bilatr_cmdstan_available()]. [bilatr_example_fit()] returns a
#' small precomputed fit, so every post-fit function can be tried without
#' them.
#'
#' @section Data input:
#' [download_gdelt()], [read_gdelt()], [gdelt_files()], [gdelt_columns()],
#' [count_gdelt_actor_types()], [recode_cameo()], [cameo_lookup],
#' [validate_bilatr_events()] ([bilatr_event_data]).
#'
#' Country codes: [gdelt_countries], [gdelt_country_name()],
#' [gdelt_country_type()], [gdelt_is_state()], [gdelt_country_iso3()],
#' [filter_gdelt_states()], [gdelt_state_spells()].
#'
#' @section Preparation and fitting:
#' [aggregate_events()], [validate_bilatr_dyad_period()]
#' ([bilatr_dyad_period_data]), [prepare_stan_data()], [fit_panel()],
#' [fit_dyad_ts()], [compile_bilatr_model()],
#' [bilatr_clear_model_cache()], [bilatr_cmdstan_available()].
#'
#' @section Extraction:
#' [extract_theta()], [extract_alpha()], [extract_mu_intercept()],
#' [extract_gamma()], [gamma_norm()], [event_class_labels()],
#' [read_bilatr_draws()], [diagnose_and_extract_bilatr()] (diagnostics and
#' extraction from CmdStan CSV files in one read).
#'
#' @section Diagnostics and checks:
#' [diagnose_convergence()], [check_compositional_residuals()],
#' [write_residual_check()], [icc_curves()], [icc_crossings()],
#' [diagnose_category_merges()], [merge_cost()], [alpha_prior_moments()].
#'
#' @section Tuning:
#' [benchmark_bilatr()], [benchmark_grid()], [insert_reduce_sum_profile()],
#' [summarize_benchmark()], [benchmark_pareto()], [plot_benchmark()].
#'
#' @section Data:
#' [theta_trajectories] (estimated bilateral relations for every country
#' pair, 1991-2025, downloaded with [download_bilatr_theta()]; see also
#' [bilatr_theta_versions()] and [bilatr_clear_data_cache()]),
#' [theta_sample] (a shipped extract), [cameo_lookup] (the CAMEO
#' recoding table), and [bilatr_example_fit()].
#'
#' @section Vignettes:
#' * `vignette("dyad_time_series", package = "bilatr")`: one dyad, from
#'   events to a fitted trajectory;
#' * `vignette("panel_model", package = "bilatr")`: a panel of dyads;
#' * `vignette("diagnostics", package = "bilatr")`: tiered convergence
#'   diagnostics;
#' * `vignette("tuning", package = "bilatr")`: threads, chunk size and
#'   chains.
#'
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom rlang %||% .data :=
"_PACKAGE"
