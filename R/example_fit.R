#' A tiny example fit shipped with bilatr
#'
#' Returns the paths of the CmdStan CSV files of a small, real fit that
#' ships with the package, so that the post-fit functions
#' ([extract_theta()], [diagnose_convergence()], [icc_curves()],
#' [check_compositional_residuals()], ...) can be tried, and their
#' examples run, without CmdStan.
#'
#' **What the fit is.** A simulated directed panel of 3 countries (`"USA"`,
#' `"CHN"`, `"RUS"`; 6 directed dyads) observed yearly over 2015-2020: 749
#' events whose CAMEO codes are more cooperative in more cooperative
#' dyad-years, recoded with [recode_cameo()] and aggregated to PentaClass
#' counts with `aggregate_events(grouping_var = "PentaClass")` (36
#' dyad-years). It was fitted with [fit_panel()] with
#' `reference_category = "1"` (verbal cooperation, so higher `theta` means
#' more cooperative relations), 2 chains, `iter_warmup = 150`,
#' `iter_sampling = 60`, `seed = 1`, and CSV output with 4 significant
#' figures. Two models are available:
#'
#' * `"stable"`, fitted with `compute_theta_filtered = 1` and
#'   `compute_log_lik = 1`, so it also has `theta_filtered`,
#'   `theta_filtered_sd` and `log_lik`;
#' * `"stable_gamma"`, the experimental model with country-level offsets
#'   `gamma` (see [extract_gamma()]).
#'
#' The script that builds it is `data-raw/make_example_fit.R` in the
#' source repository.
#'
#' **It is not a valid inference.** 60 draws per chain from a short warmup
#' on six dyads are far too few for reliable estimates or convergence
#' diagnostics, and the data are simulated. The fit exists only for
#' examples and tests.
#'
#' **The Stan data is found automatically.** Each fit's directory also
#' holds the `example_stan_data.rds` file that [fit_panel()] wrote next
#' to the CSVs, and the post-fit functions read it from there when given
#' these paths (with a message saying so). It records the model that was
#' fit, so `stan_model` never needs to be passed.
#'
#' @param stan_model Which fit: `"stable"` (default) or `"stable_gamma"`.
#' @return A character vector of the two CSV file paths (one per chain),
#'   sorted.
#' @seealso [read_bilatr_draws()] to read selected variables from the
#'   files.
#' @examples
#' csv <- bilatr_example_fit()
#' extract_theta(csv)
#'
#' basename(bilatr_example_fit("stable_gamma"))
#' @export
bilatr_example_fit <- function(stan_model = c("stable", "stable_gamma")) {
  stan_model <- rlang::arg_match(stan_model)
  dir <- system.file("extdata", "example_fit", stan_model, package = "bilatr")
  if (!nzchar(dir)) {
    stop(
      "bilatr_example_fit(): the example fit `", stan_model, "` is not installed.",
      call. = FALSE
    )
  }
  sort(list.files(dir, pattern = "\\.csv$", full.names = TRUE))
}

#' Is CmdStan available for fitting?
#'
#' `TRUE` when the 'cmdstanr' package is installed at the version bilatr
#' needs (0.8.0 or later) *and* cmdstanr can find a CmdStan installation;
#' `FALSE` otherwise. It never errors, so it can guard code that fits a
#' model, as the examples of [fit_panel()] and [compile_bilatr_model()] do.
#'
#' Compiling or fitting without CmdStan errors with a message saying what
#' is missing and how to install it.
#'
#' @return `TRUE` or `FALSE`.
#' @examples
#' bilatr_cmdstan_available()
#' @export
bilatr_cmdstan_available <- function() {
  is.null(.cmdstan_problem())
}

#' What is missing for compiling and fitting, if anything
#'
#' The logic shared by [bilatr_cmdstan_available()] and
#' [.require_cmdstan()]. Never errors.
#'
#' @return `NULL` when cmdstanr (at least [.BILATR_MIN_CMDSTANR]) and
#'   CmdStan are available; otherwise `"cmdstanr"` (cmdstanr is not
#'   installed, or too old) or `"cmdstan"` (cmdstanr cannot find CmdStan).
#' @keywords internal
.cmdstan_problem <- function() {
  if (!rlang::is_installed("cmdstanr", version = .BILATR_MIN_CMDSTANR)) {
    return("cmdstanr")
  }
  version <- tryCatch(
    suppressWarnings(suppressMessages(cmdstanr::cmdstan_version(error_on_NA = FALSE))),
    error = function(e) NULL
  )
  if (is.null(version)) "cmdstan" else NULL
}
