# Registry of the Stan models shipped under inst/stan/.
#
# fit_dyad_ts(), fit_panel() and compile_bilatr_model() accept a
# `stan_model` name and resolve it here; an unregistered name errors in
# .canonical_stan_model() rather than silently falling through.
#
# To register a new model: add an entry below (the `file` must exist
# under inst/stan/), and it becomes fittable via
# fit_dyad_ts()/fit_panel(stan_model = "<name>"). `status` is `"stable"`
# for the default model (see .BILATR_DEFAULT_MODEL) or `"experimental"`;
# fitting an experimental model emits a one-time-per-session message (see
# .bilatr_warn_if_experimental() below). `description` is shown in
# ?compile_bilatr_model's Models section (see .bilatr_models_rd()).

.bilatr_stan_models <- list(
  stable = list(
    file = "bilatr_alphanorm.stan",
    description = paste(
      "The dynamic IRT model. Each dyad's latent state `theta` follows a",
      "random walk; its event counts follow a Dirichlet-multinomial",
      "whose class probabilities depend on `theta` through the",
      "discrimination vector `alpha` and the intercepts `mu_intercept`.",
      "Process noise, the dispersion `phi` and the initial states are",
      "partially pooled across dyads."
    ),
    status = "stable"
  ),
  ou = list(
    file = "bilatr_alphanorm_ou.stan",
    description = paste(
      "As `\"stable\"`, but `theta` follows a mean-reverting",
      "Ornstein-Uhlenbeck (AR(1)) process around a dyad-specific level",
      "`mu_dyad`, with a shared persistence `rho` (prior",
      "`beta(rho_prior_a, rho_prior_b)`; see [prepare_stan_data()]).",
      "`within_between_ratio` is the ratio of within-dyad to between-dyad",
      "standard deviation."
    ),
    status = "experimental"
  ),
  stable_gamma = list(
    file = "bilatr_alphanorm_gamma.stan",
    description = paste(
      "As `\"stable\"`, plus a country-level offset `gamma` on each action",
      "class, which absorbs a country's habitual repertoire, so `theta` is",
      "measured net of each country's baseline mix of event types. Directed",
      "and undirected data estimate the same `gamma`; for an undirected",
      "pair the two countries' offsets are averaged by event share, a",
      "first-order approximation. Needs actor columns in the data; see",
      "[extract_gamma()]."
    ),
    status = "experimental"
  )
)

#' The Models section of ?compile_bilatr_model, from the registry
#'
#' Used by inline roxygen code, so the help page and the registry share
#' one description per model.
#' @return A single markdown string, one bullet per registered model.
#' @noRd
.bilatr_models_rd <- function() {
  label <- function(name, status) {
    paste0(
      "`\"", name, "\"`",
      if (identical(name, .BILATR_DEFAULT_MODEL)) " (the default)" else if (identical(status, "experimental")) " (experimental)"
    )
  }
  items <- vapply(names(.bilatr_stan_models), function(name) {
    m <- .bilatr_stan_models[[name]]
    paste0("* ", label(name, m$status), ": ", m$description)
  }, character(1))
  paste(items, collapse = "\n")
}

.BILATR_DEFAULT_MODEL <- "stable"

#' Resolve a `stan_model` argument to its canonical registered name
#'
#' Every `stan_model`-accepting entry point should call this before doing
#' anything else with the value: an unrecognized name must never be
#' allowed to silently fall through as though it were a no-op.
#'
#' @param name A `stan_model` value as given by the caller.
#' @return The canonical registered name; anything else is an error.
#' @noRd
.canonical_stan_model <- function(name) {
  if (!rlang::is_string(name)) {
    stop(
      "`stan_model` must be a single string, one of: ",
      paste(names(.bilatr_stan_models), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (name %in% names(.bilatr_stan_models)) {
    return(name)
  }
  stop(
    "Unknown stan_model '", name, "'. Registered models: ",
    paste(names(.bilatr_stan_models), collapse = ", "), ".",
    call. = FALSE
  )
}

#' Resolve a registered Stan model name to its file path
#'
#' @param name A name registered in `.bilatr_stan_models` (see
#'   [.canonical_stan_model()]).
#' @return The path to the model's `.stan` file (works both from the
#'   package source tree under `devtools::load_all()` and from an
#'   installed package, via `system.file()`).
#' @noRd
.resolve_stan_model <- function(name) {
  name <- .canonical_stan_model(name)

  stan_file <- .bilatr_stan_models[[name]]$file
  path <- system.file("stan", stan_file, package = "bilatr")
  if (!nzchar(path)) {
    stop(
      "Stan file '", stan_file, "' registered for model '", name, "' ",
      "was not found under inst/stan/.",
      call. = FALSE
    )
  }
  path
}

#' Does a registered Stan model declare a country-level `gamma` offset?
#'
#' The single source of truth for whether a `stan_model` needs `gamma`-
#' aware behavior downstream: [check_compositional_residuals()] and
#' [icc_curves()] both branch on this, via the registry rather than a
#' user-facing flag, so a future `gamma`-bearing model only needs an
#' entry here, not a change at every call site. `stable_gamma` is
#' currently the only one.
#'
#' @param stan_model Name registered in `.bilatr_stan_models`; see
#'   [.canonical_stan_model()].
#' @return `TRUE`/`FALSE`.
#' @noRd
.bilatr_model_has_gamma <- function(stan_model) {
  .canonical_stan_model(stan_model) %in% c("stable_gamma")
}

#' Session-lifetime guard for the once-per-model experimental-model message
#' @noRd
.bilatr_experimental_messaged <- new.env(parent = emptyenv())

#' Reset the once-per-session experimental-model message guard
#'
#' Internal, for tests: [.bilatr_warn_if_experimental()] only messages once
#' per `stan_model` name per session, so a test exercising that message
#' needs to reset the guard first rather than relying on being the first
#' caller in the process.
#' @noRd
.reset_bilatr_experimental_messaged <- function() {
  rm(list = ls(.bilatr_experimental_messaged), envir = .bilatr_experimental_messaged)
}

#' Message once per session when fitting a registry `status = "experimental"`
#' model
#'
#' [fit_dyad_ts()]/[fit_panel()] call this after resolving `stan_model`, so a
#' caller who picks an experimental variant (`ou`, `stable_gamma`; see
#' `.bilatr_stan_models`) is told so exactly once per model name per
#' session, not on every fit.
#'
#' @param stan_model A canonical name (already resolved via
#'   [.canonical_stan_model()]).
#' @return `NULL`, invisibly.
#' @noRd
.bilatr_warn_if_experimental <- function(stan_model) {
  if (identical(.bilatr_stan_models[[stan_model]]$status, "experimental") &&
    !exists(stan_model, envir = .bilatr_experimental_messaged, inherits = FALSE)) {
    message(
      "stan_model = \"", stan_model, "\" is experimental (not yet ",
      "prior-predictive calibrated or fit at production scale; see ",
      "?compile_bilatr_model). This message is shown once per session."
    )
    assign(stan_model, TRUE, envir = .bilatr_experimental_messaged)
  }
  invisible(NULL)
}
