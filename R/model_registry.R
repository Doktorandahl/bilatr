# Internal registry of Stan model variants shipped under inst/stan/.
#
# Exported functions (fit_dyad_ts(), fit_panel(), compile_bilatr_model())
# always resolve `.BILATR_DEFAULT_MODEL` and never expose model choice to
# callers. The `_dev` variants in R/fit_dev.R accept a `stan_model` name
# and resolve it here, for use during model development only.
#
# To register a new model variant: add an entry below (the `file` must
# exist under inst/stan/), and it becomes immediately fittable via
# fit_dyad_ts_dev()/fit_panel_dev(stan_model = "<name>"). This registry is
# intentionally general (not hardcoded to exactly two entries) so it can
# later absorb the planned consolidation of the original pre-bilatr model
# variants without a redesign.

.bilatr_stan_models <- list(
  stable = list(
    file = "bilatr_dirmult_irt.stan",
    description = "Consolidated model (base + optional dyad/period/action weighting); per-dyad constant phi.",
    status = "stable"
  ),
  phi_logn = list(
    file = "bilatr_phi_logn.stan",
    description = "Experimental: phi modeled as a function of centered log(n_dt), per-dyad-period rather than per-dyad constant.",
    status = "experimental"
  ),
  stable_ncproc = list(
    file = "bilatr_dirmult_irt_ncproc.stan",
    description = "Experimental: stable model with the per-dyad process_noise lognormal hierarchy sampled in non-centered form (to break the process_noise funnel); all other math identical to 'stable'.",
    status = "experimental"
  ),
  phi_logn_ncproc = list(
    file = "bilatr_phi_logn_ncproc.stan",
    description = "Experimental: phi_logn model with the per-dyad process_noise lognormal hierarchy also sampled in non-centered form (on top of its existing non-centered log_phi0 hierarchy); all other math identical to 'phi_logn'.",
    status = "experimental"
  )
)

.BILATR_DEFAULT_MODEL <- "stable"

#' Resolve a registered Stan model name to its file path
#'
#' @param name A name registered in `.bilatr_stan_models`.
#' @return The path to the model's `.stan` file (works both from the
#'   package source tree under `devtools::load_all()` and from an
#'   installed package, via `system.file()`).
#' @keywords internal
.resolve_stan_model <- function(name) {
  if (!name %in% names(.bilatr_stan_models)) {
    stop(
      "Unknown stan_model '", name, "'. Registered models: ",
      paste(names(.bilatr_stan_models), collapse = ", "), ".",
      call. = FALSE
    )
  }

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
