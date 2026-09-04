# Internal registry of Stan model variants shipped under inst/stan/.
#
# Exported functions (fit_dyad_ts(), fit_panel(), compile_bilatr_model())
# always resolve `.BILATR_DEFAULT_MODEL` and never expose model choice to
# callers. The `_dev` variants in R/fit_dev.R accept a `stan_model` name
# and resolve it here, for use during model development only.
#
# To register a new model variant: add an entry below (the `file` must
# exist under inst/stan/), and it becomes immediately fittable via
# fit_dyad_ts_dev()/fit_panel_dev(stan_model = "<name>"). `status` is
# `"stable"` for the single default model (see .BILATR_DEFAULT_MODEL) or
# `"experimental"` for anything else; exported functions only ever fit the
# stable one.

.bilatr_stan_models <- list(
  stable = list(
    file = "bilatr_dirmult_irt.stan",
    description = paste(
      "Consolidated model (base + optional dyad/period/action weighting);",
      "per-dyad constant phi; non-centered process_noise hierarchy;",
      "alpha[1] fixed to 1, all other alpha freely estimated."
    ),
    status = "stable"
  ),
  alphanorm = list(
    file = "bilatr_alphanorm.stan",
    description = paste(
      "Experimental: closes the affine ridge in stable's identification by",
      "hard-pinning location (mu_theta0 removed) and normalizing alpha's",
      "RMS to 1 (sum_to_zero_vector) instead of pinning alpha[1] = 1.",
      "process_noise's prior is now a ratio to sigma_theta0, not an",
      "absolute theta-unit quantity -- see header comment; not yet",
      "prior-predictive calibrated."
    ),
    status = "experimental"
  ),
  ou = list(
    file = "bilatr_ou.stan",
    description = paste(
      "Experimental: replaces stable's random-walk theta with an OU/AR(1)",
      "process with dyad-specific equilibria (mu_dyad) and a global",
      "persistence rho, giving cross-dyad ordering a restoring force.",
      "stable's identification (alpha[1] = 1, mu_intercept[1] = 0,",
      "mu_theta_bar) is left untouched so the comparison isolates the",
      "dynamics; not yet prior-predictive calibrated."
    ),
    status = "experimental"
  ),
  alphanorm_ou = list(
    file = "bilatr_alphanorm_ou.stan",
    description = paste(
      "Experimental: combines alphanorm's identification (hard location",
      "pin, RMS-1 alpha normalization) with ou's OU/AR(1) dynamics.",
      "sd_stat is relative to sigma_mu, so exp(mu_log_sd_stat) is",
      "directly the within/between-dyad SD ratio; not yet",
      "prior-predictive calibrated."
    ),
    status = "experimental"
  )
)

# Historical note: several experimental variants were retired to
# inst/stan/legacy/ (gitignored, kept for local reference only, not
# registered):
#   - the pre-0.3.0 centered "stable"/"phi_logn" models (centered
#     process_noise hierarchy, hostile-anchored alpha[A]) and the
#     transitional "stable_ncproc"/"phi_logn_ncproc" entries, in 0.3.0;
#   - "phi_logn" (phi as a function of centered log(n_dt), per-dyad-period
#     rather than per-dyad constant), in 0.3.2 -- its final source is
#     inst/stan/legacy/bilatr_phi_logn.stan, the pre-0.3.0 centered one is
#     bilatr_phi_logn_pre_0.3.0.stan.

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
