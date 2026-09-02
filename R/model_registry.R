# Internal registry of Stan model variants shipped under inst/stan/.
#
# Exported functions (fit_dyad_ts(), fit_panel(), compile_bilatr_model())
# always resolve `.BILATR_DEFAULT_MODEL` and never expose model choice to
# callers. The `_dev` variants in R/fit_dev.R accept a `stan_model` name
# and resolve it here, for use during model development only.
#
# To register a new model variant: add an entry below (the `file` must
# exist under inst/stan/), and it becomes immediately fittable via
# fit_dyad_ts_dev()/fit_panel_dev(stan_model = "<name>"). The registry is
# intentionally general (not hardcoded to a single entry) so it can later
# absorb new model variants without a redesign, even though only the
# `stable` model is currently registered.

.bilatr_stan_models <- list(
  stable = list(
    file = "bilatr_dirmult_irt.stan",
    description = paste(
      "Consolidated model (base + optional dyad/period/action weighting);",
      "per-dyad constant phi; non-centered process_noise hierarchy;",
      "alpha[1] fixed to 1, all other alpha freely estimated."
    ),
    status = "stable"
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
