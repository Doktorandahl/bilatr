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
    file = "bilatr_alphanorm.stan",
    description = paste(
      "Promoted from the experimental `alphanorm` variant in 0.4.0 (see",
      "NEWS.md): closes the affine ridge in the previous stable model's",
      "identification by hard-pinning location (no mu_theta0) and",
      "normalizing alpha's RMS to 1 (sum_to_zero_vector) instead of",
      "pinning alpha[1] = 1. This leaves an exact alpha/theta reflection",
      "symmetry (no fixed alpha element), broken by a soft sign anchor on",
      "alpha[1] (anchor_scale data field, default 0.1) oriented so",
      "positive alpha[1] means better relations. process_noise's prior is",
      "a ratio to sigma_theta0, not an absolute theta-unit quantity -- see",
      "the .stan file's header comment."
    ),
    status = "stable"
  ),
  ou = list(
    file = "bilatr_alphanorm_ou.stan",
    description = paste(
      "Promoted from the experimental `alphanorm_ou` variant in 0.4.0 (see",
      "NEWS.md): combines stable's identification (hard location pin,",
      "RMS-1 alpha normalization, soft alpha[1] sign anchor via",
      "anchor_scale) with an OU/AR(1) theta process (dyad-specific",
      "equilibria `mu_dyad` and a global persistence `rho`), giving",
      "cross-dyad ordering a restoring force in place of stable's",
      "random-walk theta. sd_stat is relative to sigma_mu, so",
      "exp(mu_log_sd_stat) is directly the within/between-dyad SD ratio.",
      "Still experimental: not yet prior-predictive calibrated."
    ),
    status = "experimental"
  )
)

# Historical note: several models were retired to inst/stan/legacy/
# (gitignored, kept for local reference only, not registered):
#   - the pre-0.3.0 centered "stable"/"phi_logn" models (centered
#     process_noise hierarchy, hostile-anchored alpha[A]) and the
#     transitional "stable_ncproc"/"phi_logn_ncproc" entries, in 0.3.0;
#   - "phi_logn" (phi as a function of centered log(n_dt), per-dyad-period
#     rather than per-dyad constant), in 0.3.2 -- its final source is
#     inst/stan/legacy/bilatr_phi_logn.stan, the pre-0.3.0 centered one is
#     bilatr_phi_logn_pre_0.3.0.stan;
#   - the original "stable" (the consolidated Dirichlet-multinomial model,
#     `alpha[1]` hard-fixed to 1, `mu_theta0`-anchored theta0 location)
#     and the original "ou" (its OU/AR(1) theta variant), both replaced in
#     0.4.0 by the promoted `alphanorm`/`alphanorm_ou` models above -- see
#     inst/stan/legacy/bilatr_dirmult_irt_pre_0.4.0.stan and
#     bilatr_ou_pre_0.4.0.stan.

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
