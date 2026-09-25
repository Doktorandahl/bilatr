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
      "Since 0.4.2b (see NEWS.md): identifies alpha[1]'s sign via an",
      "ORIENTATION FOLD -- alpha_raw is a free sum_to_zero_vector (as",
      "before 0.4.2), but alpha and theta are both multiplied by",
      "sign(alpha_raw[1]) in transformed parameters, so alpha[1] >= 0",
      "always in the REPORTED draws regardless of which raw-space basin a",
      "chain occupies (a 0.4.2 hard-positivity constraint on alpha_raw",
      "tried this first and failed on the cluster -- see NEWS.md and the",
      ".stan file's header). Also hard-pins location (no mu_theta0) and",
      "normalizes alpha's RMS to 1. positive alpha[1] means better",
      "relations, for every fit made under this program (no post-hoc",
      "orientation needed for alpha/theta; the raw alpha_raw/z_theta0/",
      "theta_raw remain sign-ambiguous and are excluded from tiered",
      "diagnostics -- see R/diagnose_convergence.R). process_noise's",
      "prior is a ratio to sigma_theta0, not an absolute theta-unit",
      "quantity -- see the .stan file's header comment."
    ),
    status = "stable"
  ),
  ou = list(
    file = "bilatr_alphanorm_ou.stan",
    description = paste(
      "Since 0.4.2b (see NEWS.md): combines stable's identification (hard",
      "location pin, RMS-1 alpha normalization, alpha[1] >= 0 via the",
      "orientation fold) with an OU/AR(1) theta process (dyad-specific",
      "equilibria `mu_dyad` and a global persistence `rho`), giving",
      "cross-dyad ordering a restoring force in place of stable's",
      "random-walk theta. sd_stat is relative to sigma_mu, so",
      "exp(mu_log_sd_stat) is directly the within/between-dyad SD ratio.",
      "Still experimental: not yet prior-predictive calibrated."
    ),
    status = "experimental"
  ),
  stable_gamma = list(
    file = "bilatr_alphanorm_gamma.stan",
    description = paste(
      "Since 0.7.0 (see NEWS.md): adds a country-level category-offset",
      "vector `gamma` on top of stable's identification (same alpha/theta/",
      "mu_intercept construction, orientation fold, and RMS-1 alpha",
      "normalization, all unchanged).",
      "eta_dt = alpha*theta_dt - mu_intercept - g_d, with g_d built from the",
      "SAME `gamma` for directed (g_d = gamma[, ctry_a[d]]) and undirected",
      "(g_d = w_send[d]*gamma[, ctry_a[d]] + (1-w_send[d])*gamma[, ctry_b[d]])",
      "dyads -- one program, one shared parameter, the difference following",
      "mechanically from what a dyad-year is in each. `gamma` is",
      "identified by three constraints enforced in transformed parameters",
      "(orthogonal to `alpha`, since theta already absorbs any component",
      "along it; orthogonal to `1` per country, softmax shift invariance;",
      "centred across countries per category, since the common part is",
      "absorbed by mu_intercept) -- see the .stan file's header for the",
      "four-step derivation and why the order matters. Averaging",
      "gamma_i/gamma_j for undirected dyads is a first-order approximation",
      "to the true two-component mixture likelihood, accurate when gamma is",
      "small (deliberate, documented trade-off -- see the .stan file's",
      "header, \"HONEST CAVEAT\"). Nests `stable` exactly as",
      "sigma_gamma -> 0 (and exactly, not just in the limit, at",
      "n_countries = 1). Motivated by check_compositional_residuals()",
      "(0.6.0), which found a dyad-level compositional residual a",
      "dyad-specific beta_d cannot identify (theta already absorbs it) but",
      "a country-level offset can (a country appears in many dyads).",
      "Needs `stan_data` assembled by assemble_stan_data() >= 0.7.0",
      "(`n_countries`/`ctry_a`/`ctry_b`/`w_send`); check_compositional_residuals()",
      "and icc_curves() are gamma-aware for this model (see NEWS.md). Still",
      "experimental: not yet fit at production scale."
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
#     0.4.0 by the promoted `alphanorm`/`alphanorm_ou` models -- see
#     inst/stan/legacy/bilatr_dirmult_irt_pre_0.4.0.stan and
#     bilatr_ou_pre_0.4.0.stan.
#   - 0.4.2: `stable`/`ou`'s free sum_to_zero_vector alpha_raw (with a
#     soft sign anchor on alpha[1]) replaced with a hard-positivity
#     constraint (real<lower=0> alpha_raw_1) intended to remove the
#     alpha/theta reflection symmetry entirely. A production cluster run
#     falsified this: the constraint removes one mode from the parameter
#     space but not the likelihood barrier between the two, so a chain
#     that would have landed in the excluded mode instead slides to the
#     constraint boundary and parks there (see NEWS.md and each .stan
#     file's header for the numbers). Reverted in 0.4.2b.
#   - 0.4.2b: `stable`/`ou`'s alpha_raw is back to a free
#     sum_to_zero_vector (as it was pre-0.4.2), and alpha[1]'s sign is
#     instead identified via an ORIENTATION FOLD: alpha and theta are
#     both multiplied by sign(alpha_raw[1]) in transformed parameters, so
#     alpha[1] >= 0 always in the REPORTED draws without excluding any
#     region of the parameter space (see NEWS.md and each current .stan
#     file's header, "IDENTIFICATION: ORIENTATION FOLD"). The pre-0.4.2
#     `stable_soft_anchor`/`ou_soft_anchor` programs (what `stable`/`ou`
#     were named 0.4.0-0.4.1) were kept registered (status = "legacy")
#     from 0.4.2 through 0.9.1 so runs made under them stayed readable.
#   - 0.10.0: `stable_soft_anchor`/`ou_soft_anchor`, the pre-0.4.0
#     `alphanorm`/`alphanorm_ou` aliases, and the entire post-hoc
#     orientation stack they needed (`bilatr_orient()`,
#     `.bilatr_flip_variables()`, `.warn_if_wrong_basin()`) were retired
#     fully -- unlike every earlier retirement above, these two are no
#     longer registered at all: [.canonical_stan_model()] now errors on
#     all four names, naming bilatr <= 0.9.1 as what's needed to read
#     fits made under them. This is safe because `stable`/`ou`'s
#     orientation fold (since 0.4.2b) never needed the soft-anchor
#     programs for anything except reading old fits, and none remain in
#     active use.

.BILATR_DEFAULT_MODEL <- "stable"

#' Names retired in 0.10.0 that once resolved to a registered model
#'
#' `alphanorm`/`alphanorm_ou` (the pre-0.4.0 names) and
#' `stable_soft_anchor`/`ou_soft_anchor` (the pre-0.4.2 soft-sign-anchor
#' programs those aliases ultimately resolved to) are no longer registered
#' at all -- every fit made under any of these four names predates 0.4.2,
#' so it needs bilatr <= 0.9.1 (the last version to still register and
#' orient them) to read, not this or any later version.
#' @keywords internal
.BILATR_RETIRED_MODEL_NAMES <- c("alphanorm", "alphanorm_ou", "stable_soft_anchor", "ou_soft_anchor")

#' Resolve a `stan_model` argument to its canonical registered name
#'
#' Every `stan_model`-accepting entry point should call this before doing
#' anything else with the value: an unrecognized name must never be
#' allowed to silently fall through as though it were a no-op.
#'
#' @param name A `stan_model` value as given by the caller: a name
#'   currently in `.bilatr_stan_models`, one of
#'   `.BILATR_RETIRED_MODEL_NAMES` (a dedicated "retired in 0.10.0" error),
#'   or anything else (a generic "unknown" error).
#' @return The canonical registered name.
#' @keywords internal
.canonical_stan_model <- function(name) {
  if (name %in% names(.bilatr_stan_models)) {
    return(name)
  }
  if (name %in% .BILATR_RETIRED_MODEL_NAMES) {
    stop(
      "stan_model '", name, "' was retired in bilatr 0.10.0 (see NEWS.md): ",
      "its sign-ambiguous soft-anchor identification was replaced by the ",
      "current models' orientation fold, and the post-hoc orientation ",
      "machinery it needed no longer exists. Install bilatr <= 0.9.1 to ",
      "read/re-derive a fit made under this name.",
      call. = FALSE
    )
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
#' @keywords internal
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
#' aware behaviour downstream (0.7.0+): [check_compositional_residuals()]
#' and [icc_curves()] both branch on this, via the registry rather than a
#' user-facing flag, so a future `gamma`-bearing variant only needs an
#' entry here, not a change at every call site. `stable_gamma` is
#' currently the only one; written as an explicit set (not, say, a
#' `has_gamma` field on `.bilatr_stan_models` entries) because only one
#' entry needs it so far and this keeps the registry's existing shape
#' unchanged for every other model.
#'
#' @param stan_model Name registered in `.bilatr_stan_models`; see
#'   [.canonical_stan_model()].
#' @return `TRUE`/`FALSE`.
#' @keywords internal
.bilatr_model_has_gamma <- function(stan_model) {
  .canonical_stan_model(stan_model) %in% c("stable_gamma")
}
