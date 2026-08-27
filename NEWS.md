# bilatr 0.2.0

## New features

* Added `diagnose_convergence()` and its `print.bilatr_diagnostics()`
  method (`R/diagnose_convergence.R`), a post-hoc MCMC diagnostic triage
  helper for the panel model. It runs `posterior::summarise_draws()` and
  splits the result into three tiers -- global/shared parameters, per-dyad
  hierarchical parameters, and per-dyad-period latent states -- so that
  Rhat/ESS problems in shared parameters are always surfaced, while the
  wide, poorly-identified posteriors expected for sparse dyads' own
  `theta` trajectories are aggregated rather than flooding the report.
  Per-dyad tiers are screened against a supplied `n_dt` (per-dyad
  observation count) table to flag dyads whose diagnostics are worse than
  their sparsity alone would predict. A `tiers` argument (default `1:3`)
  restricts computation to a subset of tiers -- e.g. `tiers = 1` checks
  only the global/shared parameters, skips computing Rhat/ESS for the
  (typically far more numerous) Tier 2/3 quantities entirely via
  `posterior::subset_draws()`, and needs no `n_dt` argument. See
  `vignette("diagnostics")`.

## Internal

* Removed the internal `relevant_actors()` helper; its value (`c("GOV",
  "MIL", "SPY")`) is now the default for a new `relevant_actors`
  argument on [extract_all_relevant_gdelt()], matching how
  [ingest_icews()] already exposes its own actor filter via
  `relevant_sectors`. Backward compatible (existing calls with no
  argument are unaffected); callers can now override which GDELT
  actor-type codes count as "relevant" without editing package code.
* Renamed all kebab-case source, test, and vignette files to snake_case
  (e.g. `R/model-registry.R` -> `R/model_registry.R`,
  `vignettes/panel-model.Rmd` -> `vignettes/panel_model.Rmd`) for
  consistent file naming across the package. No functional changes;
  `vignette()` cross-references were updated to match the new vignette
  file names.
* Declared `stats` and `utils` in `Imports` (both already used via `::`
  elsewhere in the package, e.g. `stats::quantile()` in
  `extract_theta()`, but not previously declared).

## Development notes (internal, not part of the public API)

* Added an experimental Stan model variant, `inst/stan/bilatr_phi_logn.stan`,
  in which the Dirichlet-multinomial concentration `phi` is modeled as a
  per-dyad-period quantity (`log(phi[d,t]) = log_phi0[d] + beta_logn *
  centered_log_n[d,t]`) rather than a per-dyad constant. Shares the stable
  model's data interface exactly (no changes to `assemble_stan_data()` were
  needed) and nests the stable model's dispersion structure as
  `beta_logn -> 0`.
* Added an internal Stan model registry (`R/model_registry.R`,
  `.bilatr_stan_models`) mapping short names to `inst/stan/` files, plus a
  shared, cached compilation path (`.compile_stan_model()`) used by both
  the exported fitters and the new internal dev entry points below.
  `fit_dyad_ts()` and `fit_panel()` are unaffected: their signatures and
  behavior are unchanged, and they always resolve the registry's
  `"stable"` entry internally.
* Added internal (non-exported) development entry points
  `fit_dyad_ts_dev()`/`fit_panel_dev()` in `R/fit_dev.R`, identical to
  `fit_dyad_ts()`/`fit_panel()` but with a `stan_model` argument for
  selecting a registered model during active development. Not part of the
  public API.
