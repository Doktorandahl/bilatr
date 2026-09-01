
# bilatr 0.3.0

## New features

* `recode_cameo()` now also attaches `EventRootCode2`, `BilatrClass`, and
  `BilatrClassName` columns (alongside the existing
  `QuadClass`/`PentaClass`/`PentaClass_modified`). `EventRootCode2` is a
  coarser regrouping of the CAMEO root codes (root 04 "Consult" split
  three ways; roots 09/14/15/18/20 folded into related roots), and
  `BilatrClass` is an 11-level (0-10) action-class scheme intended as the
  model's default `grouping_var`. New internal helpers
  `assign_eventrootcode2()`, `assign_bilatr_class()`, and
  `bilatr_class_name()` in `R/cameo_recode.R` encode the mapping (derived
  from `original_code/cameo_df.csv`), and the `cameo_lookup` package data
  object gained the three columns.

## Breaking changes

* **Removed the `reference_hostile` argument** from
  `assemble_stan_data()` and `grouped_events_to_dyad_period()`, and the
  `reference_hostile` parameter from the internal `order_event_classes()`.
  The Stan models no longer anchor the known-hostile action class:
  `alpha[A]` is no longer fixed to `-alpha_hostile` (with
  `alpha_hostile > 0`) but is now freely estimated like every other
  `alpha[2:A]`. `alpha[1] = 1` remains the sole discrimination-scale/sign
  reference (with `mu_intercept[1] = 0` for location), which is still
  sufficient for identification. Existing calls that passed
  `reference_hostile = ...` must drop that argument; the resulting class
  ordering changes only in that the previously-last "hostile" class now
  sorts alphabetically with the rest.
* **Promoted the non-centered-`process_noise` models to the primary
  models and retired the previous ones.** `inst/stan/bilatr_dirmult_irt.stan`
  (registry `"stable"`) and `inst/stan/bilatr_phi_logn.stan` (registry
  `"phi_logn"`) now sample the per-dyad `process_noise` lognormal
  hierarchy in non-centered form (`log_process_noise_raw ~ std_normal()`,
  with `process_noise` built in transformed parameters), which broke the
  `process_noise` funnel that was driving non-convergence. The
  transitional `"stable_ncproc"` / `"phi_logn_ncproc"` registry entries
  (added in 0.2.2) are removed. The pre-0.3.0 centered sources are kept
  under `inst/stan/legacy/` (gitignored) for reference only and are not
  registered or shipped.

## Internal

* `bilatr_init_fn()` now initializes `log_process_noise_raw` (not
  `process_noise`) and an `alpha_raw` of length `A - 1` (not `A - 2`),
  and no longer sets `alpha_hostile`, to match the new parameterization.
* Dropped `"alpha_hostile"` from `diagnose_convergence()`'s Tier 1
  parameter-name list (the parameter no longer exists).

# bilatr 0.2.2

* Added two new experimental models

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
