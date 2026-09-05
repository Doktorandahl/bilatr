
# bilatr 0.3.4

## Changes

* `alphanorm`/`alphanorm_ou` gained a soft sign anchor on `alpha[1]`
  (`target += log_inv_logit(alpha[1] * inv(anchor_scale))`, new
  `anchor_scale` data field on `assemble_stan_data()`, default `0.1`),
  breaking an exact reflection symmetry those two models' identification
  (`sum_to_zero_vector` alpha with no fixed element) otherwise leaves in
  place: negating `alpha` together with `theta` and its upstream raw
  parameters leaves the likelihood and every prior unchanged, so the
  posterior has two mirror modes of equal mass and chains could land in
  either, making Rhat on `alpha`/`theta` uninterpretable. The anchor
  orients the positive mode as canonical, so higher `theta` means better
  (less hostile) relations, matching `stable`/`ou`. `ou` is unaffected
  (its hard `alpha[1] = 1` already selects a mode) and does not get this
  data field.

# bilatr 0.3.3

## New features

* Added three EXPERIMENTAL Stan model variants (`R/model_registry.R`,
  `status = "experimental"`), fittable only via
  `fit_dyad_ts_dev()`/`fit_panel_dev(stan_model = ...)`; the stable model
  is unchanged and remains the only one the exported fitters use.
  * `alphanorm` (`inst/stan/bilatr_alphanorm.stan`) targets the residual
    affine ridge in stable's identification: with `alpha[1] = 1` and
    `mu_intercept[1] = 0`, location and scale are only softly pinned, and
    an affine reparameterization of `theta`/`alpha`/`mu_intercept` leaves
    the likelihood exactly invariant -- inflating posterior intervals and
    likely contributing to treedepth pathology. `alphanorm` pins location
    hard (`mu_theta0` removed) and pins scale on the alpha side instead
    of theta's (`sum_to_zero_vector` alpha, normalized to RMS 1).
  * `ou` (`inst/stan/bilatr_ou.stan`) targets the stable model's random
    walk having no restoring force: cross-dyad SD at t=1 is small
    relative to a panel's worth of accumulated drift, so any dyad can
    traverse the entire cross-dyad range within a decade and initial
    ordering carries little information about later ordering -- this is
    why dyads that should be clearly separated in level can invert during
    a temporary detente. `ou` replaces the random walk with an OU/AR(1)
    process with dyad-specific equilibria and a global persistence `rho`,
    giving cross-dyad ordering a permanent component.
  * `alphanorm_ou` (`inst/stan/bilatr_alphanorm_ou.stan`) combines both:
    `alphanorm`'s identification with `ou`'s dynamics.
  * All three are new and UNCALIBRATED: several priors carry over from
    stable with materially changed units/meaning (documented per-model in
    each `.stan` file's header comment), and none has been checked
    against prior predictive simulation yet.

## Changes

* `assemble_stan_data()` gained `rho_prior_a`/`rho_prior_b` (default `8`,
  `2`) and `compute_log_lik` (default `0`) arguments, threaded into the
  returned data list. Consumed only by the three new experimental
  variants (`rho_prior_a`/`rho_prior_b` by `ou`/`alphanorm_ou`;
  `compute_log_lik` gates a per-dyad-period `log_lik` in `generated
  quantities` for all three, off by default since it is `D x T x` draws);
  ignored by `stable`, so existing callers are unaffected.
* `bilatr_init_fn()` now branches on `stan_model` to build
  parameterization-appropriate inits for each registered model, instead
  of always returning stable's init list regardless of the argument.

## Internal

* The `reduce_sum` likelihood (`partial_log_lik`) shared by all four
  registered models is no longer duplicated inline in each `.stan` file.
  Its single canonical source is
  `inst/stan/include/partial_log_lik.stanfunctions`;
  `data-raw/sync_stan_functions.R` splices it into a marker-delimited
  block in each model file (`R/stan_includes.R`), checked for drift by
  `tests/testthat/test_stan_includes.R`. This is NOT a Stan `#include`:
  cmdstanr has a bug (upstream issue
  [stan-dev/cmdstanr#820](https://github.com/stan-dev/cmdstanr/issues/820),
  fixed but not yet released) where `#include` resolution breaks at
  `$sample()`-time whenever the include path contains a space, which
  bites this project's own `devtools::load_all()` working tree.

# bilatr 0.3.2

## New features

* Added `BilatrClass2`, a 9-level (0-8) coarsening of `BilatrClass` that
  merges the two adjacent pairs of hostile levels: `BilatrClass` 7
  ("Investigate, demand, reject, or reduce relations") + 8 ("Disapprove")
  -> `BilatrClass2` 7 ("Disapprove, demand, reject, or reduce relations"),
  and `BilatrClass` 9 ("Threaten or coerce") + 10 ("Assault, fight, or
  mass violence") -> `BilatrClass2` 8 ("Threaten, coerce, or use force").
  `recode_cameo()` now attaches `BilatrClass2` / `BilatrClass2Name`
  alongside the other recode columns, `cameo_lookup` gained the two
  columns, and new internal helpers `assign_bilatr_class2()` /
  `bilatr_class2_name()` (`R/cameo_recode.R`) encode the mapping.

## Changes

* `recode_cameo()` now checks whether any of the recode columns it would
  add are already present in `data`; those are left untouched (no more
  silent `.x`/`.y` suffixing from the join) and a warning names them.
* Retired the experimental `phi_logn` model. `stable` is now the only
  entry in the internal model registry; `.resolve_stan_model("phi_logn")`
  (and `fit_*_dev(stan_model = "phi_logn")`) now error like any other
  unknown name. Its final source is kept at
  `inst/stan/legacy/bilatr_phi_logn.stan` (gitignored, not shipped), with
  the pre-0.3.0 centered version alongside it as
  `bilatr_phi_logn_pre_0.3.0.stan`. `bilatr_init_fn()` no longer branches
  on model name (no `log_phi0_raw` / `beta_logn` init), and `beta_logn`
  was dropped from `diagnose_convergence()`'s Tier 1 name list. The
  exported fitters were already `stable`-only and are unaffected.

# bilatr 0.3.1

## Internal

* Added internal helpers for the `reduce_sum` threading benchmark sweep
  (`R/benchmark.R`, all `@keywords internal`, not exported):
  `build_reduce_sum_grid()` (sizes a cores × grainsize × replicate grid
  from a scenario's dyad count), `insert_reduce_sum_profile()` (wraps the
  model's `reduce_sum` likelihood call in a Stan `profile()` block),
  `summarise_benchmark()`, `pareto_front()`, and `plot_benchmark()`.
  Covered by `tests/testthat/test_benchmark.R`. The SLURM sweep that
  uses them — re-tuning `threads_per_chain` × grainsize for the `stable`
  model on the politically-relevant vs full dyad sets — lives in
  `runscripts/benchmark/` (not shipped).
* `ggplot2` added to `Suggests` (used by `plot_benchmark()`).

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
