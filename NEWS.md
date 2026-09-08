
# bilatr 0.4.0

## Breaking changes

* Promoted the experimental `alphanorm` model to `stable`: it is now the
  model `fit_dyad_ts()`/`fit_panel()`/`compile_bilatr_model()` fit by
  default, replacing the previous `stable` (the consolidated
  Dirichlet-multinomial model with `alpha[1]` hard-fixed to 1 and a
  `mu_theta0`-anchored `theta0` location), which is retired to
  `inst/stan/legacy/bilatr_dirmult_irt_pre_0.4.0.stan` (gitignored, kept
  for local reference only, no longer registered/fittable). Renamed the
  experimental `alphanorm_ou` model to `ou`, still experimental,
  replacing the previous `ou` (an OU/AR(1) variant of the old stable
  model), similarly retired to
  `inst/stan/legacy/bilatr_ou_pre_0.4.0.stan`. The promoted `stable`/`ou`
  identify `alpha`/`mu_intercept` via `sum_to_zero_vector`s with a soft
  sign anchor on `alpha[1]` (new `anchor_scale`/`compute_log_lik`
  arguments to [assemble_stan_data()], already present since 0.3.x for
  the experimental variants) rather than hard-fixing `alpha[1] = 1`; this
  leaves an exact alpha/theta reflection symmetry that
  [bilatr_orient()] -- already applied by default in
  [extract_theta()]/[extract_alpha()]/[extract_mu_intercept()] -- now
  corrects for the default model too, not just the experimental one.
  `theta0`'s population mean is pinned at exactly 0 (no separate
  `mu_theta0`). CSVs from a fit made under the pre-0.4.0 `stable`/`ou`
  are unaffected and remain fully readable by
  [diagnose_convergence()]/[extract_theta()]/etc. (nothing about reading
  raw CmdStan CSVs depends on the registry), but code assuming
  `mu_theta0`, `mu_intercept_raw`, `alpha[1] == 1`, or
  `mu_intercept[1] == 0` will need updating for any model *fit* under
  this release's `stable`/`ou`, since `stan_model = "stable"`/`"ou"` now
  names a different Stan program than before.

# bilatr 0.3.10

## Bug fixes

* Fixed a memory-model bug in every CSV-file-path read
  (`diagnose_convergence()`, `diagnose_and_extract_bilatr()`,
  `extract_theta()`, `extract_alpha()`, `extract_mu_intercept()`) that
  could OOM-kill a job however small `max_memory_mb`/`chunk_size` was
  set, or however many cores `parallel = TRUE` was given: internally,
  every one of these previously called
  `cmdstanr::read_cmdstan_csv(csv_files, variables = ...)`, which reads
  each chain's CSV via `data.table::fread(cmd = "grep -v '^#' <file>")`
  -- and per `?data.table::fread`, a `cmd=`/piped input is always
  written to a full temporary copy in `tempdir()` before being read "as
  normal," regardless of how few `variables` are requested. That
  temp-file write happened on *every* call -- once per Tier 3 chunk, not
  once per run -- so a sweep with hundreds of chunks repeated a
  near-complete copy of each (often many-GB) chain file hundreds of
  times; and if `tempdir()`/`$TMPDIR` resolves to a RAM-backed `tmpfs`
  (a common per-node HPC/SLURM default), each copy was a direct hit
  against the job's memory allocation, independent of anything
  `max_memory_mb`/`chunk_size`/`n_workers` controlled.
* All of these now read via two new internal helpers,
  `.prepare_fast_csv_read()`/`.fast_read_post_warmup_draws()`: each
  chain file has its comment lines stripped exactly once per call (via
  a portable, streaming `readr::read_lines_chunked()` pass -- no
  external `grep` dependency), written next to the source file by
  default rather than through `tempdir()`, and every subsequent
  Tier 1/2/3 read/chunk reads `variables` straight from that one cleaned
  file via `data.table::fread(file = ...)` (a real path, so `fread` can
  select columns without buffering the whole file), verified to match
  `cmdstanr::read_cmdstan_csv()`'s output exactly. This restores the
  memory bound `max_memory_mb`/`chunk_size` was always meant to provide,
  makes `parallel = TRUE` safe to use at its documented memory cost
  again, and removes the redundant per-chunk file-copying regardless. A
  new `scratch_dir` argument (default `NULL`, meaning alongside each
  source file) lets you redirect the one-time cleaned copies elsewhere
  if needed.

# bilatr 0.3.9

## New features

* Added `diagnose_and_extract_bilatr()`, the single-shared-read follow-up
  0.3.8 anticipated: given raw CmdStan CSV file paths, it reads Tier 1/2
  once and Tier 3 once (chunked, same as `diagnose_convergence()`),
  reorients whichever columns `bilatr_orient()` would flip as part of
  that one read, and derives `diagnose_convergence()`,
  `extract_theta()`, `extract_alpha()`, and `extract_mu_intercept()`'s
  outputs from those two summaries -- rather than calling the four
  public functions independently, which re-reads `theta`'s columns
  twice (once for Tier 3 diagnostics, once for `extract_theta()`) and
  pays a full separate file-scan for each of `extract_alpha()`/
  `extract_mu_intercept()`, on top of `diagnose_convergence()`'s own
  Tier 1/2 read of the same columns. Meaningful for CSV-path,
  I/O-bound reads specifically (see 0.3.8); an in-memory `fit` has no
  repeated-read cost to fuse away, so it isn't supported here.

# bilatr 0.3.8

## Changes

* A benchmark against a production-scale CmdStan CSV (~1.9M Tier 3
  columns) found `read_cmdstan_csv()` strongly I/O-bound, not
  parsing-bound: per-chunk wall-time barely depends on how many
  variables are requested (a 1000x range in chunk size changed per-call
  time by under 10%), so total sweep time scales with chunk COUNT, not
  memory saved. `diagnose_convergence()`'s chunking is unchanged, but
  its default-`max_memory_mb` message now explains this tradeoff
  explicitly, so a smaller budget isn't mistaken for a "safer" choice.
* `extract_theta()`'s `fit` argument now also accepts a character vector
  of raw CmdStan CSV file paths, reading/summarising `theta` in the same
  memory-bounded, optionally-parallel chunks as
  `diagnose_convergence()`'s CSV path (new `max_memory_mb`/`chunk_size`/
  `parallel`/`n_workers` arguments, sharing its underlying helper so a
  worker script can eventually feed both from one pass over the CSVs
  instead of two). Sign orientation (`bilatr_orient()`) is applied to
  raw draws before summarising, per chunk, matching the in-memory path
  exactly rather than adjusting already-computed quantiles post hoc.
  Only supports the default `probs`. `extract_alpha()`/
  `extract_mu_intercept()` also gained CSV-path support (no chunking --
  Tier 1/2 is small regardless of dyad-set size), for cross-chain
  roll-up scripts run after several independently-submitted SLURM jobs
  have completed, with no in-memory fit available.

# bilatr 0.3.7

## New features

* Added preliminary human-readable labels for `EventRootCode2`
  (`eventrootcode2_name()`, matching the existing `bilatr_class_name()`/
  `bilatr_class2_name()` convention), attached as a new
  `EventRootCode2Name` column in `cameo_lookup` and by `recode_cameo()`.
  The two root-04 splits (`"044"`/`"046"`) reuse `bilatr_class_name()`'s
  text for the same category; every other label is new, since
  `EventRootCode2` is finer-grained than `BilatrClass` elsewhere.

# bilatr 0.3.6

## Changes

* `diagnose_convergence()`'s `fit` argument now also accepts a character
  vector of raw CmdStan CSV file paths (one per chain), for diagnosing a
  completed run without loading it into memory first via
  `cmdstanr::as_cmdstan_fit()`. In this mode, Tier 1/2 variables are read
  in one small call as before, but Tier 3 (`theta`/`theta_raw`, typically
  the overwhelming majority of monitored quantities for production-sized
  panels) is read and summarised in memory-bounded chunks via
  `cmdstanr::read_cmdstan_csv()`'s `variables` argument, discarding each
  chunk's draws before reading the next -- the full draws array is never
  materialized at once. New arguments `max_memory_mb` (default `8192`;
  drives an automatically-derived chunk size, so callers don't have to
  guess a variable count directly), `chunk_size` (explicit override),
  `parallel` (process chunks concurrently via `furrr::future_map_dfr()`;
  trades the sequential path's memory bound for wall-clock speed, and
  falls back from `future::multicore` to `future::multisession` with a
  warning on Windows), and `n_workers`. Existing calls (an in-memory fit
  or draws object, the previous signature) are unaffected -- the new
  arguments are unused on that path.

# bilatr 0.3.5

## Changes

* `alphanorm`/`alphanorm_ou`'s `alpha[1]` sign orientation is now
  deterministic across runs: their soft sign anchor alone could not
  guarantee it (it corrects relative posterior mass between the two
  mirror modes, but the modes are separated by a likelihood barrier a
  single chain essentially never crosses, so whichever mode a chain's
  init landed in was the one it reported). `bilatr_init_fn()` now biases
  initial values toward the canonical (`alpha[1] > 0`) mode and gives the
  latent states real initial spread, `fit_bilatr()` warns after sampling
  if a fit still came back in the wrong mode, and the new
  `bilatr_orient()` (`R/orient.R`) deterministically relabels a fit's
  draws to the canonical orientation regardless of which mode the
  sampler found -- now wired into `extract_theta()`/`extract_alpha()`/
  `extract_mu_intercept()` via their new `stan_model` argument, so no
  consumer of those functions sees an unoriented fit. Fits produced
  before this change may be sign-flipped relative to fits produced after
  it; `bilatr_orient()` relabels them to match.

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
