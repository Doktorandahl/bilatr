# bilatr 1.0.0

First CRAN release.

bilatr estimates latent bilateral relations between pairs of actors over
time from dyadic event data, with a hierarchical dynamic item-response
model fitted in 'Stan'.

## Data input and validation

* `download_gdelt()`, `read_gdelt()`, `gdelt_files()`, `gdelt_columns()` and
  `count_gdelt_actor_types()` download, read and filter GDELT 1.0 event
  files, with checksum verification and memory-bounded reading.
* `gdelt_countries` classifies every GDELT 1.0 country code (state,
  former state, territory, region or other) with its Correlates of War
  membership years; `gdelt_country_name()`, `gdelt_country_type()`,
  `gdelt_is_state()`, `gdelt_country_iso3()`, `filter_gdelt_states()` and
  `gdelt_state_spells()` look codes up, keep state-to-state events in the
  year of the event, and harmonize codes to ISO3 (e.g. `ROM` to `ROU`).
  `read_gdelt()` and `count_gdelt_actor_types()` gain `states_only`,
  `harmonize_codes` and `recover_country_codes` (which fills the empty
  country codes GDELT leaves for Romanian, Taiwanese (before 2013) and
  South Sudanese actors).
* `recode_cameo()` and the `cameo_lookup` table map CAMEO event codes to
  QuadClass (as in GDELT), PentaClass and an 18-class ModifiedRootCode
  scheme.
* `validate_bilatr_events()` and `validate_bilatr_dyad_period()` check
  event tables and dyad-period tables against documented formats
  (`?bilatr_event_data`, `?bilatr_dyad_period_data`), so any categorically
  coded event data can be used.

## Models and fitting

* `aggregate_events()` counts events per dyad and period (directed or
  undirected; yearly, monthly, weekly or daily); `prepare_stan_data()`
  builds the Stan data.
* `aggregate_events(actor_periods = )` marks dyad-periods in which an
  actor does not exist (e.g. a country before it enters the state system)
  with an `active` column, which `prepare_stan_data()` validates and
  carries into `dyad_ids`; `extract_theta()` and
  `diagnose_and_extract_bilatr()` drop those periods by default
  (`active_only = TRUE`), and Tier 3 diagnostics leave them out.
* `fit_panel()` and `fit_dyad_ts()` fit three models: `"stable"` (the
  default, a random-walk latent state) and the experimental `"ou"`
  (mean-reverting) and `"stable_gamma"` (country-level offsets). The
  latent scale is oriented by `reference_category`. Fits are reproducible
  with `seed`, and record which model and data they used.
* In the `"stable_gamma"` model, the prior on the per-category country-offset
  spread `sigma_gamma` is widened from half-normal(0, 0.3) to
  half-normal(0, 0.5).
* Models compile with `compile_bilatr_model()` into the user cache
  directory, never into the package library; `bilatr_clear_model_cache()`
  removes them.

## Post-fit tools

* `extract_theta()`, `extract_alpha()`, `extract_mu_intercept()` and
  `extract_gamma()` return tidy posterior summaries, from a fit object or
  directly from CmdStan CSV files, read in memory-bounded chunks.
* `diagnose_convergence()` reports R-hat and ESS by parameter tier, so
  global parameters are never masked by sparse dyads;
  `diagnose_and_extract_bilatr()` diagnoses and extracts in one read.
* `check_compositional_residuals()`, `icc_curves()`,
  `diagnose_category_merges()` and `merge_cost()` check model fit and the
  action-class scheme.
* `bilatr_example_fit()` ships a small precomputed fit, so every post-fit
  function can be tried without CmdStan.

## Tuning

* `benchmark_bilatr()`, `benchmark_grid()`, `insert_reduce_sum_profile()`,
  `summarize_benchmark()`, `benchmark_pareto()` and `plot_benchmark()`
  measure threads per chain and chunk size; see `vignette("tuning")`.

## Data

* Estimated yearly trajectories for every country pair, 1991-2025
  (`?theta_trajectories`), are downloaded once and cached with
  `download_bilatr_theta()` from a versioned data release (`v2026.11.1`,
  doi:10.5281/zenodo.24109618, CC BY 4.0); `theta_sample` ships a small
  extract. Known gap: Slovenia, Bosnia and Herzegovina, Kosovo and
  Montenegro are absent, because GDELT 1.0 has no country code for them.
  The earlier data release `v2026.11` is superseded and is kept only for
  reproducing earlier work.
