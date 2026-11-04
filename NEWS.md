# bilatr 1.0.0

First CRAN release.

bilatr estimates latent bilateral relations between pairs of actors over
time from dyadic event data, with a hierarchical dynamic item-response
model fitted in 'Stan'.

## Data input and validation

* `download_gdelt()`, `read_gdelt()`, `gdelt_files()`, `gdelt_columns()` and
  `count_gdelt_actor_types()` download, read and filter GDELT 1.0 event
  files, with checksum verification and memory-bounded reading.
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
* `fit_panel()` and `fit_dyad_ts()` fit three models: `"stable"` (the
  default, a random-walk latent state) and the experimental `"ou"`
  (mean-reverting) and `"stable_gamma"` (country-level offsets). The
  latent scale is oriented by `reference_category`. Fits are reproducible
  with `seed`, and record which model and data they used.
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
  (`?theta_trajectories`), are downloaded with `bilatr_theta()` from a
  versioned data release (`v2026.11`); `theta_sample` ships a small
  extract. Known gap: six states (Serbia, Romania, Bosnia and Herzegovina,
  Montenegro, Slovenia and Timor-Leste) are missing because of a
  country-code crosswalk error, to be corrected in the next data release,
  expected in January 2027.
