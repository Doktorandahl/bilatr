# bilatr

`bilatr` fits a hierarchical Bayesian dynamic IRT-style model to dyadic
event data (CAMEO-coded, or any categorically coded event type) to
estimate latent bilateral relations (`theta`) between country pairs
over time, on a scale running from conflict to cooperation. Action types
discriminate between more and less cooperative states through a Dirichlet-multinomial likelihood, with dyad-specific
latent states following a random-walk process. The package supports
both single-dyad time-series estimation and multi-dyad panel estimation
with hierarchical pooling. You supply the event data (see `?bilatr_event_data`);
GDELT download/ingest helpers and a CAMEO recoding table are included.

## Installation

Fitting models needs [`cmdstanr`](https://mc-stan.org/cmdstanr/), which
in turn requires a working CmdStan installation and a C++ toolchain.
`cmdstanr` is a suggested package, so it is not installed with `bilatr`:
install it as below. Data preparation, and reading and checking
CmdStan output CSV files, work without it.
Set those up first:

```r
# 1. Install cmdstanr itself (not on CRAN)
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# 2. Check for a C++ toolchain (installs one if needed, on some platforms)
cmdstanr::check_cmdstan_toolchain()

# 3. Install CmdStan (the command-line Stan interface bilatr compiles models against)
cmdstanr::install_cmdstan(cores = 4)
```

See the [cmdstanr installation
guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) if you hit
toolchain issues (Windows in particular sometimes needs RTools
installed first).

Then install `bilatr` itself. From a local clone:

```r
# install.packages("pak")
pak::pak(".")

# or
devtools::install(".")
```

Or, once pushed to GitHub:

```r
pak::pak("doktorandahl/bilatr")
```

## Quick start

```r
library(bilatr)
library(dplyr)

# --- 1. Data prep -----------------------------------------------------
# Bring an event table in the format documented at ?bilatr_event_data:
# an actor1/actor2 pair (any code alphabet -- ISO3, ISO2, COW numeric,
# free-text labels), an event date, and an event-class column. GDELT
# helpers are included but optional -- two download modes:

# a) a quick look: download, filter, and read into memory, nothing kept
#    on disk
events <- download_gdelt("2020-01-01")

# -- OR --

# b) for real work: keep the zips as a cache (skips files already
#    present on a later call) and read them yourself
status <- download_gdelt("2020-01-01", "2020-01-31", dest_dir = "data/gdelt_raw")
events <- read_gdelt(status$path[status$status %in% c("downloaded", "cached")])

# Recode CAMEO event codes to QuadClass/PentaClass using the package's
# built-in lookup table (no external CAMEO reference package needed):
events <- recode_cameo(events)

# --- 2. Aggregate to a dyad-period table ---------------------------------
# actor1/actor2/date name your event table's columns (the defaults match
# GDELT's own column names, so they can be omitted for GDELT data).
# The result has one row per dyad and year, and one count column per
# event class (see ?bilatr_dyad_period_data). Any other source of
# already aggregated dyad-period counts can be used in its place.
agg <- aggregate_events(
  events,
  grouping_var = "PentaClass",
  resolution = "yearly",
  start = "2015-01-01",
  end = "2020-12-31"
)

# --- 3. Fit --------------------------------------------------------------
# reference_category orients theta: higher theta makes that class more
# likely relative to the other action types, so with a cooperative
# reference higher theta means more cooperative relations (see
# ?prepare_stan_data). Every other action class's discrimination
# (alpha[2:A]) is freely estimated.
# A full panel (D > 1), with hierarchical pooling on theta0/phi/process noise:
fit <- fit_panel(
  agg,
  reference_category = "1", # verbal cooperation
  min_n_events = 10,
  chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1
)

# Or one dyad at a time (D == 1):
# usa_chn <- filter(agg, dyad == "USA_CHN")
# fit <- fit_dyad_ts(usa_chn, reference_category = "1", chains = 4) # verbal cooperation

# --- 4. Extract ------------------------------------------------------------
# The fit carries the Stan data it was fit to (attr(fit, "stan_data")),
# so the extractors need nothing else.
theta <- extract_theta(fit)            # latent trajectories, with dyad IDs reattached
alpha <- extract_alpha(fit)            # labelled with the event classes
mu_intercept <- extract_mu_intercept(fit)

# --- 5. Diagnose -----------------------------------------------------------
fit$summary(variables = c("alpha", "mu_intercept", "phi"))  # includes rhat, ess_bulk, ess_tail
fit$diagnostic_summary()                                     # divergences, tree-depth saturation
```

See `vignette("dyad_time_series")` and `vignette("panel_model")` for
walkthroughs covering both estimation modes end to end, including how
to read the convergence diagnostics.

## Model overview

This describes the registered `stable` Stan program
(`inst/stan/bilatr_alphanorm.stan`), fit by `fit_dyad_ts()`/`fit_panel()`
by default.

- **Likelihood**: dyad-period event-type counts follow a
  Dirichlet-multinomial, with concentration
  `phi[d] * softmax(alpha .* theta[d,t] - mu_intercept)`.
- **Identification**: `alpha` has RMS (population SD) exactly 1 and sums
  to zero, by construction; `mu_intercept` is likewise a
  `sum_to_zero_vector`, with no fixed-first-to-0 element. An orientation
  fold (see the `.stan` file's header, "IDENTIFICATION: ORIENTATION
  FOLD") reports every fit with `alpha[1] >= 0`, without truncating or
  excluding any part of the sampled parameter space. θ is therefore
  oriented so that higher values make `reference_category` more likely
  relative to the other action types; with a cooperative reference
  category, higher θ means more cooperative relations (see the
  `reference_category` argument of `?prepare_stan_data`). There is no
  dyad-specific intercept — cross-dyad level differences are absorbed
  into the global `mu_intercept`, which is what keeps `theta` comparable
  across dyads.
- **Dynamics**: `theta` follows a random walk per dyad, starting from
  `theta0 = sigma_theta0 * z_theta0` (population mean pinned at 0, no
  separate `mu_theta0` location parameter).
- **Pooling** (panel mode only): process noise and `phi` are partially
  pooled across dyads via lognormal hyperpriors; the `process_noise`
  hierarchy is sampled non-centered.

Two other Stan programs are registered as experimental variants: `ou`
(an Ornstein-Uhlenbeck `theta` with mean reversion instead of a pure
random walk) and `stable_gamma` (adds a country-level category-offset
`gamma` on top of `stable`'s likelihood; it needs actor columns in the
dyad-period table). See `?compile_bilatr_model` and `?prepare_stan_data`
for details. Compiled models are cached in a user cache directory (see
`?compile_bilatr_model`), never in the package library.

## Performance notes

- The defaults of `fit_panel()` (`threads_per_chain = 1`,
  `parallel_chains = 1`, `chunk_size = 1`) are conservative and use one
  core. On a multi-core machine, run chains in parallel
  (`parallel_chains`) and give each chain several threads
  (`threads_per_chain`) with a larger `chunk_size` (the `reduce_sum`
  grainsize used to chunk the likelihood across dyads); tune them
  together for your data size and hardware, as described in
  `vignette("tuning")`.
- The Stan data's count array `Y` is dense, dyads x periods x action
  types (about `D * T * A * 4` bytes), so very fine time resolutions
  (e.g. daily) are only practical for few dyads or short windows.
- `compile_bilatr_model()`'s `opt_level` argument controls the C++
  compiler optimization level. The `reduce_sum`-based likelihood
  generates a lot of template code, and compiling at the default
  `opt_level = 3` can use several GB of RAM — enough to get killed by
  the OOM killer on laptops or small CI runners. Drop to `opt_level = 1`
  if compilation is failing or crashing for that reason; it trades a
  slower-to-sample model for a much cheaper compile.
