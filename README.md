# bilatr

`bilatr` fits a hierarchical Bayesian dynamic IRT-style model to dyadic
event data (CAMEO-coded, or any categorically coded event type) to
estimate latent conflict trajectories (`theta`) between country pairs
over time. Action types discriminate between high- and low-conflict
states through a Dirichlet-multinomial likelihood, with dyad-specific
latent states following a random-walk process. The package supports
both single-dyad time-series estimation and multi-dyad panel estimation
with hierarchical pooling. You supply the event data (see `?bilatr_event_data`);
GDELT download/ingest helpers and a CAMEO recoding table are included.

## Installation

`bilatr` depends on [`cmdstanr`](https://mc-stan.org/cmdstanr/), which
in turn requires a working CmdStan installation and a C++ toolchain.
Set those up first:

```r
# 1. Install cmdstanr itself (not on CRAN)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

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

# --- 2. Assemble Stan data ---------------------------------------------
# actor1/actor2/date name your event table's columns (defaults shown
# here match GDELT's own column names, so they can be omitted for GDELT
# data). reference_category anchors the model's scale/sign reference;
# every other action class's discrimination (alpha[2:A]) is freely
# estimated.
stan_data <- assemble_stan_data(
  events,
  years = 2015:2020,
  resolution = "yearly",
  grouping_var = "PentaClass",
  reference_category = 0, # verbal cooperation
  min_n_events = 10,
  actor1 = "Actor1CountryCode",
  actor2 = "Actor2CountryCode",
  date = "SQLDATE"
)

# --- 3. Fit --------------------------------------------------------------
# One dyad at a time (D == 1):
fit <- fit_dyad_ts(stan_data, chains = 4, iter_warmup = 1000, iter_sampling = 1000)

# Or a full panel (D > 1), with hierarchical pooling on theta0/phi/process noise:
# fit <- fit_panel(stan_data, chains = 4, threads_per_chain = 16, iter_warmup = 1000, iter_sampling = 1000)

# --- 4. Extract ------------------------------------------------------------
theta <- extract_theta(fit, stan_data)          # latent trajectories, with dyad IDs reattached
alpha <- extract_alpha(fit, event_classes = attr(stan_data, "event_classes"))
mu_intercept <- extract_mu_intercept(fit, event_classes = attr(stan_data, "event_classes"))

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
  FOLD") reports every fit with `alpha[1] >= 0`, so positive `alpha[1]`
  always means higher `theta` corresponds to better (less hostile)
  relations at the reference/neutral action class -- without truncating
  or excluding any part of the sampled parameter space. There is no
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
`gamma` on top of `stable`'s likelihood). See `R/model_registry.R` and
[`assemble_stan_data()`] for details.

## Performance notes

- `assemble_stan_data()`'s `chunk_size` argument sets the `reduce_sum`
  grainsize used to chunk the likelihood across dyads. 16 threads with
  `chunk_size = 600` was Pareto-optimal in this project's own
  benchmarking for panel-sized data; tune both together for your data
  size and hardware.
- `compile_bilatr_model()`'s `opt_level` argument controls the C++
  compiler optimization level. The `reduce_sum`-based likelihood
  generates a lot of template code, and compiling at the default
  `opt_level = 3` can use several GB of RAM — enough to get killed by
  the OOM killer on laptops or small CI runners. Drop to `opt_level = 1`
  if compilation is failing or crashing for that reason; it trades a
  slower-to-sample model for a much cheaper compile.
