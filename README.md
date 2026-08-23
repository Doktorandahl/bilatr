# bilatr

`bilatr` fits a hierarchical Bayesian dynamic IRT-style model to dyadic
event data (GDELT/ICEWS, CAMEO-coded) to estimate latent conflict
trajectories (`theta`) between country pairs over time. Action types
discriminate between high- and low-conflict states through a
Dirichlet-multinomial likelihood, with dyad-specific latent states
following a random-walk process. The package supports both single-dyad
time-series estimation and multi-dyad panel estimation with
hierarchical pooling.

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
pak::pak("randahl/bilatr")
```

## Quick start

```r
library(bilatr)
library(dplyr)

# --- 1. Data prep -----------------------------------------------------
# Start from a raw GDELT export (or ingest_icews() for ICEWS data):
events <- extract_all_relevant_gdelt("data/gdelt_raw/20200101.zip")

# Recode CAMEO event codes to QuadClass/PentaClass using the package's
# built-in lookup table (no external CAMEO reference package needed):
events <- recode_cameo(events)

# --- 2. Assemble Stan data ---------------------------------------------
# reference_category anchors the model's scale/neutral action (alpha[1] = 1);
# reference_hostile anchors the known-hostile end (alpha[A], forced negative).
stan_data <- assemble_stan_data(
  events,
  years = 2015:2020,
  resolution = "yearly",
  grouping_var = "PentaClass",
  reference_category = 0, # verbal cooperation
  reference_hostile = 4,  # most severe conflict class
  min_n_events = 10
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

See `vignette("dyad-time-series")` and `vignette("panel-model")` for
walkthroughs covering both estimation modes end to end, including how
to read the convergence diagnostics.

## Model overview

- **Likelihood**: dyad-period event-type counts follow a
  Dirichlet-multinomial, with concentration `phi[d] * softmax(alpha .*
  theta[d,t] - mu_intercept)` (optionally rescaled by
  `dyad_weight`/`period_weight`/`action_weight`, all default to 1s).
- **Identification**: `alpha[1] = 1` and `mu_intercept[1] = 0` are fixed
  reference points; a known hostile action class anchors the negative
  end of `alpha`. There is no dyad-specific intercept — cross-dyad level
  differences are absorbed into the global `mu_intercept`, which is what
  keeps `theta` comparable across dyads.
- **Dynamics**: `theta` follows a random walk per dyad, starting from a
  fully hierarchical `theta0` (`mu_theta0`, `sigma_theta0`), not a
  fixed or data-supplied prior.
- **Pooling** (panel mode only): process noise, `phi`, and `theta0` are
  partially pooled across dyads via lognormal/normal hyperpriors.

One Stan program (`inst/stan/bilatr_dirmult_irt.stan`) covers both
estimation modes and all reweighting configurations — see
[`assemble_stan_data()`] for how the weight vectors are constructed.

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
