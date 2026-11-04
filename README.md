# bilatr

bilatr estimates latent bilateral relations between pairs of actors
over time from dyadic event data. A hierarchical Bayesian dynamic
item-response model, fitted in Stan, places every dyad-period on a scale
from conflict to cooperation: action types (for example CAMEO-coded
GDELT events) discriminate between more and less cooperative states
through a Dirichlet-multinomial likelihood, and each dyad's latent state
follows a random walk, pooled across dyads. The package covers the whole
workflow: reading GDELT events or your own, aggregating them, fitting
single dyads or panels, extracting the estimates, diagnosing convergence
by parameter tier, checking model fit, and tuning parallel performance.

## Installation

```r
install.packages("bilatr")

# development version
# install.packages("pak")
pak::pak("doktorandahl/bilatr")
```

Fitting models needs the cmdstanr package and CmdStan. cmdstanr is not
on CRAN; install it from the Stan r-universe, then install CmdStan
(this also needs a C++ toolchain; see the [cmdstanr
guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)):

```r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()

bilatr::bilatr_cmdstan_available() # TRUE when fitting will work
```

Everything else -- data preparation, and every post-fit function on
saved CmdStan output -- works without them.

## Quick start

`bilatr_example_fit()` returns the CmdStan output of a small fit that
ships with the package, so you can try the post-fit tools without
CmdStan:

```r
library(bilatr)
csv <- bilatr_example_fit()          # 6 dyads, 2015-2020, 5 action classes

diagnose_convergence(csv)            # R-hat and ESS by parameter tier
theta <- extract_theta(csv)          # latent relations per dyad and year
head(theta)
extract_alpha(csv)                   # discrimination of each action class
curves <- icc_curves(csv, theta_range = c(-2, 2))
plot(curves)                         # needs ggplot2
```

## Workflow

```r
library(bilatr)

# 1. Events: GDELT 1.0 (or any table in the format of ?bilatr_event_data)
status <- download_gdelt("2020-01-01", "2020-12-31", dest_dir = "gdelt")
events <- read_gdelt(status$path[status$status %in% c("downloaded", "cached")])
events <- recode_cameo(events)       # QuadClass, PentaClass, ModifiedRootCode

# 2. Counts per dyad and period (see ?bilatr_dyad_period_data)
agg <- aggregate_events(events, grouping_var = "PentaClass", resolution = "yearly")

# 3. Fit. reference_category orients the scale: higher theta makes that
#    class more likely, so a cooperative class makes higher theta mean
#    more cooperative relations.
fit <- fit_panel(
  agg,
  reference_category = "1", # verbal cooperation
  chains = 4, parallel_chains = 4, seed = 1
)

# 4. Check and extract
diagnose_convergence(fit)
theta <- extract_theta(fit)
check_compositional_residuals(fit)
```

`fit_dyad_ts()` fits a single dyad. For large panels, `threads_per_chain`
and `chunk_size` control within-chain parallelism; see
`vignette("tuning")`.

## The model

For dyad *d* in period *t*, the event counts over the *A* action classes
follow a Dirichlet-multinomial with class probabilities
`softmax(alpha * theta[d, t] - mu_intercept)` and dyad-specific
concentration `phi[d]`. `theta` follows a random walk per dyad;
process noise, `phi` and the initial states are partially pooled across
dyads. `alpha` and `mu_intercept` sum to zero, and `alpha` has
root-mean-square 1, which fixes the scale of `theta`.

**Orientation.** θ is oriented so that higher values make
`reference_category` more likely relative to the other action types;
with a cooperative reference category, higher θ means more cooperative
relations (see `?prepare_stan_data`).

Two experimental variants are available through `stan_model`: `"ou"`,
with mean-reverting latent states, and `"stable_gamma"`, with
country-level offsets for each action class (see
`?compile_bilatr_model`).

## Estimated bilateral relations, 1991-2025

The package's estimates for every country pair with at least one
recorded event, 1991-2025, from GDELT 1.0 government, military and
intelligence events, are published as a versioned data release (CC BY
4.0) and downloaded on first use:

```r
theta_directed <- bilatr_theta("directed")       # 16,750 directed dyads
theta_undirected <- bilatr_theta("undirected")   # 9,139 country pairs
```

`theta_sample` is a small extract that ships with the package. Known
gap: six states (Serbia, Romania, Bosnia and Herzegovina, Montenegro,
Slovenia and Timor-Leste) are missing from the current release because
of a country-code matching error; this will be corrected in the next
data update, expected in January 2027. See `?theta_trajectories`.

## Vignettes

- `vignette("dyad_time_series")`: one dyad, from events to a fitted
  trajectory.
- `vignette("panel_model")`: a panel of dyads.
- `vignette("diagnostics")`: tiered convergence diagnostics.
- `vignette("tuning")`: threads, chunk size and chains.

## Citation

```r
citation("bilatr")
```

When you use the estimated trajectories, please also cite the data
release and GDELT (see `?theta_trajectories`).

## License

MIT. The estimated trajectories are licensed CC BY 4.0.
