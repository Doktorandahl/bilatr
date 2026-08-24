#' Cache of compiled CmdStanModel objects, keyed by "<stan file>::<opt_level>"
#' @keywords internal
.bilatr_model_cache <- new.env(parent = emptyenv())

#' Compile (or fetch from cache) a registered Stan model
#'
#' Shared compilation path for both the exported [compile_bilatr_model()]
#' (always the stable model) and the internal `_dev` fitters (any
#' registered model; see `R/model-registry.R`). Compiled `CmdStanModel`
#' objects are cached in `.bilatr_model_cache` by resolved file path and
#' `opt_level`, so repeated calls (e.g. across dev iterations) don't pay
#' `cmdstanr`'s own hash-checking/Makefile overhead every time.
#'
#' @param stan_model Name registered in `.bilatr_stan_models`.
#' @param opt_level C++ compiler optimization level; see
#'   [compile_bilatr_model()].
#' @param force_recompile If `TRUE`, bypass the cache and rebuild.
#' @return A `CmdStanModel` object.
#' @keywords internal
.compile_stan_model <- function(
  stan_model,
  opt_level,
  force_recompile = FALSE
) {
  stan_file <- .resolve_stan_model(stan_model)
  cache_key <- paste(stan_file, opt_level, sep = "::")

  if (
    !force_recompile &&
      exists(cache_key, envir = .bilatr_model_cache, inherits = FALSE)
  ) {
    return(get(cache_key, envir = .bilatr_model_cache, inherits = FALSE))
  }

  mod <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = TRUE, O = opt_level),
    force_recompile = force_recompile
  )
  assign(cache_key, mod, envir = .bilatr_model_cache)
  mod
}

#' Compile the bilatr Stan model
#'
#' Compiles the package's single collapsed Stan program
#' (`inst/stan/bilatr_dirmult_irt.stan`), used by both [fit_dyad_ts()] and
#' [fit_panel()]. Compiled executables are cached internally by
#' `opt_level`, so repeated calls with the same `opt_level` are cheap
#' after the first.
#'
#' @param opt_level C++ compiler optimization level passed to CmdStan's
#'   `O` Makefile variable (0-3). The reduce_sum-based likelihood in this
#'   model generates a lot of C++ template code; at the default `O = 3`,
#'   compilation can use several GB of RAM and occasionally exhausts
#'   memory on constrained machines (laptops, small CI runners). Dropping
#'   to `opt_level = 1` cuts compilation RAM substantially at the cost of
#'   somewhat slower sampling; `opt_level = 0` is fastest to compile but
#'   slowest to sample and rarely worth it. If compilation is being
#'   killed by the OOM killer, try `opt_level = 1` before anything else.
#' @param force_recompile Set `TRUE` to force a rebuild (e.g. after
#'   changing `opt_level`), bypassing the internal cache.
#' @return A `CmdStanModel` object.
#' @examples
#' \dontrun{
#' mod <- compile_bilatr_model()
#' mod_low_ram <- compile_bilatr_model(opt_level = 1, force_recompile = TRUE)
#' }
#' @export
compile_bilatr_model <- function(opt_level = 3, force_recompile = FALSE) {
  .compile_stan_model(.BILATR_DEFAULT_MODEL, opt_level, force_recompile)
}

#' Build an initial-value generator matching a model's parameterization
#'
#' `stan_model == "stable"` (the default, used by the exported fitters)
#' initializes the per-dyad-constant `phi`. Other registered models may
#' need different/additional init fields for parameters not present in
#' the stable model; `"phi_logn"` initializes `log_phi0_raw`/`beta_logn`
#' instead of `phi` (which isn't a free parameter in that variant).
#'
#' @param stan_data A Stan data list as returned by [assemble_stan_data()].
#' @param stan_model Name registered in `.bilatr_stan_models`.
#' @return A zero-argument function suitable for `cmdstanr`'s `init`
#'   argument.
#' @keywords internal
bilatr_init_fn <- function(stan_data, stan_model = .BILATR_DEFAULT_MODEL) {
  function() {
    init <- list(
      mu_theta0 = 0,
      sigma_theta0 = 0.5,
      z_theta0 = rep(0, stan_data$D),
      theta_raw = matrix(0, stan_data$D, stan_data$T),
      mu_log_noise = log(0.2),
      sigma_log_noise = 0.3,
      mu_log_phi = 0,
      sigma_log_phi = 0.5,
      process_noise = rep(0.2, stan_data$D),
      alpha_raw = rep(0, stan_data$A - 2),
      alpha_hostile = 0.5,
      mu_intercept_raw = rep(0, stan_data$A - 1)
    )

    if (identical(stan_model, "phi_logn")) {
      init$log_phi0_raw <- rep(0, stan_data$D)
      init$beta_logn <- 0
    } else {
      init$phi <- rep(1, stan_data$D)
    }

    init
  }
}

#' Shared sampling logic behind fit_dyad_ts()/fit_panel() and their _dev
#' counterparts in R/fit-dev.R
#' @keywords internal
fit_bilatr <- function(
  stan_data,
  chains,
  parallel_chains,
  threads_per_chain,
  iter_warmup,
  iter_sampling,
  seed,
  opt_level,
  output_dir,
  stan_model = .BILATR_DEFAULT_MODEL,
  ...
) {
  mod <- .compile_stan_model(stan_model, opt_level)
  mod$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    seed = seed,
    init = bilatr_init_fn(stan_data, stan_model = stan_model),
    output_dir = output_dir,
    ...
  )
}

#' Fit the bilatr model to a single dyad's time series
#'
#' Fits the collapsed dyadic IRT model (see the package's Stan file,
#' `bilatr_dirmult_irt.stan`) to a single dyad (`D == 1`), i.e. estimates
#' one latent conflict trajectory `theta` with no cross-dyad pooling on
#' `theta0`, `phi`, or process noise (those hierarchical parameters still
#' exist in the model but are estimated from a single unit, so their
#' pooling has no effect). Use [fit_panel()] instead when fitting more
#' than one dyad, to get the benefit of partial pooling.
#'
#' `reduce_sum` chunking across dyads has nothing to parallelize over
#' when `D == 1`, so `threads_per_chain` defaults to 1 here (raising it
#' will not speed up sampling for a single-dyad fit).
#'
#' @param stan_data A Stan data list from [assemble_stan_data()] with
#'   `D == 1`.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of chains to run in parallel.
#' @param threads_per_chain Threads per chain for within-chain
#'   parallelization. Not useful for `D == 1`; left as an argument mainly
#'   for consistency with [fit_panel()].
#' @param iter_warmup Number of warmup iterations per chain.
#' @param iter_sampling Number of post-warmup sampling iterations per
#'   chain.
#' @param seed Random seed, or `NULL` for `cmdstanr`'s default.
#' @param opt_level Compiler optimization level; see
#'   [compile_bilatr_model()].
#' @param output_dir Directory to write CmdStan's raw output CSVs to, or
#'   `NULL` for a temporary directory.
#' @param ... Additional arguments passed to `CmdStanModel$sample()`.
#' @return A `CmdStanMCMC` fit object.
#' @examples
#' \dontrun{
#' stan_data <- assemble_stan_data(
#'   dplyr::filter(events, dyad == "USA_CHN"),
#'   years = 2015:2020, resolution = "yearly", grouping_var = "PentaClass",
#'   reference_category = 0, reference_hostile = 4
#' )
#' fit <- fit_dyad_ts(stan_data, chains = 4, iter_sampling = 1000)
#' }
#' @export
fit_dyad_ts <- function(
  stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 2,
  output_dir = NULL,
  ...
) {
  if (stan_data$D != 1) {
    stop(
      "fit_dyad_ts() expects a single-dyad stan_data (D == 1), got D = ",
      stan_data$D,
      ". ",
      "Use fit_panel() for multiple dyads.",
      call. = FALSE
    )
  }
  fit_bilatr(
    stan_data,
    chains,
    parallel_chains,
    threads_per_chain,
    iter_warmup,
    iter_sampling,
    seed,
    opt_level,
    output_dir,
    stan_model = .BILATR_DEFAULT_MODEL,
    ...
  )
}

#' Fit the bilatr model to a panel of dyads with hierarchical pooling
#'
#' Fits the collapsed dyadic IRT model (see the package's Stan file,
#' `bilatr_dirmult_irt.stan`) across multiple dyads simultaneously
#' (`D > 1`), with hierarchical partial pooling on process noise (via
#' `mu_log_noise`/`sigma_log_noise`), the dispersion parameter `phi` (via
#' `mu_log_phi`/`sigma_log_phi`), and the initial latent state `theta0`
#' (via `mu_theta0`/`sigma_theta0`). Use [fit_dyad_ts()] instead for a
#' single dyad.
#'
#' `threads_per_chain` controls `reduce_sum` parallelization across
#' dyads within a chain, combined with the `chunk_size` (`C`) baked into
#' `stan_data` by [assemble_stan_data()]. 16 threads with `chunk_size =
#' 600` was the Pareto-optimal setting found in this project's own
#' benchmarking for panel-sized data; tune both together for your data
#' size and hardware rather than one in isolation.
#'
#' @inheritParams fit_dyad_ts
#' @param stan_data A Stan data list from [assemble_stan_data()] with
#'   `D > 1`.
#' @param threads_per_chain Threads per chain for `reduce_sum`
#'   parallelization across dyads. Defaults to 16, the setting found to
#'   be Pareto-optimal (together with `chunk_size = 600` at the data-
#'   assembly stage) in this project's own benchmarking; reduce it on
#'   machines with fewer cores.
#' @return A `CmdStanMCMC` fit object.
#' @examples
#' \dontrun{
#' stan_data <- assemble_stan_data(
#'   events, years = 2015:2020, resolution = "yearly", grouping_var = "PentaClass",
#'   reference_category = 0, reference_hostile = 4, chunk_size = 600
#' )
#' fit <- fit_panel(stan_data, chains = 4, threads_per_chain = 16, iter_sampling = 1000)
#' }
#' @export
fit_panel <- function(
  stan_data,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 16,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 3,
  output_dir = NULL,
  ...
) {
  if (stan_data$D < 2) {
    stop(
      "fit_panel() expects multiple dyads (D >= 2), got D = ",
      stan_data$D,
      ". ",
      "Use fit_dyad_ts() for a single dyad.",
      call. = FALSE
    )
  }
  fit_bilatr(
    stan_data,
    chains,
    parallel_chains,
    threads_per_chain,
    iter_warmup,
    iter_sampling,
    seed,
    opt_level,
    output_dir,
    stan_model = .BILATR_DEFAULT_MODEL,
    ...
  )
}
