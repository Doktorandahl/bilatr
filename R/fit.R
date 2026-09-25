#' Cache of compiled CmdStanModel objects, keyed by "<stan file>::<opt_level>"
#' @keywords internal
.bilatr_model_cache <- new.env(parent = emptyenv())

#' Compile (or fetch from cache) a registered Stan model
#'
#' Shared compilation path for both the exported [compile_bilatr_model()]
#' (always the stable model) and the internal `_dev` fitters (any
#' registered model; see `R/model_registry.R`). Compiled `CmdStanModel`
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
#' Compiles the package's registered `stable` Stan program (see
#' `R/model_registry.R`), used by both [fit_dyad_ts()] and [fit_panel()].
#' Compiled executables are cached internally by `opt_level`, so repeated
#' calls with the same `opt_level` are cheap after the first.
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

#' A sign-biased sum-to-zero init draw for `sum_to_zero_vector[A]` alpha
#' parameters
#'
#' Every registered model's `alpha_raw` is a free `sum_to_zero_vector[A]`,
#' normalized by `sqrt(A / dot_self(alpha_raw))`, so an all-zero init
#' (otherwise the natural default) would divide by zero; this draws a
#' random vector instead, centered to sum to exactly 0 (required by
#' `sum_to_zero_vector`'s constrained representation) and away from the
#' `dot_self(alpha_raw) == 0` degeneracy.
#'
#' Every registered model also folds `alpha_raw[1]`'s sign into the
#' REPORTED `alpha`/`theta` (`orientation_sign()`, each `.stan` file's
#' header, "IDENTIFICATION: ORIENTATION FOLD"), which makes the reported
#' quantities correct regardless of which raw-space basin a chain
#' occupies -- so this function's sign bias is a nicety, not a
#' correctness requirement: it just means a chain's RAW parameters
#' (`alpha_raw` itself, and the theta-side raws that share its sign)
#' usually won't label-switch either, keeping their own diagnostics
#' interpretable more often (though never guaranteed; see
#' `R/diagnose_convergence.R`'s exclusion of those raw names from the
#' tiered diagnostics, which is what actually protects a caller when this
#' bias doesn't hold).
#'
#' @param A Number of action types.
#' @return A length-`A` numeric vector summing to exactly 0, with its
#'   first element positive.
#' @keywords internal
.alpha_raw_sum0_init <- function(A) {
  v <- stats::rnorm(A, 0, 0.5)
  v <- v - mean(v)
  if (v[1] < 0) {
    v <- -v
  }
  v
}

#' Build an initial-value generator matching a model's parameterization
#'
#' Initial values are model-specific: `stan_model` selects among the
#' registered models' distinct parameter sets (a non-centered
#' `process_noise` hierarchy vs. an OU/AR(1) `sd_stat` hierarchy, etc.).
#' `stable_gamma` (0.7.0) gets its own branch: `stable`'s init list plus
#' `gamma_z`/`sigma_gamma`, both started small (not exactly 0 -- same
#' reasoning as `alpha_raw`'s non-zero start: an all-zero init is valid
#' but leaves no spread to move away from the prior mode with).
#' Unrecognised names are rejected upstream by [.resolve_stan_model()],
#' so the `stop()` below should be unreachable in practice.
#'
#' @param stan_data A Stan data list as returned by [assemble_stan_data()].
#' @param stan_model Name registered in `.bilatr_stan_models`.
#' @return A zero-argument function suitable for `cmdstanr`'s `init`
#'   argument.
#' @keywords internal
bilatr_init_fn <- function(stan_data, stan_model = .BILATR_DEFAULT_MODEL) {
  D <- stan_data$D
  Tn <- stan_data$T
  A <- stan_data$A

  init_list <- switch(
    stan_model,
    stable = list(
      theta_raw = matrix(0, D, Tn),
      mu_intercept = rep(0, A),
      alpha_raw = .alpha_raw_sum0_init(A),
      sigma_theta0 = 0.5,
      # real initial spread, not near 0: with theta near 0 the likelihood
      # is nearly flat in alpha's direction, leaving an early-warmup
      # window in which alpha could still rotate before the data locks
      # the orientation in
      z_theta0 = stats::rnorm(D, 0, 0.5),
      log_process_noise_raw = rep(0, D),
      mu_log_noise = log(0.2),
      sigma_log_noise = 0.3,
      phi = rep(1, D),
      mu_log_phi = 0,
      sigma_log_phi = 0.5
    ),
    stable_gamma = {
      n_countries <- stan_data$n_countries
      c(
        list(
          theta_raw = matrix(0, D, Tn),
          mu_intercept = rep(0, A),
          alpha_raw = .alpha_raw_sum0_init(A),
          sigma_theta0 = 0.5,
          z_theta0 = stats::rnorm(D, 0, 0.5),
          log_process_noise_raw = rep(0, D),
          mu_log_noise = log(0.2),
          sigma_log_noise = 0.3,
          phi = rep(1, D),
          mu_log_phi = 0,
          sigma_log_phi = 0.5
        ),
        # gamma_z near (but not at) 0: an all-zero init is valid (gamma_z
        # is unconstrained) but, like alpha_raw above, leaves the
        # projected-out directions exactly at their prior mode with no
        # spread to break out of -- a small non-zero draw is the more
        # conservative default. sigma_gamma starts small (its
        # half-normal(0, 0.3) prior means the data should pull it up if
        # warranted, not down from an over-confident large start).
        list(
          gamma_z = matrix(stats::rnorm(A * n_countries, 0, 0.1), A, n_countries),
          sigma_gamma = rep(0.1, A)
        )
      )
    },
    ou = list(
      theta_raw = matrix(0, D, Tn),
      mu_intercept = rep(0, A),
      alpha_raw = .alpha_raw_sum0_init(A),
      sigma_mu = 0.5,
      # real initial spread, not near 0 -- see the stable branch above
      mu_dyad_raw = stats::rnorm(D, 0, 0.5),
      rho = 0.8,
      mu_log_sd_stat = log(1),
      sigma_log_sd_stat = 0.3,
      log_sd_stat_raw = rep(0, D),
      phi = rep(1, D),
      mu_log_phi = 0,
      sigma_log_phi = 0.5
    ),
    stop(
      "bilatr_init_fn(): no init generator registered for stan_model '",
      stan_model, "'.",
      call. = FALSE
    )
  )

  function() {
    init_list
  }
}

#' Shared sampling logic behind fit_dyad_ts()/fit_panel() and their _dev
#' counterparts in R/fit_dev.R
#'
#' Canonicalises `stan_model` once, here, so `.compile_stan_model()` and
#' [bilatr_init_fn()] agree on the same resolved name.
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
  stan_model <- .canonical_stan_model(stan_model)
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
#' `bilatr_alphanorm.stan`, registered as the `stable` model -- see
#' `R/model_registry.R`) to a single dyad (`D == 1`), i.e. estimates
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
#'   reference_category = 0
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
#' `bilatr_alphanorm.stan`, registered as the `stable` model -- see
#' `R/model_registry.R`) across multiple dyads simultaneously (`D > 1`),
#' with hierarchical partial pooling on process noise (via
#' `mu_log_noise`/`sigma_log_noise`), the dispersion parameter `phi` (via
#' `mu_log_phi`/`sigma_log_phi`), and the initial latent state `theta0`
#' (via `sigma_theta0`; `theta0`'s population mean is pinned at exactly 0
#' rather than a separately pooled location). Use [fit_dyad_ts()] instead
#' for a single dyad.
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
#'   reference_category = 0, chunk_size = 600
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
