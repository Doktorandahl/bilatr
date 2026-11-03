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
#' @param stan_data A Stan data list as returned by [prepare_stan_data()]
#'   (or attached to a fit by [fit_panel()]/[fit_dyad_ts()]).
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

#' Shared sampling logic behind fit_dyad_ts()/fit_panel()
#'
#' Compiles (or fetches from cache) the model and samples. Canonicalises
#' `stan_model` again (idempotent), so [.compile_stan_model()] and
#' [bilatr_init_fn()] agree on the same resolved name, and emits the
#' once-per-session experimental-model message (see
#' [.bilatr_warn_if_experimental()]). Every check on the data has
#' already run in [.bilatr_fit_entry()].
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
  .bilatr_warn_if_experimental(stan_model)
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

#' The data-preparation arguments of fit_dyad_ts()/fit_panel()
#'
#' Every argument that [prepare_stan_data()] consumes. With a prepared
#' `stan_data` as `data`, supplying any of them is an error.
#' @keywords internal
.bilatr_prep_args <- c(
  "time", "action_cols", "dyad", "actor1", "actor2", "w_send", "directed",
  "reference_category", "min_n_events", "time_range", "chunk_size",
  "rho_prior_a", "rho_prior_b", "compute_log_lik", "prior_only",
  "compute_theta_filtered", "filter_dyads"
)

#' Which data-preparation arguments the caller supplied
#'
#' Uses `missing()` in the calling `fit_*()` frame, so an argument passed
#' with its default value still counts as supplied.
#'
#' @param env The calling `fit_*()` function's frame.
#' @return A character vector, a subset of [.bilatr_prep_args].
#' @keywords internal
.bilatr_supplied_prep_args <- function(env) {
  Filter(function(arg) !eval(call("missing", as.name(arg)), env), .bilatr_prep_args)
}

#' Does `x` look like a prepared Stan data list?
#'
#' A list that is not a data frame and has a `D` field. The remaining
#' required fields are checked by [.bilatr_fit_entry()] after the D check.
#'
#' @param x Anything.
#' @return `TRUE`/`FALSE`.
#' @keywords internal
.is_stan_data_list <- function(x) {
  is.list(x) && !is.data.frame(x) && "D" %in% names(x)
}

#' Write an `.rds` file atomically
#'
#' Writes to a temporary file in the same directory, then renames it into
#' place, so concurrent writers (e.g. one SLURM job per chain) never
#' leave a half-written file.
#'
#' @param object The object to save.
#' @param path The destination path.
#' @return `path`, invisibly.
#' @keywords internal
.write_rds_atomic <- function(object, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "_"), tmpdir = dirname(path), fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) {
    stop("Could not write `", path, "`.", call. = FALSE)
  }
  invisible(path)
}

#' Shared entry point of fit_dyad_ts()/fit_panel()
#'
#' Resolves `data` to a Stan data list, runs every check that does not
#' need the compiled model, writes the `stan_data` sidecar file, samples
#' via [fit_bilatr()], and attaches `stan_data` to the fit. In order:
#'
#' 1. canonicalises `stan_model`;
#' 2. dispatches on `data`: a prepared `stan_data` list is used as-is
#'    (error if any data-preparation argument was supplied); a data frame
#'    goes through [prepare_stan_data()]; anything else is an error;
#' 3. the D check (`D == 1` for `fit_dyad_ts()`, `D >= 2` for
#'    `fit_panel()`), then that a `stan_data` list has `T`, `A`, `Y` and
#'    `is_obs`;
#' 4. for `stable_gamma`: errors if the country fields are absent, and
#'    warns if `w_send` was imputed.
#'
#' @param caller `"fit_dyad_ts()"` or `"fit_panel()"`.
#' @param data The `data` argument.
#' @param prep Named list of the data-preparation arguments, with
#'   `action_cols` as a quosure.
#' @param supplied The data-preparation arguments the caller supplied.
#' @param sampler Named list of the sampler arguments of [fit_bilatr()].
#' @param stan_model The `stan_model` argument.
#' @param ... Passed to `CmdStanModel$sample()`.
#' @return The `CmdStanMCMC` fit, with `attr(fit, "stan_data")` set.
#' @keywords internal
.bilatr_fit_entry <- function(caller, data, prep, supplied, sampler, stan_model, ...) {
  # 1. stan_model
  stan_model <- .canonical_stan_model(stan_model)

  # 2. data
  n_dyads_in <- NULL
  if (.is_stan_data_list(data)) {
    if (length(supplied) > 0) {
      stop(
        caller, ": `data` is an already prepared stan_data list, so the ",
        "data-preparation argument(s) ", paste0("`", supplied, "`", collapse = ", "),
        " cannot be applied. Pass them to prepare_stan_data() instead, or pass ",
        "the dyad-period table itself as `data`.",
        call. = FALSE
      )
    }
    stan_data <- data
  } else if (is.data.frame(data)) {
    n_dyads_in <- .n_dyads_in_table(data, prep$dyad, prep$actor1, prep$actor2)
    stan_data <- prepare_stan_data(
      data,
      time = prep$time,
      action_cols = !!prep$action_cols,
      dyad = prep$dyad,
      actor1 = prep$actor1,
      actor2 = prep$actor2,
      w_send = prep$w_send,
      directed = prep$directed,
      reference_category = prep$reference_category,
      min_n_events = prep$min_n_events,
      time_range = prep$time_range,
      chunk_size = prep$chunk_size,
      rho_prior_a = prep$rho_prior_a,
      rho_prior_b = prep$rho_prior_b,
      compute_log_lik = prep$compute_log_lik,
      prior_only = prep$prior_only,
      compute_theta_filtered = prep$compute_theta_filtered,
      filter_dyads = prep$filter_dyads
    )
  } else {
    stop(
      caller, ": `data` must be a dyad-period data frame (see ?bilatr_dyad_period_data) ",
      "or a stan_data list from prepare_stan_data(); got an object of class ",
      paste0("<", class(data), ">", collapse = "/"), ".",
      call. = FALSE
    )
  }

  # 3. D, then the remaining required fields
  D <- stan_data$D
  filtered_note <- if (!is.null(n_dyads_in) && n_dyads_in != D) {
    sprintf(
      " after `min_n_events`/`time_range` filtering (`data` has %d dyad(s); `min_n_events` = %s)",
      n_dyads_in, format(prep$min_n_events)
    )
  } else {
    ""
  }
  if (caller == "fit_dyad_ts()" && D != 1) {
    stop(
      "fit_dyad_ts() expects a single-dyad stan_data (D == 1), got D = ", D,
      filtered_note, ". Use fit_panel() for multiple dyads.",
      call. = FALSE
    )
  }
  if (caller == "fit_panel()" && D < 2) {
    stop(
      "fit_panel() expects multiple dyads (D >= 2), got D = ", D,
      filtered_note, ". Use fit_dyad_ts() for a single dyad.",
      call. = FALSE
    )
  }
  missing_fields <- setdiff(c("D", "T", "A", "Y", "is_obs"), names(stan_data))
  if (length(missing_fields) > 0) {
    stop(
      caller, ": `data` looks like a stan_data list but is missing field(s) ",
      paste0("`", missing_fields, "`", collapse = ", "),
      "; build it with prepare_stan_data().",
      call. = FALSE
    )
  }

  # 4. stable_gamma
  if (.bilatr_model_has_gamma(stan_model)) {
    gamma_fields <- c("n_countries", "ctry_a", "ctry_b", "w_send")
    absent <- setdiff(gamma_fields, names(stan_data))
    if (length(absent) > 0) {
      stop(
        caller, ": stan_model = \"stable_gamma\" needs the country fields ",
        paste0("`", absent, "`", collapse = ", "),
        ", which prepare_stan_data() builds only from actor columns. Supply ",
        "the `actor1`/`actor2` columns in `data` (or name them with the ",
        "`actor1`/`actor2` arguments).",
        call. = FALSE
      )
    }
    if (isTRUE(attr(stan_data, "w_send_imputed"))) {
      warning(
        caller, ": `w_send` was not supplied, so 0.5 is assumed for every ",
        "undirected dyad in stable_gamma's country offset (each actor sends ",
        "half of the dyad's events). aggregate_events() computes the actual ",
        "share from event data.",
        call. = FALSE
      )
    }
  }

  # sidecar, before sampling, so a crashed run still leaves it
  output_dir <- sampler$output_dir
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    basename <- list(...)$output_basename %||% "bilatr"
    .write_rds_atomic(stan_data, file.path(output_dir, paste0(basename, "_stan_data.rds")))
  }

  fit <- fit_bilatr(
    stan_data,
    chains = sampler$chains,
    parallel_chains = sampler$parallel_chains,
    threads_per_chain = sampler$threads_per_chain,
    iter_warmup = sampler$iter_warmup,
    iter_sampling = sampler$iter_sampling,
    seed = sampler$seed,
    opt_level = sampler$opt_level,
    output_dir = output_dir,
    stan_model = stan_model,
    ...
  )
  attr(fit, "stan_data") <- stan_data
  fit
}

#' Number of distinct dyads in a dyad-period table
#'
#' For the "after filtering" note in the D-check error. Mirrors how
#' [prepare_stan_data()] builds the dyad key.
#'
#' @param data,dyad,actor1,actor2 As in [prepare_stan_data()].
#' @return An integer, or `NULL` if the dyad key cannot be built.
#' @keywords internal
.n_dyads_in_table <- function(data, dyad, actor1, actor2) {
  if (rlang::is_string(dyad) && dyad %in% names(data)) {
    return(length(unique(data[[dyad]])))
  }
  if (rlang::is_string(actor1) && rlang::is_string(actor2) && all(c(actor1, actor2) %in% names(data))) {
    return(nrow(unique(data[c(actor1, actor2)])))
  }
  NULL
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
#' `data` is either a dyad-period table (see `?bilatr_dyad_period_data`;
#' e.g. the output of [aggregate_events()]), which is passed to
#' [prepare_stan_data()] together with the data-preparation arguments
#' (`time` through `filter_dyads`), or a Stan data list already built by
#' [prepare_stan_data()] (or by `assemble_stan_data()` in bilatr 0.11.0
#' and earlier), which is used as-is. With a prepared list, supplying any
#' data-preparation argument is an error.
#'
#' Before compiling anything, `fit_dyad_ts()` checks that the data has
#' exactly one dyad (after `min_n_events` filtering), and, for
#' `stan_model = "stable_gamma"`, that the country fields exist (the
#' table needs actor columns) and warns if `w_send` was imputed as 0.5
#' for undirected dyads.
#'
#' `reduce_sum` chunking across dyads has nothing to parallelize over
#' when `D == 1`, so `threads_per_chain` defaults to 1 here (raising it
#' will not speed up sampling for a single-dyad fit).
#'
#' @param data A dyad-period data frame (see `?bilatr_dyad_period_data`)
#'   or a Stan data list from [prepare_stan_data()]; see Details.
#' @inheritParams prepare_stan_data
#' @param directed Whether the dyads are directed; see
#'   [prepare_stan_data()]. Defaults to `attr(data, "directed")` if set
#'   (by [aggregate_events()]), else `TRUE`.
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
#'   `NULL` for a temporary directory. When given, the Stan data list is
#'   also saved there (see Value).
#' @param stan_model (0.10.1) Name registered in `.bilatr_stan_models`
#'   (see [.canonical_stan_model()]); an unregistered name errors,
#'   including the four names retired in 0.10.0. Defaults to
#'   `"stable"` (`.BILATR_DEFAULT_MODEL`). Fitting a model whose registry
#'   `status` is `"experimental"` (currently `ou`/`stable_gamma`) emits a
#'   `message()` the first time it is fit in a session (see
#'   [.bilatr_warn_if_experimental()]).
#' @param ... Additional arguments passed to `CmdStanModel$sample()`
#'   (e.g. `refresh`, `chain_ids`, `output_basename`).
#' @return A `CmdStanMCMC` fit object, with the Stan data list it was fit
#'   to attached as `attr(fit, "stan_data")` (with its own attributes:
#'   `dyad_ids`, `event_classes`, ...). The extractors ([extract_theta()],
#'   [extract_gamma()], [diagnose_category_merges()],
#'   [check_compositional_residuals()], ...) use it when their own
#'   `stan_data` argument is not given. `CmdStanMCMC` objects are
#'   environments, so the attribute travels with every reference to the
#'   fit, and survives `fit$save_object()`/`saveRDS()` and `readRDS()`.
#'
#'   When `output_dir` is given, the Stan data list is also written,
#'   before sampling starts, to
#'   `file.path(output_dir, paste0(<basename>, "_stan_data.rds"))`, where
#'   `<basename>` is the `output_basename` passed through `...`, else
#'   `"bilatr"`. The write is atomic (a temporary file, then a rename), so
#'   concurrent per-chain jobs sharing an `output_dir` are safe. The
#'   extractors find this file next to the CmdStan CSVs when given CSV
#'   paths.
#' @examples
#' \dontrun{
#' agg <- aggregate_events(
#'   dplyr::filter(events, Actor1CountryCode == "USA", Actor2CountryCode == "CHN"),
#'   grouping_var = "PentaClass", start = "2015-01-01", end = "2020-12-31"
#' )
#' fit <- fit_dyad_ts(agg, reference_category = "0", chains = 4, iter_sampling = 1000)
#' theta <- extract_theta(fit)
#' }
#' @export
fit_dyad_ts <- function(
  data,
  time = NULL,
  action_cols = NULL,
  dyad = "dyad",
  actor1 = "actor1",
  actor2 = "actor2",
  w_send = "w_send",
  directed = attr(data, "directed") %||% TRUE,
  reference_category = NULL,
  min_n_events = 1,
  time_range = NULL,
  chunk_size = 100,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 3,
  output_dir = NULL,
  stan_model = .BILATR_DEFAULT_MODEL,
  ...
) {
  supplied <- .bilatr_supplied_prep_args(environment())
  .bilatr_fit_entry(
    "fit_dyad_ts()", data,
    prep = .bilatr_collect_prep(environment(), rlang::enquo(action_cols)),
    supplied = supplied,
    sampler = list(
      chains = chains, parallel_chains = parallel_chains,
      threads_per_chain = threads_per_chain, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, seed = seed, opt_level = opt_level,
      output_dir = output_dir
    ),
    stan_model = stan_model,
    ...
  )
}

#' Collect the data-preparation arguments of a fit_*() call
#'
#' Evaluates them in the `fit_*()` frame, except `action_cols`, which is
#' kept as a quosure so tidyselect expressions reach
#' [prepare_stan_data()] unevaluated. Nothing is forced when `data` is a
#' prepared Stan data list and no argument was supplied, so the
#' `directed` default never looks at a list's attributes.
#'
#' @param env The `fit_*()` frame.
#' @param action_quo The `action_cols` quosure.
#' @return A named list.
#' @keywords internal
.bilatr_collect_prep <- function(env, action_quo) {
  if (.is_stan_data_list(get("data", envir = env))) {
    return(list())
  }
  args <- setdiff(.bilatr_prep_args, "action_cols")
  prep <- rlang::set_names(lapply(args, function(a) get(a, envir = env)), args)
  prep$action_cols <- action_quo
  prep
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
#' `data` is either a dyad-period table or a prepared Stan data list,
#' exactly as for [fit_dyad_ts()]; see its Details for the dispatch and
#' the checks run before compiling. `fit_panel()` requires `D >= 2`.
#'
#' `threads_per_chain` controls `reduce_sum` parallelization across
#' dyads within a chain, combined with `chunk_size` (`C` in the Stan
#' data). 16 threads with `chunk_size = 600` was the Pareto-optimal
#' setting found in this project's own benchmarking for panel-sized data;
#' tune both together for your data size and hardware rather than one in
#' isolation.
#'
#' @inheritParams fit_dyad_ts
#' @param threads_per_chain Threads per chain for `reduce_sum`
#'   parallelization across dyads. Defaults to 16, the setting found to
#'   be Pareto-optimal (together with `chunk_size = 600`) in this
#'   project's own benchmarking; reduce it on machines with fewer cores.
#' @return As for [fit_dyad_ts()]: a `CmdStanMCMC` fit object with
#'   `attr(fit, "stan_data")` set, and the Stan data list saved to
#'   `output_dir` when one is given.
#' @examples
#' \dontrun{
#' agg <- aggregate_events(
#'   events, grouping_var = "PentaClass", resolution = "yearly",
#'   start = "2015-01-01", end = "2020-12-31"
#' )
#' fit <- fit_panel(
#'   agg, reference_category = "0", chunk_size = 600,
#'   chains = 4, threads_per_chain = 16, iter_sampling = 1000
#' )
#' theta <- extract_theta(fit)
#' }
#' @export
fit_panel <- function(
  data,
  time = NULL,
  action_cols = NULL,
  dyad = "dyad",
  actor1 = "actor1",
  actor2 = "actor2",
  w_send = "w_send",
  directed = attr(data, "directed") %||% TRUE,
  reference_category = NULL,
  min_n_events = 1,
  time_range = NULL,
  chunk_size = 100,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 16,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 3,
  output_dir = NULL,
  stan_model = .BILATR_DEFAULT_MODEL,
  ...
) {
  supplied <- .bilatr_supplied_prep_args(environment())
  .bilatr_fit_entry(
    "fit_panel()", data,
    prep = .bilatr_collect_prep(environment(), rlang::enquo(action_cols)),
    supplied = supplied,
    sampler = list(
      chains = chains, parallel_chains = parallel_chains,
      threads_per_chain = threads_per_chain, iter_warmup = iter_warmup,
      iter_sampling = iter_sampling, seed = seed, opt_level = opt_level,
      output_dir = output_dir
    ),
    stan_model = stan_model,
    ...
  )
}
