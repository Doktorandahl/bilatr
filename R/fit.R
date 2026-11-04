#' Cache of compiled CmdStanModel objects
#'
#' Keyed by "<copied .stan path>::<opt_level>::<stan_threads>"; see
#' [.compile_stan_model()]. Emptied by [bilatr_clear_model_cache()].
#' @noRd
.bilatr_model_cache <- new.env(parent = emptyenv())

#' Oldest supported cmdstanr version
#'
#' Kept equal to the `cmdstanr (>= x.y.z)` entry in `DESCRIPTION`'s
#' Suggests. bilatr's own calls (`cmdstan_version(error_on_NA = )`,
#' `cmdstan_model(cpp_options = , force_recompile = )`, `$sample()`,
#' `$draws()`, `$metadata()`) need 0.4.0, the tests' `$init_model_methods()`
#' needs 0.6.0, and 0.8.0 is the first release compatible with CmdStan
#' 2.34 and 2.35, and the oldest bilatr has been checked with.
#' @noRd
.BILATR_MIN_CMDSTANR <- "0.8.0"

#' Check that cmdstanr and CmdStan are available
#'
#' `cmdstanr` is in Suggests, so every code path that compiles or samples
#' calls this first. Requires at least [.BILATR_MIN_CMDSTANR]. Shares its
#' checks with [bilatr_cmdstan_available()] through [.cmdstan_problem()].
#'
#' @param caller The calling function's name, for the error message.
#' @return `TRUE`, invisibly; errors otherwise.
#' @noRd
.require_cmdstan <- function(caller) {
  problem <- .cmdstan_problem()
  if (identical(problem, "cmdstanr")) {
    # errors (or, interactively, offers to install)
    rlang::check_installed(
      "cmdstanr",
      version = .BILATR_MIN_CMDSTANR,
      reason = "to compile and fit bilatr models."
    )
    problem <- .cmdstan_problem()
  }
  if (identical(problem, "cmdstan")) {
    stop(
      caller, ": CmdStan is not installed, or cmdstanr cannot find it. ",
      "Install it with cmdstanr::install_cmdstan(), or point cmdstanr to an ",
      "existing installation with cmdstanr::set_cmdstan_path().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' The user cache directory of bilatr
#'
#' `tools::R_user_dir("bilatr", "cache")`; a function so tests can point
#' it elsewhere.
#' @return A directory path.
#' @noRd
.bilatr_cache_root <- function() {
  tools::R_user_dir("bilatr", "cache")
}

#' Directory that compiled bilatr models are cached in
#'
#' `getOption("bilatr.model_dir")` if set, otherwise a per-version
#' directory under the user cache directory,
#' `file.path(tools::R_user_dir("bilatr", "cache"), "models",
#' <package version>)`, so an upgrade never reuses an executable compiled
#' by an older version. Created if it does not exist. Nothing is ever
#' written into the installed package.
#'
#' @param create If `TRUE` (default), create the directory.
#' @return A directory path.
#' @noRd
.bilatr_model_dir <- function(create = TRUE) {
  dir <- getOption("bilatr.model_dir")
  if (is.null(dir)) {
    dir <- file.path(.bilatr_cache_root(), "models", as.character(utils::packageVersion("bilatr")))
  } else if (!rlang::is_string(dir) || !nzchar(dir)) {
    stop("`options(bilatr.model_dir = )` must be a single directory path.", call. = FALSE)
  }
  if (create && !dir.exists(dir)) {
    ok <- dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(dir)) {
      stop(
        "Could not create the model cache directory `", dir, "`. Set ",
        "`options(bilatr.model_dir = )` to a writable directory.",
        call. = FALSE
      )
    }
  }
  dir
}

#' Copy a registered model's `.stan` file into the model cache directory
#'
#' Copies only when the target is missing or its MD5 sum differs from the
#' installed source, so an unchanged file keeps its modification time and
#' `cmdstanr` keeps reusing the compiled executable next to it.
#'
#' With `profile = TRUE`, writes instead the program returned by
#' [insert_reduce_sum_profile()] for that model, under its own name
#' (`<file>_profile.stan`, e.g. `bilatr_alphanorm_profile.stan`), so the
#' profiled and unprofiled models are cached side by side. It is
#' rewritten only when its content changes.
#'
#' @param stan_model A canonical registered name.
#' @param dir The cache directory, from [.bilatr_model_dir()].
#' @param profile Whether to stage the profiled copy.
#' @return The path of the copy.
#' @noRd
.bilatr_stage_stan_file <- function(stan_model, dir, profile = FALSE) {
  src <- .resolve_stan_model(stan_model)
  if (profile) {
    target <- file.path(dir, paste0(tools::file_path_sans_ext(basename(src)), "_profile.stan"))
    code <- insert_reduce_sum_profile(paste(readLines(src, warn = FALSE), collapse = "\n"))
    current <- if (file.exists(target)) paste(readLines(target, warn = FALSE), collapse = "\n") else NULL
    if (!identical(current, code)) {
      writeLines(code, target)
      stale <- startsWith(ls(.bilatr_model_cache), paste0(target, "::"))
      rm(list = ls(.bilatr_model_cache)[stale], envir = .bilatr_model_cache)
    }
    return(target)
  }
  target <- file.path(dir, basename(src))
  if (!file.exists(target) || !identical(unname(tools::md5sum(src)), unname(tools::md5sum(target)))) {
    if (!file.copy(src, target, overwrite = TRUE)) {
      stop("Could not copy `", src, "` to the model cache directory `", dir, "`.", call. = FALSE)
    }
    # an in-session object compiled from the old file is stale now
    stale <- startsWith(ls(.bilatr_model_cache), paste0(target, "::"))
    rm(list = ls(.bilatr_model_cache)[stale], envir = .bilatr_model_cache)
  }
  target
}

#' Compile (or fetch from cache) a registered Stan model
#'
#' Shared compilation path behind [compile_bilatr_model()] and
#' [fit_panel()]/[fit_dyad_ts()]. Copies the model's `.stan` file into
#' [.bilatr_model_dir()] and compiles it there, so the executable is
#' written next to the copy, never into the installed package. Compiled
#' `CmdStanModel` objects are also kept in `.bilatr_model_cache` for the
#' session, keyed by the copied path, `opt_level` and `stan_threads`.
#'
#' @param stan_model A registered model name.
#' @param opt_level C++ compiler optimization level; see
#'   [compile_bilatr_model()].
#' @param force_recompile If `TRUE`, bypass both caches and rebuild.
#' @param stan_threads Whether to compile with `STAN_THREADS`. Always
#'   `TRUE` for the package's own fits.
#' @param profile Whether to compile the copy whose `reduce_sum`
#'   likelihood is profiled (see [.bilatr_stage_stan_file()]), as
#'   [benchmark_bilatr()] does.
#' @return A `CmdStanModel` object.
#' @noRd
.compile_stan_model <- function(stan_model, opt_level, force_recompile = FALSE, stan_threads = TRUE,
                                profile = FALSE) {
  .require_cmdstan("compile_bilatr_model()")
  stan_model <- .canonical_stan_model(stan_model)
  stan_file <- .bilatr_stage_stan_file(stan_model, .bilatr_model_dir(), profile = profile)
  cache_key <- paste(stan_file, opt_level, stan_threads, sep = "::")

  if (!force_recompile && exists(cache_key, envir = .bilatr_model_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .bilatr_model_cache, inherits = FALSE))
  }

  mod <- cmdstanr::cmdstan_model(
    stan_file,
    cpp_options = list(stan_threads = stan_threads, O = opt_level),
    force_recompile = force_recompile
  )
  assign(cache_key, mod, envir = .bilatr_model_cache)
  mod
}

#' Compile a bilatr Stan model
#'
#' Compiles one of the package's Stan models with CmdStan, or returns the
#' already compiled model. [fit_panel()] and [fit_dyad_ts()] call this
#' themselves, so you only need it to compile ahead of time (for example,
#' before submitting batch jobs) or with a non-default `opt_level`.
#'
#' @section Models:
#' Three models are available:
#'
#' `r .bilatr_models_rd()`
#'
#' The experimental models have not yet been calibrated or tested at
#' scale as thoroughly as `"stable"`; fitting one shows a message once per
#' session.
#'
#' @section Where compiled models are stored:
#' The `.stan` file is copied into a cache directory and compiled there:
#' by default `file.path(tools::R_user_dir("bilatr", "cache"), "models",
#' <bilatr version>)`, so each package version compiles its own
#' executables. Set `options(bilatr.model_dir = "<directory>")` to use
#' another directory, for example one shared by several batch jobs. Compiled models are reused across sessions until the `.stan` file
#' changes. [bilatr_clear_model_cache()] deletes them.
#'
#' Compiling needs the 'cmdstanr' package and a working CmdStan
#' installation (see `cmdstanr::install_cmdstan()`).
#'
#' @param stan_model Which model to compile: `"stable"` (default), `"ou"`
#'   or `"stable_gamma"`. See the Models section.
#' @param opt_level C++ compiler optimization level passed to CmdStan's
#'   `O` Makefile variable (0-3). The reduce_sum-based likelihood in this
#'   model generates a lot of C++ template code; at the default `O = 3`,
#'   compilation can use several GB of RAM and occasionally exhausts
#'   memory on constrained machines (laptops, small CI runners). Dropping
#'   to `opt_level = 1` cuts compilation RAM substantially at the cost of
#'   somewhat slower sampling; `opt_level = 0` is fastest to compile but
#'   slowest to sample and rarely worth it. If compilation is being
#'   killed by the OOM killer, try `opt_level = 1` before anything else.
#' @param force_recompile Set `TRUE` to force a rebuild, bypassing both
#'   the in-session cache and the compiled executable on disk.
#' @return A `CmdStanModel` object.
#' @seealso [bilatr_clear_model_cache()]
#' @examplesIf bilatr_cmdstan_available()
#' \donttest{
#' # compiling takes a minute or more the first time; later calls reuse
#' # the compiled model
#' mod <- compile_bilatr_model()
#' mod_gamma <- compile_bilatr_model("stable_gamma")
#'
#' }
#' # less memory for compiling, somewhat slower sampling (rebuilds, so it
#' # takes as long as the first compile):
#' # mod_low_ram <- compile_bilatr_model(opt_level = 1, force_recompile = TRUE)
#' @export
compile_bilatr_model <- function(stan_model = "stable", opt_level = 3, force_recompile = FALSE) {
  stan_model <- .canonical_stan_model(stan_model)
  .compile_stan_model(stan_model, opt_level, force_recompile)
}

#' Delete compiled bilatr models
#'
#' Deletes the compiled models that [compile_bilatr_model()],
#' [fit_panel()] and [fit_dyad_ts()] cache on disk (see
#' [compile_bilatr_model()], "Where compiled models are stored"), and
#' forgets the models compiled in this session. The next fit compiles
#' again.
#'
#' By default the cache directory of the installed bilatr version is
#' deleted. If `options(bilatr.model_dir = )` is set, only the files
#' bilatr put in that directory (the `.stan` copies, the profiled copies
#' that [benchmark_bilatr()] compiles, and their compiled executables) are
#' deleted, not the directory itself.
#'
#' @param all If `TRUE`, also delete the cached models of every other
#'   bilatr version under `tools::R_user_dir("bilatr", "cache")`.
#' @return The deleted paths, invisibly.
#' @seealso [compile_bilatr_model()]
#' @examples
#' old <- options(bilatr.model_dir = tempfile("bilatr_models_"))
#' bilatr_clear_model_cache()
#' options(old)
#' @export
bilatr_clear_model_cache <- function(all = FALSE) {
  if (!rlang::is_bool(all)) {
    stop("bilatr_clear_model_cache(): `all` must be TRUE or FALSE.", call. = FALSE)
  }
  rm(list = ls(.bilatr_model_cache), envir = .bilatr_model_cache)

  dir <- .bilatr_model_dir(create = FALSE)
  targets <- if (is.null(getOption("bilatr.model_dir"))) {
    dir
  } else {
    stems <- tools::file_path_sans_ext(vapply(.bilatr_stan_models, `[[`, "", "file"))
    pattern <- paste0("^(", paste(stems, collapse = "|"), ")(_profile)?(\\.(stan|exe|hpp|o))?$")
    list.files(dir, pattern = pattern, full.names = TRUE)
  }
  if (all) {
    targets <- c(targets, file.path(.bilatr_cache_root(), "models"))
  }
  targets <- unique(targets[file.exists(targets)])
  unlink(targets, recursive = TRUE)
  invisible(targets)
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
#' header, "Sign: the ORIENTATION FOLD"), which makes the reported
#' quantities correct regardless of which raw-space basin a chain
#' occupies -- so this function's sign bias is a nicety, not a
#' correctness requirement: it just means a chain's RAW parameters
#' (`alpha_raw` itself, and the theta-side raws that share its sign)
#' usually won't label-switch either, keeping their own diagnostics
#' interpretable more often (though never guaranteed; the exclusion of
#' those raw names from the tiered diagnostics is what actually protects
#' a caller when this bias doesn't hold).
#'
#' @param A Number of action types.
#' @return A length-`A` numeric vector summing to exactly 0, with its
#'   first element positive.
#' @noRd
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
#' `stable_gamma` gets `stable`'s init list plus `gamma_z`/`sigma_gamma`,
#' both started small (not exactly 0 -- same reasoning as `alpha_raw`'s
#' non-zero start: an all-zero init is valid but leaves no spread to move
#' away from the prior mode with).
#'
#' The returned function draws fresh random values (`alpha_raw`,
#' `z_theta0`/`mu_dyad_raw`, `gamma_z`) on every call, so calling it once
#' per chain gives each chain its own starting point. It draws from the
#' global random-number stream; [.bilatr_draw_inits()] seeds each call
#' when the fit has a `seed`.
#'
#' @param stan_data A Stan data list as returned by [prepare_stan_data()]
#'   (or attached to a fit by [fit_panel()]/[fit_dyad_ts()]).
#' @param stan_model A registered model name.
#' @return A zero-argument function returning a named list of initial
#'   values.
#' @noRd
bilatr_init_fn <- function(stan_data, stan_model = .BILATR_DEFAULT_MODEL) {
  stan_model <- .canonical_stan_model(stan_model)
  D <- stan_data$D
  Tn <- stan_data$T
  A <- stan_data$A
  n_countries <- stan_data$n_countries

  # stable's parameters, shared by stable_gamma
  stable_init <- function() {
    list(
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
    )
  }

  switch(
    stan_model,
    stable = stable_init,
    stable_gamma = function() {
      c(
        stable_init(),
        # gamma_z near (but not at) 0: an all-zero init is valid (gamma_z
        # is unconstrained) but, like alpha_raw above, leaves the
        # projected-out directions exactly at their prior mode with no
        # spread to break out of. sigma_gamma starts small (its
        # half-normal(0, 0.3) prior means the data should pull it up if
        # warranted, not down from an over-confident large start).
        list(
          gamma_z = matrix(stats::rnorm(A * n_countries, 0, 0.1), A, n_countries),
          sigma_gamma = rep(0.1, A)
        )
      )
    },
    ou = function() {
      list(
        theta_raw = matrix(0, D, Tn),
        mu_intercept = rep(0, A),
        alpha_raw = .alpha_raw_sum0_init(A),
        sigma_mu = 0.5,
        # real initial spread, not near 0 -- see stable_init() above
        mu_dyad_raw = stats::rnorm(D, 0, 0.5),
        rho = 0.8,
        mu_log_sd_stat = log(1),
        sigma_log_sd_stat = 0.3,
        log_sd_stat_raw = rep(0, D),
        phi = rep(1, D),
        mu_log_phi = 0,
        sigma_log_phi = 0.5
      )
    },
    stop(
      "bilatr_init_fn(): no init generator registered for stan_model '",
      stan_model, "'.",
      call. = FALSE
    )
  )
}

#' Draw one list of initial values per chain
#'
#' With a `seed`, chain `k`'s values are drawn under
#' `.with_seed(seed + chain_ids[k], ...)`, so they are reproducible and
#' the caller's own random-number stream is left untouched. With
#' `seed = NULL`, they are drawn from the global stream.
#'
#' @param stan_data,stan_model As for [bilatr_init_fn()].
#' @param chains Number of chains.
#' @param seed The fit's `seed`, or `NULL`.
#' @param chain_ids The `chain_ids` passed to `CmdStanModel$sample()`, or
#'   `NULL` for `seq_len(chains)`.
#' @return A list of `chains` named lists, as `cmdstanr`'s `init` accepts.
#' @noRd
.bilatr_draw_inits <- function(stan_data, stan_model, chains, seed, chain_ids = NULL) {
  init_fn <- bilatr_init_fn(stan_data, stan_model)
  ids <- chain_ids %||% seq_len(chains)
  if (length(ids) != chains) {
    stop(
      "`chain_ids` has length ", length(ids), " but `chains` is ", chains, ".",
      call. = FALSE
    )
  }
  lapply(ids, function(id) {
    if (is.null(seed)) init_fn() else .with_seed(seed + id, init_fn())
  })
}

#' Shared sampling logic behind fit_dyad_ts()/fit_panel()
#'
#' Compiles (or fetches from cache) the model and samples, after emitting
#' the once-per-session experimental-model message (see
#' [.bilatr_warn_if_experimental()]). Every check on the data, and the
#' initial values, are handled by [.bilatr_fit_entry()].
#'
#' @param init A list of per-chain initial-value lists, or whatever the
#'   user passed as `init`.
#' @noRd
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
  init = NULL,
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
    init = init,
    output_dir = output_dir,
    ...
  )
}

#' The data-preparation arguments of fit_dyad_ts()/fit_panel()
#'
#' Every argument that [prepare_stan_data()] consumes. With a prepared
#' `stan_data` as `data`, supplying any of them is an error.
#' @noRd
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
#' @noRd
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
#' @noRd
.is_stan_data_list <- function(x) {
  is.list(x) && !is.data.frame(x) && "D" %in% names(x)
}

#' Write an `.rds` file atomically
#'
#' Writes to a temporary file in the same directory, then renames it into
#' place, so concurrent writers (e.g. one batch job per chain) never
#' leave a half-written file.
#'
#' @param object The object to save.
#' @param path The destination path.
#' @return `path`, invisibly.
#' @noRd
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
#' 1. canonicalizes `stan_model`;
#' 2. dispatches on `data`: a prepared `stan_data` list is used as-is
#'    (error if any data-preparation argument was supplied); a data frame
#'    goes through [prepare_stan_data()]; anything else is an error;
#' 3. the D check (`D == 1` for `fit_dyad_ts()`, `D >= 2` for
#'    `fit_panel()`), then that a `stan_data` list has `T`, `A`, `Y` and
#'    `is_obs`;
#' 4. for `stable_gamma`: errors if the country fields are absent, and
#'    warns if `w_send` was imputed;
#' 5. records the model as `attr(stan_data, "stan_model")` (replacing,
#'    with a `message()`, a different value already on a prepared list);
#' 6. draws one list of initial values per chain with
#'    [.bilatr_draw_inits()], unless `init` was passed through `...`.
#'
#' @param caller `"fit_dyad_ts()"` or `"fit_panel()"`.
#' @param data The `data` argument.
#' @param prep Named list of the data-preparation arguments, with
#'   `action_cols` as a quosure.
#' @param supplied The data-preparation arguments the caller supplied.
#' @param sampler Named list of the sampler arguments of [fit_bilatr()].
#' @param stan_model The `stan_model` argument.
#' @param ... Passed to `CmdStanModel$sample()`.
#' @return The `CmdStanMCMC` fit, with `attr(fit, "stan_data")` and
#'   `attr(fit, "stan_model")` set.
#' @noRd
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

  # 5. record the model, so the sidecar and the fit both carry it
  recorded <- attr(stan_data, "stan_model", exact = TRUE)
  if (!is.null(recorded) && !identical(recorded, stan_model)) {
    message(
      caller, ": `data` records stan_model = \"", recorded, "\"; replacing it with \"",
      stan_model, "\", the model being fit."
    )
  }
  attr(stan_data, "stan_model") <- stan_model

  # 6. initial values: one list per chain, seeded when `seed` is given.
  # A user-supplied `init` (through `...`) is used as given.
  dots <- list(...)
  if ("init" %in% names(dots)) {
    init <- dots$init
    dots$init <- NULL
  } else {
    init <- .bilatr_draw_inits(stan_data, stan_model, sampler$chains, sampler$seed, dots$chain_ids)
  }

  # sidecar, before sampling, so a crashed run still leaves it
  output_dir <- sampler$output_dir
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    basename <- dots$output_basename %||% "bilatr"
    .write_rds_atomic(stan_data, file.path(output_dir, paste0(basename, "_stan_data.rds")))
  }

  fit <- rlang::exec(
    fit_bilatr,
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
    init = init,
    !!!dots
  )
  attr(fit, "stan_data") <- stan_data
  attr(fit, "stan_model") <- stan_model
  fit
}

#' Number of distinct dyads in a dyad-period table
#'
#' For the "after filtering" note in the D-check error. Mirrors how
#' [prepare_stan_data()] builds the dyad key.
#'
#' @param data,dyad,actor1,actor2 As in [prepare_stan_data()].
#' @return An integer, or `NULL` if the dyad key cannot be built.
#' @noRd
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
#' Fits the dynamic IRT model (by default the `"stable"` model; see
#' [compile_bilatr_model()] for the available models) to a single dyad
#' (`D == 1`), i.e. estimates
#' one latent bilateral-relations trajectory `theta` with no cross-dyad pooling on
#' `theta0`, `phi`, or process noise (those hierarchical parameters still
#' exist in the model but are estimated from a single unit, so their
#' pooling has no effect). Use [fit_panel()] instead when fitting more
#' than one dyad, to get the benefit of partial pooling.
#'
#' `data` is either a dyad-period table (see `?bilatr_dyad_period_data`;
#' e.g. the output of [aggregate_events()]), which is passed to
#' [prepare_stan_data()] together with the data-preparation arguments
#' (`time` through `filter_dyads`), or a Stan data list already built by
#' [prepare_stan_data()], which is used as-is. With a prepared list,
#' supplying any data-preparation argument is an error.
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
#' **Reproducibility.** Each chain starts from its own random initial
#' values. With a `seed`, chain `k`'s initial values are drawn with seed
#' `seed + k` (`seed + chain_ids[k]` if `chain_ids` is passed through
#' `...`), and `seed` is also passed to CmdStan, so two calls with the
#' same `seed` and data give identical draws (on the same machine and
#' CmdStan version). With `seed = NULL`, the initial values come from R's
#' current random-number stream. Pass `init` through `...` to supply your
#' own initial values instead.
#'
#' The fitted model's name is recorded as `attr(fit, "stan_model")` and
#' `attr(attr(fit, "stan_data"), "stan_model")`, so the post-fit functions
#' ([diagnose_convergence()], [check_compositional_residuals()],
#' [icc_curves()], ...) know which model was fit without being told.
#'
#' @param data A dyad-period data frame (see `?bilatr_dyad_period_data`)
#'   or a Stan data list from [prepare_stan_data()]; see Details.
#' @inheritParams prepare_stan_data
#' @param directed Whether the dyads are directed; see
#'   [prepare_stan_data()]. Defaults to `attr(data, "directed")` if set
#'   (by [aggregate_events()]), else `TRUE`.
#' @param chains Number of MCMC chains.
#' @param parallel_chains Number of chains to run in parallel. Defaults
#'   to 1 (chains run one after another); set it to `chains` to run them
#'   all at once if you have the cores.
#' @param threads_per_chain Threads per chain for within-chain
#'   parallelization. Not useful for `D == 1`; left as an argument mainly
#'   for consistency with [fit_panel()].
#' @param iter_warmup Number of warmup iterations per chain.
#' @param iter_sampling Number of post-warmup sampling iterations per
#'   chain.
#' @param seed Random seed for CmdStan and for the initial values (see
#'   Details), or `NULL` for no seed.
#' @param opt_level Compiler optimization level; see
#'   [compile_bilatr_model()].
#' @param output_dir Directory to write CmdStan's raw output CSVs to, or
#'   `NULL` for a temporary directory. When given, the Stan data list is
#'   also saved there (see Value).
#' @param stan_model Which model to fit: `"stable"` (default), or the
#'   experimental `"ou"` or `"stable_gamma"`; see [compile_bilatr_model()].
#'   Fitting an experimental model shows a `message()` once per session.
#' @param ... Additional arguments passed to `CmdStanModel$sample()`
#'   (e.g. `refresh`, `chain_ids`, `output_basename`, `init`).
#' @return A `CmdStanMCMC` fit object, with the Stan data list it was fit
#'   to attached as `attr(fit, "stan_data")` (with its own attributes:
#'   `dyad_ids`, `event_classes`, `stan_model`, ...) and the model's name
#'   as `attr(fit, "stan_model")`. The extractors ([extract_theta()],
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
#' @examplesIf bilatr_cmdstan_available()
#' # 300 simulated USA -> CHN events over 2015-2020
#' set.seed(1)
#' days <- seq(as.Date("2015-01-01"), as.Date("2020-12-31"), by = "day")
#' events <- recode_cameo(tibble::tibble(
#'   Actor1CountryCode = "USA",
#'   Actor2CountryCode = "CHN",
#'   SQLDATE = as.integer(format(sample(days, 300, replace = TRUE), "%Y%m%d")),
#'   EventCode = sample(cameo_lookup$CAMEOEVENTCODE, 300, replace = TRUE)
#' ))
#' agg <- aggregate_events(
#'   events,
#'   grouping_var = "PentaClass", start = "2015-01-01", end = "2020-12-31"
#' )
#' \donttest{
#' # reference "1" (verbal cooperation): higher theta = more cooperative.
#' # Far too short for real use; the defaults are 4 chains of 1000 + 1000.
#' fit <- fit_dyad_ts(
#'   agg,
#'   reference_category = "1",
#'   chains = 2, iter_warmup = 200, iter_sampling = 200, seed = 1, refresh = 0
#' )
#' extract_theta(fit)
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
  chunk_size = 1,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  chains = 4,
  parallel_chains = 1,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 3,
  output_dir = NULL,
  stan_model = "stable",
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
#' @noRd
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
#' Fits the dynamic IRT model (by default the `"stable"` model; see
#' [compile_bilatr_model()] for the available models) across multiple
#' dyads simultaneously (`D > 1`),
#' with hierarchical partial pooling on process noise (via
#' `mu_log_noise`/`sigma_log_noise`), the dispersion parameter `phi` (via
#' `mu_log_phi`/`sigma_log_phi`), and the initial latent state `theta0`
#' (via `sigma_theta0`; `theta0`'s population mean is pinned at exactly 0
#' rather than a separately pooled location). Use [fit_dyad_ts()] instead
#' for a single dyad.
#'
#' `data` is either a dyad-period table or a prepared Stan data list,
#' exactly as for [fit_dyad_ts()]; see its Details for the dispatch, the
#' checks run before compiling, and reproducibility with `seed`.
#' `fit_panel()` requires `D >= 2`.
#'
#' `threads_per_chain` controls `reduce_sum` parallelization across
#' dyads within a chain, combined with `chunk_size` (`C` in the Stan
#' data). The defaults (`threads_per_chain = 1`, `parallel_chains = 1`,
#' `chunk_size = 1`) are conservative: they run on any machine, but use
#' one core. On a multi-core machine, more threads per chain and a larger
#' `chunk_size` can be much faster for large panels; tune them together
#' for your data size and hardware, as described in `vignette("tuning", package = "bilatr")`.
#'
#' @inheritParams fit_dyad_ts
#' @param threads_per_chain Threads per chain for `reduce_sum`
#'   parallelization across dyads. Defaults to 1; see Details.
#' @return As for [fit_dyad_ts()]: a `CmdStanMCMC` fit object with
#'   `attr(fit, "stan_data")` set, and the Stan data list saved to
#'   `output_dir` when one is given.
#' @examplesIf bilatr_cmdstan_available()
#' # 600 simulated events among three countries over 2015-2020
#' set.seed(1)
#' countries <- c("USA", "CHN", "RUS")
#' days <- seq(as.Date("2015-01-01"), as.Date("2020-12-31"), by = "day")
#' pairs <- expand.grid(a1 = countries, a2 = countries, stringsAsFactors = FALSE)
#' pairs <- pairs[pairs$a1 != pairs$a2, ]
#' idx <- sample(nrow(pairs), 600, replace = TRUE)
#' events <- recode_cameo(tibble::tibble(
#'   Actor1CountryCode = pairs$a1[idx],
#'   Actor2CountryCode = pairs$a2[idx],
#'   SQLDATE = as.integer(format(sample(days, 600, replace = TRUE), "%Y%m%d")),
#'   EventCode = sample(cameo_lookup$CAMEOEVENTCODE, 600, replace = TRUE)
#' ))
#' agg <- aggregate_events(
#'   events, grouping_var = "PentaClass", resolution = "yearly",
#'   start = "2015-01-01", end = "2020-12-31"
#' )
#' \donttest{
#' # reference "1" (verbal cooperation): higher theta = more cooperative.
#' # Far too short for real use; the defaults are 4 chains of 1000 + 1000.
#' fit <- fit_panel(
#'   agg,
#'   reference_category = "1",
#'   chains = 2, iter_warmup = 200, iter_sampling = 200, seed = 1, refresh = 0
#' )
#' extract_theta(fit)
#'
#' # on a machine with 16 cores and a large panel: 4 chains in parallel,
#' # 4 threads each (see vignette("tuning", package = "bilatr"))
#' # fit <- fit_panel(
#' #   agg, reference_category = "1",
#' #   chains = 4, parallel_chains = 4,
#' #   threads_per_chain = 4, chunk_size = 50
#' # )
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
  chunk_size = 1,
  rho_prior_a = 8,
  rho_prior_b = 2,
  compute_log_lik = 0,
  prior_only = 0,
  compute_theta_filtered = 0,
  filter_dyads = NULL,
  chains = 4,
  parallel_chains = 1,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  opt_level = 3,
  output_dir = NULL,
  stan_model = "stable",
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
