skip_if_no_cmdstan <- function() {
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- tryCatch({
    cmdstanr::cmdstan_path()
    TRUE
  }, error = function(e) FALSE)
  skip_if_not(has_cmdstan, "CmdStan is not installed")
}

make_fake_events <- function(n = 400, seed = 42, years = 2015:2019) {
  set.seed(seed)
  countries <- c("USA", "CHN", "RUS")
  dates <- format(
    sample(
      seq(as.Date(paste0(min(years), "-01-01")), as.Date(paste0(max(years), "-12-31")), by = "day"),
      n,
      replace = TRUE
    ),
    "%Y%m%d"
  )
  events <- tibble::tibble(
    Actor1CountryCode = sample(countries, n, replace = TRUE),
    Actor2CountryCode = sample(countries, n, replace = TRUE),
    SQLDATE = as.integer(dates),
    EventCode = sample(bilatr::cameo_lookup$CAMEOEVENTCODE, n, replace = TRUE)
  )
  events[events$Actor1CountryCode != events$Actor2CountryCode, ]
}

# --- fixture for CSV-file-path branches (diagnose_convergence(),
# extract_theta()/extract_alpha()/extract_mu_intercept()): a real, tiny
# multi-chain CmdStan run, since read_cmdstan_csv()/
# cmdstanr:::read_csv_metadata() need real CmdStan CSV files, not a
# hand-built posterior::draws_array like make_fake_draws() (in
# test_diagnose_convergence.R) or the synthetic posterior::draws_df
# fixtures in test_orient.R. D/T/A are kept small (a few dozen variables
# spanning all three tiers) purely to keep the test fast; the chunking
# logic itself is exercised via a deliberately tiny chunk_size/
# max_memory_mb, not by the fixture's own size.
#
# stan_model defaults to "stable", which (since 0.4.0's promotion of
# alphanorm -- see NEWS.md) DOES have a reflection symmetry; the default
# call exercises the ordinary (right-basin) path since bilatr_init_fn()'s
# actual init biases toward it. Pass a deliberately wrong-basin init (see
# test_orient.R's pattern) to exercise bilatr_orient()'s sign flip through
# the CSV-chunked path instead -- also pass ... = adapt_engaged = FALSE,
# step_size = <tiny>, max_treedepth = <small> (test_orient.R's pinning
# trick) if the point is to keep the chain from adapting its way out of
# that basin during ordinary warmup, since this fixture's default 30
# warmup iterations are otherwise enough to escape a wrong-basin init on
# a dataset this small. `compute_log_lik`/`anchor_scale` (needed by both
# registered models) default to the same values
# [assemble_stan_data()] does; pass `extra_data` to override them.
make_csv_diagnostics_fixture <- function(stan_model = "stable", init = NULL, extra_data = list(), ...) {
  set.seed(1)
  D <- 6
  Tn <- 4
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- utils::modifyList(
    list(
      T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
      dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
      compute_log_lik = 0, anchor_scale = 0.1
    ),
    extra_data
  )
  mod <- .compile_stan_model(stan_model, opt_level = 1)
  outdir <- tempfile()
  dir.create(outdir)
  sample_args <- utils::modifyList(
    list(
      data = data_list, chains = 2, parallel_chains = 2,
      iter_warmup = 30, iter_sampling = 15, seed = 1, refresh = 0,
      output_dir = outdir, show_messages = FALSE, threads_per_chain = 1
    ),
    list(...)
  )
  if (!is.null(init)) sample_args$init <- init
  fit <- suppressWarnings(do.call(mod$sample, sample_args))

  # a stan_data-shaped list with a manually-built dyad_ids attribute,
  # matching assemble_stan_data()'s output shape, for extract_theta()'s
  # dyad/dyad2/year join -- built directly rather than via
  # assemble_stan_data() so this fixture's D/T/A stay exactly controlled
  # and independent of make_fake_events()'s randomness.
  stan_data <- data_list
  attr(stan_data, "dyad_ids") <- tibble::tibble(
    dyad_id = rep(seq_len(D), each = Tn),
    time_index = rep(seq_len(Tn), D),
    dyad = paste0("dyad", rep(seq_len(D), each = Tn)),
    dyad2 = paste0("dyad", rep(seq_len(D), each = Tn)),
    year = 2000L + rep(seq_len(Tn), D)
  )

  list(
    fit = fit,
    stan_data = stan_data,
    csv_files = list.files(outdir, pattern = "\\.csv$", full.names = TRUE),
    n_dt = tibble::tibble(dyad_id = seq_len(D), n_dt = apply(Y, 1, sum))
  )
}
