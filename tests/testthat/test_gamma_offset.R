# Tests for the 0.7.0 country-level gamma offset (stable_gamma). See
# dev/claude_code_prompt_0.7.0_country_offsets.md for the full design
# this file implements against, and inst/stan/bilatr_alphanorm_gamma.stan's
# header for the model itself.

#' Pure-R replication of bilatr_alphanorm_gamma.stan's transformed-parameters
#' construction of `gamma` from `gamma_z`/`sigma_gamma`/`alpha`
#'
#' Mirrors the four Stan steps exactly (diag_pre_multiply, row-centre
#' across countries, column-centre across categories, project out alpha
#' via dot_self(alpha)) so the constraint algebra and orientation
#' invariance can be tested without a CmdStan fit.
#'
#' @keywords internal
.gamma_construction_r <- function(gamma_z, sigma_gamma, alpha) {
  raw <- sigma_gamma * gamma_z # diag_pre_multiply(sigma_gamma, gamma_z): row k scaled by sigma_gamma[k]
  row_centred <- raw - rowMeans(raw)
  col_centred <- sweep(row_centred, 2, colMeans(row_centred), "-")
  n_countries <- ncol(gamma_z)
  A <- nrow(gamma_z)
  gamma <- matrix(0, A, n_countries)
  for (c in seq_len(n_countries)) {
    v <- col_centred[, c]
    gamma[, c] <- v - (sum(v * alpha) / sum(alpha^2)) * alpha
  }
  gamma
}

.make_random_alpha <- function(A, seed) {
  set.seed(seed)
  alpha_raw <- stats::rnorm(A)
  alpha_raw <- alpha_raw - mean(alpha_raw)
  alpha <- alpha_raw * sqrt(A / sum(alpha_raw^2))
  if (alpha[1] < 0) alpha <- -alpha
  alpha
}

test_that("the four-step gamma construction satisfies all three identification constraints simultaneously", {
  A <- 6
  n_countries <- 9
  alpha <- .make_random_alpha(A, seed = 1)

  for (seed in c(2, 3, 4)) {
    set.seed(seed)
    gamma_z <- matrix(stats::rnorm(A * n_countries), A, n_countries)
    sigma_gamma <- abs(stats::rnorm(A)) + 0.05

    gamma <- .gamma_construction_r(gamma_z, sigma_gamma, alpha)

    # constraint 3: centred across countries, per category
    expect_equal(rowSums(gamma), rep(0, A), tolerance = 1e-10)
    # constraint 2: orthogonal to 1, per country
    expect_equal(colSums(gamma), rep(0, n_countries), tolerance = 1e-10)
    # constraint 1: orthogonal to alpha, per country
    expect_equal(as.vector(alpha %*% gamma), rep(0, n_countries), tolerance = 1e-8)
  }
})

test_that("gamma is identically 0 at n_countries = 1, for any gamma_z/sigma_gamma/alpha", {
  A <- 5
  alpha <- .make_random_alpha(A, seed = 11)

  for (seed in c(12, 13, 14)) {
    set.seed(seed)
    gamma_z <- matrix(stats::rnorm(A), A, 1)
    sigma_gamma <- abs(stats::rnorm(A)) + 0.05
    gamma <- .gamma_construction_r(gamma_z, sigma_gamma, alpha)
    expect_equal(gamma, matrix(0, A, 1))
  }
})

test_that("gamma is unchanged by alpha -> -alpha (orientation invariance)", {
  A <- 6
  n_countries <- 4
  alpha <- .make_random_alpha(A, seed = 21)

  set.seed(22)
  gamma_z <- matrix(stats::rnorm(A * n_countries), A, n_countries)
  sigma_gamma <- abs(stats::rnorm(A)) + 0.05

  gamma_pos <- .gamma_construction_r(gamma_z, sigma_gamma, alpha)
  gamma_neg <- .gamma_construction_r(gamma_z, sigma_gamma, -alpha)

  expect_equal(gamma_pos, gamma_neg, tolerance = 1e-10)
})

test_that("stable_gamma nests stable exactly at n_countries = 1 (CmdStan log_prob check)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  D <- 5
  Tn <- 3
  A <- 4
  Y <- array(sample(0:5, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  base_data <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0,
    n_filter_dyads = 0, filter_dyads = integer(0)
  )
  gamma_data <- utils::modifyList(base_data, list(
    n_countries = 1L, ctry_a = rep(1L, D), ctry_b = rep(1L, D), w_send = rep(1, D)
  ))

  # force_recompile = TRUE: cmdstanr's $init_model_methods() (needed for
  # $log_prob()/$unconstrain_variables()/$constrain_variables() below)
  # refuses to attach to a reused, previously-compiled executable.
  mod_stable <- .compile_stan_model("stable", opt_level = 1, force_recompile = TRUE)
  mod_gamma <- .compile_stan_model("stable_gamma", opt_level = 1, force_recompile = TRUE)

  outdir1 <- tempfile()
  dir.create(outdir1)
  outdir2 <- tempfile()
  dir.create(outdir2)
  fit_stable <- mod_stable$sample(
    data = base_data, chains = 1, iter_warmup = 1, iter_sampling = 1, seed = 1,
    refresh = 0, show_messages = FALSE, output_dir = outdir1, threads_per_chain = 1
  )
  fit_gamma <- mod_gamma$sample(
    data = gamma_data, chains = 1, iter_warmup = 1, iter_sampling = 1, seed = 1,
    refresh = 0, show_messages = FALSE, output_dir = outdir2, threads_per_chain = 1
  )
  suppressMessages(fit_stable$init_model_methods())
  suppressMessages(fit_gamma$init_model_methods())

  make_shared_pars <- function(seed) {
    set.seed(seed)
    mu <- stats::rnorm(A)
    mu <- mu - mean(mu)
    al <- stats::rnorm(A)
    al <- al - mean(al)
    list(
      theta_raw = matrix(stats::rnorm(D * Tn), D, Tn),
      mu_intercept = mu, alpha_raw = al,
      sigma_theta0 = 0.7, z_theta0 = stats::rnorm(D),
      log_process_noise_raw = stats::rnorm(D), mu_log_noise = -1.5, sigma_log_noise = 0.3,
      phi = rep(5, D), mu_log_phi = 1.0, sigma_log_phi = 0.4
    )
  }
  gamma_ref <- list(gamma_z = matrix(0, A, 1), sigma_gamma = rep(0.4, A))

  # (a) the likelihood + every shared prior is bit-identical between the
  # two programs at n_countries = 1: holding gamma_z/sigma_gamma fixed,
  # lp_gamma - lp_stable must be the SAME constant regardless of the
  # shared parameters' own values -- if the two programs' treatment of
  # theta/alpha/mu_intercept/phi/etc. diverged in any way (e.g. a bug in
  # partial_log_lik_offset's g construction), this constant would move
  # as those values move. This constant is NOT asserted to equal a
  # hand-derived closed-form value: cmdstanr's log_prob() bridge folds in
  # Stan's own sum_to_zero_vector/<lower=0> Jacobian bookkeeping, which
  # is a private implementation detail this test does not reimplement --
  # what matters, and what this checks, is that it doesn't move.
  diffs <- vapply(c(11, 22, 33, 44), function(s) {
    shared <- make_shared_pars(s)
    up_s <- fit_stable$unconstrain_variables(variables = shared)
    up_g <- fit_gamma$unconstrain_variables(variables = c(shared, gamma_ref))
    fit_gamma$log_prob(up_g) - fit_stable$log_prob(up_s)
  }, numeric(1))
  expect_equal(diffs, rep(diffs[1], length(diffs)), tolerance = 1e-6)

  # (b) holding the shared parameters fixed, the SENSITIVITY of lp_gamma
  # to gamma_z/sigma_gamma matches the closed-form prior + Jacobian
  # contribution exactly: gamma_z ~ std_normal() (identity transform, no
  # Jacobian) and sigma_gamma ~ normal(0, 0.3) (<lower=0>, Jacobian
  # log(sigma_gamma)).
  shared_fixed <- make_shared_pars(99)
  up_shared_fixed <- fit_stable$unconstrain_variables(variables = shared_fixed)
  lp_shared_fixed <- fit_stable$log_prob(up_shared_fixed)

  eval_gamma_extra <- function(gamma_z, sigma_gamma) {
    up <- fit_gamma$unconstrain_variables(variables = c(shared_fixed, list(gamma_z = gamma_z, sigma_gamma = sigma_gamma)))
    fit_gamma$log_prob(up) - lp_shared_fixed
  }
  gamma_prior_contrib <- function(gamma_z, sigma_gamma) {
    sum(stats::dnorm(as.vector(gamma_z), 0, 1, log = TRUE)) +
      sum(stats::dnorm(sigma_gamma, 0, 0.3, log = TRUE)) +
      sum(log(sigma_gamma))
  }

  set.seed(7)
  gz1 <- matrix(stats::rnorm(A), A, 1)
  sg1 <- abs(stats::rnorm(A)) + 0.1
  gz2 <- matrix(stats::rnorm(A), A, 1)
  sg2 <- abs(stats::rnorm(A)) + 0.1

  actual_delta <- eval_gamma_extra(gz1, sg1) - eval_gamma_extra(gz2, sg2)
  predicted_delta <- gamma_prior_contrib(gz1, sg1) - gamma_prior_contrib(gz2, sg2)
  expect_equal(actual_delta, predicted_delta, tolerance = 1e-6)

  # (c) gamma itself is exactly 0 in the ACTUAL compiled program's
  # transformed parameters (not just the R reimplementation above).
  up_ref <- fit_gamma$unconstrain_variables(variables = c(shared_fixed, gamma_ref))
  gamma_tp <- fit_gamma$constrain_variables(up_ref)$gamma
  expect_equal(as.vector(gamma_tp), rep(0, A))
})

test_that("w_send is exactly 1 for directed dyads, and the observed sender share for undirected dyads", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN"),
    Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA"),
    SQLDATE = 20180101L,
    PentaClass = c(0, 1, 0, 1, 0)
  )

  directed <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = TRUE
  )
  w_send_directed <- attr(directed, "w_send")
  expect_true(all(w_send_directed$w_send == 1))

  undirected <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE
  )
  w_send_undirected <- attr(undirected, "w_send")
  # USA/RUS pair (rows 1-3): pmin("USA","RUS") = "RUS" is side A. 2 of the
  # 3 events have USA (side B) as Actor1, 1 has RUS (side A) as Actor1 --
  # an asymmetric pair, share strictly between 0 and 1.
  usa_rus <- dplyr::filter(w_send_undirected, dyad == "RUS_USA")
  expect_equal(usa_rus$ctry_a_code, "RUS")
  expect_equal(usa_rus$ctry_b_code, "USA")
  expect_equal(usa_rus$w_send, 1 / 3)

  # CHN/USA pair (rows 4-5): pmin("CHN","USA") = "CHN" is side A. BOTH
  # events have CHN as Actor1 -- a one-direction-only pair, giving
  # exactly 1 (not 0 or 1 by coincidence: every event in this pair goes
  # the same way).
  chn_usa <- dplyr::filter(w_send_undirected, dyad == "CHN_USA")
  expect_equal(chn_usa$ctry_a_code, "CHN")
  expect_equal(chn_usa$w_send, 1)
})

test_that("gamma/gamma_z are classified correctly by .classify_bilatr_tier() (Tier 1, and prior-only-excluded respectively)", {
  variable <- c("gamma[1,1]", "gamma[4,12]", "gamma_z[1,1]", "gamma_z[6,3]", "sigma_gamma[2]", "alpha[1]", "theta[3,7]")
  tiers <- .classify_bilatr_tier(variable)

  expect_identical(tiers$tier[tiers$variable %in% c("gamma[1,1]", "gamma[4,12]")], c(1L, 1L))
  expect_identical(tiers$tier[tiers$variable == "sigma_gamma[2]"], 1L)
  expect_true(all(is.na(tiers$tier[tiers$variable %in% c("gamma_z[1,1]", "gamma_z[6,3]")])))
  # unaffected pre-existing classifications
  expect_identical(tiers$tier[tiers$variable == "alpha[1]"], 1L)
  expect_identical(tiers$tier[tiers$variable == "theta[3,7]"], 3L)

  # gamma/gamma_z never carry a dyad_id/time_index
  expect_true(all(is.na(tiers$dyad_id[grepl("^gamma", tiers$variable)])))
})

test_that("assemble_stan_data() attaches n_countries/ctry_a/ctry_b/w_send unconditionally, and is deterministic", {
  events <- make_fake_events()
  events <- recode_cameo(events, code_col = "EventCode")

  sd1 <- assemble_stan_data(
    events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )
  sd2 <- assemble_stan_data(
    events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )

  expect_true(all(c("n_countries", "ctry_a", "ctry_b", "w_send") %in% names(sd1)))
  expect_length(sd1$ctry_a, sd1$D)
  expect_length(sd1$ctry_b, sd1$D)
  expect_length(sd1$w_send, sd1$D)
  expect_true(all(sd1$ctry_a >= 1 & sd1$ctry_a <= sd1$n_countries))
  expect_true(all(sd1$ctry_b >= 1 & sd1$ctry_b <= sd1$n_countries))
  expect_true(all(sd1$w_send >= 0 & sd1$w_send <= 1))
  expect_identical(attr(sd1, "country_codes"), sort(unique(attr(sd1, "country_codes"))))
  expect_length(attr(sd1, "country_codes"), sd1$n_countries)

  expect_identical(sd1$n_countries, sd2$n_countries)
  expect_identical(sd1$ctry_a, sd2$ctry_a)
  expect_identical(sd1$ctry_b, sd2$ctry_b)
  expect_identical(sd1$w_send, sd2$w_send)
  expect_identical(attr(sd1, "country_codes"), attr(sd2, "country_codes"))
})

test_that("assemble_stan_data()'s n_countries excludes countries only present in dropped dyads", {
  # A country ("ISO") that appears in exactly one low-volume dyad should
  # not consume a country index once min_n_events drops that dyad.
  events <- tibble::tibble(
    Actor1CountryCode = c(rep("USA", 20), rep("RUS", 20), "ISO"),
    Actor2CountryCode = c(rep("RUS", 20), rep("USA", 20), "USA"),
    SQLDATE = 20180101L,
    EventCode = sample(bilatr::cameo_lookup$CAMEOEVENTCODE, 41, replace = TRUE)
  )
  events <- recode_cameo(events, code_col = "EventCode")

  sd <- assemble_stan_data(
    events, years = 2018, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 5
  )
  expect_false("ISO" %in% attr(sd, "country_codes"))
})

test_that("stable_gamma recovers a known country-level gamma from simulated data", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  A <- 4
  n_countries <- 4
  Tn <- 6
  sigma_gamma_true <- 0.5
  alpha <- .make_random_alpha(A, seed = 201)

  set.seed(202)
  mu_raw <- stats::rnorm(A)
  mu_intercept <- mu_raw - mean(mu_raw)
  gamma_z <- matrix(stats::rnorm(A * n_countries), A, n_countries)
  gamma_true <- .gamma_construction_r(gamma_z, rep(sigma_gamma_true, A), alpha)

  pairs <- expand.grid(a = seq_len(n_countries), b = seq_len(n_countries))
  pairs <- pairs[pairs$a != pairs$b, ]
  D <- nrow(pairs)
  ctry_a <- pairs$a
  ctry_b <- pairs$b
  w_send <- rep(1, D) # directed

  process_noise <- 0.2
  theta0 <- stats::rnorm(D)
  theta <- matrix(NA_real_, D, Tn)
  theta[, 1] <- theta0 + process_noise * stats::rnorm(D)
  for (t in 2:Tn) theta[, t] <- theta[, t - 1] + process_noise * stats::rnorm(D)
  phi <- rep(15, D)
  is_obs <- matrix(1L, D, Tn)
  n_dt <- 100

  Y <- array(0L, dim = c(D, Tn, A))
  for (d in seq_len(D)) {
    g_d <- gamma_true[, ctry_a[d]]
    for (t in seq_len(Tn)) {
      eta <- alpha * theta[d, t] - mu_intercept - g_d
      p <- exp(eta - max(eta))
      p <- p / sum(p)
      conc <- phi[d] * p
      q <- stats::rgamma(A, shape = conc, rate = 1)
      q <- q / sum(q)
      Y[d, t, ] <- stats::rmultinom(1, n_dt, q)[, 1]
    }
  }

  stan_data <- list(
    D = D, T = Tn, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0,
    n_filter_dyads = 0, filter_dyads = integer(0),
    n_countries = n_countries, ctry_a = ctry_a, ctry_b = ctry_b, w_send = w_send
  )
  attr(stan_data, "dyad_ids") <- tibble::tibble(
    dyad_id = rep(seq_len(D), each = Tn), time_index = rep(seq_len(Tn), D),
    dyad = paste0("dyad", rep(seq_len(D), each = Tn)), dyad2 = paste0("dyad", rep(seq_len(D), each = Tn))
  )
  attr(stan_data, "event_classes") <- as.character(seq_len(A))
  attr(stan_data, "country_codes") <- paste0("C", seq_len(n_countries))

  fit <- suppressWarnings(fit_panel_dev(
    stan_data,
    chains = 2, parallel_chains = 2, threads_per_chain = 1,
    iter_warmup = 300, iter_sampling = 300, seed = 1, opt_level = 1,
    stan_model = "stable_gamma", refresh = 0, show_messages = FALSE
  ))

  gamma_est <- extract_gamma(fit, stan_data)
  truth_tbl <- tibble::tibble(
    country_index = rep(seq_len(n_countries), each = A),
    action_index = rep(seq_len(A), n_countries),
    truth = as.vector(gamma_true)
  )
  cmp <- dplyr::inner_join(gamma_est, truth_tbl, by = c("country_index", "action_index"))
  expect_gt(stats::cor(cmp$mean, cmp$truth), 0.6)

  sigma_gamma_draws <- posterior::as_draws_matrix(.get_draws(fit, "sigma_gamma"))
  sigma_gamma_mean <- mean(sigma_gamma_draws)
  expect_lt(abs(sigma_gamma_mean - sigma_gamma_true), 0.35)
})

test_that("check_compositional_residuals()'s implied_beta_rms is materially smaller under stable_gamma when the simulated truth has country structure, and not when it doesn't", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  simulate_and_fit <- function(has_country_structure, seed) {
    A <- 4
    n_countries <- 4
    Tn <- 6
    alpha <- .make_random_alpha(A, seed = seed)

    set.seed(seed + 1)
    mu_raw <- stats::rnorm(A)
    mu_intercept <- mu_raw - mean(mu_raw)
    if (has_country_structure) {
      gamma_z <- matrix(stats::rnorm(A * n_countries), A, n_countries)
      gamma_true <- .gamma_construction_r(gamma_z, rep(0.6, A), alpha)
    } else {
      gamma_true <- matrix(0, A, n_countries)
    }

    pairs <- expand.grid(a = seq_len(n_countries), b = seq_len(n_countries))
    pairs <- pairs[pairs$a != pairs$b, ]
    D <- nrow(pairs)
    ctry_a <- pairs$a
    ctry_b <- pairs$b
    w_send <- rep(1, D)

    process_noise <- 0.2
    theta0 <- stats::rnorm(D)
    theta <- matrix(NA_real_, D, Tn)
    theta[, 1] <- theta0 + process_noise * stats::rnorm(D)
    for (t in 2:Tn) theta[, t] <- theta[, t - 1] + process_noise * stats::rnorm(D)
    phi <- rep(15, D)
    is_obs <- matrix(1L, D, Tn)
    n_dt <- 100

    Y <- array(0L, dim = c(D, Tn, A))
    for (d in seq_len(D)) {
      g_d <- gamma_true[, ctry_a[d]]
      for (t in seq_len(Tn)) {
        eta <- alpha * theta[d, t] - mu_intercept - g_d
        p <- exp(eta - max(eta))
        p <- p / sum(p)
        conc <- phi[d] * p
        q <- stats::rgamma(A, shape = conc, rate = 1)
        q <- q / sum(q)
        Y[d, t, ] <- stats::rmultinom(1, n_dt, q)[, 1]
      }
    }

    stan_data <- list(
      D = D, T = Tn, A = A, C = 1, is_obs = is_obs, Y = Y,
      compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0,
      n_filter_dyads = 0, filter_dyads = integer(0),
      n_countries = n_countries, ctry_a = ctry_a, ctry_b = ctry_b, w_send = w_send
    )
    attr(stan_data, "dyad_ids") <- tibble::tibble(
      dyad_id = rep(seq_len(D), each = Tn), time_index = rep(seq_len(Tn), D),
      dyad = paste0("dyad", rep(seq_len(D), each = Tn)), dyad2 = paste0("dyad", rep(seq_len(D), each = Tn))
    )
    attr(stan_data, "event_classes") <- as.character(seq_len(A))
    attr(stan_data, "country_codes") <- paste0("C", seq_len(n_countries))

    fit_stable <- suppressWarnings(fit_panel_dev(
      stan_data,
      chains = 2, parallel_chains = 2, threads_per_chain = 1,
      iter_warmup = 250, iter_sampling = 250, seed = 1, opt_level = 1,
      stan_model = "stable", refresh = 0, show_messages = FALSE
    ))
    fit_gamma <- suppressWarnings(fit_panel_dev(
      stan_data,
      chains = 2, parallel_chains = 2, threads_per_chain = 1,
      iter_warmup = 250, iter_sampling = 250, seed = 1, opt_level = 1,
      stan_model = "stable_gamma", refresh = 0, show_messages = FALSE
    ))

    rc_stable <- check_compositional_residuals(
      fit_stable, stan_data,
      stan_model = "stable", n_dyads = D, n_draws = 100, seed = 1
    )
    rc_gamma <- check_compositional_residuals(
      fit_gamma, stan_data,
      stan_model = "stable_gamma", n_dyads = D, n_draws = 100, seed = 1
    )

    list(stable = rc_stable$global$implied_beta_rms, gamma = rc_gamma$global$implied_beta_rms)
  }

  with_structure <- simulate_and_fit(TRUE, seed = 301)
  without_structure <- simulate_and_fit(FALSE, seed = 401)

  expect_lt(with_structure$gamma, with_structure$stable * 0.8)
  expect_gte(without_structure$gamma, without_structure$stable * 0.7)
})
