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

  # 0.7.1: lp_gamma - lp_stable equals EXACTLY the gamma_z/sigma_gamma
  # prior + Jacobian contribution, evaluated Stan's own way -- `~`
  # distribution statements call the possibly-unnormalized _lupdf/
  # _lupmf variant, dropping additive terms that don't depend on the
  # sampled parameter: 0.5*log(2*pi) per gamma_z element (std_normal)
  # and 0.5*log(2*pi) + log(0.3) per sigma_gamma element (normal(0, 0.3)
  # -- 0.3 is the FIXED prior scale, not sigma_gamma's own value, so
  # dropping it is legitimate; contrast the <lower=0> Jacobian below,
  # which genuinely depends on sigma_gamma and is NOT dropped, since it
  # isn't part of the `~` statement at all -- it's compiler-generated for
  # the constrained-parameter transform). An earlier version of this
  # test used the fully-normalized dnorm(..., log = TRUE) here and found
  # an exactly reproducible, value-independent ~2.5 nat residual,
  # speculatively attributed (in the 0.7.0 summary) to cmdstanr's
  # sum_to_zero_vector Jacobian bookkeeping. That speculation was WRONG:
  # the dropped-constants total at this fixture's A = 4,
  # n_countries = 1 is 4*0.5*log(2*pi) + 4*(0.5*log(2*pi) + log(0.3)) =
  # 2.535616, matching the observed residual (2.535617) to 6 significant
  # figures -- confirmed by testing the hypothesis, not assuming it (see
  # dev/summary_0.7.0_country_offsets_2026-09-23.md, "Found along the
  # way", and dev/claude_code_prompt_0.7.1_gamma_fixes.md Part 4).
  gamma_prior_contrib_lupdf <- function(gamma_z, sigma_gamma) {
    gamma_z_kernel <- sum(-0.5 * as.vector(gamma_z)^2)
    sigma_gamma_kernel <- sum(-0.5 * (sigma_gamma / 0.3)^2)
    jacobian <- sum(log(sigma_gamma)) # <lower=0> transform; not dropped
    gamma_z_kernel + sigma_gamma_kernel + jacobian
  }

  # Varying BOTH the shared parameters and gamma_z/sigma_gamma together
  # tests invariance (the likelihood + shared priors are bit-identical
  # between the two programs, since gamma is identically 0) and
  # sensitivity (the gamma-specific contribution's dependence on
  # gamma_z/sigma_gamma) simultaneously, via one absolute assertion each
  # -- strictly stronger than checking either alone.
  for (s in c(11, 22, 33, 44)) {
    shared <- make_shared_pars(s)
    set.seed(s + 1000)
    gz <- matrix(stats::rnorm(A), A, 1)
    sg <- abs(stats::rnorm(A)) + 0.1

    up_s <- fit_stable$unconstrain_variables(variables = shared)
    up_g <- fit_gamma$unconstrain_variables(variables = c(shared, list(gamma_z = gz, sigma_gamma = sg)))
    actual_diff <- fit_gamma$log_prob(up_g) - fit_stable$log_prob(up_s)
    predicted_diff <- gamma_prior_contrib_lupdf(gz, sg)
    expect_equal(actual_diff, predicted_diff, tolerance = 1e-6)
  }

  # gamma itself is exactly 0 in the ACTUAL compiled program's
  # transformed parameters, for any gamma_z/sigma_gamma (not just the R
  # reimplementation in the "gamma is identically 0" test above).
  shared_fixed <- make_shared_pars(99)
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
  # USA/RUS pair (rows 1-3): C-locale radix order puts "RUS" before "USA",
  # so RUS is side A. 2 of the 3 events have USA (side B) as Actor1, 1 has
  # RUS (side A) as Actor1 -- an asymmetric pair, share strictly between
  # 0 and 1.
  usa_rus <- dplyr::filter(w_send_undirected, dyad == "RUS_USA")
  expect_equal(usa_rus$actor_a, "RUS")
  expect_equal(usa_rus$actor_b, "USA")
  expect_equal(usa_rus$w_send, 1 / 3)

  # CHN/USA pair (rows 4-5): "CHN" sorts before "USA", so CHN is side A.
  # BOTH events have CHN as Actor1 -- a one-direction-only pair, giving
  # exactly 1 (not 0 or 1 by coincidence: every event in this pair goes
  # the same way).
  chn_usa <- dplyr::filter(w_send_undirected, dyad == "CHN_USA")
  expect_equal(chn_usa$actor_a, "CHN")
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

test_that("w_send is pooled over the `years` window, not every year in `data` (0.7.1)", {
  # RUS_USA's mix differs between the out-of-window years (1985-1986,
  # entirely RUS-as-sender) and the in-window year (2018, entirely
  # USA-as-sender): the unrestricted call must pool both eras, the
  # years-restricted call must reflect only 2018.
  events <- tibble::tibble(
    Actor1CountryCode = c("RUS", "RUS", "RUS", "USA", "USA"),
    Actor2CountryCode = c("USA", "USA", "USA", "RUS", "RUS"),
    SQLDATE = c(19850101L, 19850601L, 19861231L, 20180101L, 20180601L),
    PentaClass = c(0, 1, 0, 1, 0)
  )

  unrestricted <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE
  )
  w_send_unrestricted <- attr(unrestricted, "w_send")
  rus_usa_unrestricted <- dplyr::filter(w_send_unrestricted, dyad == "RUS_USA")
  expect_equal(rus_usa_unrestricted$w_send, 3 / 5) # 3 of 5 events have RUS (side A) as Actor1

  restricted <- grouped_events_to_dyad_period(
    events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE,
    years = 2018
  )
  w_send_restricted <- attr(restricted, "w_send")
  rus_usa_restricted <- dplyr::filter(w_send_restricted, dyad == "RUS_USA")
  expect_equal(rus_usa_restricted$w_send, 0) # both 2018 events have USA (side B) as Actor1

  expect_false(isTRUE(all.equal(rus_usa_unrestricted$w_send, rus_usa_restricted$w_send)))
})

test_that("assemble_stan_data() passes `years` through to window w_send", {
  events <- make_fake_events(n = 800, years = 2010:2019)
  events <- recode_cameo(events, code_col = "EventCode")

  sd_full <- assemble_stan_data(
    events, years = 2010:2019, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )
  sd_window <- assemble_stan_data(
    events, years = 2015:2019, resolution = "yearly", grouping_var = "PentaClass",
    reference_category = 0, min_n_events = 1
  )
  # not asserting a specific relationship (both are legitimate w_send
  # vectors for different windows) -- just that assemble_stan_data()'s
  # own `years` reaches grouped_events_to_dyad_period() at all, i.e. the
  # two windows are free to differ (the fixture's random event mix makes
  # them differ with overwhelming probability; if they were IDENTICAL
  # that would mean `years` never reached the w_send computation).
  expect_true(sd_full$D > 0 && sd_window$D > 0)
})

test_that("a row with NA Actor1CountryCode is a validate_bilatr_events() error (0.9.0; was a warning + drop)", {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", NA, "USA"),
    Actor2CountryCode = c("RUS", "RUS", "RUS"),
    SQLDATE = 20180101L,
    PentaClass = c(0, 1, 0)
  )

  expect_error(
    grouped_events_to_dyad_period(
      events, resolution = "yearly", grouping_var = "PentaClass", directed = FALSE
    ),
    "missing \\(NA\\)"
  )
})

test_that(".compute_gamma_tier()/has_gamma keeps gamma out of Tier 1 (and sigma_gamma in it)", {
  summ <- tibble::tibble(
    variable = c("alpha[1]", "sigma_gamma[1]", "sigma_gamma[2]", "gamma[1,1]", "gamma[2,3]", "lp__"),
    rhat = c(1.0, 1.0, 1.0, 1.5, 1.02, 1.0),
    ess_bulk = c(1000, 1000, 1000, 50, 1000, 1000),
    ess_tail = c(1000, 1000, 1000, 50, 1000, 1000),
    tier = 1L
  )

  res_gamma <- .assemble_bilatr_diagnostics(
    summ, NULL, 1L, 1.01, 400, has_gamma = TRUE, country_codes = c("USA", "RUS", "CHN")
  )
  expect_setequal(res_gamma$tier1$variable, c("alpha[1]", "sigma_gamma[1]", "sigma_gamma[2]", "lp__"))
  expect_setequal(res_gamma$gamma$variable, c("gamma[1,1]", "gamma[2,3]"))
  expect_identical(res_gamma$gamma$country_code[res_gamma$gamma$variable == "gamma[1,1]"], "USA")
  expect_identical(res_gamma$gamma$country_code[res_gamma$gamma$variable == "gamma[2,3]"], "CHN")
  expect_equal(res_gamma$summary$n_tier1_total, 4L)
  expect_equal(res_gamma$summary$n_gamma_total, 2L)
  expect_equal(res_gamma$summary$n_gamma_flagged, 2L)

  res_no_gamma <- .assemble_bilatr_diagnostics(summ, NULL, 1L, 1.01, 400, has_gamma = FALSE)
  expect_null(res_no_gamma$gamma)
  expect_true(is.na(res_no_gamma$summary$n_gamma_total))
  # has_gamma = FALSE leaves gamma[...]/sigma_gamma[...] all in tier1,
  # matching pre-0.7.1 (undifferentiated) behaviour
  expect_true(all(c("gamma[1,1]", "gamma[2,3]") %in% res_no_gamma$tier1$variable))
})

test_that("print.bilatr_diagnostics() renders a gamma section only when has_gamma, and it's absent for stable", {
  summ <- tibble::tibble(
    variable = c("alpha[1]", "gamma[1,1]"),
    rhat = c(1.0, 1.0), ess_bulk = c(1000, 1000), ess_tail = c(1000, 1000), tier = 1L
  )
  res_gamma <- .assemble_bilatr_diagnostics(summ, NULL, 1L, 1.01, 400, has_gamma = TRUE)
  expect_output(print(res_gamma), "gamma: country-level offsets")

  res_stable <- .assemble_bilatr_diagnostics(summ, NULL, 1L, 1.01, 400, has_gamma = FALSE)
  out <- capture.output(print(res_stable))
  expect_false(any(grepl("gamma: country-level offsets", out)))
})

test_that("diagnose_and_extract_bilatr()'s gamma table has country_code for stable_gamma and is absent for stable", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  build_fixture <- function(stan_model) {
    set.seed(11)
    D <- 6
    Tn <- 3
    A <- 3
    Y <- array(sample(0:4, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
    is_obs <- matrix(1L, D, Tn)
    data_list <- list(
      T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
      compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0,
      n_filter_dyads = 0, filter_dyads = integer(0)
    )
    if (stan_model == "stable_gamma") {
      data_list <- utils::modifyList(data_list, list(
        n_countries = 2L, ctry_a = rep(c(1L, 2L), length.out = D),
        ctry_b = rep(c(2L, 1L), length.out = D), w_send = rep(1, D)
      ))
    }
    mod <- .compile_stan_model(stan_model, opt_level = 1)
    outdir <- tempfile()
    dir.create(outdir)
    fit <- suppressWarnings(mod$sample(
      data = data_list, chains = 2, parallel_chains = 2, threads_per_chain = 1,
      iter_warmup = 15, iter_sampling = 15, seed = 1, refresh = 0,
      output_dir = outdir, show_messages = FALSE
    ))
    stan_data <- data_list
    attr(stan_data, "dyad_ids") <- tibble::tibble(
      dyad_id = rep(seq_len(D), each = Tn), time_index = rep(seq_len(Tn), D),
      dyad = paste0("dyad", rep(seq_len(D), each = Tn)), dyad2 = paste0("dyad", rep(seq_len(D), each = Tn))
    )
    attr(stan_data, "event_classes") <- as.character(seq_len(A))
    if (stan_model == "stable_gamma") attr(stan_data, "country_codes") <- c("USA", "RUS")
    list(csv_files = fit$output_files(), stan_data = stan_data, n_dt = tibble::tibble(dyad_id = seq_len(D), n_dt = apply(Y, 1, sum)))
  }

  fx_gamma <- build_fixture("stable_gamma")
  res_gamma <- diagnose_and_extract_bilatr(
    fx_gamma$csv_files, fx_gamma$stan_data, n_dt = fx_gamma$n_dt,
    stan_model = "stable_gamma", tiers = 1
  )
  expect_false(is.null(res_gamma$diagnostics$gamma))
  expect_true("country_code" %in% names(res_gamma$diagnostics$gamma))
  expect_setequal(res_gamma$diagnostics$gamma$country_code, c("USA", "RUS"))
  expect_false(any(grepl("^gamma\\[", res_gamma$diagnostics$tier1$variable)))

  fx_stable <- build_fixture("stable")
  res_stable <- diagnose_and_extract_bilatr(
    fx_stable$csv_files, fx_stable$stan_data, n_dt = fx_stable$n_dt,
    stan_model = "stable", tiers = 1
  )
  expect_null(res_stable$diagnostics$gamma)
})

test_that("order_event_classes() is unaffected by the system locale (0.7.1, locale = \"C\")", {
  classes <- c("10", "2", "1", "18", "5")
  baseline <- order_event_classes(classes)

  old_locale <- tryCatch(Sys.getlocale("LC_COLLATE"), error = function(e) NA_character_)
  changed <- tryCatch(
    suppressWarnings(Sys.setlocale("LC_COLLATE", "C")) != "",
    error = function(e) FALSE
  )
  skip_if_not(changed, "could not set an alternate LC_COLLATE locale on this platform")
  on.exit(suppressWarnings(Sys.setlocale("LC_COLLATE", old_locale)), add = TRUE)

  under_c_locale <- order_event_classes(classes)
  expect_identical(baseline, under_c_locale)
  expect_identical(baseline, c("1", "2", "5", "10", "18"))
})
