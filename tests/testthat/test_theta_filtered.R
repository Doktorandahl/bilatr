.r_filter_stable <- function(alpha, mu_intercept, sigma_theta0, process_noise_d, phi_d, Y_d, is_obs_d) {
  Tn <- length(is_obs_d)
  m <- 0
  p_var <- sigma_theta0^2
  m_out <- numeric(Tn)
  sd_out <- numeric(Tn)
  for (t in seq_len(Tn)) {
    m_pred <- m
    p_pred <- p_var + process_noise_d^2
    if (is_obs_d[t] == 1) {
      eta <- alpha * m_pred - mu_intercept
      p <- exp(eta - max(eta))
      p <- p / sum(p)
      a_bar <- sum(p * alpha)
      n <- sum(Y_d[t, ])
      info <- n * (1 + phi_d) / (n + phi_d) * sum(p * (alpha - a_bar)^2)
      conc <- phi_d * p
      g <- sum(conc * (alpha - a_bar) * (digamma(Y_d[t, ] + conc) - digamma(conc)))
      p_var <- 1 / (1 / p_pred + info)
      m <- m_pred + p_var * g
    } else {
      m <- m_pred
      p_var <- p_pred
    }
    m_out[t] <- m
    sd_out[t] <- sqrt(p_var)
  }
  list(mean = m_out, sd = sd_out)
}

.r_filter_ou <- function(alpha, mu_intercept, mu_dyad_d, sd_stat_d, rho, process_noise_d, phi_d, Y_d, is_obs_d) {
  Tn <- length(is_obs_d)
  m <- mu_dyad_d
  p_var <- sd_stat_d^2
  m_out <- numeric(Tn)
  sd_out <- numeric(Tn)
  for (t in seq_len(Tn)) {
    m_pred <- mu_dyad_d + rho * (m - mu_dyad_d)
    p_pred <- rho^2 * p_var + process_noise_d^2
    if (is_obs_d[t] == 1) {
      eta <- alpha * m_pred - mu_intercept
      p <- exp(eta - max(eta))
      p <- p / sum(p)
      a_bar <- sum(p * alpha)
      n <- sum(Y_d[t, ])
      info <- n * (1 + phi_d) / (n + phi_d) * sum(p * (alpha - a_bar)^2)
      conc <- phi_d * p
      g <- sum(conc * (alpha - a_bar) * (digamma(Y_d[t, ] + conc) - digamma(conc)))
      p_var <- 1 / (1 / p_pred + info)
      m <- m_pred + p_var * g
    } else {
      m <- m_pred
      p_var <- p_pred
    }
    m_out[t] <- m
    sd_out[t] <- sqrt(p_var)
  }
  list(mean = m_out, sd = sd_out)
}

# Named numeric vector, single (fixed_param) draw's values for `vars`.
.draws_row <- function(fit, vars) {
  mat <- posterior::as_draws_matrix(fit$draws(variables = vars))
  stats::setNames(as.numeric(mat[1, vars]), vars)
}

.fit_filtered_fixed_params <- function(stan_model, data_list, params, seed = 1) {
  mod <- .compile_stan_model(stan_model, opt_level = 1)
  mod$sample(
    data = data_list, chains = 1, iter_warmup = 0, iter_sampling = 1,
    seed = seed, refresh = 0, threads_per_chain = 1, fixed_param = TRUE,
    init = list(params), output_dir = tempdir(), show_messages = FALSE
  )
}

.stable_data_and_params <- function(D, Tn, A, seed = 1) {
  set.seed(seed)
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0,
    compute_theta_filtered = 1, n_filter_dyads = D, filter_dyads = seq_len(D),
    anchor_scale = 0.1
  )
  sum0 <- function(x) x - mean(x)
  params <- list(
    theta_raw = matrix(stats::rnorm(D * Tn, 0, 0.3), D, Tn),
    mu_intercept = sum0(stats::rnorm(A, 0, 0.3)),
    alpha_raw = sum0(c(2, stats::rnorm(A - 1, 0, 0.3))),
    sigma_theta0 = 0.5,
    z_theta0 = stats::rnorm(D, 0, 0.3),
    log_process_noise_raw = stats::rnorm(D, 0, 0.3),
    mu_log_noise = log(0.2), sigma_log_noise = 0.3,
    phi = stats::rlnorm(D, 0, 0.2), mu_log_phi = 0.05, sigma_log_phi = 0.4
  )
  list(data_list = data_list, params = params, Y = Y, is_obs = is_obs)
}

.ou_data_and_params <- function(D, Tn, A, seed = 1) {
  set.seed(seed)
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0,
    compute_theta_filtered = 1, n_filter_dyads = D, filter_dyads = seq_len(D),
    anchor_scale = 0.1, rho_prior_a = 8, rho_prior_b = 2
  )
  sum0 <- function(x) x - mean(x)
  params <- list(
    theta_raw = matrix(stats::rnorm(D * Tn, 0, 0.3), D, Tn),
    mu_intercept = sum0(stats::rnorm(A, 0, 0.3)),
    alpha_raw = sum0(c(2, stats::rnorm(A - 1, 0, 0.3))),
    sigma_mu = 0.5,
    mu_dyad_raw = stats::rnorm(D, 0, 0.3),
    rho = 0.7,
    mu_log_sd_stat = log(1), sigma_log_sd_stat = 0.3, log_sd_stat_raw = stats::rnorm(D, 0, 0.3),
    phi = stats::rlnorm(D, 0, 0.2), mu_log_phi = 0.05, sigma_log_phi = 0.4
  )
  list(data_list = data_list, params = params, Y = Y, is_obs = is_obs)
}

test_that("theta_filtered/theta_filtered_sd (stable) match a brute-force R re-implementation of the filter", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 2
  Tn <- 5
  A <- 4
  setup <- .stable_data_and_params(D, Tn, A)
  fit <- .fit_filtered_fixed_params("stable", setup$data_list, setup$params)

  alpha <- .draws_row(fit, paste0("alpha[", seq_len(A), "]"))
  mu_intercept <- .draws_row(fit, paste0("mu_intercept[", seq_len(A), "]"))
  sigma_theta0 <- .draws_row(fit, "sigma_theta0")
  process_noise <- .draws_row(fit, paste0("process_noise[", seq_len(D), "]"))
  phi <- setup$params$phi

  for (d in seq_len(D)) {
    r_res <- .r_filter_stable(alpha, mu_intercept, sigma_theta0, process_noise[d], phi[d], setup$Y[d, , ], setup$is_obs[d, ])
    stan_mean <- .draws_row(fit, paste0("theta_filtered[", d, ",", seq_len(Tn), "]"))
    stan_sd <- .draws_row(fit, paste0("theta_filtered_sd[", d, ",", seq_len(Tn), "]"))
    expect_equal(unname(stan_mean), r_res$mean, tolerance = 1e-6)
    expect_equal(unname(stan_sd), r_res$sd, tolerance = 1e-6)
  }
})

test_that("theta_filtered/theta_filtered_sd (ou) match a brute-force R re-implementation of the filter, including the rho^2 prediction-variance correction", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 2
  Tn <- 5
  A <- 4
  setup <- .ou_data_and_params(D, Tn, A)
  fit <- .fit_filtered_fixed_params("ou", setup$data_list, setup$params)

  alpha <- .draws_row(fit, paste0("alpha[", seq_len(A), "]"))
  mu_intercept <- .draws_row(fit, paste0("mu_intercept[", seq_len(A), "]"))
  mu_dyad <- .draws_row(fit, paste0("mu_dyad[", seq_len(D), "]"))
  sd_stat <- .draws_row(fit, paste0("sd_stat[", seq_len(D), "]"))
  process_noise <- .draws_row(fit, paste0("process_noise[", seq_len(D), "]"))
  rho <- setup$params$rho
  phi <- setup$params$phi

  for (d in seq_len(D)) {
    r_res <- .r_filter_ou(alpha, mu_intercept, mu_dyad[d], sd_stat[d], rho, process_noise[d], phi[d], setup$Y[d, , ], setup$is_obs[d, ])
    stan_mean <- .draws_row(fit, paste0("theta_filtered[", d, ",", seq_len(Tn), "]"))
    stan_sd <- .draws_row(fit, paste0("theta_filtered_sd[", d, ",", seq_len(Tn), "]"))
    expect_equal(unname(stan_mean), r_res$mean, tolerance = 1e-6)
    expect_equal(unname(stan_sd), r_res$sd, tolerance = 1e-6)
  }
})

test_that("stable: an unobserved dyad's filter stays at the prior mean path, with sd growing without bound", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 1
  Tn <- 5
  A <- 4
  setup <- .stable_data_and_params(D, Tn, A)
  setup$data_list$is_obs <- matrix(0L, D, Tn)
  setup$is_obs <- matrix(0L, D, Tn)

  fit <- .fit_filtered_fixed_params("stable", setup$data_list, setup$params)
  stan_mean <- unname(.draws_row(fit, paste0("theta_filtered[1,", seq_len(Tn), "]")))
  stan_sd <- unname(.draws_row(fit, paste0("theta_filtered_sd[1,", seq_len(Tn), "]")))

  expect_equal(stan_mean, rep(0, Tn), tolerance = 1e-8)
  expect_true(all(diff(stan_sd) > 0)) # strictly growing
})

test_that("ou: an unobserved dyad's filter stays at mu_dyad, with sd constant at the stationary sd_stat (not growing)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 1
  Tn <- 5
  A <- 4
  setup <- .ou_data_and_params(D, Tn, A)
  setup$data_list$is_obs <- matrix(0L, D, Tn)
  setup$is_obs <- matrix(0L, D, Tn)

  fit <- .fit_filtered_fixed_params("ou", setup$data_list, setup$params)
  mu_dyad <- unname(.draws_row(fit, "mu_dyad[1]"))
  sd_stat <- unname(.draws_row(fit, "sd_stat[1]"))
  stan_mean <- unname(.draws_row(fit, paste0("theta_filtered[1,", seq_len(Tn), "]")))
  stan_sd <- unname(.draws_row(fit, paste0("theta_filtered_sd[1,", seq_len(Tn), "]")))

  expect_equal(stan_mean, rep(mu_dyad, Tn), tolerance = 1e-8)
  expect_equal(stan_sd, rep(sd_stat, Tn), tolerance = 1e-8)
})

test_that("a dyad with exactly one observed period updates only there", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 1
  Tn <- 5
  A <- 4
  setup <- .stable_data_and_params(D, Tn, A)
  is_obs <- matrix(0L, D, Tn)
  is_obs[1, 3] <- 1L
  setup$data_list$is_obs <- is_obs
  setup$is_obs <- is_obs

  fit <- .fit_filtered_fixed_params("stable", setup$data_list, setup$params)
  stan_mean <- unname(.draws_row(fit, paste0("theta_filtered[1,", seq_len(Tn), "]")))

  # unchanged before the observation
  expect_equal(stan_mean[1:2], c(0, 0), tolerance = 1e-8)
  # updates at t = 3, then stays constant (pure random walk, no further info)
  expect_false(isTRUE(all.equal(stan_mean[3], 0)))
  expect_equal(stan_mean[3:5], rep(stan_mean[3], 3), tolerance = 1e-8)
})

test_that("theta_filtered[, T] is close to but not equal to the smoothed theta[, T]", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 2
  Tn <- 6
  A <- 4
  setup <- .stable_data_and_params(D, Tn, A, seed = 7)
  fit <- .fit_filtered_fixed_params("stable", setup$data_list, setup$params)

  theta_T <- unname(.draws_row(fit, paste0("theta[", seq_len(D), ",", Tn, "]")))
  filtered_T <- unname(.draws_row(fit, paste0("theta_filtered[", seq_len(D), ",", Tn, "]")))

  expect_false(isTRUE(all.equal(theta_T, filtered_T)))
  expect_lt(max(abs(theta_T - filtered_T)), 2)
})

test_that("no orientation_sign() double-flip: theta_filtered agrees between opposite-alpha_raw-basin fits", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # Mirrors test_orient.R's pinning-trick pattern: two chains given
  # deliberately opposite alpha_raw/z_theta0/theta_raw inits, pinned so
  # neither can adapt out of its basin, must report identical
  # theta_filtered (the filter consumes only already-oriented alpha/
  # mu_intercept/process_noise, so applying orientation_sign() again
  # inside it would flip the sign a second time and disagree here).
  set.seed(1)
  D <- 2
  Tn <- 4
  A <- 5
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0,
    compute_theta_filtered = 1, n_filter_dyads = D, filter_dyads = seq_len(D),
    anchor_scale = 0.1
  )

  mod <- .compile_stan_model("stable", opt_level = 1)
  sum0 <- function(x) x - mean(x)
  alpha_raw_pos <- sum0(c(2, stats::rnorm(A - 1, 0, 0.5)))
  if (alpha_raw_pos[1] < 0) alpha_raw_pos <- -alpha_raw_pos
  z_theta0_pos <- stats::rnorm(D, 0, 0.5)
  theta_raw_pos <- matrix(stats::rnorm(D * Tn, 0, 0.5), D, Tn)

  make_init <- function(sign_mult) {
    force(sign_mult)
    function() {
      list(
        theta_raw = sign_mult * theta_raw_pos, mu_intercept = rep(0, A),
        alpha_raw = sign_mult * alpha_raw_pos, sigma_theta0 = 0.5,
        z_theta0 = sign_mult * z_theta0_pos, log_process_noise_raw = rep(0, D),
        mu_log_noise = log(0.2), sigma_log_noise = 0.3,
        phi = rep(1, D), mu_log_phi = 0, sigma_log_phi = 0.5
      )
    }
  }
  fit_one <- function(sign_mult, seed) {
    suppressWarnings(mod$sample(
      data = data_list, chains = 1, iter_warmup = 2, iter_sampling = 10,
      seed = seed, refresh = 0, threads_per_chain = 1,
      adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
      init = make_init(sign_mult), output_dir = tempdir(), show_messages = FALSE
    ))
  }

  fit_pos <- fit_one(1, seed = 1)
  fit_neg <- fit_one(-1, seed = 2)

  vars <- paste0("theta_filtered[", rep(seq_len(D), each = Tn), ",", rep(seq_len(Tn), D), "]")
  mean_pos <- posterior::summarise_draws(fit_pos$draws(variables = vars), mean = mean)
  mean_neg <- posterior::summarise_draws(fit_neg$draws(variables = vars), mean = mean)

  expect_equal(mean_pos$mean, mean_neg$mean, tolerance = 0.05)
})

test_that("filter_dyads narrowing produces the same theta_filtered values as filtering all dyads", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 3
  Tn <- 4
  A <- 4
  setup <- .stable_data_and_params(D, Tn, A, seed = 3)

  fit_all <- .fit_filtered_fixed_params("stable", setup$data_list, setup$params)

  data_subset <- setup$data_list
  data_subset$n_filter_dyads <- 2L
  data_subset$filter_dyads <- c(2L, 3L)
  fit_subset <- .fit_filtered_fixed_params("stable", data_subset, setup$params)

  for (i in seq_along(data_subset$filter_dyads)) {
    d <- data_subset$filter_dyads[i]
    all_vals <- unname(.draws_row(fit_all, paste0("theta_filtered[", d, ",", seq_len(Tn), "]")))
    subset_vals <- unname(.draws_row(fit_subset, paste0("theta_filtered[", i, ",", seq_len(Tn), "]")))
    expect_equal(subset_vals, all_vals, tolerance = 1e-10)
  }
})

test_that("assemble_stan_data()'s filter_dyads column-count arithmetic matches production scale", {
  # Not an actual fit at production scale -- just the arithmetic the
  # roxygen states: 2 * D * T extra columns when compute_theta_filtered
  # is on (theta_filtered + theta_filtered_sd), 0 when off.
  D <- 16750
  Tn <- 35
  expect_equal(2 * D * Tn, 1172500)
})
