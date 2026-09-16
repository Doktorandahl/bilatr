.make_synthetic_draws <- function(extra) {
  n <- length(extra[[1]])
  df <- data.frame(
    .chain = 1L, .iteration = seq_len(n), .draw = seq_len(n),
    as.data.frame(extra, check.names = FALSE),
    check.names = FALSE
  )
  posterior::as_draws_df(df)
}

test_that("bilatr_orient() flips exactly the legacy stable_soft_anchor flip-list and leaves the rest unchanged", {
  set.seed(1)
  n <- 20
  extra <- list(
    `alpha[1]` = stats::rnorm(n, -2, 0.1), # deliberately negative median
    `alpha[2]` = stats::rnorm(n, 1, 0.1),
    `alpha_raw[1]` = stats::rnorm(n, -1, 0.1),
    `theta[1,1]` = stats::rnorm(n, -0.5, 0.1),
    `theta[1,2]` = stats::rnorm(n, -0.6, 0.1),
    `theta0[1]` = stats::rnorm(n, -0.3, 0.1),
    `z_theta0[1]` = stats::rnorm(n, -0.4, 0.1),
    `theta_raw[1,1]` = stats::rnorm(n, -0.1, 0.1),
    `mu_intercept[1]` = stats::rnorm(n, 0.2, 0.1),
    `phi[1]` = stats::rnorm(n, 1, 0.1),
    `sigma_theta0` = stats::rnorm(n, 0.5, 0.1)
  )
  draws <- .make_synthetic_draws(extra)

  flip_vars <- c("alpha", "alpha_raw", "theta", "theta0", "z_theta0", "theta_raw")
  unchanged_vars <- c("mu_intercept", "phi", "sigma_theta0")
  oriented <- bilatr_orient(draws, stan_model = "stable_soft_anchor", variables = c(flip_vars, unchanged_vars))

  for (v in c("alpha[1]", "alpha[2]", "alpha_raw[1]", "theta[1,1]", "theta[1,2]", "theta0[1]", "z_theta0[1]", "theta_raw[1,1]")) {
    expect_equal(
      posterior::extract_variable(oriented, v), -extra[[v]],
      info = paste("expected", v, "to be negated")
    )
  }
  for (v in c("mu_intercept[1]", "phi[1]", "sigma_theta0")) {
    expect_equal(
      posterior::extract_variable(oriented, v), extra[[v]],
      info = paste("expected", v, "to be unchanged")
    )
  }

  # the invariant the whole symmetry rests on: alpha .* theta unaffected
  alpha1_theta11_before <- extra[["alpha[1]"]] * extra[["theta[1,1]"]]
  alpha1_theta11_after <- posterior::extract_variable(oriented, "alpha[1]") *
    posterior::extract_variable(oriented, "theta[1,1]")
  expect_equal(alpha1_theta11_after, alpha1_theta11_before)
})

test_that("bilatr_orient() flips exactly the legacy ou_soft_anchor flip-list and leaves the rest unchanged", {
  set.seed(2)
  n <- 20
  extra <- list(
    `alpha[1]` = stats::rnorm(n, -1.5, 0.1),
    `alpha_raw[1]` = stats::rnorm(n, -0.8, 0.1),
    `theta[1,1]` = stats::rnorm(n, -0.4, 0.1),
    `mu_dyad[1]` = stats::rnorm(n, -0.2, 0.1),
    `mu_dyad_raw[1]` = stats::rnorm(n, -0.3, 0.1),
    `theta_raw[1,1]` = stats::rnorm(n, -0.1, 0.1),
    `mu_intercept[1]` = stats::rnorm(n, 0.2, 0.1),
    `sigma_mu` = stats::rnorm(n, 0.5, 0.1),
    `sd_stat[1]` = stats::rnorm(n, 0.4, 0.1),
    `rho` = stats::runif(n, 0.5, 0.9),
    `within_between_ratio` = stats::rnorm(n, 1, 0.1)
  )
  draws <- .make_synthetic_draws(extra)

  flip_vars <- c("alpha", "alpha_raw", "theta", "mu_dyad", "mu_dyad_raw", "theta_raw")
  unchanged_vars <- c("mu_intercept", "sigma_mu", "sd_stat", "rho", "within_between_ratio")
  oriented <- bilatr_orient(draws, stan_model = "ou_soft_anchor", variables = c(flip_vars, unchanged_vars))

  for (v in c("alpha[1]", "alpha_raw[1]", "theta[1,1]", "mu_dyad[1]", "mu_dyad_raw[1]", "theta_raw[1,1]")) {
    expect_equal(posterior::extract_variable(oriented, v), -extra[[v]], info = paste("expected", v, "to be negated"))
  }
  for (v in c("mu_intercept[1]", "sigma_mu", "sd_stat[1]", "rho", "within_between_ratio")) {
    expect_equal(posterior::extract_variable(oriented, v), extra[[v]], info = paste("expected", v, "to be unchanged"))
  }

  alpha1_theta11_before <- extra[["alpha[1]"]] * extra[["theta[1,1]"]]
  alpha1_theta11_after <- posterior::extract_variable(oriented, "alpha[1]") *
    posterior::extract_variable(oriented, "theta[1,1]")
  expect_equal(alpha1_theta11_after, alpha1_theta11_before)
})

test_that("bilatr_orient() does not flip when alpha[1]'s median is already positive", {
  set.seed(3)
  n <- 20
  extra <- list(
    `alpha[1]` = stats::rnorm(n, 2, 0.1),
    `theta[1,1]` = stats::rnorm(n, 0.5, 0.1)
  )
  draws <- .make_synthetic_draws(extra)
  oriented <- bilatr_orient(draws, stan_model = "stable_soft_anchor", variables = c("alpha", "theta"))
  expect_equal(posterior::extract_variable(oriented, "alpha[1]"), extra[["alpha[1]"]])
  expect_equal(posterior::extract_variable(oriented, "theta[1,1]"), extra[["theta[1,1]"]])
})

test_that("stable/ou (0.4.2+) need no orientation: empty flip list, and extract_theta() is unchanged whether or not draws are routed through bilatr_orient()", {
  for (name in c("stable", "ou")) {
    expect_identical(.bilatr_flip_variables(name), character(0))
  }

  # Simulate what extract_theta()'s in-memory branch does: a "fit"-like
  # object whose $draws() errors if bilatr_orient() ever asked it for
  # alpha[1] alongside theta (which the modern branch must not do).
  set.seed(5)
  n <- 15
  theta_vals <- stats::rnorm(n, 0.3, 0.1)
  fake_fit <- list(
    draws = function(variables) {
      if (!identical(variables, "theta")) {
        stop("extract_theta() requested more than 'theta' for a model with no reflection symmetry")
      }
      .make_synthetic_draws(list(`theta[1,1]` = theta_vals))
    }
  )
  stan_data <- structure(list(), dyad_ids = data.frame(dyad_id = 1L, time_index = 1L, dyad = "A_B"))

  theta <- extract_theta(fake_fit, stan_data, stan_model = "stable")
  expect_equal(theta$mean, mean(theta_vals))

  # bilatr_orient() itself, called directly on the same draws, must be a
  # true no-op (not just unreached) for a model with an empty flip list
  draws <- .make_synthetic_draws(list(`alpha[1]` = stats::rnorm(n, -5, 0.1), `theta[1,1]` = theta_vals))
  oriented <- bilatr_orient(draws, stan_model = "stable", variables = "theta")
  expect_equal(posterior::extract_variable(oriented, "theta[1,1]"), theta_vals)
})

test_that(".bilatr_flip_variables()/bilatr_orient() error on an unrecognized stan_model rather than silently no-op (B1)", {
  # B1: an unrecognized stan_model used to make .bilatr_flip_variables()
  # return character(0) -- indistinguishable from "this model has no
  # reflection symmetry" -- silently disabling orientation for a
  # wrong-basin fit with no error or warning. Both now go through
  # .canonical_stan_model() first, which errors instead.
  expect_error(.bilatr_flip_variables("some_unregistered_model"), "Unknown stan_model")

  set.seed(4)
  n <- 10
  extra <- list(`alpha[1]` = stats::rnorm(n, -3, 0.1), `alpha[2]` = stats::rnorm(n, 1, 0.1))
  draws <- .make_synthetic_draws(extra)
  expect_error(
    bilatr_orient(draws, stan_model = "some_unregistered_model", variables = "alpha"),
    "Unknown stan_model"
  )
})

test_that(".bilatr_flip_variables() accepts the pre-0.4.0 alphanorm/alphanorm_ou aliases and matches the legacy soft-anchor entries exactly", {
  .reset_bilatr_alias_messaged()
  expect_message(
    alphanorm_flip <- .bilatr_flip_variables("alphanorm"),
    "pre-0.4.0 name of 'stable_soft_anchor'"
  )
  expect_identical(alphanorm_flip, .bilatr_flip_variables("stable_soft_anchor"))

  expect_message(
    alphanorm_ou_flip <- .bilatr_flip_variables("alphanorm_ou"),
    "pre-0.4.0 name of 'ou_soft_anchor'"
  )
  expect_identical(alphanorm_ou_flip, .bilatr_flip_variables("ou_soft_anchor"))
})

test_that("bilatr_orient() errors informatively if alpha[1] is not present in draws", {
  n <- 5
  draws <- .make_synthetic_draws(list(`theta[1,1]` = stats::rnorm(n)))
  expect_error(
    bilatr_orient(draws, stan_model = "stable_soft_anchor", variables = "theta"),
    "alpha\\[1\\]"
  )
})

test_that(".warn_if_wrong_basin() fires when alpha[1] is negative, and bilatr_orient() recovers a positive orientation from that fit (legacy stable_soft_anchor)", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # Repointed to the legacy stable_soft_anchor program (0.4.2): that is
  # where a WRONG basin (as opposed to a merely sign-ambiguous one) still
  # exists -- the current `stable` folds alpha[1]'s sign into the
  # REPORTED alpha/theta (see inst/stan/bilatr_alphanorm.stan's header,
  # "IDENTIFICATION: ORIENTATION FOLD"), so there is no wrong orientation
  # for a caller to land in, even though its RAW alpha_raw remains
  # sign-ambiguous (see the dedicated fold-mechanism test below).
  set.seed(1)
  D <- 3
  Tn <- 5
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
    compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0, n_filter_dyads = 0, filter_dyads = integer(0), anchor_scale = 0.1
  )

  mod <- .compile_stan_model("stable_soft_anchor", opt_level = 1)

  # deliberately seed into the WRONG basin (opposite sign convention from
  # .legacy_alpha_raw_sum0_init()), and pin the sampler near its init
  # (adapt_engaged = FALSE, tiny step_size, capped max_treedepth) so a
  # short run can't cross the (data-scale-dependent) likelihood barrier
  # between the two modes regardless of how large it happens to be for
  # this particular synthetic dataset -- this is what makes the test
  # deterministic rather than dependent on the barrier actually holding
  # for a tiny dataset.
  bad_init <- function() {
    a <- stats::rnorm(A, 0, 0.5)
    a <- a - mean(a)
    if (a[1] > 0) a <- -a
    list(
      theta_raw = matrix(0, D, Tn), mu_intercept = rep(0, A), alpha_raw = a,
      sigma_theta0 = 0.5, z_theta0 = stats::rnorm(D, 0, 0.5),
      log_process_noise_raw = rep(0, D), mu_log_noise = log(0.2), sigma_log_noise = 0.3,
      phi = rep(1, D), mu_log_phi = 0, sigma_log_phi = 0.5
    )
  }

  fit_wrong <- suppressWarnings(mod$sample(
    data = data_list, chains = 1, iter_warmup = 2, iter_sampling = 5,
    seed = 1, refresh = 0, threads_per_chain = 1,
    adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
    init = bad_init, output_dir = tempdir(), show_messages = FALSE
  ))
  expect_lt(stats::median(posterior::extract_variable(fit_wrong$draws("alpha[1]"), "alpha[1]")), 0)

  expect_warning(
    expect_message(.warn_if_wrong_basin(fit_wrong, "stable_soft_anchor"), "Posterior median of alpha\\[1\\]"),
    "wrong-sign basin"
  )

  oriented <- bilatr_orient(fit_wrong$draws(variables = "alpha"), stan_model = "stable_soft_anchor", variables = "alpha")
  expect_gt(stats::median(posterior::extract_variable(oriented, "alpha[1]")), 0)

  # and the right-basin case (bilatr_init_fn()'s actual, anchored init)
  # should report but not warn
  good_init <- bilatr_init_fn(list(D = D, T = Tn, A = A), stan_model = "stable_soft_anchor")
  fit_right <- suppressWarnings(mod$sample(
    data = data_list, chains = 1, iter_warmup = 20, iter_sampling = 5,
    seed = 1, refresh = 0, threads_per_chain = 1,
    init = good_init, output_dir = tempdir(), show_messages = FALSE
  ))
  expect_no_warning(expect_message(.warn_if_wrong_basin(fit_right, "stable_soft_anchor"), "Posterior median of alpha\\[1\\]"))
})

test_that("the orientation fold reports alpha[1] > 0 and agreeing alpha/theta from two chains pinned in opposite alpha_raw basins, while alpha_raw/z_theta0/theta_raw's own Rhat stays large", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # This is the actual proof the orientation fold -- not coincidence --
  # is doing the work (0.4.2b; see inst/stan/bilatr_alphanorm.stan's
  # header, "IDENTIFICATION: ORIENTATION FOLD"). Two chains are given
  # EXACT mirror-image raw-parameter inits (alpha_raw, z_theta0, theta_raw
  # all negated between the two -- same exact-mirror-image construction
  # test_stan_smoke.R's log-prob-gap tests use, just as an init here
  # instead of a log_prob() comparison), each pinned via
  # adapt_engaged = FALSE / tiny step_size / small max_treedepth so
  # neither can move far from its init during a short run (this
  # file's established pattern -- see the wrong-basin test above). If the
  # fold is doing its job, the REPORTED alpha/theta must come out
  # (near-)identical between the two chains regardless -- while the RAW
  # alpha_raw/z_theta0/theta_raw, which the fold does NOT correct, must
  # show the opposite: large Rhat, since one chain's raw values are the
  # exact negation of the other's.
  set.seed(1)
  D <- 3
  Tn <- 4
  A <- 5
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0, n_filter_dyads = 0, filter_dyads = integer(0)
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
        theta_raw = sign_mult * theta_raw_pos,
        mu_intercept = rep(0, A),
        alpha_raw = sign_mult * alpha_raw_pos,
        sigma_theta0 = 0.5,
        z_theta0 = sign_mult * z_theta0_pos,
        log_process_noise_raw = rep(0, D),
        mu_log_noise = log(0.2), sigma_log_noise = 0.3,
        phi = rep(1, D), mu_log_phi = 0, sigma_log_phi = 0.5
      )
    }
  }

  fit_one <- function(sign_mult, seed) {
    suppressWarnings(mod$sample(
      data = data_list, chains = 1, iter_warmup = 2, iter_sampling = 50,
      seed = seed, refresh = 0, threads_per_chain = 1,
      adapt_engaged = FALSE, step_size = 0.001, max_treedepth = 2,
      init = make_init(sign_mult), output_dir = tempdir(), show_messages = FALSE
    ))
  }

  fit_pos <- fit_one(1, seed = 1)
  fit_neg <- fit_one(-1, seed = 2)

  # confirm the fixture actually landed in opposite RAW basins
  expect_true(all(posterior::extract_variable(fit_pos$draws("alpha_raw[1]"), "alpha_raw[1]") > 0))
  expect_true(all(posterior::extract_variable(fit_neg$draws("alpha_raw[1]"), "alpha_raw[1]") < 0))

  vars <- c("alpha", "theta", "alpha_raw", "z_theta0", "theta_raw")
  combined <- posterior::bind_draws(
    fit_pos$draws(variables = vars), fit_neg$draws(variables = vars),
    along = "chain"
  )

  # assertion 1: alpha[1] > 0 in BOTH chains
  expect_true(all(posterior::extract_variable_matrix(combined, "alpha[1]") > 0))

  summ <- posterior::summarise_draws(combined, mean = mean, rhat = posterior::rhat)
  reported <- summ[grepl("^(alpha|theta)\\[", summ$variable), ]
  raw <- summ[grepl("^(alpha_raw|z_theta0|theta_raw)\\[", summ$variable), ]
  expect_gt(nrow(reported), 0)
  expect_gt(nrow(raw), 0)

  # assertion 2: alpha/theta summaries (means, across the pooled 2-chain
  # draws) agree closely between the mirror-image inits -- checked
  # directly here, not just implied by assertion 3's Rhat
  mean_pos <- posterior::summarise_draws(fit_pos$draws(variables = vars), mean = mean)
  mean_neg <- posterior::summarise_draws(fit_neg$draws(variables = vars), mean = mean)
  reported_pos <- mean_pos$mean[grepl("^(alpha|theta)\\[", mean_pos$variable)]
  reported_neg <- mean_neg$mean[grepl("^(alpha|theta)\\[", mean_neg$variable)]
  expect_equal(reported_pos, reported_neg, tolerance = 0.05)

  # assertion 3: Rhat on alpha/theta is small -- RELATIVE to the raw
  # parameters' own Rhat, not an absolute "converged" threshold. With both
  # chains deliberately pinned near their (mirror-image) inits via
  # adapt_engaged = FALSE and a tiny step_size, within-chain variance is
  # tiny for EVERY variable, raw and reported alike (that is the whole
  # point of the pinning trick -- see the wrong-basin test above), which
  # makes split-Rhat's absolute scale hypersensitive to residual
  # floating-point-scale noise regardless of variable (confirmed
  # empirically: an absolute Rhat < 1.1 threshold is not achievable here
  # for ANY variable, reported or raw, without abandoning the pinning that
  # keeps the two chains deterministically separated). The qualitative
  # claim that actually matters -- and that holds robustly across the
  # configurations checked while writing this test -- is that the
  # REPORTED quantities agree substantially better, on average, than the
  # RAW ones do.
  expect_lt(mean(reported$rhat), mean(raw$rhat))
  expect_lt(mean(reported$rhat), 1.8)

  # assertion 4 -- THE assertion that proves the fold, not coincidence, is
  # doing the work: Rhat on the raw parameters is large. If this were
  # small too, the "agreement" above could just mean the two chains
  # happened to converge to the same point, not that the fold is actively
  # reconciling two genuinely different raw-space basins.
  expect_true(all(raw$rhat > 1.5))
})

test_that("legacy stable_soft_anchor and the new stable agree on the identified quantities from the same data (post-orientation), confirming this is a re-parameterization, not a model change", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  set.seed(1)
  D <- 6
  Tn <- 5
  A <- 5
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    dyad_weight = rep(1, D), period_weight = rep(1, Tn), action_weight = rep(1, A),
    compute_log_lik = 0, prior_only = 0, compute_theta_filtered = 0, n_filter_dyads = 0, filter_dyads = integer(0), anchor_scale = 0.1
  )

  fit_one <- function(stan_model) {
    mod <- .compile_stan_model(stan_model, opt_level = 1)
    suppressWarnings(suppressMessages(mod$sample(
      data = data_list, chains = 2, parallel_chains = 2,
      iter_warmup = 200, iter_sampling = 200, seed = 1, refresh = 0,
      threads_per_chain = 1, init = bilatr_init_fn(list(D = D, T = Tn, A = A), stan_model = stan_model),
      output_dir = tempdir(), show_messages = FALSE
    )))
  }

  fit_legacy <- fit_one("stable_soft_anchor")
  fit_new <- fit_one("stable")

  # alpha/theta: orient the legacy fit (last legitimate use of
  # bilatr_orient() -- see NEWS.md), compare directly to the new fit
  oriented_legacy <- bilatr_orient(
    fit_legacy$draws(variables = c("alpha[1]", "alpha", "theta")),
    stan_model = "stable_soft_anchor", variables = c("alpha", "theta")
  )
  new_draws <- fit_new$draws(variables = c("alpha", "theta"))

  legacy_mean <- posterior::summarise_draws(oriented_legacy, mean = mean)
  new_mean <- posterior::summarise_draws(new_draws, mean = mean)
  cmp <- merge(legacy_mean, new_mean, by = "variable", suffixes = c("_legacy", "_new"))
  expect_gt(nrow(cmp), 0)
  expect_gt(stats::cor(cmp$mean_legacy, cmp$mean_new), 0.9)
  expect_lt(max(abs(cmp$mean_legacy - cmp$mean_new)), 0.5)

  # mu_intercept/phi/sigma_theta0: never flip, compare directly, no
  # orientation needed
  other_legacy <- posterior::summarise_draws(
    fit_legacy$draws(variables = c("mu_intercept", "phi", "sigma_theta0")), mean = mean
  )
  other_new <- posterior::summarise_draws(
    fit_new$draws(variables = c("mu_intercept", "phi", "sigma_theta0")), mean = mean
  )
  cmp_other <- merge(other_legacy, other_new, by = "variable", suffixes = c("_legacy", "_new"))
  expect_gt(nrow(cmp_other), 0)
  expect_gt(stats::cor(cmp_other$mean_legacy, cmp_other$mean_new), 0.8)
})

test_that(".warn_if_wrong_basin() is a true no-op for stable/ou (0.4.2+): no fit$draws() call at all", {
  # Unlike the legacy stable_soft_anchor case above, stable/ou have no
  # reflection symmetry left to check -- confirmed here by a fake fit
  # whose $draws() errors if ever called.
  fake_fit <- list(draws = function(...) stop("should not be called"))
  for (name in c("stable", "ou")) {
    expect_no_message(expect_no_warning(.warn_if_wrong_basin(fake_fit, name)))
  }
})

test_that(".warn_if_wrong_basin() errors on an unrecognized stan_model rather than silently no-op (B1)", {
  # .warn_if_wrong_basin() decides via .bilatr_flip_variables(), which
  # itself validates stan_model through .canonical_stan_model() -- an
  # unrecognized name errors rather than being silently treated as "no
  # reflection symmetry, nothing to check" (B1). In practice this
  # stan_model has already been validated by .compile_stan_model()
  # earlier in fit_bilatr(), so this case shouldn't arise from a real
  # call, but the guard must not paper over it if it somehow did.
  fake_fit <- list(draws = function(...) stop("should not be called"))
  expect_error(.warn_if_wrong_basin(fake_fit, "some_unregistered_model"), "Unknown stan_model")
})

test_that(".warn_if_wrong_basin() accepts the pre-0.4.0 alphanorm/alphanorm_ou aliases", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # alphanorm/alphanorm_ou resolve to the legacy stable_soft_anchor/
  # ou_soft_anchor entries (0.4.2+), which DO still have a reflection
  # symmetry -- a fake fit whose $draws() would error if ever called
  # confirms .warn_if_wrong_basin() treats the alias exactly like its
  # canonical name (this does NOT return early; it must reach $draws())
  fake_fit <- list(draws = function(...) stop("should not be called"))
  expect_error(
    suppressMessages(.warn_if_wrong_basin(fake_fit, "alphanorm")),
    "should not be called"
  )
})
