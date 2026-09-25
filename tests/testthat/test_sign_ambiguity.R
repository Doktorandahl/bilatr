test_that(".bilatr_sign_ambiguous_raw_names()'s output is pinned: the runscripts call it directly by name", {
  # This function's name and output must not change across versions (see
  # its own docs in R/diagnose_convergence.R) -- the runscripts call
  # bilatr:::.bilatr_sign_ambiguous_raw_names() directly for job sizing.
  expect_setequal(
    .bilatr_sign_ambiguous_raw_names(),
    c("alpha_raw", "z_theta0", "theta_raw", "gamma_z", "mu_dyad_raw")
  )
})

test_that(".bilatr_sign_tied_names() returns the raw sign-tied names per model, character(0) for an unregistered symmetry", {
  expect_setequal(.bilatr_sign_tied_names("stable")$raw, c("alpha_raw", "z_theta0", "theta_raw"))
  expect_setequal(.bilatr_sign_tied_names("ou")$raw, c("alpha_raw", "mu_dyad_raw", "theta_raw"))
  expect_setequal(
    .bilatr_sign_tied_names("stable_gamma")$raw,
    c("alpha_raw", "z_theta0", "theta_raw", "gamma_z")
  )
})

test_that("the orientation fold reports alpha[1] > 0 and agreeing alpha/theta from two chains pinned in opposite alpha_raw basins, while alpha_raw/z_theta0/theta_raw's own Rhat stays large", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  # This is the actual proof the orientation fold -- not coincidence --
  # is doing the work (0.4.2b; see inst/stan/bilatr_alphanorm.stan's
  # header, "IDENTIFICATION: ORIENTATION FOLD"). Two chains are given
  # EXACT mirror-image raw-parameter inits (alpha_raw, z_theta0, theta_raw
  # all negated between the two), each pinned via
  # adapt_engaged = FALSE / tiny step_size / small max_treedepth so
  # neither can move far from its init during a short run. If the
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
  # point of the pinning trick), which makes split-Rhat's absolute scale
  # hypersensitive to residual floating-point-scale noise regardless of
  # variable (confirmed empirically: an absolute Rhat < 1.1 threshold is
  # not achievable here for ANY variable, reported or raw, without
  # abandoning the pinning that keeps the two chains deterministically
  # separated). The qualitative claim that actually matters -- and that
  # holds robustly across the configurations checked while writing this
  # test -- is that the REPORTED quantities agree substantially better,
  # on average, than the RAW ones do.
  expect_lt(mean(reported$rhat), mean(raw$rhat))
  expect_lt(mean(reported$rhat), 1.8)

  # assertion 4 -- THE assertion that proves the fold, not coincidence, is
  # doing the work: Rhat on the raw parameters is large. If this were
  # small too, the "agreement" above could just mean the two chains
  # happened to converge to the same point, not that the fold is actively
  # reconciling two genuinely different raw-space basins.
  expect_true(all(raw$rhat > 1.5))
})
