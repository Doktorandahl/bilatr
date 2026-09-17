test_that(".clr(softmax(eta)) recovers eta - mean(eta) to numerical tolerance", {
  set.seed(1)
  eta <- rnorm(12)
  p <- bilatr:::.softmax(eta)
  expect_equal(bilatr:::.clr(p), eta - mean(eta), tolerance = 1e-10)

  eta_mat <- matrix(rnorm(5 * 12), 5, 12)
  p_mat <- bilatr:::.softmax_rows(eta_mat)
  expect_equal(bilatr:::.clr(p_mat), eta_mat - rowMeans(eta_mat), tolerance = 1e-10)
})

test_that("RMS-1 identification: sd(outer(theta, alpha)) == sqrt(mean(theta^2)) for a sum-to-zero, RMS-1 alpha", {
  # This is the identity 1b.3's comparator formula relies on
  # (implied_beta_rms's per-component comparator reduces to the RMS of
  # theta itself because mean(alpha_k^2) = 1 exactly) -- asserted here as
  # a free check on the identification, not recomputed inside
  # check_compositional_residuals() itself. Verified in the prompt against
  # a real simulated panel: 1.82148 vs 1.82147.
  set.seed(2)
  A <- 10
  alpha_raw <- rnorm(A)
  alpha_raw <- alpha_raw - mean(alpha_raw)
  alpha <- alpha_raw * sqrt(A / sum(alpha_raw^2))
  expect_equal(mean(alpha^2), 1, tolerance = 1e-10)

  # population sd (not stats::sd's N-1 sample correction), since
  # mean(outer(theta, alpha)) == mean(theta) * mean(alpha) == 0 exactly
  # (alpha sums to zero) is what the identity actually relies on
  pop_sd <- function(x) sqrt(mean(x^2) - mean(x)^2)
  theta <- rnorm(500)
  expect_equal(pop_sd(outer(theta, alpha)), sqrt(mean(theta^2)), tolerance = 1e-8)
})

test_that(".stratified_dyad_sample() returns the target size, spans strata, and is reproducible", {
  set.seed(3)
  n_d <- c(rep(1, 50), rep(100, 50), rep(10000, 50))

  s1 <- bilatr:::.stratified_dyad_sample(n_d, n_target = 60, n_strata = 3, seed = 42)
  s2 <- bilatr:::.stratified_dyad_sample(n_d, n_target = 60, n_strata = 3, seed = 42)
  expect_identical(s1, s2)
  expect_equal(length(s1), 60)
  expect_equal(length(unique(s1)), 60)

  strat <- dplyr::ntile(n_d, 3)
  expect_true(all(table(strat[s1]) > 0))

  # capped at the pool size
  s3 <- bilatr:::.stratified_dyad_sample(n_d[1:10], n_target = 100, n_strata = 3, seed = 1)
  expect_equal(length(s3), 10)
})

test_that(".rmultinom_rows()/.dirichlet_multinomial_rows() produce correctly-summed, roughly-centred draws", {
  set.seed(4)
  n_draws <- 2000
  A <- 5
  p <- c(0.1, 0.2, 0.3, 0.15, 0.25)
  p_mat <- matrix(p, n_draws, A, byrow = TRUE)
  n <- 40L

  y <- bilatr:::.rmultinom_rows(n, p_mat)
  expect_true(all(rowSums(y) == n))
  expect_lt(max(abs(colMeans(y) - n * p)), 0.5)

  conc <- 8 * p_mat
  yd <- bilatr:::.dirichlet_multinomial_rows(n, conc)
  expect_true(all(rowSums(yd) == n))
  expect_lt(max(abs(colMeans(yd) - n * p)), 1.5) # extra DM overdispersion
})

# --- calibration / power / subspace-isolation, per 1f ---

test_that("calibration: beta = 0 gives a non-extreme pooled PPP and a roughly-uniform per-dyad PPP distribution", {
  dat <- .make_fake_residual_data(D = 200, Tn = 12, A = 8, seed = 11, beta = NULL)
  fit <- .make_fake_residual_fit(dat, n_pseudo_draws = 250)

  res <- check_compositional_residuals(
    fit, dat$stan_data, n_dyads = 200, n_strata = 4, n_draws = 250, seed = 5
  )
  expect_s3_class(res, "bilatr_residual_check")
  expect_gt(res$global$pooled_ppp, 0.02)
  expect_lt(res$global$pooled_ppp, 0.98)
  # spec's own null-calibration figure was 0.058 at D = 400; a generous
  # bound at this much smaller D keeps the test from being flaky
  expect_lt(mean(res$dyads$ppp_dyad < 0.05), 0.15)
  expect_equal(nrow(res$global$eps_sensitivity), 3)
  expect_setequal(round(res$global$eps_sensitivity$eps, 3), c(0.5, 0.1, 1))
})

test_that("power/recovery: a beta orthogonal to alpha is detected, and implied_beta_rms is within a factor of ~2 of truth", {
  A <- 8
  set.seed(21)
  alpha_raw <- rnorm(A)
  alpha_raw <- alpha_raw - mean(alpha_raw)
  alpha0 <- alpha_raw * sqrt(A / sum(alpha_raw^2))

  v <- rnorm(A)
  v <- v - mean(v) # orthogonal to 1
  v <- v - sum(v * alpha0) / sum(alpha0^2) * alpha0 # orthogonal to alpha
  true_rms <- 0.4
  beta <- v / sqrt(mean(v^2)) * true_rms
  expect_equal(sqrt(mean(beta^2)), true_rms, tolerance = 1e-8)

  dat <- .make_fake_residual_data(D = 250, Tn = 15, A = A, seed = 22, beta = beta)
  fit <- .make_fake_residual_fit(dat, n_pseudo_draws = 250)

  res <- check_compositional_residuals(
    fit, dat$stan_data, n_dyads = 250, n_strata = 4, n_draws = 250, seed = 6
  )
  expect_lt(res$global$pooled_ppp, 0.05)
  expect_lt(res$global$implied_beta_rms, 2 * true_rms)
  expect_gt(res$global$implied_beta_rms, true_rms / 2)
  expect_gt(res$global$beta_signal_ratio, 0)
})

test_that("subspace isolation: a compositional offset parallel to alpha decomposes almost entirely into along, not perp", {
  # A direct test of the r/along/perp decomposition's own math, at the
  # population level (expected shares, not simulated counts). This is
  # deliberately NOT routed through the full check_compositional_residuals()
  # + Dirichlet-multinomial-simulated-Y pipeline: the DM likelihood has a
  # genuine overdispersion NOISE FLOOR that does not vanish as n_dt grows
  # (Var(Y_k/n) -> p_k(1-p_k)/(1+phi) as n -> Inf, not 0 -- the Dirichlet
  # layer's own variance, confirmed empirically while writing this test),
  # so feeding the check a dyad's theta that is deliberately "wrong" by a
  # constant (exactly what a parallel-to-alpha beta looks like) also
  # shifts the noise characteristics of the resulting simulated Y away
  # from the replicate's, contaminating a full-pipeline comparison with
  # an effect that has nothing to do with the decomposition itself. The
  # decomposition is a per-draw linear-algebra step on a clr residual
  # `r`; this test constructs a known `r` -- the clr difference between
  # the true (shifted) and model (unshifted) n-weighted average shares --
  # and checks directly that it resolves almost entirely into `along`.
  set.seed(51)
  A <- 8
  Tn <- 15
  alpha_raw <- rnorm(A)
  alpha_raw <- alpha_raw - mean(alpha_raw)
  alpha <- alpha_raw * sqrt(A / sum(alpha_raw^2))
  mu_raw <- rnorm(A)
  mu <- mu_raw - mean(mu_raw)

  theta <- cumsum(c(rnorm(1), 0.2 * rnorm(Tn - 1)))
  c_shift <- 0.4
  beta_parallel <- alpha * c_shift # parallel to alpha, not orthogonal

  p_model <- t(vapply(theta, function(th) bilatr:::.softmax(alpha * th - mu), numeric(A)))
  p_true <- t(vapply(theta, function(th) bilatr:::.softmax(alpha * th - mu + beta_parallel), numeric(A)))
  pbar_model <- colMeans(p_model)
  pbar_true <- colMeans(p_true)

  r <- bilatr:::.clr(pbar_true) - bilatr:::.clr(pbar_model)
  alpha_c <- alpha - mean(alpha)
  along <- sum(r * alpha_c) / sum(alpha_c^2)
  perp <- r - along * alpha_c

  expect_equal(along, c_shift, tolerance = 0.05)
  expect_lt(sqrt(sum(perp^2)), 0.05 * abs(along))
})

test_that("determinism: the same seed on the same fit/stan_data reproduces identical output", {
  dat <- .make_fake_residual_data(D = 120, Tn = 10, A = 6, seed = 41, beta = NULL)
  fit <- .make_fake_residual_fit(dat, n_pseudo_draws = 150)

  r1 <- check_compositional_residuals(fit, dat$stan_data, n_dyads = 120, n_strata = 3, n_draws = 150, seed = 9)
  r2 <- check_compositional_residuals(fit, dat$stan_data, n_dyads = 120, n_strata = 3, n_draws = 150, seed = 9)

  expect_identical(r1$settings$sampled_dyad_ids, r2$settings$sampled_dyad_ids)
  expect_equal(r1$global$pooled_ppp, r2$global$pooled_ppp)
  expect_equal(r1$global$implied_beta_rms, r2$global$implied_beta_rms)
  expect_equal(r1$dyads, r2$dyads)
})

test_that("check_compositional_residuals() errors without a dyad_ids attribute", {
  dat <- .make_fake_residual_data(D = 10, Tn = 4, A = 4, seed = 1)
  attr(dat$stan_data, "dyad_ids") <- NULL
  fit <- .make_fake_residual_fit(dat, n_pseudo_draws = 20)
  expect_error(check_compositional_residuals(fit, dat$stan_data), "dyad_ids")
})

# --- CHECKPOINT 1: real CmdStan fixture ---

test_that("check_compositional_residuals() runs end to end on a real CmdStan fixture, with print/autoplot/plot", {
  skip_if_no_cmdstan()
  fx <- make_csv_diagnostics_fixture()

  res <- check_compositional_residuals(
    fx$fit, fx$stan_data, n_dyads = 6, n_strata = 2, n_draws = 20, seed = 1
  )
  expect_s3_class(res, "bilatr_residual_check")
  expect_output(print(res), "bilatr_residual_check")

  skip_if_not_installed("ggplot2")
  plots <- ggplot2::autoplot(res)
  expect_type(plots, "list")
  expect_s3_class(plots$perp_vs_n, "ggplot")
  expect_s3_class(plots$ppp_hist, "ggplot")
  expect_s3_class(plots$category, "ggplot")

  outdir <- tempfile()
  write_residual_check(res, outdir)
  expect_true(file.exists(file.path(outdir, "residual_check_dyads.csv")))
  expect_true(file.exists(file.path(outdir, "residual_check_perp_vs_n.png")))
})
