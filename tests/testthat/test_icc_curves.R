.make_icc_fixture <- function(A = 6, D = 5, Tn = 4, seed = 61, n_pseudo_draws = 50) {
  truth <- .make_fake_residual_data(D = D, Tn = Tn, A = A, seed = seed, beta = NULL)
  fit <- .make_fake_residual_fit(truth, n_pseudo_draws = n_pseudo_draws)
  list(truth = truth, fit = fit)
}

test_that("icc_curves(type = 'probability') matches a hand-computed softmax at known alpha/mu/theta", {
  fx <- .make_icc_fixture()
  res <- icc_curves(fx$fit, theta_range = c(-1.5, 1.5), n_grid = 7, probs = c(0.05, 0.95))
  expect_s3_class(res, "bilatr_icc_curves")
  expect_equal(nrow(res), 7 * length(fx$truth$alpha))

  theta_star <- res$theta[4] # a grid point, exact since evenly spaced
  p_expected <- bilatr:::.softmax(fx$truth$alpha * theta_star - fx$truth$mu_intercept)
  got <- dplyr::filter(res, abs(.data$theta - theta_star) < 1e-8) %>%
    dplyr::arrange(.data$action_index)
  expect_equal(got$median, p_expected, tolerance = 1e-6)
  # a point-mass "posterior" (see .make_fake_residual_fit()) has zero
  # posterior spread, so the interval collapses onto the median exactly
  expect_equal(got$lower, p_expected, tolerance = 1e-6)
  expect_equal(got$upper, p_expected, tolerance = 1e-6)
})

test_that("icc_curves(type = 'information') matches .var_pi_alpha() directly at a spot-check theta", {
  fx <- .make_icc_fixture()
  res <- icc_curves(fx$fit, theta_range = c(-1.5, 1.5), n_grid = 9, type = "information")
  expect_true(!is.null(attr(res, "total_information")))

  theta_star <- res$theta[3]
  p <- bilatr:::.softmax(fx$truth$alpha * theta_star - fx$truth$mu_intercept)
  vp <- bilatr:::.var_pi_alpha(matrix(fx$truth$alpha, nrow = 1), matrix(p, nrow = 1))
  expected_contrib <- as.vector(vp$info_share * vp$var_pi)

  got <- dplyr::filter(res, abs(.data$theta - theta_star) < 1e-8) %>%
    dplyr::arrange(.data$action_index)
  expect_equal(got$median, expected_contrib, tolerance = 1e-6)

  total_at_star <- dplyr::filter(attr(res, "total_information"), abs(.data$theta - theta_star) < 1e-8)
  expect_equal(total_at_star$median, unname(vp$var_pi), tolerance = 1e-6)
})

test_that("icc_crossings() matches the closed-form formula and where the curves actually cross", {
  fx <- .make_icc_fixture(A = 6)
  res <- icc_curves(fx$fit, theta_range = c(-3, 3), n_grid = 4001)
  cross <- icc_crossings(res)
  expect_equal(nrow(cross), choose(6, 2))

  # pick a pair with a well-separated alpha and check the formula directly
  alpha <- fx$truth$alpha
  mu <- fx$truth$mu_intercept
  j <- 1L
  k <- which.max(abs(alpha - alpha[1]))
  theta_star <- (mu[j] - mu[k]) / (alpha[j] - alpha[k])

  row <- dplyr::filter(cross, .data$index_a == min(j, k), .data$index_b == max(j, k))
  expect_equal(row$crossing_median, theta_star, tolerance = 1e-6)

  # and confirm the two categories' probabilities are actually equal there
  p_star <- bilatr:::.softmax(alpha * theta_star - mu)
  expect_equal(p_star[j], p_star[k], tolerance = 1e-8)

  if (row$in_range) {
    expect_true(theta_star >= -3 && theta_star <= 3)
  }
})

test_that("theta_range resolution: explicit range takes precedence, then theta_summary, then a messaged fallback", {
  fx <- .make_icc_fixture(D = 8, Tn = 3)

  r1 <- icc_curves(fx$fit, theta_range = c(-2, 2), n_grid = 3)
  expect_equal(attr(r1, "theta_range"), c(-2, 2))

  ts <- tibble::tibble(mean = seq(-5, 5, length.out = 100))
  r2 <- icc_curves(fx$fit, theta_summary = ts, quantile_probs = c(0.1, 0.9), n_grid = 3)
  expect_equal(attr(r2, "theta_range"), unname(stats::quantile(ts$mean, c(0.1, 0.9))))

  expect_message(
    r3 <- icc_curves(fx$fit, stan_data = fx$truth$stan_data, n_grid = 3, n_dyads_fallback = 5, seed = 2),
    "last resort"
  )
  expect_true(all(is.finite(attr(r3, "theta_range"))))

  expect_error(
    icc_curves(fx$fit, n_grid = 3),
    "theta_range"
  )
})

test_that("categories subsets the output to just the requested action classes", {
  fx <- .make_icc_fixture(A = 6)
  res <- icc_curves(fx$fit, theta_range = c(-1, 1), n_grid = 5, categories = c("2", "4"))
  expect_setequal(unique(res$action_index), c(2L, 4L))
  expect_equal(nrow(res), 5 * 2)
})

# --- CHECKPOINT 2: real CmdStan fixture, both types, plots, crossings ---

test_that("icc_curves()/icc_crossings() run end to end on a real CmdStan fixture, both types, with print/autoplot/plot", {
  skip_if_no_cmdstan()
  fx <- make_csv_diagnostics_fixture()

  res_prob <- icc_curves(fx$fit, stan_data = fx$stan_data, n_grid = 21, type = "probability")
  expect_s3_class(res_prob, "bilatr_icc_curves")
  expect_output(print(res_prob), "bilatr_icc_curves")

  res_info <- icc_curves(fx$fit, theta_range = attr(res_prob, "theta_range"), n_grid = 21, type = "information")
  expect_true(!is.null(attr(res_info, "total_information")))

  cross <- icc_crossings(res_info)
  expect_true(all(c("crossing_median", "in_range") %in% names(cross)))

  skip_if_not_installed("ggplot2")
  p1 <- ggplot2::autoplot(res_prob)
  p2 <- ggplot2::autoplot(res_info)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})
