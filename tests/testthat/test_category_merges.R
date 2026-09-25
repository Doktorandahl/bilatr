test_that("alpha_prior_moments() matches the closed-form analytic result at A = 4 and A = 18", {
  # Closed form (derived from alpha = sqrt(A) * U, U uniform on the unit
  # sphere in the sum-zero hyperplane; U_1's marginal is a scaled Beta):
  #   E[alpha_1] = sqrt(A-1) * Gamma((A-1)/2) / (Gamma(A/2) * sqrt(pi))
  #   E[alpha_k], k>1 = -E[alpha_1] / (A-1)
  #   Var(alpha_1) = 1 - E[alpha_1]^2 ; Var(alpha_k) = 1 - E[alpha_k]^2
  # Checked by hand against the task's own reported A=10 numbers (mean
  # 0.820/sd 0.572 for the reference class, mean -0.091/sd 0.996 for the
  # rest) before writing this test.
  analytic <- function(A) {
    e1 <- sqrt(A - 1) * gamma((A - 1) / 2) / (gamma(A / 2) * sqrt(pi))
    erest <- -e1 / (A - 1)
    list(
      mean = c(e1, rep(erest, A - 1)),
      sd = sqrt(1 - c(e1, rep(erest, A - 1))^2)
    )
  }

  for (A in c(4, 18)) {
    got <- alpha_prior_moments(A, n_sim = 2e6, seed = 1)
    want <- analytic(A)
    # absolute, not expect_equal()'s default mean-relative-difference --
    # several components are near zero, where a relative comparison is
    # dominated by noise even when the absolute deviation is tiny
    expect_lt(max(abs(got$prior_mean - want$mean)), 0.01)
    expect_lt(max(abs(got$prior_sd - want$sd)), 0.01)
  }
})

test_that("alpha_prior_moments() is memoised: repeated calls with the same arguments return identical, cached results", {
  a <- alpha_prior_moments(5, n_sim = 1e4, seed = 2)
  b <- alpha_prior_moments(5, n_sim = 1e4, seed = 2)
  expect_identical(a, b)
})

test_that(".simulate_alpha_prior() satisfies the exact per-draw sum-to-zero identity, not just in expectation", {
  # alpha sums to exactly 0 in every draw (fold multiplies ALL coordinates
  # by the same sign, so it never breaks this), hence
  # mean(alpha[, -1]) == -alpha[, 1] / (A - 1) exactly, no Monte Carlo
  # error -- a much stronger check than comparing summary statistics.
  A <- 7
  alpha <- .simulate_alpha_prior(A, n_sim = 500, seed = 3)
  expect_equal(rowMeans(alpha[, -1, drop = FALSE]), -alpha[, 1] / (A - 1), tolerance = 1e-10)
})

test_that("loss(j, k) matches a brute-force recomputation of Var_pi(alpha) before and after merging", {
  set.seed(4)
  A <- 5
  alpha <- rnorm(A)
  shares <- runif(A)
  shares <- shares / sum(shares)

  alpha_bar <- sum(shares * alpha)
  var_before <- sum(shares * (alpha - alpha_bar)^2)

  j <- 2
  k <- 4
  pi_jk <- shares[j] + shares[k]
  alpha_jk <- (shares[j] * alpha[j] + shares[k] * alpha[k]) / pi_jk
  new_shares <- shares[-c(j, k)]
  new_alpha <- alpha[-c(j, k)]
  new_shares <- c(new_shares, pi_jk)
  new_alpha <- c(new_alpha, alpha_jk)
  new_alpha_bar <- sum(new_shares * new_alpha) # should equal alpha_bar
  var_after <- sum(new_shares * (new_alpha - new_alpha_bar)^2)

  brute_force_loss <- var_before - var_after
  formula_loss <- shares[j] * shares[k] * (alpha[j] - alpha[k])^2 / pi_jk

  expect_equal(new_alpha_bar, alpha_bar, tolerance = 1e-10)
  expect_equal(brute_force_loss, formula_loss, tolerance = 1e-10)

  # and against the package's own pairwise helper
  pw <- bilatr:::.pairwise_merge_loss(matrix(alpha, nrow = 1), shares)
  col <- which(pw$pairs[1, ] == j & pw$pairs[2, ] == k)
  expect_equal(pw$loss[1, col], formula_loss, tolerance = 1e-10)
})

test_that("loss(S) for a two-element S equals the pairwise formula", {
  set.seed(5)
  A <- 6
  alpha_mat <- matrix(rnorm(A), nrow = 1)
  shares <- runif(A)
  shares <- shares / sum(shares)
  vp <- bilatr:::.var_pi_alpha(alpha_mat, shares)

  group_loss <- bilatr:::.group_merge_loss(alpha_mat, shares, vp$alpha_bar, c(3, 5))
  pairwise_loss <- shares[3] * shares[5] * (alpha_mat[1, 3] - alpha_mat[1, 5])^2 / (shares[3] + shares[5])

  expect_equal(group_loss, pairwise_loss, tolerance = 1e-10)
})

test_that("the greedy ladder reproduces stats::hclust(method = 'ward.D')'s merge order on a small synthetic case", {
  # hclust's `members` argument does NOT reweight the leaf-level merge
  # heights the way one might expect for weighted Ward -- verified this
  # empirically while writing this test. The construction that DOES
  # reproduce weighted Ward exactly is replication: represent each
  # category k by round(pi_k * Nrep) unit-weight copies of alpha_k, and
  # run plain unweighted hclust on the replicated points.
  set.seed(42)
  A <- 6
  alpha <- rnorm(A)
  shares <- runif(A)
  shares <- shares / sum(shares)

  ladder <- bilatr:::.greedy_merge_ladder(
    matrix(alpha, nrow = 1), shares, as.character(seq_len(A)), probs = c(0.05, 0.5, 0.95)
  )

  n_rep <- 20000
  counts <- round(shares * n_rep)
  rep_vals <- rep(alpha, counts)
  rep_cat <- rep(seq_len(A), counts)

  hc <- stats::hclust(stats::dist(rep_vals)^2, method = "ward.D")
  n <- length(rep_vals)
  cluster_cats <- vector("list", n - 1)
  get_cats <- function(idx) if (idx < 0) rep_cat[-idx] else cluster_cats[[idx]]

  # only the LAST (A - 1) merges combine genuinely distinct categories
  # (all earlier merges just combine same-category replicate copies,
  # height 0) -- collect those, in order, as this construction's answer.
  real_merges <- list()
  for (i in seq_len(nrow(hc$merge))) {
    a <- hc$merge[i, 1]
    b <- hc$merge[i, 2]
    cats_a <- unique(get_cats(a))
    cats_b <- unique(get_cats(b))
    cats <- unique(c(cats_a, cats_b))
    cluster_cats[[i]] <- cats
    if (!setequal(cats_a, cats_b) || length(cats_a) > 1 || length(cats_b) > 1 || length(cats) > length(cats_a)) {
      if (!(length(cats_a) == 1 && length(cats_b) == 1 && cats_a == cats_b)) {
        real_merges[[length(real_merges) + 1]] <- list(a = sort(cats_a), b = sort(cats_b))
      }
    }
  }
  # keep only genuine cross-category merges (drop any same-category
  # replicate merges that slipped through)
  real_merges <- Filter(function(m) !setequal(m$a, m$b), real_merges)

  expect_equal(length(real_merges), A - 1)

  ladder_pairs <- Map(function(a, b) list(a = sort(as.integer(strsplit(a, "\\+")[[1]])), b = sort(as.integer(strsplit(b, "\\+")[[1]]))),
    ladder$merged_a, ladder$merged_b)

  for (i in seq_len(A - 1)) {
    got <- unique(c(ladder_pairs[[i]]$a, ladder_pairs[[i]]$b))
    want <- unique(c(real_merges[[i]]$a, real_merges[[i]]$b))
    expect_setequal(got, want)
  }
})

test_that("the per-draw mean and the plug-in-at-posterior-mean value of loss(j,k) actually differ", {
  set.seed(6)
  n_draws <- 2000
  A <- 4
  # induce correlation the way the sum-to-zero constraint does: draw raw
  # normals, subtract row means (mirrors .simulate_alpha_prior(), but no
  # fold, since the point here is just "nonlinearity + correlation makes
  # E[f(alpha)] != f(E[alpha])", not orientation)
  raw <- matrix(rnorm(n_draws * A, sd = 0.5), n_draws, A)
  raw <- raw - rowMeans(raw)
  alpha_mat <- raw + matrix(c(1, -0.5, 0.2, -0.3), n_draws, A, byrow = TRUE)
  shares <- c(0.4, 0.3, 0.2, 0.1)

  pw <- bilatr:::.pairwise_merge_loss(alpha_mat, shares)
  per_draw_mean <- colMeans(pw$loss)

  alpha_at_mean <- matrix(colMeans(alpha_mat), nrow = 1)
  plug_in <- bilatr:::.pairwise_merge_loss(alpha_at_mean, shares)$loss[1, ]

  expect_false(isTRUE(all.equal(per_draw_mean, plug_in)))
})

test_that("the merge ranking (pairwise and ladder order) is invariant to phi", {
  # phi only ever scales the DM-corrected `effective_info` summary
  # figure, never the pct_info_lost ranking used for pairwise ordering or
  # the greedy ladder's merge choices.
  set.seed(7)
  D <- 4
  Tn <- 3
  A <- 5
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  stan_data <- list(D = D, T = Tn, A = A, Y = Y, is_obs = is_obs)
  attr(stan_data, "event_classes") <- as.character(seq_len(A))

  alpha_draws <- matrix(rnorm(50 * A), 50, A)
  colnames(alpha_draws) <- paste0("alpha[", seq_len(A), "]")
  phi_low <- matrix(0.5, 50, D)
  colnames(phi_low) <- paste0("phi[", seq_len(D), "]")
  phi_high <- matrix(20, 50, D)
  colnames(phi_high) <- paste0("phi[", seq_len(D), "]")

  make_draws <- function(phi_mat) {
    posterior::as_draws_array(array(cbind(alpha_draws, phi_mat), dim = c(50, 1, A + D),
      dimnames = list(NULL, NULL, c(colnames(alpha_draws), colnames(phi_mat)))))
  }

  fit_low <- list(draws = function(variables) make_draws(phi_low))
  fit_high <- list(draws = function(variables) make_draws(phi_high))

  res_low <- diagnose_category_merges(fit_low, stan_data)
  res_high <- diagnose_category_merges(fit_high, stan_data)

  expect_equal(res_low$pairwise$index_a, res_high$pairwise$index_a)
  expect_equal(res_low$pairwise$index_b, res_high$pairwise$index_b)
  expect_equal(res_low$ladder$merged_a, res_high$ladder$merged_a)
  expect_equal(res_low$ladder$merged_b, res_high$ladder$merged_b)
  expect_false(isTRUE(all.equal(res_low$summary$effective_info_mean, res_high$summary$effective_info_mean)))
})

test_that("the empirical-share default matches shares computed by hand from a small stan_data$Y", {
  set.seed(8)
  D <- 3
  Tn <- 2
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  stan_data <- list(D = D, T = Tn, A = A, Y = Y, is_obs = is_obs)
  attr(stan_data, "event_classes") <- as.character(seq_len(A))

  alpha_draws <- matrix(rnorm(20 * A), 20, A)
  colnames(alpha_draws) <- paste0("alpha[", seq_len(A), "]")
  fake_fit <- list(draws = function(variables) {
    posterior::as_draws_array(array(alpha_draws, dim = c(20, 1, A),
      dimnames = list(NULL, NULL, colnames(alpha_draws))))
  })

  res <- diagnose_category_merges(fake_fit, stan_data, phi = 1)

  hand_shares <- apply(Y, 3, sum)
  hand_shares <- hand_shares / sum(hand_shares)
  expect_equal(res$categories$share[order(res$categories$action_index)], hand_shares, tolerance = 1e-10)
})

test_that("merge_cost() on an explicit grouping matches the ladder's cumulative loss at the matching step", {
  set.seed(9)
  D <- 4
  Tn <- 3
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  stan_data <- list(D = D, T = Tn, A = A, Y = Y, is_obs = is_obs)
  attr(stan_data, "event_classes") <- as.character(seq_len(A))

  alpha_draws <- matrix(rnorm(500 * A), 500, A)
  colnames(alpha_draws) <- paste0("alpha[", seq_len(A), "]")
  fake_fit <- list(draws = function(variables) {
    posterior::as_draws_array(array(alpha_draws, dim = c(500, 1, A),
      dimnames = list(NULL, NULL, colnames(alpha_draws))))
  })

  res <- diagnose_category_merges(fake_fit, stan_data, phi = 1)

  # price the exact grouping the ladder's step 1 represents
  step1 <- res$ladder[1, ]
  groups <- list(c(strsplit(step1$merged_a, "\\+")[[1]], strsplit(step1$merged_b, "\\+")[[1]]))
  priced <- merge_cost(res, groups)

  expect_equal(priced$pct_info_lost, step1$loss_mean, tolerance = 1e-8)
})

test_that("a prior_only = 1 short fit reproduces alpha_prior_moments() within MCMC error", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  D <- 3
  Tn <- 3
  A <- 4
  Y <- array(sample(0:6, D * Tn * A, replace = TRUE), dim = c(D, Tn, A))
  is_obs <- matrix(1L, D, Tn)
  data_list <- list(
    T = Tn, D = D, A = A, C = 1, is_obs = is_obs, Y = Y,
    compute_log_lik = 0, prior_only = 1, compute_theta_filtered = 0, n_filter_dyads = 0, filter_dyads = integer(0)
  )

  # A prior-only fit's alpha has no likelihood pulling on it at all, so
  # its only geometry is the radial-degeneracy funnel the .stan header
  # already documents (some divergences here are expected/harmless) --
  # needs enough chains/draws for adequate ESS on the exchangeable
  # non-reference classes, or this reproduces alpha_prior_moments() only
  # noisily. Confirmed 4 chains x 2000 draws gets comfortably within
  # 0.05 of both mean and sd in exploratory runs.
  mod <- .compile_stan_model("stable", opt_level = 1)
  fit <- suppressWarnings(mod$sample(
    data = data_list, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000, seed = 1, refresh = 0,
    threads_per_chain = 1, init = bilatr_init_fn(list(D = D, T = Tn, A = A), stan_model = "stable"),
    output_dir = tempdir(), show_messages = FALSE
  ))

  alpha_summ <- posterior::summarise_draws(fit$draws("alpha"), mean = mean, sd = stats::sd)
  alpha_summ <- dplyr::mutate(
    alpha_summ, action_index = as.integer(stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])"))
  )
  alpha_summ <- dplyr::arrange(alpha_summ, action_index)

  prior <- alpha_prior_moments(A)

  expect_lt(max(abs(alpha_summ$mean - prior$prior_mean)), 0.05)
  expect_lt(max(abs(alpha_summ$sd - prior$prior_sd)), 0.05)
})

test_that("print.bilatr_category_merges() runs without error on the CmdStan fixture", {
  skip_if_no_cmdstan()
  skip_on_cran()
  skip_on_ci()

  fx <- make_csv_diagnostics_fixture(stan_model = "stable")
  res <- diagnose_category_merges(fx$fit, fx$stan_data)
  expect_s3_class(res, "bilatr_category_merges")
  expect_output(print(res), "bilatr_category_merges")
  expect_output(print(res), "not a severity")

  priced <- merge_cost(res, list(res$event_classes[1:2]))
  expect_true(is.data.frame(priced))
  expect_true(all(c("pct_info_lost", "pct_info_lost_lower", "pct_info_lost_upper") %in% names(priced)))
})
