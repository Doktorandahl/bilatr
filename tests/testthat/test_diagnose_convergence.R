make_fake_draws <- function() {
  set.seed(1)
  n_iter <- 400
  n_chain <- 4
  variables <- c(
    "lp__", "alpha[1]", "alpha[2]", "mu_theta0", "sigma_theta0",
    "phi[1]", "phi[2]", "phi[3]", "phi[4]",
    "process_noise[1]", "process_noise[2]", "process_noise[3]", "process_noise[4]",
    "theta0[1]", "theta0[2]", "theta0[3]", "theta0[4]",
    "theta[1,1]", "theta[1,2]", "theta[2,1]", "theta[2,2]",
    "theta[3,1]", "theta[3,2]", "theta[4,1]", "theta[4,2]"
  )
  arr <- array(
    stats::rnorm(n_iter * n_chain * length(variables)),
    dim = c(n_iter, n_chain, length(variables)),
    dimnames = list(NULL, NULL, variables)
  )

  # mu_theta0: chain-shifted mean -> bad Rhat, the one Tier 1 issue.
  for (ch in seq_len(n_chain)) {
    arr[, ch, "mu_theta0"] <- arr[, ch, "mu_theta0"] + ch * 5
  }

  # phi[1]: strong within-chain autocorrelation (random walk) -> low ESS,
  # even though its dyad (dyad 1) is the most data-rich (see n_dt fixture
  # below), so it should be flagged as worse-than-expected for its n_dt.
  for (ch in seq_len(n_chain)) {
    arr[, ch, "phi[1]"] <- cumsum(arr[, ch, "phi[1]"]) / 20
  }

  posterior::as_draws_array(arr)
}

make_fake_n_dt <- function() {
  # dyad 1 is the most data-rich but (per make_fake_draws()) has an
  # anomalously noisy phi; dyad 5 exists in n_dt but not in the draws.
  tibble::tibble(
    dyad_id = c(1, 2, 3, 4, 5),
    n_dt = c(500, 400, 300, 200, 999)
  )
}

# --- tier classification -----------------------------------------------

test_that(".classify_bilatr_tier assigns Tier 1 by fixed name, regardless of index shape", {
  out <- .classify_bilatr_tier(c("alpha[1]", "alpha[2]", "mu_theta0", "lp__", "mu_log_phi"))
  expect_equal(out$tier, rep(1L, 5))
  expect_true(all(is.na(out$dyad_id)))
})

test_that(".classify_bilatr_tier assigns Tier 2 to single-index, non-Tier-1 parameters", {
  out <- .classify_bilatr_tier(c("phi[3]", "process_noise[1]", "theta0[2]", "z_theta0[5]"))
  expect_equal(out$tier, rep(2L, 4))
  expect_equal(out$dyad_id, c(3L, 1L, 2L, 5L))
  expect_true(all(is.na(out$time_index)))
})

test_that(".classify_bilatr_tier assigns Tier 3 to double-index parameters", {
  out <- .classify_bilatr_tier(c("theta[3,12]", "theta_raw[2,1]"))
  expect_equal(out$tier, c(3L, 3L))
  expect_equal(out$dyad_id, c(3L, 2L))
  expect_equal(out$time_index, c(12L, 1L))
})

test_that(".classify_bilatr_tier classifies phi structurally by index shape", {
  # stable model: phi is per-dyad (Tier 2).
  per_dyad <- .classify_bilatr_tier("phi[4]")
  expect_equal(per_dyad$tier, 2L)
  expect_equal(per_dyad$dyad_id, 4L)

  # a per-dyad-period phi[d, t] (as a future variant might use) falls to
  # Tier 3 on the same base name, no special-casing needed.
  per_dyad_period <- .classify_bilatr_tier("phi[4,7]")
  expect_equal(per_dyad_period$tier, 3L)
  expect_equal(per_dyad_period$dyad_id, 4L)
  expect_equal(per_dyad_period$time_index, 7L)
})

test_that(".classify_bilatr_tier keeps bracketed non-centered globals in Tier 1, not Tier 2", {
  # alpha_raw[i] / mu_intercept_raw[i] are single-indexed by *action type*,
  # not dyad -- the structural single-index rule alone would wrongly route
  # them into Tier 2's dyad join, so they must be caught by name first.
  out <- .classify_bilatr_tier(c("alpha_raw[1]", "mu_intercept_raw[2]"))
  expect_equal(out$tier, c(1L, 1L))
  expect_true(all(is.na(out$dyad_id)))
})

test_that(".classify_bilatr_tier folds unclassified scalars into Tier 1 rather than dropping them", {
  out <- .classify_bilatr_tier("some_new_global_param")
  expect_equal(out$tier, 1L)
})

# --- n_dt normalization / join behavior ---------------------------------

test_that(".normalize_n_dt accepts a data frame with a dyad-id-like and n_dt-like column", {
  out <- .normalize_n_dt(tibble::tibble(dyad = c(1, 2), n_dt = c(10, 20)))
  expect_equal(out$dyad_id, c(1L, 2L))
  expect_equal(out$n_dt, c(10, 20))
})

test_that(".normalize_n_dt accepts a named numeric vector keyed by dyad id", {
  out <- .normalize_n_dt(c(`1` = 10, `2` = 20))
  expect_equal(out$dyad_id, c(1L, 2L))
  expect_equal(out$n_dt, c(10, 20))
})

test_that(".normalize_n_dt errors when it can't find the expected columns", {
  expect_error(.normalize_n_dt(tibble::tibble(foo = 1, bar = 2)), "dyad-id-like")
})

test_that("diagnose_convergence() joins Tier 2 against n_dt and surfaces both-sided mismatches", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())

  # dyad 5 is in n_dt but not in the draws.
  expect_equal(diag$summary$n_dyads_missing_from_draws, 1L)
  expect_true(5L %in% diag$tier2$dyad_id)
  expect_true(is.na(diag$tier2$phi_rhat[diag$tier2$dyad_id == 5L]))

  # all 4 drawn dyads are matched in n_dt, so nothing missing from n_dt.
  expect_equal(diag$summary$n_dyads_missing_from_n_dt, 0L)
})

test_that("diagnose_convergence() joins Tier 2 against n_dt with an extra draws-side dyad", {
  extra_draws <- make_fake_draws()
  n_dt <- dplyr::filter(make_fake_n_dt(), dyad_id != 4) # drop dyad 4 from n_dt

  diag <- diagnose_convergence(extra_draws, n_dt = n_dt)

  expect_equal(diag$summary$n_dyads_missing_from_n_dt, 1L)
  expect_true(4L %in% diag$tier2$dyad_id)
  expect_true(is.na(diag$tier2$n_dt[diag$tier2$dyad_id == 4L]))
})

# --- threshold flagging ---------------------------------------------------

test_that(".flag_diagnostic() flags on Rhat, ESS_bulk, or ESS_tail breaching threshold", {
  expect_true(.flag_diagnostic(1.02, 1000, 1000, 1.01, 400))
  expect_true(.flag_diagnostic(1.0, 300, 1000, 1.01, 400))
  expect_true(.flag_diagnostic(1.0, 1000, 300, 1.01, 400))
  expect_false(.flag_diagnostic(1.0, 1000, 1000, 1.01, 400))
})

test_that("diagnose_convergence() flags Tier 1 issues and never drops them regardless of count", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())

  expect_true("mu_theta0" %in% diag$tier1$variable[diag$tier1$flagged])
  expect_gte(diag$summary$n_tier1_flagged, 1L)
  # a well-mixed global parameter should not be flagged.
  expect_false(diag$tier1$flagged[diag$tier1$variable == "sigma_theta0"])
})

test_that("diagnose_convergence() threshold arguments change what gets flagged", {
  lenient <- diagnose_convergence(
    make_fake_draws(), n_dt = make_fake_n_dt(),
    rhat_threshold = 100, ess_threshold = 0
  )
  expect_equal(lenient$summary$n_tier1_flagged, 0L)
})

# --- Tier 2/3 shape and outlier flag -------------------------------------

test_that(".flag_worse_than_expected flags within-bracket outliers, not just sparse dyads", {
  n_dt <- c(10, 12, 11, 500, 480, 510)
  # last dyad in the high-n_dt bracket is an ESS outlier relative to its peers.
  ess <- c(50, 55, 48, 3000, 2900, 40)

  flags <- .flag_worse_than_expected(n_dt, ess, n_bins = 2)

  expect_false(flags[1]) # sparse dyad with low ESS: expected, not flagged
  expect_true(flags[6]) # data-rich dyad with low ESS: outlier, flagged
  expect_false(flags[4]) # data-rich dyad with high (expected) ESS
})

test_that("diagnose_convergence() Tier 2 output has one row per dyad with a worse_than_expected column", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())

  expect_true(all(c("phi_ess_bulk", "process_noise_ess_bulk", "theta0_ess_bulk", "worse_than_expected") %in% names(diag$tier2)))
  expect_equal(nrow(diag$tier2), length(unique(c(diag$tier2$dyad_id))))
})

test_that("diagnose_convergence() Tier 3 is aggregated per dyad, not per dyad-period", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())

  expect_equal(nrow(diag$tier3), 4L) # 4 dyads with theta[d,t], not 8 dyad-periods
  expect_true(all(c("min_ess_bulk", "min_ess_tail", "share_ess_below_threshold") %in% names(diag$tier3)))
})

# --- class and print method ------------------------------------------------

test_that("diagnose_convergence() returns a bilatr_diagnostics object", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())
  expect_s3_class(diag, "bilatr_diagnostics")
  expect_named(diag, c("tier1", "tier2", "tier3", "summary"))
})

# --- tiers argument ------------------------------------------------------

test_that(".validate_tiers accepts subsets of 1:3 and rejects everything else", {
  expect_equal(.validate_tiers(1), 1L)
  expect_equal(.validate_tiers(c(3, 1, 1)), c(1L, 3L))
  expect_equal(.validate_tiers(1:3), 1:3)
  expect_error(.validate_tiers(4), "subset of 1:3")
  expect_error(.validate_tiers(0), "subset of 1:3")
  expect_error(.validate_tiers(character(0)), "subset of 1:3")
})

test_that("diagnose_convergence(tiers = 1) computes only Tier 1 and needs no n_dt", {
  diag <- diagnose_convergence(make_fake_draws(), tiers = 1)

  expect_false(is.null(diag$tier1))
  expect_null(diag$tier2)
  expect_null(diag$tier3)
  expect_equal(diag$summary$tiers_computed, 1L)
  expect_true("mu_theta0" %in% diag$tier1$variable[diag$tier1$flagged])
})

test_that("diagnose_convergence() errors if n_dt is missing but Tier 2/3 were requested", {
  expect_error(
    diagnose_convergence(make_fake_draws(), tiers = 2),
    "n_dt.*required"
  )
  expect_error(
    diagnose_convergence(make_fake_draws(), tiers = c(1, 3)),
    "n_dt.*required"
  )
})

test_that("diagnose_convergence(tiers = 2) computes only Tier 2", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt(), tiers = 2)

  expect_null(diag$tier1)
  expect_false(is.null(diag$tier2))
  expect_null(diag$tier3)
  expect_equal(diag$summary$tiers_computed, 2L)
  expect_true(all(c("phi_ess_bulk", "worse_than_expected") %in% names(diag$tier2)))
})

test_that("diagnose_convergence(tiers = c(2, 3)) matches the Tier 2/3 output of tiers = 1:3", {
  full <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())
  partial <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt(), tiers = c(2, 3))

  expect_null(partial$tier1)
  expect_equal(partial$tier2, full$tier2)
  expect_equal(partial$tier3, full$tier3)
})

test_that("print.bilatr_diagnostics() only prints sections for computed tiers", {
  diag <- diagnose_convergence(make_fake_draws(), tiers = 1)
  out <- capture.output(print(diag))

  expect_true(any(grepl("Tier 1", out)))
  expect_false(any(grepl("Tier 2", out)))
  expect_false(any(grepl("Tier 3", out)))
})

test_that("print.bilatr_diagnostics() always prints flagged Tier 1 rows in full", {
  diag <- diagnose_convergence(make_fake_draws(), n_dt = make_fake_n_dt())
  out <- capture.output(print(diag))
  expect_true(any(grepl("mu_theta0", out)))
  expect_true(any(grepl("Tier 1", out)))
  expect_true(any(grepl("Tier 2", out)))
  expect_true(any(grepl("Tier 3", out)))
})
