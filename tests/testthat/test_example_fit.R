# Every post-fit function on the example fits shipped in inst/extdata. No
# CmdStan and no skip_on_cran(): these are the post-fit tests that run on
# CRAN. The assertions are light (shape, names, the recorded model, no
# error); the functions' behavior is tested in their own files.

example_stan_data <- function(model) {
  readRDS(file.path(dirname(bilatr_example_fit(model)[1]), "example_stan_data.rds"))
}

test_that("bilatr_example_fit() returns the two sorted CSVs of each fit", {
  for (model in c("stable", "stable_gamma")) {
    csv <- bilatr_example_fit(model)
    expect_identical(basename(csv), c("example-1.csv", "example-2.csv"))
    expect_true(all(file.exists(csv)))
    expect_identical(csv, sort(csv))
    expect_true(file.exists(file.path(dirname(csv[1]), "example_stan_data.rds")))
  }
  expect_identical(bilatr_example_fit(), bilatr_example_fit("stable"))
  expect_error(bilatr_example_fit("ou"), "stan_model")
})

test_that("the example fits record their model, in the sidecar and the CSV header", {
  for (model in c("stable", "stable_gamma")) {
    csv <- bilatr_example_fit(model)
    sd <- example_stan_data(model)
    expect_identical(attr(sd, "stan_model"), model)
    expect_identical(.stan_model_from_csv_header(csv[1]), model)
    expect_identical(
      suppressMessages(.resolve_stan_model_for_fit(csv, sd, NULL, "test")),
      model
    )
    expect_identical(c(sd$D, sd$T, sd$A), c(6L, 6L, 5L))
    expect_identical(attr(sd, "grouping_var"), "PentaClass")
    expect_identical(attr(sd, "event_classes")[1], "1")
  }
})

test_that("the example CSV headers carry no absolute paths", {
  for (csv in c(bilatr_example_fit("stable"), bilatr_example_fit("stable_gamma"))) {
    header <- grep("^#", readLines(csv), value = TRUE)
    expect_false(any(grepl("= /", header, fixed = TRUE)))
    expect_false(any(grepl("=/", header, fixed = TRUE)))
    parsed <- .scan_stan_csv_header_and_skip(csv)
    expect_identical(parsed$num_samples, 60L)
  }
})

test_that("read_bilatr_draws() reads the example fits", {
  draws <- read_bilatr_draws(bilatr_example_fit(), c("lp__", "alpha", "divergent__"))
  expect_s3_class(draws, "draws_array")
  expect_identical(dim(draws), c(60L, 2L, 7L))
  expect_identical(posterior::variables(draws)[1:2], c("lp__", "alpha[1]"))
})

test_that("extract_*() run on the example fits", {
  csv <- bilatr_example_fit()
  theta <- suppressMessages(extract_theta(csv, max_memory_mb = 4096))
  expect_equal(nrow(theta), 36L)
  expect_true(all(c("dyad_id", "time_index", "mean", "5%", "50%", "95%", "dyad", "actor1", "actor2", "year") %in% names(theta)))
  expect_setequal(unique(theta$actor1), c("CHN", "RUS", "USA"))

  alpha <- extract_alpha(csv, event_classes = attr(example_stan_data("stable"), "event_classes"))
  expect_equal(nrow(alpha), 5L)
  expect_identical(alpha$event_class[1], "1")
  # alpha is folded so that alpha[1] >= 0
  expect_gt(alpha$mean[1], 0)

  mu <- extract_mu_intercept(csv)
  expect_equal(nrow(mu), 5L)
  expect_true(all(c("variable", "mean", "action_index") %in% names(mu)))

  gamma <- suppressMessages(extract_gamma(bilatr_example_fit("stable_gamma")))
  expect_s3_class(gamma, "bilatr_gamma")
  expect_equal(nrow(gamma), 3L * 5L)
  expect_setequal(unique(gamma$country_code), c("CHN", "RUS", "USA"))
  norms <- gamma_norm(gamma)
  expect_equal(nrow(norms), 3L)
  expect_false(is.unsorted(rev(norms$gamma_norm)))
})

test_that("diagnose_convergence() and diagnose_and_extract_bilatr() run on both fits", {
  for (model in c("stable", "stable_gamma")) {
    csv <- bilatr_example_fit(model)
    diag <- suppressMessages(diagnose_convergence(csv, max_memory_mb = 4096))
    expect_s3_class(diag, "bilatr_diagnostics")
    expect_equal(nrow(diag$tier2), 6L)
    expect_equal(nrow(diag$tier3), 6L)
    if (model == "stable_gamma") {
      expect_equal(nrow(diag$gamma), 15L)
      expect_false(any(grepl("^gamma", diag$tier1$variable)))
    } else {
      expect_null(diag$gamma)
      expect_setequal(diag$tier3_other$base_name %||% diag$tier3_other[[1]],
        c("log_lik", "theta_filtered", "theta_filtered_sd"))
    }
    expect_output(print(diag), "<bilatr_diagnostics>")
  }

  out <- suppressMessages(diagnose_and_extract_bilatr(bilatr_example_fit(), max_memory_mb = 4096))
  expect_named(out, c("diagnostics", "theta", "theta_filtered", "theta_filtered_sd", "alpha", "mu_intercept"))
  expect_equal(nrow(out$theta), 36L)
  expect_equal(nrow(out$theta_filtered), 36L)
})

test_that("icc_curves() and icc_crossings() run on both fits", {
  curves <- suppressMessages(icc_curves(bilatr_example_fit(), theta_range = c(-2, 2), n_grid = 11))
  expect_s3_class(curves, "bilatr_icc_curves")
  expect_equal(nrow(curves), 5L * 11L)
  expect_identical(attr(curves, "settings")$stan_model, "stable")
  expect_output(print(curves), "<bilatr_icc_curves> type = probability")
  # the probabilities at each theta sum to (roughly) one
  sums <- tapply(curves$median, curves$theta, sum)
  expect_true(all(abs(sums - 1) < 0.2))

  crossings <- icc_crossings(curves)
  expect_equal(nrow(crossings), choose(5, 2))
  expect_true(all(c("index_a", "index_b", "crossing_median", "in_range") %in% names(crossings)))

  gamma_curves <- suppressMessages(icc_curves(
    bilatr_example_fit("stable_gamma"), theta_range = c(-2, 2), n_grid = 11, country = "USA"
  ))
  expect_identical(attr(gamma_curves, "settings")$stan_model, "stable_gamma")
  expect_identical(attr(gamma_curves, "settings")$country, "USA")

  info <- suppressMessages(icc_curves(bilatr_example_fit(), theta_range = c(-2, 2), n_grid = 11, type = "information"))
  expect_equal(nrow(attr(info, "total_information")), 11L)
})

test_that("check_compositional_residuals() and write_residual_check() run on both fits", {
  for (model in c("stable", "stable_gamma")) {
    rc <- suppressMessages(check_compositional_residuals(bilatr_example_fit(model), n_dyads = 6, n_draws = 20))
    expect_s3_class(rc, "bilatr_residual_check")
    expect_equal(nrow(rc$dyads), 6L)
    expect_equal(nrow(rc$categories), 5L)
    expect_identical(rc$settings$stan_model, model)
    expect_output(print(rc), "<bilatr_residual_check>")
  }

  skip_if_not_installed("ggplot2")
  dir <- withr::local_tempdir()
  write_residual_check(rc, dir)
  expect_setequal(list.files(dir), c(
    "residual_check_dyads.csv", "residual_check_categories.csv",
    "residual_check_eps_sensitivity.csv", "residual_check_global.csv",
    "residual_check_perp_vs_n.png", "residual_check_ppp_hist.png",
    "residual_check_category_contrib.png"
  ))
})

test_that("diagnose_category_merges() and merge_cost() run on both fits", {
  for (model in c("stable", "stable_gamma")) {
    merges <- suppressMessages(diagnose_category_merges(bilatr_example_fit(model)))
    expect_s3_class(merges, "bilatr_category_merges")
    expect_equal(nrow(merges$categories), 5L)
    expect_equal(nrow(merges$pairwise), choose(5, 2))
    expect_equal(nrow(merges$ladder), 4L)
    expect_output(print(merges), "<bilatr_category_merges>")
  }
  cost <- merge_cost(merges, groups = list(c("1", "2"), c("3", "4")))
  expect_named(cost, c("pct_info_lost", "pct_info_lost_lower", "pct_info_lost_upper"))
  expect_true(cost$pct_info_lost >= 0 && cost$pct_info_lost <= 1)
})

test_that("event_class_labels() labels the example fit's PentaClass classes", {
  labels <- event_class_labels(example_stan_data("stable"))
  expect_length(labels, 5L)
  expect_identical(unname(labels[1]), "Verbal cooperation")
})

test_that("the plot methods run on the example fits", {
  skip_if_not_installed("ggplot2")
  curves <- suppressMessages(icc_curves(bilatr_example_fit(), theta_range = c(-2, 2), n_grid = 11))
  gamma <- suppressMessages(extract_gamma(bilatr_example_fit("stable_gamma")))
  rc <- suppressMessages(check_compositional_residuals(bilatr_example_fit(), n_dyads = 6, n_draws = 20))

  expect_s3_class(ggplot2::autoplot(curves), "ggplot")
  expect_s3_class(ggplot2::autoplot(gamma), "ggplot")
  expect_length(ggplot2::autoplot(rc), 3L)
  expect_s3_class(ggplot2::autoplot(rc, which = "category"), "ggplot")

  withr::local_pdf(NULL)
  expect_identical(plot(curves), curves)
  expect_identical(plot(gamma), gamma)
  expect_identical(plot(rc, which = "ppp_hist"), rc)
})
