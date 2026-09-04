test_that("the registry contains stable plus the three experimental variants", {
  expect_setequal(
    names(.bilatr_stan_models),
    c("stable", "alphanorm", "ou", "alphanorm_ou")
  )
  expect_identical(.BILATR_DEFAULT_MODEL, "stable")
  # retired variants must not be registered
  expect_false(any(
    c("phi_logn", "stable_ncproc", "phi_logn_ncproc") %in% names(.bilatr_stan_models)
  ))
})

test_that("only stable is status 'stable'; the rest are 'experimental'", {
  statuses <- vapply(.bilatr_stan_models, function(x) x$status, character(1))
  expect_identical(statuses[["stable"]], "stable")
  expect_true(all(
    statuses[setdiff(names(statuses), "stable")] == "experimental"
  ))
})

test_that(".resolve_stan_model() resolves valid names to existing files", {
  stable_path <- .resolve_stan_model("stable")
  expect_true(file.exists(stable_path))
  expect_match(stable_path, "bilatr_dirmult_irt\\.stan$")

  expect_match(.resolve_stan_model("alphanorm"), "bilatr_alphanorm\\.stan$")
  expect_match(.resolve_stan_model("ou"), "bilatr_ou\\.stan$")
  expect_match(.resolve_stan_model("alphanorm_ou"), "bilatr_alphanorm_ou\\.stan$")
})

test_that(".resolve_stan_model() errors informatively on an unknown name", {
  expect_error(.resolve_stan_model("not_a_real_model"), "Unknown stan_model")
  expect_error(.resolve_stan_model("not_a_real_model"), "stable")
  # phi_logn is retired -> also an unknown name now
  expect_error(.resolve_stan_model("phi_logn"), "Unknown stan_model")
})

test_that("every registered model's file actually exists under inst/stan/", {
  for (name in names(.bilatr_stan_models)) {
    path <- .resolve_stan_model(name)
    expect_true(file.exists(path), info = paste("missing Stan file for model:", name))
  }
})
