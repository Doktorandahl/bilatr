test_that("the registry contains exactly the stable model", {
  expect_setequal(names(.bilatr_stan_models), "stable")
  expect_identical(.BILATR_DEFAULT_MODEL, "stable")
  # retired variants must not be registered
  expect_false(any(
    c("phi_logn", "stable_ncproc", "phi_logn_ncproc") %in% names(.bilatr_stan_models)
  ))
})

test_that(".resolve_stan_model() resolves valid names to existing files", {
  stable_path <- .resolve_stan_model("stable")
  expect_true(file.exists(stable_path))
  expect_match(stable_path, "bilatr_dirmult_irt\\.stan$")
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
