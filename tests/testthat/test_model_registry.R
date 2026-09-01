test_that("the registry contains exactly the stable and phi_logn models", {
  expect_setequal(names(.bilatr_stan_models), c("stable", "phi_logn"))
  expect_identical(.BILATR_DEFAULT_MODEL, "stable")
  # the transitional "_ncproc" entries were retired in 0.3.0
  expect_false(any(c("stable_ncproc", "phi_logn_ncproc") %in% names(.bilatr_stan_models)))
})

test_that(".resolve_stan_model() resolves valid names to existing files", {
  stable_path <- .resolve_stan_model("stable")
  expect_true(file.exists(stable_path))
  expect_match(stable_path, "bilatr_dirmult_irt\\.stan$")

  phi_logn_path <- .resolve_stan_model("phi_logn")
  expect_true(file.exists(phi_logn_path))
  expect_match(phi_logn_path, "bilatr_phi_logn\\.stan$")
})

test_that(".resolve_stan_model() errors informatively on an unknown name", {
  expect_error(.resolve_stan_model("not_a_real_model"), "Unknown stan_model")
  expect_error(.resolve_stan_model("not_a_real_model"), "stable")
  expect_error(.resolve_stan_model("not_a_real_model"), "phi_logn")
})

test_that("every registered model's file actually exists under inst/stan/", {
  for (name in names(.bilatr_stan_models)) {
    path <- .resolve_stan_model(name)
    expect_true(file.exists(path), info = paste("missing Stan file for model:", name))
  }
})
