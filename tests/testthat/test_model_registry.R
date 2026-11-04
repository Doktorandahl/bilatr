test_that("the registry contains exactly stable/ou/stable_gamma", {
  expect_setequal(
    names(.bilatr_stan_models),
    c("stable", "ou", "stable_gamma")
  )
  expect_identical(.BILATR_DEFAULT_MODEL, "stable")
  # names of unsupported development models must not be registered
  expect_false(any(
    c(
      "phi_logn", "stable_ncproc", "phi_logn_ncproc",
      "alphanorm", "alphanorm_ou", "stable_soft_anchor", "ou_soft_anchor"
    ) %in%
      names(.bilatr_stan_models)
  ))
})

test_that("stable is 'stable'; ou/stable_gamma are 'experimental'", {
  statuses <- vapply(.bilatr_stan_models, function(x) x$status, character(1))
  expect_identical(statuses[["stable"]], "stable")
  expect_identical(statuses[["ou"]], "experimental")
  expect_identical(statuses[["stable_gamma"]], "experimental")
})

test_that(".resolve_stan_model() resolves valid names to existing files", {
  stable_path <- .resolve_stan_model("stable")
  expect_true(file.exists(stable_path))
  expect_match(stable_path, "bilatr_alphanorm\\.stan$")

  expect_match(.resolve_stan_model("ou"), "bilatr_alphanorm_ou\\.stan$")
  expect_match(.resolve_stan_model("stable_gamma"), "bilatr_alphanorm_gamma\\.stan$")
})

test_that(".bilatr_model_has_gamma() is TRUE only for stable_gamma", {
  expect_true(.bilatr_model_has_gamma("stable_gamma"))
  expect_false(.bilatr_model_has_gamma("stable"))
  expect_false(.bilatr_model_has_gamma("ou"))
})

test_that(".resolve_stan_model() errors informatively on an unknown name", {
  expect_error(.resolve_stan_model("not_a_real_model"), "Unknown stan_model")
  expect_error(.resolve_stan_model("not_a_real_model"), "stable")
  # an unsupported development model name is also unknown
  expect_error(.resolve_stan_model("phi_logn"), "Unknown stan_model")
})

test_that("names that are no longer registered get the generic 'Unknown stan_model' error", {
  for (name in c("alphanorm", "alphanorm_ou", "stable_soft_anchor", "ou_soft_anchor")) {
    expect_error(.canonical_stan_model(name), "Unknown stan_model")
    expect_error(.resolve_stan_model(name), "Unknown stan_model")
  }
})

test_that(".canonical_stan_model() rejects a non-string", {
  expect_error(.canonical_stan_model(c("stable", "ou")), "single string")
  expect_error(.canonical_stan_model(NULL), "single string")
  expect_error(.canonical_stan_model(NA_character_), "single string")
})

test_that("every registered model's file actually exists under inst/stan/", {
  for (name in names(.bilatr_stan_models)) {
    path <- .resolve_stan_model(name)
    expect_true(file.exists(path), info = paste("missing Stan file for model:", name))
  }
})

# The following two blocks check that fit_dyad_ts()/fit_panel() validate
# `stan_model` before doing anything else.

test_that("fit_dyad_ts()/fit_panel() reject an unknown stan_model before touching Stan", {
  fake_data_d1 <- list(D = 1, T = 3, A = 4)
  fake_data_d3 <- list(D = 3, T = 3, A = 4)
  expect_error(
    fit_dyad_ts(fake_data_d1, stan_model = "not_a_real_model"),
    "Unknown stan_model"
  )
  expect_error(
    fit_panel(fake_data_d3, stan_model = "not_a_real_model"),
    "Unknown stan_model"
  )
})

test_that("fit_dyad_ts()/fit_panel() still validate D == 1 / D >= 2 for a registered stan_model", {
  fake_data_d3 <- list(D = 3, T = 3, A = 4)
  fake_data_d1 <- list(D = 1, T = 3, A = 4)
  expect_error(
    fit_dyad_ts(fake_data_d3, stan_model = "stable"),
    "expects a single-dyad"
  )
  expect_error(
    fit_panel(fake_data_d1, stan_model = "stable"),
    "expects multiple dyads"
  )
})

test_that(".bilatr_warn_if_experimental() messages once per stan_model name per session, and never for 'stable'", {
  .reset_bilatr_experimental_messaged()

  expect_message(.bilatr_warn_if_experimental("ou"), "experimental.*\\?compile_bilatr_model")
  expect_no_message(.bilatr_warn_if_experimental("ou"))
  expect_message(.bilatr_warn_if_experimental("stable_gamma"), "experimental")
  expect_no_message(.bilatr_warn_if_experimental("stable_gamma"))
  expect_no_message(.bilatr_warn_if_experimental("stable"))

  .reset_bilatr_experimental_messaged()
  expect_message(.bilatr_warn_if_experimental("ou"), "experimental")
})
