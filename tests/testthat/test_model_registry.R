test_that("the registry contains stable plus the ou experimental variant", {
  expect_setequal(
    names(.bilatr_stan_models),
    c("stable", "ou")
  )
  expect_identical(.BILATR_DEFAULT_MODEL, "stable")
  # retired variants (including the pre-0.4.0 stable/ou) must not be
  # registered
  expect_false(any(
    c("phi_logn", "stable_ncproc", "phi_logn_ncproc", "alphanorm", "alphanorm_ou") %in%
      names(.bilatr_stan_models)
  ))
})

test_that("only stable is status 'stable'; ou is 'experimental'", {
  statuses <- vapply(.bilatr_stan_models, function(x) x$status, character(1))
  expect_identical(statuses[["stable"]], "stable")
  expect_true(all(
    statuses[setdiff(names(statuses), "stable")] == "experimental"
  ))
})

test_that(".resolve_stan_model() resolves valid names to existing files", {
  stable_path <- .resolve_stan_model("stable")
  expect_true(file.exists(stable_path))
  expect_match(stable_path, "bilatr_alphanorm\\.stan$")

  expect_match(.resolve_stan_model("ou"), "bilatr_alphanorm_ou\\.stan$")
})

test_that(".resolve_stan_model() errors informatively on an unknown name", {
  expect_error(.resolve_stan_model("not_a_real_model"), "Unknown stan_model")
  expect_error(.resolve_stan_model("not_a_real_model"), "stable")
  # phi_logn is retired -> also an unknown name now (never registered,
  # and not a recognized alias either)
  expect_error(.resolve_stan_model("phi_logn"), "Unknown stan_model")
})

test_that(".canonical_stan_model()/.resolve_stan_model() accept the pre-0.4.0 alphanorm/alphanorm_ou aliases (B1)", {
  # the alias message is once-per-session (per name), not once-per-call
  # (see .reset_bilatr_alias_messaged()'s docs) -- reset before each
  # call this test expects to message again
  .reset_bilatr_alias_messaged()

  # alphanorm/alphanorm_ou were never registered names themselves --
  # .canonical_stan_model() maps them to stable/ou with a message, it
  # doesn't add them to .bilatr_stan_models
  expect_message(canonical <- .canonical_stan_model("alphanorm"), "pre-0.4.0 name of 'stable'")
  expect_identical(canonical, "stable")
  expect_message(canonical_ou <- .canonical_stan_model("alphanorm_ou"), "pre-0.4.0 name of 'ou'")
  expect_identical(canonical_ou, "ou")

  .reset_bilatr_alias_messaged()
  expect_message(alphanorm_path <- .resolve_stan_model("alphanorm"), "pre-0.4.0 name")
  expect_identical(alphanorm_path, .resolve_stan_model("stable"))
  expect_message(alphanorm_ou_path <- .resolve_stan_model("alphanorm_ou"), "pre-0.4.0 name")
  expect_identical(alphanorm_ou_path, .resolve_stan_model("ou"))
})

test_that(".canonical_stan_model() messages only once per session per alias name", {
  .reset_bilatr_alias_messaged()
  expect_message(.canonical_stan_model("alphanorm"), "pre-0.4.0 name of 'stable'")
  expect_no_message(.canonical_stan_model("alphanorm"))
  # a different alias still messages independently
  expect_message(.canonical_stan_model("alphanorm_ou"), "pre-0.4.0 name of 'ou'")
})

test_that("every registered model's file actually exists under inst/stan/", {
  for (name in names(.bilatr_stan_models)) {
    path <- .resolve_stan_model(name)
    expect_true(file.exists(path), info = paste("missing Stan file for model:", name))
  }
})
