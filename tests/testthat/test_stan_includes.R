test_that("every registered model's GENERATED partial_log_lik block matches the canonical source", {
  fragment_path <- system.file(
    "stan", "include", "partial_log_lik.stanfunctions",
    package = "bilatr"
  )
  fragment_lines <- readLines(fragment_path)

  for (name in names(.bilatr_stan_models)) {
    stan_path <- .resolve_stan_model(name)
    current_lines <- readLines(stan_path)
    expected_lines <- .splice_partial_log_lik_block(current_lines, fragment_lines)
    expect_identical(
      current_lines, expected_lines,
      info = paste0(
        "inst/stan/", basename(stan_path), " has drifted from ",
        "inst/stan/include/partial_log_lik.stanfunctions -- rerun ",
        "`Rscript data-raw/sync_stan_functions.R`"
      )
    )
  }
})

test_that(".splice_partial_log_lik_block() errors on a malformed or missing marker pair", {
  fragment_lines <- c("real f(real x) { return x; }")

  expect_error(
    .splice_partial_log_lik_block(c("functions {", "}"), fragment_lines),
    "well-formed GENERATED partial_log_lik block"
  )

  duplicated <- c(
    "functions {",
    "  // <<< BEGIN GENERATED partial_log_lik (source: inst/stan/include/partial_log_lik.stanfunctions) >>>",
    "  // <<< END GENERATED partial_log_lik >>>",
    "  // <<< BEGIN GENERATED partial_log_lik (source: inst/stan/include/partial_log_lik.stanfunctions) >>>",
    "  // <<< END GENERATED partial_log_lik >>>",
    "}"
  )
  expect_error(
    .splice_partial_log_lik_block(duplicated, fragment_lines),
    "well-formed GENERATED partial_log_lik block"
  )
})
