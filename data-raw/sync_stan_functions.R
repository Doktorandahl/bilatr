# Splices inst/stan/include/partial_log_lik.stanfunctions (the single
# canonical source for the shared reduce_sum likelihood) verbatim into the
# marker-delimited GENERATED block in every registered model's .stan file
# under inst/stan/. Rerun this whenever partial_log_lik.stanfunctions
# changes; tests/testthat/test_stan_includes.R fails if a model file's
# block has drifted from the canonical source, i.e. if this wasn't rerun.
#
# See R/stan_includes.R for why this exists (a cmdstanr bug makes Stan's
# own #include unusable for this project's working directory) and for the
# splice/render helpers used here.

devtools::load_all(".", quiet = TRUE)

fragment_path <- "inst/stan/include/partial_log_lik.stanfunctions"
fragment_lines <- readLines(fragment_path)

target_files <- c(
  "inst/stan/bilatr_dirmult_irt.stan",
  "inst/stan/bilatr_alphanorm.stan",
  "inst/stan/bilatr_ou.stan",
  "inst/stan/bilatr_alphanorm_ou.stan"
)

for (f in target_files) {
  if (!file.exists(f)) {
    message("Skipping ", f, " (does not exist yet)")
    next
  }
  stan_lines <- readLines(f)
  updated <- .splice_partial_log_lik_block(stan_lines, fragment_lines)
  writeLines(updated, f)
  message("Synced: ", f)
}
