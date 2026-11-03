# Builds tests/testthat/fixtures/assemble_stan_data_0.11.0_reference.rds:
# the frozen output of the pre-0.12.0 data path, assemble_stan_data(), for
# every case of tests/testthat/test_prepare_stan_data_equivalence.R.
#
# MUST BE RUN AGAINST A bilatr TREE AT OR BEFORE 0.12.0b: it calls
# assemble_stan_data(), which 0.12.0c removed. For example
#
#   git worktree add ../bilatr-0.12.0b 393f710
#   cd ../bilatr-0.12.0b
#   cp <this repo>/tests/testthat/helper_equivalence_cases.R tests/testthat/
#   Rscript <this repo>/data-raw/make_equivalence_reference.R <this repo>
#
# In-tree assemble_stan_data() at 0.12.0b was shown identical to 0.10.2 in
# 0.12.0a (data-raw/equivalence/equivalence_0.12.0a.md), so the stored objects are the 0.11.0
# output. The case inputs come from tests/testthat/helper_equivalence_cases.R,
# which the test sources too; change a case there only together with
# regenerating this file.
#
# Output: a named list, one element per case (names as in
# equivalence_cases()), each the full stan_data list with its attributes.
#
# Usage: Rscript data-raw/make_equivalence_reference.R [output repo root]
# (default: the current directory, which must be the <= 0.12.0b tree).

args <- commandArgs(trailingOnly = TRUE)
out_root <- if (length(args) >= 1) args[[1]] else "."

devtools::load_all(".", quiet = TRUE, helpers = FALSE)
if (!exists("assemble_stan_data", mode = "function")) {
  stop("assemble_stan_data() is not in this tree; run against bilatr <= 0.12.0b.")
}
# Sourced into a child of the package namespace, as testthat does with
# helpers, so they see the package's imports (%>%, %||%) and internals.
helpers <- new.env(parent = asNamespace("bilatr"))
sys.source("tests/testthat/helper_fixtures.R", envir = helpers)
sys.source(file.path(out_root, "tests/testthat/helper_equivalence_cases.R"), envir = helpers)
equivalence_cases <- helpers$equivalence_cases

old_path_stan_data <- function(case) {
  suppressMessages(do.call(assemble_stan_data, c(
    list(
      case$events,
      years = case$years, resolution = case$resolution, grouping_var = case$grouping_var,
      directed = case$directed, actor1 = case$actor1, actor2 = case$actor2, date = case$date
    ),
    case$prep
  )))
}

cases <- equivalence_cases()
reference <- lapply(cases, old_path_stan_data)

out <- file.path(out_root, "tests/testthat/fixtures/assemble_stan_data_0.11.0_reference.rds")
saveRDS(reference, out, compress = "xz")
message(sprintf(
  "wrote %d cases to %s (%.1f KB), bilatr %s",
  length(reference), out, file.size(out) / 1024, as.character(utils::packageVersion("bilatr"))
))
