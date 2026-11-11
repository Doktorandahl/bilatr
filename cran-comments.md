## Submission

This is a new submission.

## Resubmission

The pre-test on r-devel-linux-x86_64-debian-gcc failed two tests whose memory budget depended on the test process’s own resident memory on Linux. These tests now fix that baseline; the package code is unchanged.

## Test environments

* local: macOS 26.5 (aarch64), R 4.5.2, with and without 'cmdstanr' and
  'CmdStan'
* win-builder: R-devel, R-release (RESULTS PENDING)
* mac-builder: R-release (RESULT PENDING)
* R-hub: (PENDING, or remove)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  - New submission.
  - Suggests or Enhances not in mainstream repositories: cmdstanr.
    Availability using Additional_repositories specification: cmdstanr yes
    https://stan-dev.r-universe.dev

## 'cmdstanr' in Suggests

'cmdstanr' is in Suggests and is available from
<https://stan-dev.r-universe.dev> (declared in `Additional_repositories`).
The model-fitting functions need it and an external 'CmdStan'
installation (declared in `SystemRequirements`).

* Every function that needs 'cmdstanr' or 'CmdStan' checks for them first
  and stops with an informative error explaining how to install them.
  `bilatr_cmdstan_available()` reports whether they are usable.
* Examples that fit or compile models run under
  `@examplesIf bilatr_cmdstan_available()` and inside `\donttest{}`.
* The post-fit functions (extraction, diagnostics, model checks) are
  exemplified and tested on a small precomputed fit shipped in
  `inst/extdata` (`bilatr_example_fit()`), so they run without
  'cmdstanr'.
* The vignettes are precomputed and need no 'cmdstanr' to build.

## Internet use

Two functions use the internet, only when called:
`download_bilatr_theta()` downloads the package's estimated trajectories
from a versioned data release (GitHub, archived on Zenodo with a DOI), and
`download_gdelt()` downloads GDELT event files (and `bilatr_theta_versions()`
lists the data releases). Downloads are
checksum-verified and cached in `tools::R_user_dir("bilatr", "data")` or a
user-chosen directory. When offline, or when a file is missing or fails
verification, they stop with an informative error naming the URL; they
never return partial data.

No example or test run on CRAN needs the internet: network examples run
only under `@examplesIf interactive()`, and the tests use local `file://`
fixtures. Two tests that download the published data release use
`skip_on_cran()` and `skip_if_offline()`.

## Files written

Nothing is written to the user's home directory or the package library.
Compiled models go to `tools::R_user_dir("bilatr", "cache")` (or
`options(bilatr.model_dir = )`) and are removed by
`bilatr_clear_model_cache()`; downloaded data go to
`tools::R_user_dir("bilatr", "data")` and are removed by
`bilatr_clear_data_cache()`. Examples and tests run on CRAN write only to
`tempdir()`.
