## Submission

This is a new submission.

## Test environments

<!-- Filled in after the final checks. -->

* local: macOS ..., R ...
* win-builder: R-devel, R-release
* mac-builder: R-release
* R-hub: ...

## R CMD check results

0 errors | 0 warnings | N notes

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
* The vignettes are precomputed: their model output is baked into the
  `.Rmd` files, and building them needs no 'cmdstanr'.

## Internet use

`bilatr_theta()` downloads the package's estimated trajectories from a
versioned release of the public repository
<https://github.com/doktorandahl/bilatr-data> (archived on Zenodo with a
DOI) into `tools::R_user_dir("bilatr", "data")`. It verifies each file's
checksum before caching it, and when offline, or when a file is missing
or fails verification, it stops with an informative error naming the URL
and the cache path; it never returns partial data.

Examples and tests do not need the internet: examples use the shipped
`theta_sample` data, network examples run only under
`@examplesIf interactive()`, and the tests use local `file://` fixtures.
The GDELT functions (`download_gdelt()`, `read_gdelt()`) follow the same
rules.

## Tests

Tests that need 'CmdStan' use `skip_on_cran()`. The CRAN-visible suite
runs in about 40 seconds.

## Files written

Nothing is written to the user's home directory or the package library.
Compiled models go to `tools::R_user_dir("bilatr", "cache")` and are
removed by `bilatr_clear_model_cache()`; downloaded data go to
`tools::R_user_dir("bilatr", "data")` and are removed by
`bilatr_clear_data_cache()`. Examples and tests that run without
'CmdStan' (all of those run on CRAN) write only to `tempdir()`.
