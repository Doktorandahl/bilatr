# Re-knits the precomputed vignettes.
#
# The vignettes that fit a model need CmdStan, which CRAN's machines do
# not have, so their sources live in vignettes/<name>.Rmd.orig and are
# knitted here, on a machine with CmdStan, to vignettes/<name>.Rmd with
# every output baked in. Figures go to vignettes/figures/<name>-*.png.
# The .Rmd.orig files and this script are .Rbuildignore'd; only the
# knitted .Rmd files and the figures ship.
#
# Run from the package root, with cmdstanr and CmdStan installed:
#
#   Rscript vignettes/precompile.R               # all four
#   Rscript vignettes/precompile.R tuning        # one or more by name
#
# The development version in this directory is loaded with
# pkgload::load_all(), so the knitted output matches the source tree.
# tuning.Rmd.orig runs a benchmark sweep of several minutes and reports
# this machine's CPU; re-knit it on the machine you want to describe.

vignettes <- c("dyad_time_series", "panel_model", "diagnostics", "tuning")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  unknown <- setdiff(args, vignettes)
  if (length(unknown) > 0) stop("Unknown vignette(s): ", paste(unknown, collapse = ", "))
  vignettes <- args
}

if (!file.exists("DESCRIPTION") || !dir.exists("vignettes")) {
  stop("Run this script from the package root.")
}
pkgload::load_all(".", quiet = TRUE)
if (!bilatr::bilatr_cmdstan_available()) {
  stop("The precomputed vignettes need cmdstanr and CmdStan.")
}
# figures are drawn with ragg (smaller PNGs), and the plots need ggplot2
for (pkg in c("knitr", "ragg", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install ", pkg, " to re-knit the vignettes.")
}

old <- setwd("vignettes")
on.exit(setwd(old), add = TRUE)
for (name in vignettes) {
  message("== knitting ", name)
  unlink(Sys.glob(file.path("figures", paste0(name, "-*.png"))))
  knitr::knit(paste0(name, ".Rmd.orig"), paste0(name, ".Rmd"), envir = new.env(), quiet = TRUE)
}
