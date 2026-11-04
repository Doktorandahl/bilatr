# theta_trajectories.R
#
# NOTE: since bilatr 1.0.0 the theta tables are no longer in data/. They are
# published as a data release (data-raw/make_data_release.R, fetched with
# bilatr_theta()); data/theta_sample.rda is built by
# data-raw/make_theta_sample.R. This script is kept as the record of how
# the v2026.11 tables were first saved.
#
# Built data/theta_directed.rda and data/theta_undirected.rda (documented
# in R/data.R, ?theta_trajectories) from the processed theta series written
# by runscripts/process_theta_16cat.R (same code as
# runscripts/bilatr_replication/04_process_theta.R), copied to
# runscripts/thetatests/. The objects are saved unchanged; this script only
# re-saves them as maximally xz-compressed .rda files.
#
# Run from the package root:
#   Rscript data-raw/theta_trajectories.R

src <- "runscripts/thetatests"

for (nm in c("theta_directed", "theta_undirected")) {
  x <- readRDS(file.path(src, paste0(nm, ".rds")))
  stopifnot(
    !anyDuplicated(x[c("dyad", "year")]),
    !anyNA(x)
  )
  assign(nm, x)
  save(list = nm, file = file.path("data", paste0(nm, ".rda")),
       compress = "xz", compression_level = 9)
  message("saved data/", nm, ".rda (", format(nrow(x), big.mark = ","), " rows, ",
          length(unique(x$dyad)), " dyads)")
}
