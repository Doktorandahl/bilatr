# make_theta_sample.R
#
# Builds data/theta_sample.rda (documented in R/data.R, ?theta_sample): a
# small extract of the bilatr-data release for examples and vignettes.
# Twenty directed dyads (ten country pairs, both directions) and their ten
# undirected counterparts, stacked in one tibble with a `type` column.
#
# Input: the release files written by data-raw/make_data_release.R.
#
# Run from the package root:
#   Rscript data-raw/make_theta_sample.R

library(tibble)

release <- file.path("dev", "data_release", "v2026.11.1")

pairs <- list(
  c("USA", "CHN"), c("USA", "RUS"), c("RUS", "UKR"), c("CHN", "JPN"),
  c("KOR", "PRK"), c("IND", "PAK"), c("ISR", "IRN"), c("GBR", "FRA"),
  c("DEU", "FRA"), c("NOR", "SWE")
)
directed <- unlist(lapply(pairs, function(p) c(paste(p, collapse = "_"), paste(rev(p), collapse = "_"))))
undirected <- vapply(pairs, function(p) paste(sort(p, method = "radix"), collapse = "_"), "")

d <- readRDS(file.path(release, "theta_directed.rds"))
u <- readRDS(file.path(release, "theta_undirected.rds"))
stopifnot(all(directed %in% d$dyad), all(undirected %in% u$dyad))

theta_sample <- rbind(
  tibble(type = "directed", d[d$dyad %in% directed, ]),
  tibble(type = "undirected", u[u$dyad %in% undirected, ])
)
stopifnot(
  identical(names(theta_sample), c("type", names(d))),
  length(unique(theta_sample$dyad[theta_sample$type == "directed"])) == 20,
  length(unique(theta_sample$dyad[theta_sample$type == "undirected"])) == 10
)

save(theta_sample, file = file.path("data", "theta_sample.rda"),
     compress = "xz", compression_level = 9)
message("saved data/theta_sample.rda: ", nrow(theta_sample), " rows, ",
        file.size(file.path("data", "theta_sample.rda")), " bytes")
