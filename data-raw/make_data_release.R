# make_data_release.R
#
# Builds the files of one bilatr-data release (the estimated theta
# trajectories, published outside the package in the public repository
# doktorandahl/bilatr-data and archived on Zenodo; fetched in R with
# download_bilatr_theta()).
#
# Writes to dev/data_release/<version>/:
#   theta_directed.rds, theta_undirected.rds        the tibbles, xz-compressed
#   theta_directed.csv.gz, theta_undirected.csv.gz   the same, for non-R users
#   manifest.json                                    version, files, counts
#   README.md, LICENSE, CITATION.cff, .zenodo.json   repository files
#   SHA256SUMS, MD5SUMS                              checksums of every file above
#
# Input: the two tibbles written by the replication package's
# 04_process_theta.R (data/processed/theta_directed.rds and
# theta_undirected.rds there), copied to dev/data_release/source/<version>/.
# Set `version`, `release_date` and `n_observed` below for each release.
#
# Run from the package root (needs the digest and jsonlite packages):
#   Rscript data-raw/make_data_release.R
#
# The release is checked at the end: every .rds is identical() to its
# source object, every CSV reads back to the same values, and every
# checksum is recomputed.

library(tibble)

version <- "v2026.11.1" # e.g. "v2026.12"
release_date <- "2026-11-10" # the day you publish the GitHub release
bilatr_version <- "1.0.0"
src_dir <- file.path("dev", "data_release", "source", version)
concept_doi <- "10.5281/zenodo.23951818"
version_doi <- "10.5281/zenodo.24109618" # minted by Zenodo on release (NA before); fill in afterwards
badge <- "https://zenodo.org/badge/1445654273.svg"
# Dyad-years with at least one observed event, per table. The theta tables
# do not carry event counts, so these come from the fit (sum of is_obs in
# the Stan data, within the released dyad-years); update them with each refit.
n_observed <- c(directed = 130848L, undirected = 76821L)
out <- file.path("dev", "data_release", version)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

load_source <- function(nm) {
  rds <- file.path(src_dir, paste0(nm, ".rds"))
  if (!file.exists(rds)) {
    stop("missing ", rds)
  }
  readRDS(rds)
}

src <- list(
  directed = load_source("theta_directed"),
  undirected = load_source("theta_undirected")
)

# --- 2. Data files ------------------------------------------------------------

for (type in names(src)) {
  x <- src[[type]]
  stopifnot(
    inherits(x, "tbl_df"),
    !anyNA(x),
    !anyDuplicated(x[c("dyad", "year")])
  )
  rds <- file.path(out, sprintf("theta_%s.rds", type))
  saveRDS(x, rds, compress = "xz")
  stopifnot(identical(readRDS(rds), x))

  csv <- file.path(out, sprintf("theta_%s.csv.gz", type))
  con <- gzcon(file(csv, "wb"))
  utils::write.csv(as.data.frame(x), con, row.names = FALSE)
  close(con)
  back <- utils::read.csv(csv, check.names = FALSE)
  stopifnot(
    identical(names(back), names(x)),
    identical(back$dyad, x$dyad),
    identical(as.integer(back$year), x$year),
    isTRUE(all.equal(
      as.matrix(back[-(1:2)]),
      as.matrix(x[-(1:2)]),
      check.attributes = FALSE,
      tolerance = 1e-12
    ))
  )
}
data_files <- c(
  sprintf("theta_%s.rds", names(src)),
  sprintf("theta_%s.csv.gz", names(src))
)

sha256 <- function(f) digest::digest(file = f, algo = "sha256")
md5 <- function(f) unname(tools::md5sum(f))

# --- 3. manifest.json ---------------------------------------------------------

countries <- function(x) {
  sort(unique(c(sub("_.*", "", x$dyad), sub(".*_", "", x$dyad))))
}
manifest <- list(
  name = "bilatr theta trajectories",
  data_version = version,
  release_date = release_date,
  built_for_bilatr = bilatr_version,
  description = paste(
    "Yearly latent bilateral-relations trajectories (theta) for every country",
    "pair with at least one recorded event, 1991-2025, estimated with the R",
    "package bilatr from GDELT 1.0 event data. Higher theta means more",
    "cooperative relations."
  ),
  license = "CC-BY-4.0",
  years = c(1991L, 2025L),
  n_countries = length(countries(src$directed)),
  tables = lapply(src, function(x) {
    list(rows = nrow(x), dyads = length(unique(x$dyad)))
  }),
  columns = names(src$directed),
  files = lapply(data_files, function(f) {
    p <- file.path(out, f)
    list(name = f, bytes = file.size(p), sha256 = sha256(p), md5 = md5(p))
  })
)
stopifnot(identical(countries(src$directed), countries(src$undirected)))
jsonlite::write_json(
  manifest,
  file.path(out, "manifest.json"),
  auto_unbox = TRUE,
  pretty = TRUE,
  digits = NA
)

# --- 4. Repository files --------------------------------------------------------

authors <- list(
  list(given = "David", family = "Randahl", orcid = "0000-0003-1069-6067"),
  list(given = "Magnus", family = "Lundgren", orcid = "0000-0002-9961-3645"),
  list(given = "Carl", family = "Vikberg", orcid = "0000-0001-7523-9196")
)
title <- "bilatr theta trajectories: latent bilateral relations, 1991-2025"
n_fmt <- function(x) format(x, big.mark = ",")

known_gaps <- c(
  "## Known gaps",
  "",
  "- **Four absent states.** GDELT 1.0 contains no events in which the",
  "  government actors of Slovenia, Bosnia and Herzegovina, Kosovo or",
  "  Montenegro are the source or target, once the actor filter is applied:",
  "  GDELT 1.0 has no country code for them. These four Correlates of War",
  "  states are therefore absent from both tables. Their absence reflects the",
  "  coverage of the data source, not an absence of interaction, and the",
  "  tables say nothing about their bilateral relations.",
  "- **Missing and incomplete days.** For 24 days, most of them between",
  "  14 June and 1 July 2025, GDELT provides no daily event files. For a",
  "  further 13 days in 1991-1995, the files contain events, but none that pass",
  "  the filters. Because events are aggregated to dyad-years, these gaps",
  "  reduce the number of events in the affected years but do not create",
  "  missing dyad-years. The model conditions on the total number of events and",
  "  measures relations from their composition, so a gap reduces the precision",
  "  of the estimates for the affected years rather than shifting their level,",
  "  provided the missing days are not systematically different from the rest",
  "  of the year. The longest gap, 14 June to 1 July 2025, is the exception:",
  "  it coincides with the war between Israel and Iran in June 2025, so the",
  "  2025 estimates for the dyads involved in that conflict are based on a",
  "  composition of events from which the war itself is largely missing."
)

changes <- c(
  "## Changes from v2026.11",
  "",
  "This release replaces `v2026.11`, which should no longer be used.",
  "",
  "- **Country coverage.** `v2026.11` lacked six Correlates of War states",
  "  (Serbia, Romania, Bosnia and Herzegovina, Montenegro, Slovenia and",
  "  Timor-Leste) because of an error in the country-code crosswalk. The",
  "  event table is now rebuilt from raw GDELT with `bilatr::read_gdelt()` and",
  "  `bilatr::filter_gdelt_states()`, which recovers Serbia, Romania and",
  "  Timor-Leste. Only the four states listed under Known gaps remain absent.",
  sprintf(
    "  The tables now cover %s countries.",
    length(countries(src$directed))
  ),
  "- **Country codes** are harmonized to ISO 3166 alpha-3, and a country is",
  "  included only in the years it is a member of the Correlates of War state",
  "  system (see `dyad` under Columns).",
  "- **Event codes.** Events with a Goldstein score between -1 and 1 and events",
  "  with non-CAMEO codes were excluded in `v2026.11` too, although its README",
  "  said there was no filter on event codes; the filter is now part of the",
  "  scripted pipeline and documented (see Events).",
  "- **Columns.** New `median` and `filtered_median` columns; all numeric",
  "  columns are rounded to 4 decimal places (previously 3 significant",
  "  digits).",
  "- **Sampling.** The models were refitted with 8 chains (seeds 101-108",
  "  directed, 201-208 undirected) instead of 4. Correction: the README",
  "  archived with the `v2026.11.1` tag says 4 chains, but the models used 8."
)

doi_lines <- c(
  if (!is.na(version_doi)) {
    sprintf(
      "- `%s`: [doi:%s](https://doi.org/%s)",
      version,
      version_doi,
      version_doi
    )
  } else {
    sprintf(
      "- `%s`: see the Zenodo record (the DOI is minted on release)",
      version
    )
  },
  sprintf(
    "- all versions: [doi:%s](https://doi.org/%s)",
    concept_doi,
    concept_doi
  )
)

readme <- c(
  paste("#", title),
  "",
  sprintf("[![DOI](%s)](https://doi.org/%s)", badge, concept_doi),
  "",
  sprintf(
    "Data version **%s** (released %s), built for the R package",
    version,
    release_date
  ),
  sprintf(
    "[bilatr](https://github.com/doktorandahl/bilatr) %s.",
    bilatr_version
  ),
  "",
  "Yearly latent trajectories (theta) of the level of bilateral relations for",
  "every country pair with at least one recorded event, 1991-2025, estimated",
  "with bilatr's hierarchical dynamic item-response model from GDELT 1.0",
  "event data. **Higher values of theta mean more cooperative relations.**",
  "",
  "| file | contents |",
  "|---|---|",
  sprintf(
    "| `theta_directed.rds`, `theta_directed.csv.gz` | directed model: %s rows, %s sender-receiver dyads |",
    n_fmt(nrow(src$directed)),
    n_fmt(length(unique(src$directed$dyad)))
  ),
  sprintf(
    "| `theta_undirected.rds`, `theta_undirected.csv.gz` | undirected model: %s rows, %s unordered pairs |",
    n_fmt(nrow(src$undirected)),
    n_fmt(length(unique(src$undirected$dyad)))
  ),
  "| `manifest.json` | version, files, sizes, checksums, row and dyad counts |",
  "| `SHA256SUMS`, `MD5SUMS` | checksums of every file |",
  "",
  sprintf(
    "Both tables cover the same %s countries (see Known gaps).",
    length(countries(src$directed))
  ),
  "",
  "## Loading the data",
  "",
  "With bilatr (downloads once, verifies the checksum, and caches the file):",
  "",
  "```r",
  "library(bilatr)",
  sprintf(
    "theta <- download_bilatr_theta(\"directed\", version = \"%s\")",
    version
  ),
  "```",
  "",
  "Without bilatr:",
  "",
  "```r",
  sprintf(
    "base <- \"https://github.com/doktorandahl/bilatr-data/releases/download/%s/\"",
    version
  ),
  "theta <- readRDS(url(paste0(base, \"theta_directed.rds\")))   # a tibble",
  "theta <- read.csv(paste0(base, \"theta_directed.csv.gz\"), check.names = FALSE)",
  "```",
  "",
  "Outside R, any CSV reader that handles gzip will do. The interval columns",
  "have non-syntactic names (`5%`, `95%`, ...); in R, quote them with",
  "backticks.",
  "",
  "## What the data are",
  "",
  "**Events.** GDELT 1.0 events (the historical backfiles and the daily",
  "updates) dated 1991-01-01 through 2025-12-31, accessed 2026-01-31. An event",
  "is kept when both actors carry `GOV`, `MIL` or `SPY` in any of their",
  "`Type1Code`/`Type2Code`/`Type3Code` fields and the two actors' country codes",
  "differ (`bilatr::read_gdelt()` with its defaults; an empty country code is",
  "filled from the actor code's prefix with `recover_country_codes = TRUE`).",
  "All events are used, root events or not, except those with a Goldstein",
  "score between -1 and 1 (CAMEO codes 010-015, 017, 024x, 025x, 040 and",
  "041) and those with non-CAMEO event codes. Events are further restricted",
  "to those in which both countries are members of the Correlates of War",
  "state system (v2024) in the year of the event",
  "(`bilatr::filter_gdelt_states()`), which removes GDELT's regional and other",
  "non-state codes.",
  "",
  "**Model.** bilatr's `stable_gamma` model (a random-walk latent state per",
  "dyad, a Dirichlet-multinomial likelihood over event classes, and",
  "country-level category offsets), at yearly resolution, with events grouped",
  "by bilatr's `ModifiedRootCode` scheme (16 of its 18 classes occur in the",
  "data, so A = 16), reference category 5 (\"Consult: meet, discuss, or",
  "visit\"), `min_n_events = 1`, and 8 chains of 1,000 warmup and 1,000",
  "sampling iterations each. Because the model includes country-level",
  "category offsets, theta is measured net of each country's baseline mix of",
  "event types.",
  "",
  "**Orientation and scale.** Higher theta makes the reference category,",
  "a cooperative class, more likely relative to the other action types, so",
  "higher values mean more cooperative relations. The sign is fixed by the",
  "model's orientation constraint (the reference class's discrimination is",
  "non-negative), the location by the prior on the initial states (mean",
  "pinned at 0 across dyads), and the scale by normalizing the discrimination",
  "parameters to root-mean-square 1.",
  "",
  "**Coverage.** The models are fitted on each dyad's full 1991-2025 series,",
  "with the likelihood evaluated only in dyad-years with at least one event.",
  "Dyad-years before a dyad's first recorded event are then dropped, since",
  "theta before the first event carries no information from the dyad's own",
  "data, and so are dyad-years in which either country was not a member of",
  "the state system (for example, dyads with Czechoslovakia end in 1992, and",
  "those with Czechia or Slovakia start in 1993 at the earliest). From the",
  "first event on, every year is kept, including years without events (where",
  "theta is interpolated or, at the end of a series, extrapolated by the",
  "random walk), so each dyad's series has no gaps.",
  sprintf(
    paste(
      "The tables have %s directed and %s undirected dyad-years, of which",
      "%s and %s have at least one observed event."
    ),
    n_fmt(nrow(src$directed)),
    n_fmt(nrow(src$undirected)),
    n_fmt(n_observed[["directed"]]),
    n_fmt(n_observed[["undirected"]])
  ),
  "All numeric columns are rounded to 4 decimal places.",
  "",
  "**Columns.** One row per dyad-year, sorted by `dyad` then `year`.",
  "",
  "| column | description |",
  "|---|---|",
  "| `dyad` | Two three-letter country codes joined by `_`. GDELT's `Actor1CountryCode`/`Actor2CountryCode` (CAMEO country codes), harmonized to ISO 3166 alpha-3 (`bilatr::gdelt_countries`): `ROU` for Romania, `SRB` for Serbia (Correlates of War state 345, including the Federal Republic of Yugoslavia and Serbia and Montenegro, in every year), `TLS` for Timor-Leste, and `CSK` (the former ISO code) for Czechoslovakia in 1991-1992. Directed: `SENDER_RECEIVER`, so `SWE_NOR` is Sweden toward Norway. Undirected: the two codes in alphabetical (C-locale) order, so `NOR_SWE` only. |",
  "| `year` | Integer, 1991-2025. |",
  "| `mean`, `median` | Posterior mean and median of the smoothed theta (given all years). |",
  "| `5%`, `95%` | 5th and 95th percentiles of the smoothed theta posterior (a 90% credible interval). |",
  "| `filtered_mean`, `filtered_median` | Mean and median, across posterior draws, of the forward-filtered theta, which uses only events up to and including that year: for each posterior draw a linear-Bayes filter is run forward from 1991, conditional on that draw's hyperparameters, and the draw-wise filtered state means are summarized. Use these when no look-ahead is allowed (for example as a predictor in forecasting). |",
  "| `filtered_5%`, `filtered_95%` | 5th and 95th percentiles, across posterior draws, of the draw-wise filtered state mean. They reflect hyperparameter uncertainty only, not the filter's own state uncertainty, so they are much narrower than the smoothed interval and stay constant through years without events. They are not a credible interval for the filtered theta. |",
  "",
  known_gaps,
  "",
  changes,
  "",
  "## License",
  "",
  "[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You may share and",
  "adapt the data for any purpose, provided you give appropriate credit.",
  "",
  "## How to cite",
  "",
  "Please cite the data release by the DOI of the version you use:",
  "",
  doi_lines,
  "",
  "Randahl, D., Lundgren, M. and Vikberg, C. (2026). bilatr theta",
  "trajectories: latent bilateral relations, 1991-2025. Data release",
  sprintf(
    "%s. Zenodo. <https://doi.org/%s>",
    version,
    if (is.na(version_doi)) concept_doi else version_doi
  ),
  "",
  "Please also cite the bilatr package (`citation(\"bilatr\")` in R). The",
  "underlying events are GDELT's; please also cite:",
  "",
  "- Leetaru, K. and Schrodt, P. A. (2013). \"GDELT: Global Data on Events,",
  "  Location and Tone, 1979-2012.\" Paper presented at the International",
  "  Studies Association annual convention, San Francisco.",
  "- The GDELT Project. Global Database of Events, Language, and Tone (GDELT).",
  "  Accessed 31 January 2026, <https://gdeltproject.org/>."
)
writeLines(readme, file.path(out, "README.md"))

writeLines(
  c(
    sprintf("bilatr theta trajectories, data version %s", version),
    "Copyright (c) 2026 David Randahl, Magnus Lundgren and Carl Vikberg",
    "",
    "The data files in this repository are licensed under the Creative Commons",
    "Attribution 4.0 International License (CC BY 4.0).",
    "",
    "You are free to share and adapt the material for any purpose, even",
    "commercially, provided you give appropriate credit, provide a link to the",
    "license, and indicate if changes were made.",
    "",
    "Full license text: https://creativecommons.org/licenses/by/4.0/legalcode"
  ),
  file.path(out, "LICENSE")
)

zenodo <- list(
  upload_type = "dataset",
  title = title,
  version = version,
  publication_date = release_date,
  description = paste0(
    "<p>Yearly latent trajectories (theta) of the level of bilateral relations ",
    "for every country pair with at least one recorded event, 1991-2025, ",
    "estimated with the R package bilatr from GDELT 1.0 event data. Higher ",
    "theta means more cooperative relations. Directed (",
    n_fmt(length(unique(src$directed$dyad))),
    " sender-receiver dyads) and undirected (",
    n_fmt(length(unique(src$undirected$dyad))),
    " pairs) tables, as R .rds and gzipped CSV files.</p>",
    "<p>The tables cover ", length(countries(src$directed)), " Correlates of ",
    "War state system members, each in the years of its membership. ",
    "Slovenia, Bosnia and Herzegovina, Kosovo and Montenegro are absent ",
    "because GDELT 1.0 records no government events for them. GDELT provides ",
    "no daily files for 24 days, most of them between 14 June and 1 July 2025 ",
    "(the war between Israel and Iran), and 13 days in 1991-1995 have no ",
    "events that pass the filters. See README.md.</p>",
    "<p>This version replaces v2026.11, which lacked six states (Serbia, ",
    "Romania, Timor-Leste, Bosnia and Herzegovina, Montenegro and Slovenia) ",
    "because of a country-code crosswalk error; the first three are now ",
    "included. It adds median columns and rounds all values to four ",
    "decimal places.</p>"
  ),
  creators = lapply(authors, function(a) {
    list(name = paste0(a$family, ", ", a$given), orcid = a$orcid)
  }),
  license = "cc-by-4.0",
  access_right = "open",
  keywords = c(
    "international relations",
    "event data",
    "GDELT",
    "CAMEO",
    "bilateral relations",
    "cooperation",
    "conflict",
    "item response theory",
    "Bayesian measurement model",
    "Stan"
  ),
  related_identifiers = list(
    list(
      identifier = "https://github.com/doktorandahl/bilatr",
      relation = "isSupplementTo",
      resource_type = "software"
    ),
    list(
      identifier = "https://gdeltproject.org/",
      relation = "isDerivedFrom",
      resource_type = "dataset"
    )
  )
)
jsonlite::write_json(
  zenodo,
  file.path(out, ".zenodo.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cff <- c(
  "cff-version: 1.2.0",
  "message: \"If you use these data, please cite them as below.\"",
  "type: dataset",
  sprintf("title: \"%s\"", title),
  sprintf("version: \"%s\"", version),
  sprintf("date-released: \"%s\"", release_date),
  "license: CC-BY-4.0",
  "repository-code: \"https://github.com/doktorandahl/bilatr-data\"",
  "url: \"https://github.com/doktorandahl/bilatr\"",
  sprintf("doi: %s", concept_doi),
  "identifiers:",
  "  - type: doi",
  sprintf("    value: %s", concept_doi),
  "    description: \"Concept DOI: all versions of the data release\"",
  if (!is.na(version_doi)) {
    c(
      "  - type: doi",
      sprintf("    value: %s", version_doi),
      sprintf("    description: \"Version DOI of data release %s\"", version)
    )
  },
  "authors:",
  unlist(lapply(authors, function(a) {
    c(
      sprintf("  - family-names: \"%s\"", a$family),
      sprintf("    given-names: \"%s\"", a$given),
      sprintf("    orcid: \"https://orcid.org/%s\"", a$orcid)
    )
  })),
  "keywords:",
  "  - international relations",
  "  - event data",
  "  - GDELT",
  "  - bilateral relations",
  "references:",
  "  - type: conference-paper",
  "    title: \"GDELT: Global Data on Events, Location and Tone, 1979-2012\"",
  "    authors:",
  "      - family-names: \"Leetaru\"",
  "        given-names: \"Kalev\"",
  "      - family-names: \"Schrodt\"",
  "        given-names: \"Philip A.\"",
  "    year: 2013",
  "    conference:",
  "      name: \"International Studies Association Annual Convention\"",
  "      city: \"San Francisco\""
)
writeLines(cff, file.path(out, "CITATION.cff"))

# --- 5. Checksums (every file except the checksum files themselves) ----------

all_files <- sort(
  setdiff(
    list.files(out, all.files = TRUE, no.. = TRUE),
    c("SHA256SUMS", "MD5SUMS")
  ),
  method = "radix"
)
writeLines(
  sprintf("%s  %s", vapply(file.path(out, all_files), sha256, ""), all_files),
  file.path(out, "SHA256SUMS")
)
writeLines(
  sprintf("%s  %s", vapply(file.path(out, all_files), md5, ""), all_files),
  file.path(out, "MD5SUMS")
)

# --- 6. Verify ------------------------------------------------------------------

sums <- read.table(
  file.path(out, "SHA256SUMS"),
  col.names = c("sum", "file"),
  colClasses = "character"
)
stopifnot(identical(
  sums$sum,
  unname(vapply(file.path(out, sums$file), sha256, ""))
))
if (exists("sha256sum", envir = asNamespace("tools"))) {
  stopifnot(identical(
    sums$sum,
    unname(tools::sha256sum(file.path(out, sums$file)))
  ))
}
for (type in names(src)) {
  stopifnot(identical(
    readRDS(file.path(out, sprintf("theta_%s.rds", type))),
    src[[type]]
  ))
}
message("release ", version, " written to ", out, ":")
print(data.frame(
  file = all_files,
  bytes = file.size(file.path(out, all_files))
))
