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
# Input: the two tibbles exactly as shipped in data/ up to bilatr 1.0.0d
# (data/theta_directed.rda, data/theta_undirected.rda). They are no longer
# in data/, so the script reads them from git history (`src_commit`) unless
# `src_dir` holds them.
#
# Run from the package root (needs git and the digest package):
#   Rscript data-raw/make_data_release.R
#
# The release is checked at the end: every .rds is identical() to its
# source object, every CSV reads back to the same values, and every
# checksum is recomputed.

library(tibble)

version <- "v2026.11"
release_date <- "2026-11-04"
bilatr_version <- "1.0.0"
src_dir <- "data"
src_commit <- "12733c7" # last commit with data/theta_*.rda
# Zenodo DOIs. The concept DOI (all versions) is fixed. A version's own DOI
# is minted only when its GitHub release is archived, so `version_doi` is NA
# while building a release; fill it in afterwards (dev/data_release/
# POST_RELEASE_EDITS.md) and in R/theta_data.R.
concept_doi <- "10.5281/zenodo.23951818"
version_doi <- c(v2026.11 = "10.5281/zenodo.23951819")[version]
badge <- "https://zenodo.org/badge/1445654273.svg"
out <- file.path("dev", "data_release", version)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load the source tibbles ---------------------------------------------

load_source <- function(nm) {
  rda <- file.path(src_dir, paste0(nm, ".rda"))
  if (!file.exists(rda)) {
    rda <- tempfile(fileext = ".rda")
    status <- system2("git", c("show", sprintf("%s:data/%s.rda", src_commit, nm)),
                      stdout = rda)
    stopifnot(status == 0)
  }
  env <- new.env()
  load(rda, envir = env)
  env[[nm]]
}

src <- list(
  directed = load_source("theta_directed"),
  undirected = load_source("theta_undirected")
)

# --- 2. Data files ------------------------------------------------------------

for (type in names(src)) {
  x <- src[[type]]
  stopifnot(inherits(x, "tbl_df"), !anyNA(x), !anyDuplicated(x[c("dyad", "year")]))
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
    isTRUE(all.equal(as.matrix(back[-(1:2)]), as.matrix(x[-(1:2)]),
                     check.attributes = FALSE, tolerance = 1e-12))
  )
}
data_files <- c(sprintf("theta_%s.rds", names(src)), sprintf("theta_%s.csv.gz", names(src)))

sha256 <- function(f) digest::digest(file = f, algo = "sha256")
md5 <- function(f) unname(tools::md5sum(f))

# --- 3. manifest.json ---------------------------------------------------------

countries <- function(x) sort(unique(c(sub("_.*", "", x$dyad), sub(".*_", "", x$dyad))))
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
  tables = lapply(src, function(x) list(rows = nrow(x), dyads = length(unique(x$dyad)))),
  columns = names(src$directed),
  files = lapply(data_files, function(f) {
    p <- file.path(out, f)
    list(name = f, bytes = file.size(p), sha256 = sha256(p), md5 = md5(p))
  })
)
stopifnot(identical(countries(src$directed), countries(src$undirected)))
jsonlite::write_json(manifest, file.path(out, "manifest.json"),
                     auto_unbox = TRUE, pretty = TRUE, digits = NA)

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
  "- **Correlates of War restriction.** Events were further restricted to",
  "  those where both actors' country codes map to a Correlates of War state",
  "  system member, so GDELT's regional and other non-state codes do not",
  "  occur. This restriction was applied when the event table was built;",
  "  `bilatr::read_gdelt()` does not apply it on its own.",
  "- **Six missing states.** Serbia, Romania, Bosnia and Herzegovina,",
  "  Montenegro, Slovenia and Timor-Leste (`SRB`, `ROU`, `BIH`, `MNE`, `SVN`,",
  "  `TLS`) are absent from both tables, although they are Correlates of War",
  "  system members. Their CAMEO country codes were not matched in the",
  "  country-code crosswalk used to build the event table. So their dyads are",
  "  missing, and their neighbors' other dyads include no events with them.",
  "",
  "This will be corrected in the next data update (`v2027.01`), expected in",
  "January 2027."
)

doi_lines <- c(
  if (!is.na(version_doi)) {
    sprintf("- `%s`: [doi:%s](https://doi.org/%s)", version, version_doi, version_doi)
  } else {
    sprintf("- `%s`: see the Zenodo record (the DOI is minted on release)", version)
  },
  sprintf("- all versions: [doi:%s](https://doi.org/%s)", concept_doi, concept_doi)
)

readme <- c(
  paste("#", title),
  "",
  sprintf("[![DOI](%s)](https://doi.org/%s)", badge, concept_doi),
  "",
  sprintf("Data version **%s** (released %s), built for the R package", version, release_date),
  sprintf("[bilatr](https://github.com/doktorandahl/bilatr) %s.", bilatr_version),
  "",
  "Yearly latent trajectories (theta) of the level of bilateral relations for",
  "every country pair with at least one recorded event, 1991-2025, estimated",
  "with bilatr's hierarchical dynamic item-response model from GDELT 1.0",
  "event data. **Higher values of theta mean more cooperative relations.**",
  "",
  "| file | contents |",
  "|---|---|",
  sprintf("| `theta_directed.rds`, `theta_directed.csv.gz` | directed model: %s rows, %s sender-receiver dyads |",
          n_fmt(nrow(src$directed)), n_fmt(length(unique(src$directed$dyad)))),
  sprintf("| `theta_undirected.rds`, `theta_undirected.csv.gz` | undirected model: %s rows, %s unordered pairs |",
          n_fmt(nrow(src$undirected)), n_fmt(length(unique(src$undirected$dyad)))),
  "| `manifest.json` | version, files, sizes, checksums, row and dyad counts |",
  "| `SHA256SUMS`, `MD5SUMS` | checksums of every file |",
  "",
  "Both tables cover the same 189 countries (see Known gaps).",
  "",
  "## Loading the data",
  "",
  "With bilatr (downloads once, verifies the checksum, and caches the file):",
  "",
  "```r",
  "library(bilatr)",
  sprintf("theta <- download_bilatr_theta(\"directed\", version = \"%s\")", version),
  "```",
  "",
  "Without bilatr:",
  "",
  "```r",
  sprintf("base <- \"https://github.com/doktorandahl/bilatr-data/releases/download/%s/\"", version),
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
  "differ (`bilatr::read_gdelt()` with its defaults). All events are used, root",
  "events or not, with no filter on event codes. Events were further",
  "restricted to Correlates of War states (see Known gaps).",
  "",
  "**Model.** bilatr's `stable_gamma` model (a random-walk latent state per",
  "dyad, a Dirichlet-multinomial likelihood over event classes, and",
  "country-level category offsets), at yearly resolution, with events grouped",
  "by bilatr's `ModifiedRootCode` scheme (16 of its 18 classes occur in the",
  "data, so A = 16), reference category 5 (\"Consult: meet, discuss, or",
  "visit\"), `min_n_events = 1`, and 4 chains of 1,000 warmup and 1,000",
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
  "**Coverage.** Dyad-years before a dyad's first recorded event are dropped:",
  "the likelihood is evaluated only in years with events, so theta before the",
  "first event carries no information from the dyad's own data. After the",
  "first event every year through 2025 is kept, including years without",
  "events (where theta is interpolated or, at the end of a series,",
  "extrapolated by the random walk), so each dyad's series has no gaps and",
  "ends in 2025. All numeric columns are rounded to 3 significant digits.",
  "",
  "**Columns.** One row per dyad-year, sorted by `dyad` then `year`.",
  "",
  "| column | description |",
  "|---|---|",
  "| `dyad` | Two three-letter country codes joined by `_`. The codes are GDELT's `Actor1CountryCode`/`Actor2CountryCode` (CAMEO country codes); every code in these tables is also the ISO 3166-1 alpha-3 code. Directed: `SENDER_RECEIVER`, so `SWE_NOR` is Sweden toward Norway. Undirected: the two codes in alphabetical (C-locale) order, so `NOR_SWE` only. |",
  "| `year` | Integer, 1991-2025. |",
  "| `mean` | Posterior mean of the smoothed theta (given all years). |",
  "| `5%`, `95%` | 5th and 95th percentiles of the smoothed theta posterior (a 90% credible interval). |",
  "| `filtered_mean` | Posterior mean of the forward-filtered theta, which uses only events up to and including that year: for each posterior draw a linear-Bayes filter is run forward from 1991, conditional on that draw's hyperparameters, and the draw-wise filtered state means are averaged. Use this when no look-ahead is allowed (for example as a predictor in forecasting). |",
  "| `filtered_5%`, `filtered_95%` | 5th and 95th percentiles, across posterior draws, of the draw-wise filtered state mean. They reflect hyperparameter uncertainty only, not the filter's own state uncertainty, so they are much narrower than the smoothed interval and stay constant through years without events. They are not a credible interval for the filtered theta. |",
  "",
  known_gaps,
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
  sprintf("%s. Zenodo. <https://doi.org/%s>", version,
          if (is.na(version_doi)) concept_doi else version_doi),
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

writeLines(c(
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
), file.path(out, "LICENSE"))

zenodo <- list(
  upload_type = "dataset",
  title = title,
  version = version,
  publication_date = release_date,
  description = paste0(
    "<p>Yearly latent trajectories (theta) of the level of bilateral relations ",
    "for every country pair with at least one recorded event, 1991-2025, ",
    "estimated with the R package bilatr from GDELT 1.0 event data. Higher ",
    "theta means more cooperative relations. Directed (", n_fmt(length(unique(src$directed$dyad))),
    " sender-receiver dyads) and undirected (", n_fmt(length(unique(src$undirected$dyad))),
    " pairs) tables, as R .rds and gzipped CSV files.</p>",
    "<p>Known gaps: events are restricted to Correlates of War state system ",
    "members, and six states (Serbia, Romania, Bosnia and Herzegovina, ",
    "Montenegro, Slovenia, Timor-Leste) are missing because of a country-code ",
    "crosswalk error; corrected in the next data update (v2027.01), expected ",
    "in January 2027. See README.md.</p>"
  ),
  creators = lapply(authors, function(a) {
    list(name = paste0(a$family, ", ", a$given), orcid = a$orcid)
  }),
  license = "cc-by-4.0",
  access_right = "open",
  keywords = c("international relations", "event data", "GDELT", "CAMEO",
               "bilateral relations", "cooperation", "conflict",
               "item response theory", "Bayesian measurement model", "Stan"),
  related_identifiers = list(
    list(identifier = "https://github.com/doktorandahl/bilatr",
         relation = "isSupplementTo", resource_type = "software"),
    list(identifier = "https://gdeltproject.org/",
         relation = "isDerivedFrom", resource_type = "dataset")
  )
)
jsonlite::write_json(zenodo, file.path(out, ".zenodo.json"), auto_unbox = TRUE, pretty = TRUE)

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
  if (!is.na(version_doi)) c(
    "  - type: doi",
    sprintf("    value: %s", version_doi),
    sprintf("    description: \"Version DOI of data release %s\"", version)
  ),
  "authors:",
  unlist(lapply(authors, function(a) c(
    sprintf("  - family-names: \"%s\"", a$family),
    sprintf("    given-names: \"%s\"", a$given),
    sprintf("    orcid: \"https://orcid.org/%s\"", a$orcid)
  ))),
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

all_files <- sort(setdiff(list.files(out, all.files = TRUE, no.. = TRUE),
                          c("SHA256SUMS", "MD5SUMS")), method = "radix")
writeLines(sprintf("%s  %s", vapply(file.path(out, all_files), sha256, ""), all_files),
           file.path(out, "SHA256SUMS"))
writeLines(sprintf("%s  %s", vapply(file.path(out, all_files), md5, ""), all_files),
           file.path(out, "MD5SUMS"))

# --- 6. Verify ------------------------------------------------------------------

sums <- read.table(file.path(out, "SHA256SUMS"), col.names = c("sum", "file"),
                   colClasses = "character")
stopifnot(identical(sums$sum, unname(vapply(file.path(out, sums$file), sha256, ""))))
if (exists("sha256sum", envir = asNamespace("tools"))) {
  stopifnot(identical(sums$sum, unname(tools::sha256sum(file.path(out, sums$file)))))
}
for (type in names(src)) {
  stopifnot(identical(readRDS(file.path(out, sprintf("theta_%s.rds", type))), src[[type]]))
}
message("release ", version, " written to ", out, ":")
print(data.frame(file = all_files, bytes = file.size(file.path(out, all_files))))
