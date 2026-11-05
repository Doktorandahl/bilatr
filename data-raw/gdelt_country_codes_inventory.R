# Inventory of the country codes GDELT 1.0 actually uses (1.0.0f, Part 1a).
#
# Not shipped. Scans raw GDELT 1.0 export files and writes compact count
# tables to data-raw/gdelt_country_codes/:
#
#   * code_counts.csv  -- every Actor1/Actor2CountryCode value, with row
#     counts by event year (SQLDATE), side and file era, over ALL rows
#     (no row filters: the equivalent of read_gdelt(actor_types = NULL,
#     cross_border = FALSE));
#   * dyad_counts.csv  -- directed code-pair counts by event year, for rows
#     passing the production filter (both country codes present and
#     different; both actors GOV/MIL/SPY in Type1-3Code), i.e.
#     read_gdelt(actor_types = c("GOV", "MIL", "SPY"),
#     actor_type_match = "both", cross_border = TRUE), restricted to
#     SQLDATE 1991-01-01 .. 2025-12-31. Column `prod_event` flags rows
#     whose EventCode the production event table (original_code/
#     gdelt_bilatr.rds) would have kept: in bilatr's cameo_lookup, and not
#     one of the 21 "BilatrClass 0" codes (010-015, 017, 024x, 025x, 040,
#     041) that the production table excludes (see the 1.0.0f review);
#   * prefix_counts.csv -- for rows whose Actor1/Actor2CountryCode is EMPTY
#     but whose Actor1/Actor2Code is not, the first three characters of the
#     actor code, by event year, side and era. GDELT derives the country
#     code from the actor code's leading country segment, and leaves it
#     empty when its country table lacks that segment (e.g. Romania's
#     actor codes are "ROU...", but the table has "ROM"; see 1c);
#   * recovered_dyad_counts.csv -- like dyad_counts.csv, but with an empty
#     country code filled from the actor-code prefix when that prefix is an
#     ISO3 or CAMEO country code; only pairs where at least one side was
#     filled. This measures what the country-code field alone loses.
#
# Inputs: a local copy of the GDELT 1.0 events archive
# (https://data.gdeltproject.org/events/), with the yearly YYYY.zip,
# monthly YYYYMM.zip and daily YYYYMMDD[.export.CSV].zip files in one
# directory, given by GDELT_RAW_DIR. Every file from 1991 onward is
# scanned, which takes about an hour on 8 cores; set GDELT_INVENTORY_FILES
# to a comma-separated list of basenames to scan a subset.
#
# Uses data.table::fread(cmd = "unzip -p ...") with `select`, rather than
# read_gdelt(), because only 11 of the 57/58 columns are needed and the scan
# covers the whole archive; the filters are the same as read_gdelt()'s
# (.gdelt_apply_filters()).

library(data.table)

raw_dir <- Sys.getenv(
  "GDELT_RAW_DIR",
  "/Users/david/ViEWS Dropbox/David Randahl/Mac/Documents/active/gdelt_raw"
)
out_dir <- file.path("data-raw", "gdelt_country_codes")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tmp_dir <- Sys.getenv("GDELT_INVENTORY_TMP", tempdir())
n_cores <- as.integer(Sys.getenv("GDELT_INVENTORY_CORES", "6"))

files <- list.files(raw_dir, pattern = "^[0-9]{4,8}(\\.export\\.CSV)?\\.zip$")
stem <- sub("\\..*$", "", files)
files <- files[as.integer(substr(stem, 1, 4)) >= 1991]
subset <- Sys.getenv("GDELT_INVENTORY_FILES", "")
if (nzchar(subset)) {
  files <- intersect(files, strsplit(subset, ",", fixed = TRUE)[[1]])
}

era_of <- function(f) {
  n <- nchar(sub("\\..*$", "", f))
  c(`4` = "yearly", `6` = "monthly", `8` = "daily")[[as.character(n)]]
}

# 1-based positions in the GDELT 1.0 schema (checked against .gdelt_schema())
cols <- c(SQLDATE = 2L, A1C = 6L, A1 = 8L, A1T1 = 13L, A1T2 = 14L, A1T3 = 15L,
          A2C = 16L, A2 = 18L, A2T1 = 23L, A2T2 = 24L, A2T3 = 25L,
          EventCode = 27L)
types <- c("GOV", "MIL", "SPY")

# Prefixes that count as a country segment: every ISO3 code plus every
# CAMEO lookup code (CAMEO.country.txt, as of the scan date).
lookup_codes <- utils::read.delim(
  "https://www.gdeltproject.org/data/lookups/CAMEO.country.txt",
  colClasses = "character", strip.white = TRUE
)$CODE
# Event codes the production event table kept (see header)
class0_codes <- c("010", "011", "012", "013", "014", "015", "017",
                  "024", "0241", "0242", "0243", "0244",
                  "025", "0251", "0252", "0253", "0254", "0255", "0256",
                  "040", "041")
prod_event_codes <- setdiff(bilatr::cameo_lookup$CAMEOEVENTCODE, class0_codes)

country_prefixes <- union(stats::na.omit(countrycode::codelist$iso3c), lookup_codes)

scan_one <- function(f) {
  path <- file.path(raw_dir, f)
  dt <- fread(
    cmd = sprintf("unzip -p %s", shQuote(path)),
    sep = "\t", quote = "", header = FALSE, na.strings = "",
    select = unname(cols), col.names = names(cols),
    colClasses = list(character = unname(cols)),
    tmpdir = tmp_dir, showProgress = FALSE, nThread = 2L
  )
  dt[, year := as.integer(substr(SQLDATE, 1, 4))]
  era <- era_of(f)

  codes <- rbind(
    dt[, .N, by = .(year, code = A1)][, side := 1L],
    dt[, .N, by = .(year, code = A2)][, side := 2L]
  )
  codes[, `:=`(era = era, file = f)]

  keep <- dt[
    !is.na(A1) & !is.na(A2) & A1 != A2 &
      (A1T1 %chin% types | A1T2 %chin% types | A1T3 %chin% types) &
      (A2T1 %chin% types | A2T2 %chin% types | A2T3 %chin% types) &
      SQLDATE >= "19910101" & SQLDATE <= "20251231"
  ]
  dyads <- keep[, .N, by = .(year, a1 = A1, a2 = A2,
                             prod_event = EventCode %chin% prod_event_codes)]
  dyads[, `:=`(era = era, file = f)]

  dt[, `:=`(P1 = substr(A1C, 1, 3), P2 = substr(A2C, 1, 3))]
  prefixes <- rbind(
    dt[is.na(A1) & !is.na(A1C), .N, by = .(year, prefix = P1)][, side := 1L],
    dt[is.na(A2) & !is.na(A2C), .N, by = .(year, prefix = P2)][, side := 2L]
  )
  prefixes[, era := era]

  dt[, `:=`(
    R1 = fifelse(is.na(A1) & P1 %chin% country_prefixes, P1, A1),
    R2 = fifelse(is.na(A2) & P2 %chin% country_prefixes, P2, A2)
  )]
  rec <- dt[
    (is.na(A1) | is.na(A2)) &
      !is.na(R1) & !is.na(R2) & R1 != R2 &
      (A1T1 %chin% types | A1T2 %chin% types | A1T3 %chin% types) &
      (A2T1 %chin% types | A2T2 %chin% types | A2T3 %chin% types) &
      SQLDATE >= "19910101" & SQLDATE <= "20251231"
  ]
  recovered <- rec[, .N, by = .(year, a1 = R1, a2 = R2,
                                filled1 = is.na(A1), filled2 = is.na(A2),
                                prod_event = EventCode %chin% prod_event_codes)]
  recovered[, era := era]

  list(file = f, codes = codes, dyads = dyads, prefixes = prefixes,
       recovered = recovered, n = nrow(dt))
}

t0 <- Sys.time()
res <- parallel::mclapply(files, function(f) {
  tryCatch(scan_one(f), error = function(e) {
    message("FAILED ", f, ": ", conditionMessage(e))
    NULL
  })
}, mc.cores = n_cores, mc.preschedule = FALSE)
failed <- files[vapply(res, is.null, logical(1))]
res <- Filter(Negate(is.null), res)
message(sprintf("Scanned %d files (%d failed) in %.1f min",
                length(res), length(failed),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

codes <- rbindlist(lapply(res, `[[`, "codes"))
codes <- codes[, .(n = sum(N)), by = .(code, year, side, era)]
dyads <- rbindlist(lapply(res, `[[`, "dyads"))
dyads <- dyads[, .(n = sum(N)), by = .(year, a1, a2, prod_event, era)]

prefixes <- rbindlist(lapply(res, `[[`, "prefixes"))
prefixes <- prefixes[, .(n = sum(N)), by = .(prefix, year, side, era)]
recovered <- rbindlist(lapply(res, `[[`, "recovered"))
recovered <- recovered[, .(n = sum(N)), by = .(year, a1, a2, filled1, filled2, prod_event, era)]
fwrite(prefixes[order(prefix, year, side, era)], file.path(out_dir, "prefix_counts.csv"))
fwrite(recovered[order(year, a1, a2)], file.path(out_dir, "recovered_dyad_counts.csv"))

fwrite(codes[order(code, year, side, era)], file.path(out_dir, "code_counts.csv"))
fwrite(dyads[order(year, a1, a2, era)], file.path(out_dir, "dyad_counts.csv"))
scanned <- rbind(
  data.table(file = vapply(res, `[[`, character(1), "file"),
             rows = vapply(res, function(r) as.numeric(r$n), numeric(1)),
             failed = FALSE),
  data.table(file = failed, rows = NA_real_, failed = rep(TRUE, length(failed)))
)
fwrite(scanned[order(file)], file.path(out_dir, "files_scanned.csv"))
