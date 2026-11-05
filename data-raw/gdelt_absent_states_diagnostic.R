# Count-only diagnostic: events that might involve Bosnia and Herzegovina,
# Slovenia, Montenegro or Kosovo (1.0.0f review, "Still absent").
#
# Not shipped; no package code depends on it. GDELT 1.0 never writes a
# country code (or an actor-code prefix) for these four states, so their
# events cannot be assigned from the country codes. This script measures
# whether a later, validated recovery from other fields could be worth
# building. It does NOT recover anything.
#
# Scope: the production filter with the empty-country-code actor treated as
# unknown: both actors carry GOV, MIL or SPY in a Type1-3Code field, one
# actor's country code is empty and the other's is not, SQLDATE 1991-2025.
# Such a row is a candidate for a country when the empty-code actor has
#   (a) Actor*Name matching BOSNIA|SARAJEVO / SLOVEN|LJUBLJANA /
#       MONTENEG|PODGORICA / KOSOV|PRISTINA, or
#   (b) Actor*Geo_CountryCode (FIPS 10-4) BK / SI / MJ / KV.
#
# Writes data-raw/gdelt_country_codes/absent_states_counts.csv (country,
# decade, criterion: name_only / geo_only / both, n) and
# absent_states_examples.csv (20 random candidate rows per country, drawn
# from a reservoir of up to 200 per file, weighted so that every candidate
# row is about equally likely).
#
# Inputs as for gdelt_country_codes_inventory.R (GDELT_RAW_DIR,
# GDELT_INVENTORY_CORES, GDELT_INVENTORY_TMP); about 7 minutes on 6 cores.

library(data.table)

raw_dir <- Sys.getenv(
  "GDELT_RAW_DIR",
  "/Users/david/ViEWS Dropbox/David Randahl/Mac/Documents/active/gdelt_raw"
)
out_dir <- file.path("data-raw", "gdelt_country_codes")
tmp_dir <- Sys.getenv("GDELT_INVENTORY_TMP", tempdir())
n_cores <- as.integer(Sys.getenv("GDELT_INVENTORY_CORES", "6"))

files <- list.files(raw_dir, pattern = "^[0-9]{4,8}(\\.export\\.CSV)?\\.zip$")
files <- files[as.integer(substr(files, 1, 4)) >= 1991]

# 1-based positions in the GDELT 1.0 schema (checked against .gdelt_schema())
cols <- c(GLOBALEVENTID = 1L, SQLDATE = 2L,
          A1C = 6L, A1Name = 7L, A1 = 8L, A1T1 = 13L, A1T2 = 14L, A1T3 = 15L,
          A2C = 16L, A2Name = 17L, A2 = 18L, A2T1 = 23L, A2T2 = 24L, A2T3 = 25L,
          EventCode = 27L, A1Geo = 38L, A2Geo = 45L)
types <- c("GOV", "MIL", "SPY")

name_patterns <- c(BIH = "BOSNIA|SARAJEVO", SVN = "SLOVEN|LJUBLJANA",
                   MNE = "MONTENEG|PODGORICA", KOS = "KOSOV|PRISTINA")
fips <- c(BIH = "BK", SVN = "SI", MNE = "MJ", KOS = "KV")

# One candidate row per (event, side, country).
candidates <- function(dt, side) {
  other <- if (side == 1L) 2L else 1L
  empty <- is.na(dt[[paste0("A", side)]]) & !is.na(dt[[paste0("A", other)]])
  sub <- dt[empty]
  if (nrow(sub) == 0) return(NULL)
  name <- toupper(sub[[paste0("A", side, "Name")]])
  geo <- sub[[paste0("A", side, "Geo")]]
  out <- lapply(names(name_patterns), function(cty) {
    by_name <- !is.na(name) & grepl(name_patterns[[cty]], name)
    by_geo <- !is.na(geo) & geo == fips[[cty]]
    hit <- by_name | by_geo
    if (!any(hit)) return(NULL)
    data.table(
      country = cty, side = side,
      criterion = fifelse(by_name[hit] & by_geo[hit], "both", fifelse(by_name[hit], "name_only", "geo_only")),
      sub[hit, .(GLOBALEVENTID, SQLDATE, EventCode,
                 actor_code = get(paste0("A", side, "C")), actor_name = get(paste0("A", side, "Name")),
                 actor_geo = get(paste0("A", side, "Geo")),
                 other_country = get(paste0("A", other)))]
    )
  })
  rbindlist(out)
}

scan_one <- function(f) {
  dt <- fread(
    cmd = sprintf("unzip -p %s", shQuote(file.path(raw_dir, f))),
    sep = "\t", quote = "", header = FALSE, na.strings = "",
    select = unname(cols), col.names = names(cols),
    colClasses = list(character = unname(cols)),
    tmpdir = tmp_dir, showProgress = FALSE, nThread = 2L
  )
  dt <- dt[
    (A1T1 %chin% types | A1T2 %chin% types | A1T3 %chin% types) &
      (A2T1 %chin% types | A2T2 %chin% types | A2T3 %chin% types) &
      SQLDATE >= "19910101" & SQLDATE <= "20251231" &
      (is.na(A1) != is.na(A2))
  ]
  cand <- rbind(candidates(dt, 1L), candidates(dt, 2L))
  if (is.null(cand) || nrow(cand) == 0) return(NULL)
  cand[, decade := paste0(substr(SQLDATE, 1, 3), "0s")]
  counts <- cand[, .N, by = .(country, decade, criterion)]
  # a per-file reservoir for the random examples
  set.seed(as.integer(substr(gsub("\\D", "", f), 1, 8)) %% 100000L)
  # weight = rows represented by each kept row, so the final draw is a
  # uniform sample over all candidate rows, not over files
  cand[, n_grp := .N, by = country]
  examples <- cand[, .SD[sample(.N, min(.N, 200L))], by = country]
  examples[, weight := n_grp / pmin(n_grp, 200L)][, n_grp := NULL]
  examples[, file := f]
  list(counts = counts, examples = examples)
}

t0 <- Sys.time()
res <- parallel::mclapply(files, function(f) {
  tryCatch(scan_one(f), error = function(e) {
    message("FAILED ", f, ": ", conditionMessage(e))
    structure(list(), class = "scan_failed")
  })
}, mc.cores = n_cores, mc.preschedule = FALSE)
failed <- vapply(res, inherits, logical(1), "scan_failed")
message(sprintf("Scanned %d files (%d failed) in %.1f min", length(files), sum(failed),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
res <- Filter(function(r) !is.null(r) && !inherits(r, "scan_failed"), res)

counts <- rbindlist(lapply(res, `[[`, "counts"))[, .(n = sum(N)), by = .(country, decade, criterion)]
setorder(counts, country, decade, criterion)
fwrite(counts, file.path(out_dir, "absent_states_counts.csv"))

pool <- rbindlist(lapply(res, `[[`, "examples"))
set.seed(20261105)
examples <- pool[, .SD[sample(.N, min(.N, 20L), prob = weight)], by = country][, weight := NULL]
setorder(examples, country, SQLDATE)
fwrite(examples, file.path(out_dir, "absent_states_examples.csv"))
