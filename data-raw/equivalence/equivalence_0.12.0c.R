# Equivalence of the PORTED replication pipeline (0.12.0c) with 0.10.2.
#
# New side: the ported runscripts/bilatr_replication/config.R and
# R/functions.R, sourced as they are, in a callr subprocess with this tree
# installed into a temporary library (so library(bilatr) and the
# min_bilatr_version check see 0.12.0, not the system library's copy).
# output_root is pointed at a temporary directory and overwrite = TRUE;
# build_stan_data() is called for both production models in FULL mode
# (BILATR_TEST_MODE unset). Both its return value and the file it caches
# are compared.
#
# Old side: assemble_stan_data() from the ../bilatr-0.10.2 worktree
# (pkgload::load_all()), with exactly the preprocessing of the pre-0.12.0
# functions.R (runscripts/legacy/pre_0.12.0/bilatr_replication/R/functions.R)
# and the config.R settings, as in data-raw/equivalence/equivalence_0.12.0a.R.
#
# Checks:
#   * every list field, the list names, and the event_classes/country_codes/
#     grouping_var/dyad_ids attributes identical() (dyad_ids after renaming
#     actor_a/actor_b);
#   * Y/is_obs md5 equal to the hashes recorded in data-raw/equivalence/equivalence_0.12.0a.md;
#   * the cached <model>_stan_data.rds identical() to the returned list;
#   * nothing under runscripts/bilatr_replication/ changed (mtime + size of
#     every file, before vs after; output* directories listed).
# Writes data-raw/equivalence/equivalence_0.12.0c.md.
#
# Run from the package root, with the worktree in place:
#   Rscript data-raw/equivalence/equivalence_0.12.0c.R
# Pure R, no Stan.

new_tree <- normalizePath(".", mustWork = TRUE)
old_tree <- normalizePath("../bilatr-0.10.2", mustWork = TRUE)
repl_dir <- normalizePath("runscripts/bilatr_replication", mustWork = TRUE)
report_path <- "data-raw/equivalence/equivalence_0.12.0c.md"

# Y/is_obs md5 from data-raw/equivalence/equivalence_0.12.0a.md (identical old and new there)
expected_md5 <- list(
  directed = c(Y = "4ed0ffdf392325ea2334744f9ae49f8a", is_obs = "f8a7113a838ce67b0c00fcee99aa216f"),
  undirected = c(Y = "0950349b6be363670e63fdd93291ed4d", is_obs = "1ed92de73b3c72858e41b2d5a9bbd716")
)
models <- c(directed = "bilatr_16cat_directed", undirected = "bilatr_16cat_undirected")

# --- 0. snapshot the replication folder ---------------------------------------
snapshot <- function(dir) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  info <- file.info(f)
  data.frame(file = f, mtime = as.numeric(info$mtime), size = info$size, stringsAsFactors = FALSE)
}
before <- snapshot(repl_dir)
output_dirs <- list.files(repl_dir, pattern = "^output", full.names = TRUE, include.dirs = TRUE)

# --- 1. install this tree into a temporary library ----------------------------
lib <- file.path(tempdir(), "bilatr_lib")
dir.create(lib, showWarnings = FALSE)
message("installing ", new_tree, " into ", lib, " ...")
inst <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--no-docs", "--no-multiarch", "--no-test-load", "-l", shQuote(lib), shQuote(new_tree)),
  stdout = TRUE, stderr = TRUE
)
if (!dir.exists(file.path(lib, "bilatr"))) stop("install failed:\n", paste(inst, collapse = "\n"))

# --- 2. new side: the ported functions.R ---------------------------------------
new_side <- function(repl_dir, out_root, models) {
  Sys.unsetenv("BILATR_TEST_MODE")
  setwd(repl_dir)
  source("config.R", local = FALSE)
  stopifnot(!test_mode)
  output_root <<- out_root
  overwrite <<- TRUE
  source("R/functions.R", local = FALSE)
  check_bilatr_version()
  res <- list(version = as.character(utils::packageVersion("bilatr")),
              lib = dirname(find.package("bilatr")),
              stan_data_path = list(), stan_data = list(), cached = list(), seconds = list())
  for (s in names(models)) {
    m <- models[[s]]
    t0 <- proc.time()[["elapsed"]]
    sd <- suppressMessages(build_stan_data(m))
    res$seconds[[s]] <- proc.time()[["elapsed"]] - t0
    res$stan_data[[s]] <- sd
    res$stan_data_path[[s]] <- stan_data_path(m)
    res$cached[[s]] <- readRDS(stan_data_path(m))
  }
  res
}
out_root <- file.path(tempdir(), "repl_output")
message("new side: build_stan_data() from the ported functions.R ...")
new <- callr::r(
  new_side,
  args = list(repl_dir = repl_dir, out_root = out_root, models = models),
  libpath = c(lib, .libPaths()),
  show = TRUE
)

# --- 3. old side: assemble_stan_data() from 0.10.2 -----------------------------
old_side <- function(old_tree, repl_dir, models) {
  suppressMessages(pkgload::load_all(old_tree, quiet = TRUE, export_all = FALSE))
  Sys.unsetenv("BILATR_TEST_MODE")
  cfg <- new.env()
  sys.source(file.path(repl_dir, "config.R"), envir = cfg)
  stopifnot(!cfg$test_mode)
  # preprocessing of the pre-0.12.0 functions.R build_stan_data()
  cols <- c("Actor1CountryCode", "Actor2CountryCode", "SQLDATE")
  events <- readRDS(file.path(repl_dir, cfg$data_path))[, c(cols, "EventCode")]
  events <- suppressMessages(recode_cameo(events))[, c(cols, cfg$grouping_var)]
  res <- list(version = as.character(pkgload::pkg_version(old_tree)), stan_data = list())
  for (s in names(models)) {
    directed <- cfg$models$directed[cfg$models$model == models[[s]]]
    res$stan_data[[s]] <- suppressMessages(assemble_stan_data(
      events,
      years = cfg$start_year:cfg$end_year,
      resolution = cfg$resolution,
      grouping_var = cfg$grouping_var,
      directed = directed,
      reference_category = cfg$reference_category,
      min_n_events = cfg$min_n_events,
      chunk_size = cfg$grainsize,
      compute_theta_filtered = 1
    ))
  }
  res
}
message("old side: assemble_stan_data() from ", old_tree, " ...")
old <- callr::r(old_side, args = list(old_tree = old_tree, repl_dir = repl_dir, models = models), show = TRUE)

# --- 4. compare -----------------------------------------------------------------
compare_one <- function(old_sd, new_sd) {
  fields <- union(names(old_sd), names(new_sd))
  field_ok <- vapply(fields, function(f) identical(old_sd[[f]], new_sd[[f]]), logical(1))
  attr_ok <- vapply(
    c("event_classes", "country_codes", "grouping_var"),
    function(a) identical(attr(old_sd, a), attr(new_sd, a)),
    logical(1)
  )
  old_ids <- dplyr::rename(attr(old_sd, "dyad_ids"), actor1 = "actor_a", actor2 = "actor_b")
  list(
    names_ok = identical(names(old_sd), names(new_sd)),
    fields = field_ok,
    attrs = c(attr_ok, dyad_ids = identical(old_ids, attr(new_sd, "dyad_ids")))
  )
}
specs <- names(models)
results <- lapply(specs, function(s) compare_one(old$stan_data[[s]], new$stan_data[[s]]))
names(results) <- specs
md5 <- lapply(specs, function(s) {
  sd <- new$stan_data[[s]]
  c(Y = digest::digest(sd$Y, algo = "md5"), is_obs = digest::digest(sd$is_obs, algo = "md5"))
})
names(md5) <- specs
md5_ok <- vapply(specs, function(s) identical(md5[[s]], expected_md5[[s]]), logical(1))
cached_ok <- vapply(specs, function(s) identical(new$cached[[s]], new$stan_data[[s]]), logical(1))

after <- snapshot(repl_dir)
untouched <- identical(before, after)

all_ok <- all(vapply(results, function(r) r$names_ok && all(r$fields) && all(r$attrs), logical(1))) &&
  all(md5_ok) && all(cached_ok) && untouched

# --- 5. report ------------------------------------------------------------------
fmt_lgl <- function(x) ifelse(x, "TRUE", "**FALSE**")
lines <- c(
  "# Equivalence of the ported replication pipeline (0.12.0c) with 0.10.2",
  "",
  sprintf("Generated by `data-raw/equivalence/equivalence_0.12.0c.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  sprintf(paste0(
    "- new: `build_stan_data()` from the ported `runscripts/bilatr_replication/R/functions.R` ",
    "(`aggregate_events()` -> `prepare_stan_data()`), bilatr %s installed from this tree into a ",
    "temporary library; `config.R` sourced unchanged, full mode, `output_root` = a temporary ",
    "directory, `overwrite = TRUE`"
  ), new$version),
  sprintf(paste0(
    "- old: `assemble_stan_data()`, bilatr %s (`../bilatr-0.10.2` worktree at `4a4ffd5`), with the ",
    "preprocessing of the pre-0.12.0 `functions.R` and the same `config.R` settings"
  ), old$version),
  "- spec (`config.R`): `ModifiedRootCode`, reference `\"5\"`, yearly 1991-2025, `min_n_events = 1`, `grainsize = 1200L`, `compute_theta_filtered = 1`",
  "",
  sprintf("**Overall: %s**", if (all_ok) "every field and attribute identical for both models" else "NOT IDENTICAL"),
  "",
  "## Dimensions (new)",
  "",
  "| model | D | T | A | n_countries | `build_stan_data()` seconds |",
  "| --- | --- | --- | --- | --- | --- |"
)
for (s in specs) {
  sd <- new$stan_data[[s]]
  lines <- c(lines, sprintf("| %s | %d | %s | %d | %d | %.1f |", models[[s]], sd$D, format(sd$T), sd$A,
                            sd$n_countries, new$seconds[[s]]))
}
lines <- c(lines, "", "## Per-field `identical()`, new vs old", "",
           "| field | directed | undirected |", "| --- | --- | --- |")
for (f in names(results$directed$fields)) {
  lines <- c(lines, sprintf("| `%s` | %s | %s |", f,
                            fmt_lgl(results$directed$fields[[f]]), fmt_lgl(results$undirected$fields[[f]])))
}
lines <- c(lines, sprintf("| (list names, in order) | %s | %s |",
                          fmt_lgl(results$directed$names_ok), fmt_lgl(results$undirected$names_ok)))
for (a in names(results$directed$attrs)) {
  label <- if (a == "dyad_ids") "attr `dyad_ids` (actor_a/actor_b renamed)" else sprintf("attr `%s`", a)
  lines <- c(lines, sprintf("| %s | %s | %s |", label,
                            fmt_lgl(results$directed$attrs[[a]]), fmt_lgl(results$undirected$attrs[[a]])))
}
lines <- c(lines, sprintf(
  "| cached `<model>_stan_data.rds` `identical()` to the returned list | %s | %s |",
  fmt_lgl(cached_ok[["directed"]]), fmt_lgl(cached_ok[["undirected"]])
))
lines <- c(lines, "", "## md5 (`digest::digest()`) vs `data-raw/equivalence/equivalence_0.12.0a.md`", "",
           "| model | `Y` | expected | `is_obs` | expected | match |", "| --- | --- | --- | --- | --- | --- |")
for (s in specs) {
  lines <- c(lines, sprintf("| %s | `%s` | `%s` | `%s` | `%s` | %s |", s,
                            md5[[s]][["Y"]], expected_md5[[s]][["Y"]],
                            md5[[s]][["is_obs"]], expected_md5[[s]][["is_obs"]], fmt_lgl(md5_ok[[s]])))
}
lines <- c(lines, "", "## Real outputs untouched", "",
           sprintf("- `runscripts/bilatr_replication/output*` directories present: %s",
                   if (length(output_dirs)) paste0("`", basename(output_dirs), "`", collapse = ", ") else "none"),
           sprintf("- every file under `runscripts/bilatr_replication/` (%d files) has the same mtime and size before and after: %s",
                   nrow(before), fmt_lgl(untouched)),
           sprintf("- the run wrote only to `%s`", out_root))

writeLines(lines, report_path)
message("wrote ", report_path, if (all_ok) " (all identical)" else " (NOT IDENTICAL)")
