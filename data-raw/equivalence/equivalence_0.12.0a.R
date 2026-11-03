# Equivalence of the 0.12.0a data path with 0.10.2 on the production inputs.
#
# Builds the Stan data for both production 16cat specs (ModifiedRootCode,
# reference "5", yearly 1991-2025, min_n_events = 1, chunk_size = 1200L,
# compute_theta_filtered = 1; directed and undirected) twice:
#   old: assemble_stan_data() from a 0.10.2 worktree (`4a4ffd5`; 0.11.0
#        only added datasets, so this is the data-interface baseline);
#   new: prepare_stan_data(aggregate_events()) from this tree.
# Each runs in its own callr subprocess (pkgload::load_all()), applies the
# exact preprocessing of runscripts/submit_bilatr_16cat.R with its own
# package's recode_cameo(), and times the Stan-data step alone.
#
# Compares every list field and the event_classes/country_codes/
# grouping_var/dyad_ids attributes with identical(), and writes
# data-raw/equivalence/equivalence_0.12.0a.md.
#
# Run from the package root, with the worktree in place:
#   git worktree add ../bilatr-0.10.2 4a4ffd5
#   Rscript data-raw/equivalence/equivalence_0.12.0a.R
# Pure R, no Stan.

old_tree <- normalizePath("../bilatr-0.10.2", mustWork = TRUE)
new_tree <- normalizePath(".", mustWork = TRUE)
data_path <- normalizePath("runscripts/bilatr_replication/data/gdelt_bilatr.rds", mustWork = TRUE)
report_path <- "data-raw/equivalence/equivalence_0.12.0a.md"
n_reps <- 2L # timing repetitions per path and spec; the minimum is reported

spec <- list(
  grouping_var = "ModifiedRootCode",
  reference_category = "5",
  resolution = "yearly",
  start_year = 1991L,
  end_year = 2025L,
  min_n_events = 1L,
  chunk_size = 1200L,
  compute_theta_filtered = 1
)

# Runs in a fresh R process: loads `tree`, preprocesses exactly as
# submit_bilatr_16cat.R does, then builds the Stan data for directed and
# undirected with `path` ("old" or "new"), n_reps times each.
build_in_subprocess <- function(tree, path, data_path, spec, n_reps) {
  suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = FALSE))

  events <- readRDS(data_path)
  events <- events[, c("Actor1CountryCode", "Actor2CountryCode", "SQLDATE", "EventCode")]
  events <- suppressMessages(recode_cameo(events))[
    , c("Actor1CountryCode", "Actor2CountryCode", "SQLDATE", spec$grouping_var)
  ]
  input_hash <- rlang::hash(events)

  build <- function(directed) {
    if (path == "old") {
      assemble_stan_data(
        events,
        years = spec$start_year:spec$end_year,
        resolution = spec$resolution,
        grouping_var = spec$grouping_var,
        directed = directed,
        reference_category = spec$reference_category,
        min_n_events = spec$min_n_events,
        chunk_size = spec$chunk_size,
        compute_theta_filtered = spec$compute_theta_filtered
      )
    } else {
      agg <- aggregate_events(
        events,
        grouping_var = spec$grouping_var,
        resolution = spec$resolution,
        directed = directed,
        start = paste0(spec$start_year, "-01-01"),
        end = paste0(spec$end_year, "-12-31")
      )
      prepare_stan_data(
        agg,
        action_cols = -c("dyad", "actor1", "actor2", "year", "w_send"),
        reference_category = spec$reference_category,
        min_n_events = spec$min_n_events,
        chunk_size = spec$chunk_size,
        compute_theta_filtered = spec$compute_theta_filtered
      )
    }
  }

  out <- list()
  for (directed in c(TRUE, FALSE)) {
    times <- numeric(n_reps)
    for (r in seq_len(n_reps)) {
      gc()
      t0 <- proc.time()[["elapsed"]]
      sd <- suppressMessages(build(directed))
      times[r] <- proc.time()[["elapsed"]] - t0
    }
    out[[if (directed) "directed" else "undirected"]] <- list(stan_data = sd, seconds = times)
  }
  out$input_hash <- input_hash
  out$version <- as.character(pkgload::pkg_version(tree))
  gc_tab <- gc()
  out$peak_mb <- sum(gc_tab[, ncol(gc_tab)]) # last column: "max used" in Mb
  out
}

run_path <- function(tree, path) {
  message("building with the ", path, " path (", tree, ") ...")
  callr::r(
    build_in_subprocess,
    args = list(tree = tree, path = path, data_path = data_path, spec = spec, n_reps = n_reps),
    show = TRUE
  )
}

old <- run_path(old_tree, "old")
new <- run_path(new_tree, "new")

compare_one <- function(old_sd, new_sd) {
  fields <- union(names(old_sd), names(new_sd))
  field_ok <- vapply(fields, function(f) identical(old_sd[[f]], new_sd[[f]]), logical(1))
  attr_ok <- vapply(
    c("event_classes", "country_codes", "grouping_var"),
    function(a) identical(attr(old_sd, a), attr(new_sd, a)),
    logical(1)
  )
  old_ids <- dplyr::rename(attr(old_sd, "dyad_ids"), actor1 = "actor_a", actor2 = "actor_b")
  ids_ok <- identical(old_ids, attr(new_sd, "dyad_ids"))
  list(
    names_ok = identical(names(old_sd), names(new_sd)),
    fields = field_ok,
    attrs = c(attr_ok, dyad_ids = ids_ok)
  )
}

specs <- c("directed", "undirected")
results <- lapply(specs, function(s) compare_one(old[[s]]$stan_data, new[[s]]$stan_data))
names(results) <- specs
all_ok <- identical(old$input_hash, new$input_hash) &&
  all(vapply(results, function(r) r$names_ok && all(r$fields) && all(r$attrs), logical(1)))

fmt_lgl <- function(x) ifelse(x, "TRUE", "**FALSE**")
lines <- c(
  "# Equivalence of the 0.12.0a data path with 0.10.2 (production 16cat inputs)",
  "",
  sprintf("Generated by `data-raw/equivalence/equivalence_0.12.0a.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  sprintf("- old: `assemble_stan_data()`, bilatr %s, worktree at `4a4ffd5`", old$version),
  sprintf("- new: `prepare_stan_data(aggregate_events())`, bilatr %s, this tree", new$version),
  sprintf(
    "- spec: `%s`, reference `\"%s\"`, %s %d-%d (`start = \"%d-01-01\"`, `end = \"%d-12-31\"`), `min_n_events = %d`, `chunk_size = %dL`, `compute_theta_filtered = %d`",
    spec$grouping_var, spec$reference_category, spec$resolution, spec$start_year, spec$end_year,
    spec$start_year, spec$end_year, spec$min_n_events, spec$chunk_size, spec$compute_theta_filtered
  ),
  sprintf(
    "- preprocessed event table identical in both processes (`rlang::hash()`): %s (`%s`)",
    fmt_lgl(identical(old$input_hash, new$input_hash)), new$input_hash
  ),
  "",
  sprintf("**Overall: %s**", if (all_ok) "every field and attribute identical for both specs" else "NOT IDENTICAL"),
  "",
  "## Dimensions",
  "",
  "| spec | D | T | A | n_countries |",
  "| --- | --- | --- | --- | --- |"
)
for (s in specs) {
  sd <- new[[s]]$stan_data
  lines <- c(lines, sprintf("| %s | %d | %s | %d | %d |", s, sd$D, format(sd$T), sd$A, sd$n_countries))
}

lines <- c(lines, "", "## Per-field `identical()`", "", "| field | directed | undirected |", "| --- | --- | --- |")
for (f in names(results$directed$fields)) {
  lines <- c(lines, sprintf(
    "| `%s` | %s | %s |", f, fmt_lgl(results$directed$fields[[f]]), fmt_lgl(results$undirected$fields[[f]])
  ))
}
lines <- c(lines, sprintf(
  "| (list names, in order) | %s | %s |", fmt_lgl(results$directed$names_ok), fmt_lgl(results$undirected$names_ok)
))
for (a in names(results$directed$attrs)) {
  label <- if (a == "dyad_ids") "attr `dyad_ids` (actor_a/actor_b renamed)" else sprintf("attr `%s`", a)
  lines <- c(lines, sprintf(
    "| %s | %s | %s |", label, fmt_lgl(results$directed$attrs[[a]]), fmt_lgl(results$undirected$attrs[[a]])
  ))
}

lines <- c(lines, "", "## Hashes (`digest::digest()`, md5)", "",
           "| spec | path | `Y` | `is_obs` |", "| --- | --- | --- | --- |")
for (s in specs) {
  for (p in c("old", "new")) {
    sd <- get(p)[[s]]$stan_data
    lines <- c(lines, sprintf(
      "| %s | %s | `%s` | `%s` |", s, p,
      digest::digest(sd$Y, algo = "md5"), digest::digest(sd$is_obs, algo = "md5")
    ))
  }
}

lines <- c(lines, "", sprintf("## Timing (Stan-data step only, elapsed seconds, min of %d runs)", n_reps), "",
           "| spec | old | new | new / old |", "| --- | --- | --- | --- |")
ratios <- numeric(0)
for (s in specs) {
  t_old <- min(old[[s]]$seconds)
  t_new <- min(new[[s]]$seconds)
  ratios[s] <- t_new / t_old
  lines <- c(lines, sprintf("| %s | %.2f | %.2f | %.2f |", s, t_old, t_new, t_new / t_old))
}
lines <- c(lines, "", if (any(ratios > 1.5)) {
  "**FLAG: the new path is more than 1.5x slower than the old one.**"
} else {
  "The new path is not materially slower (no ratio above 1.5)."
}, "", sprintf(
  "Peak R memory (`gc()` max used, Ncells + Vcells) per subprocess: old %.0f MB, new %.0f MB.",
  old$peak_mb, new$peak_mb
))

writeLines(lines, report_path)
message("wrote ", report_path, if (all_ok) " (all identical)" else " (NOT IDENTICAL)")
