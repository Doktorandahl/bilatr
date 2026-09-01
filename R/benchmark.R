# Helpers for the reduce_sum threading benchmark sweep. These are internal
# (not exported) utilities shared by the one-off SLURM scripts under
# runscripts/benchmark/ and covered by tests/testthat/test_benchmark.R.
# The scripts themselves (grid prep, worker, submission wrapper,
# aggregation) live outside the package; only the reusable pure functions
# are here.

#' Build a (cores x grainsize x replicate) benchmark grid
#'
#' Constructs the parameter grid for a `reduce_sum` threading sweep on a
#' panel of `n_dyads` dyads. The Stan model slices `dyad_seq` (length
#' `D == n_dyads`) inside `reduce_sum`, so `grainsize` is a count of
#' dyads per task and the useful range scales with `n_dyads`: when
#' `grainsize_levels` is not supplied it is generated as a geometric
#' sequence spanning `grainsize_span * n_dyads` (default: `n_dyads / 200`
#' up to `n_dyads / 6`), which brackets the previously Pareto-optimal
#' grainsize of ~600 for `n_dyads` ~ 16750.
#'
#' @param n_dyads Integer, number of dyads (`D`) in the assembled Stan
#'   data for the scenario being benchmarked.
#' @param thread_levels Integer vector of `threads_per_chain` values to
#'   test. Default `c(4, 8, 12, 16, 24, 32, 48)`.
#' @param grainsize_levels Optional integer vector of `reduce_sum`
#'   grainsizes. If `NULL` (default), generated from `n_dyads` and
#'   `grainsize_span` / `n_grainsize`.
#' @param replicates Number of repeat runs per (cores, grainsize) cell,
#'   to average over cluster wall-clock noise. Default 3.
#' @param grainsize_span Length-2 numeric, the low/high multipliers of
#'   `n_dyads` bounding the generated grainsize sequence. Ignored if
#'   `grainsize_levels` is supplied.
#' @param n_grainsize Number of grainsize levels to generate. Ignored if
#'   `grainsize_levels` is supplied.
#' @return A tibble with one row per job: `job_id`, `cores`, `grainsize`,
#'   `replicate`, and the derived `n_chunks` (`ceiling(n_dyads /
#'   grainsize)`) and `chunks_per_core`. Rows are ordered by `cores`,
#'   then `grainsize`, then `replicate`, so a SLURM array filtered to one
#'   `cores` level indexes rows `1..(n_grainsize * replicates)` in that
#'   order.
#' @keywords internal
build_reduce_sum_grid <- function(n_dyads,
                                  thread_levels = c(4, 8, 12, 16, 24, 32, 48),
                                  grainsize_levels = NULL,
                                  replicates = 3L,
                                  grainsize_span = c(1 / 200, 1 / 6),
                                  n_grainsize = 9L) {
  stopifnot(
    length(n_dyads) == 1L, is.finite(n_dyads), n_dyads >= 2,
    length(replicates) == 1L, replicates >= 1,
    length(thread_levels) >= 1L, all(thread_levels >= 1)
  )

  if (is.null(grainsize_levels)) {
    stopifnot(length(grainsize_span) == 2L, grainsize_span[1] < grainsize_span[2])
    lo <- max(1, round(n_dyads * grainsize_span[1]))
    hi <- max(lo + 1, round(n_dyads * grainsize_span[2]))
    grainsize_levels <- unique(round(
      exp(seq(log(lo), log(hi), length.out = n_grainsize))
    ))
  }
  grainsize_levels <- sort(unique(as.integer(grainsize_levels)))
  stopifnot(all(grainsize_levels >= 1))

  tidyr::expand_grid(
    cores = sort(as.integer(thread_levels)),
    grainsize = grainsize_levels,
    replicate = seq_len(replicates)
  ) %>%
    dplyr::mutate(
      n_chunks = as.integer(ceiling(n_dyads / .data$grainsize)),
      chunks_per_core = .data$n_chunks / .data$cores,
      job_id = dplyr::row_number()
    ) %>%
    dplyr::relocate("job_id")
}

#' Wrap the reduce_sum likelihood in a Stan `profile()` block
#'
#' Takes the source of `inst/stan/bilatr_dirmult_irt.stan` (or an
#' equivalent model whose `model` block contains the exact
#' `array[D] int dyad_seq = linspaced_int_array(D, 1, D);` +
#' `target += reduce_sum(...);` pair) and returns a copy with that pair
#' wrapped in `profile("<profile_name>") { ... }`, so a fit of the
#' returned model reports per-call timing for the parallel region via
#' `CmdStanMCMC$profiles()`. Errors if the statement pair is not found
#' exactly once, so a drifted upstream model fails loudly rather than
#' benchmarking an unprofiled program.
#'
#' @param stan_code Character scalar, the full Stan program source.
#' @param profile_name Name for the `profile()` block. Default
#'   `"reduce_sum_likelihood"`.
#' @return Character scalar, the modified Stan program.
#' @keywords internal
insert_reduce_sum_profile <- function(stan_code, profile_name = "reduce_sum_likelihood") {
  stopifnot(length(stan_code) == 1L, is.character(stan_code))
  pattern <- "array\\[D\\] int dyad_seq = linspaced_int_array\\(D, 1, D\\);\\s*\\n\\s*target \\+= reduce_sum\\([^;]*\\);"

  loc <- stringr::str_locate_all(stan_code, pattern)[[1]]
  if (nrow(loc) != 1L) {
    stop(
      "insert_reduce_sum_profile(): expected exactly one reduce_sum ",
      "likelihood statement to wrap, found ", nrow(loc),
      ". The upstream model may have changed; update the pattern.",
      call. = FALSE
    )
  }

  stmt <- substr(stan_code, loc[1, "start"], loc[1, "end"])
  wrapped <- paste0(
    "profile(\"", profile_name, "\") {\n    ",
    gsub("\n", "\n  ", stmt, fixed = TRUE),
    "\n  }"
  )
  paste0(
    substr(stan_code, 1L, loc[1, "start"] - 1L),
    wrapped,
    substr(stan_code, loc[1, "end"] + 1L, nchar(stan_code))
  )
}

#' Summarise raw per-job benchmark results to one row per grid cell
#'
#' @param results A data frame of per-job benchmark rows, as written by
#'   `runscripts/benchmark/bench_worker.R` (one row per job). Must contain
#'   the grouping columns and the timing columns `sampling_sec`,
#'   `iter_per_sec`, `core_sec_per_1k_iter`, and optionally
#'   `reduce_sum_wall_sec`.
#' @param group_vars Columns identifying a grid cell. Defaults to
#'   `c("scenario", "cores", "grainsize")`; any not present in `results`
#'   are dropped.
#' @return A tibble with one row per cell: replicate count and the median
#'   (plus MAD for the headline metric) of each timing column.
#' @keywords internal
summarise_benchmark <- function(results,
                                group_vars = c("scenario", "cores", "grainsize")) {
  group_vars <- intersect(group_vars, names(results))
  has_profile <- "reduce_sum_wall_sec" %in% names(results)

  results %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      sampling_sec_median = stats::median(.data$sampling_sec, na.rm = TRUE),
      sampling_sec_mad = stats::mad(.data$sampling_sec, na.rm = TRUE),
      iter_per_sec_median = stats::median(.data$iter_per_sec, na.rm = TRUE),
      iter_per_sec_mad = stats::mad(.data$iter_per_sec, na.rm = TRUE),
      core_sec_per_1k_iter_median = stats::median(.data$core_sec_per_1k_iter, na.rm = TRUE),
      reduce_sum_wall_sec_median = if (has_profile) {
        stats::median(.data$reduce_sum_wall_sec, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

#' Flag the Pareto-optimal cells of a benchmark summary
#'
#' Marks each row as on the Pareto frontier of (higher `throughput`,
#' lower `cost`) or dominated by another row. With `by` supplied, the
#' frontier is computed within each group (e.g. per scenario).
#'
#' @param df A data frame (e.g. from [summarise_benchmark()]).
#' @param throughput Name of the throughput column (higher is better),
#'   e.g. `"iter_per_sec_median"`.
#' @param cost Name of the cost column (lower is better), e.g.
#'   `"core_sec_per_1k_iter_median"`.
#' @param by Optional character vector of grouping columns; the frontier
#'   is computed separately within each group.
#' @return `df` with an added logical `pareto` column.
#' @keywords internal
pareto_front <- function(df, throughput, cost, by = NULL) {
  stopifnot(throughput %in% names(df), cost %in% names(df))

  mark <- function(sub) {
    tp <- sub[[throughput]]
    co <- sub[[cost]]
    ok <- is.finite(tp) & is.finite(co)
    dominated <- rep(NA, nrow(sub))
    dominated[ok] <- vapply(which(ok), function(i) {
      any(tp[ok] >= tp[i] & co[ok] <= co[i] & (tp[ok] > tp[i] | co[ok] < co[i]))
    }, logical(1))
    dplyr::mutate(sub, pareto = !dominated)
  }

  if (is.null(by)) {
    return(mark(df))
  }
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::group_modify(~ mark(.x)) %>%
    dplyr::ungroup()
}

#' Heatmap of benchmark throughput over the cores x grainsize grid
#'
#' @param summary_df A benchmark summary (from [summarise_benchmark()]),
#'   optionally already passed through [pareto_front()] (if it has a
#'   `pareto` column, frontier cells are outlined).
#' @param metric Column to map to fill. Default `"iter_per_sec_median"`.
#' @param facet Column to facet by, or `NULL`. Default `"scenario"` when
#'   present.
#' @return A `ggplot` object.
#' @keywords internal
plot_benchmark <- function(summary_df,
                           metric = "iter_per_sec_median",
                           facet = "scenario") {
  rlang::check_installed("ggplot2", "for plot_benchmark()")
  stopifnot(metric %in% names(summary_df))

  p <- ggplot2::ggplot(summary_df, ggplot2::aes(
    x = factor(.data$grainsize),
    y = factor(.data$cores),
    fill = .data[[metric]]
  )) +
    ggplot2::geom_tile(colour = "grey92") +
    ggplot2::scale_fill_viridis_c(name = metric) +
    ggplot2::labs(
      x = "reduce_sum grainsize (dyads/task)",
      y = "threads per chain (cores)",
      title = "reduce_sum throughput sweep"
    ) +
    ggplot2::theme_minimal()

  if ("pareto" %in% names(summary_df)) {
    p <- p + ggplot2::geom_tile(
      data = function(d) d[which(d$pareto), , drop = FALSE],
      colour = "black", linewidth = 0.9, fill = NA
    )
  }
  if (!is.null(facet) && facet %in% names(summary_df)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[facet]]))
  }
  p
}
