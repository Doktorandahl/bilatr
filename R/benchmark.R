# Tools for tuning `threads_per_chain` and `chunk_size` (the reduce_sum
# grainsize) by measuring them: a local sweep (benchmark_bilatr()), the
# grid and profiled model for a sweep run as batch jobs (benchmark_grid(),
# insert_reduce_sum_profile()), and the summaries of either
# (summarize_benchmark(), benchmark_pareto(), plot_benchmark()). See
# vignette("tuning").

#' Build a benchmark grid of `threads_per_chain` x `chunk_size`
#'
#' Constructs the grid of a `reduce_sum` tuning sweep on a panel of
#' `n_dyads` dyads: one row per run, crossing `threads_per_chain`,
#' `chunk_size` and `replicate`. [benchmark_bilatr()] runs such a grid on
#' this machine; for a sweep run as batch jobs, each row can be one job
#' of a job array
#' (see `vignette("tuning", package = "bilatr")`).
#'
#' The models slice the dyads inside `reduce_sum`, so `chunk_size` is a
#' number of dyads per task and its useful range scales with `n_dyads`.
#' When `chunk_size` is `NULL`, it is generated as a geometric sequence of
#' `n_chunk_size` values from `n_dyads * chunk_size_span[1]` to
#' `n_dyads * chunk_size_span[2]` (by default `n_dyads / 200` to
#' `n_dyads / 6`).
#'
#' @param n_dyads Number of dyads (`D` in the Stan data) of the panel
#'   being tuned.
#' @param threads_per_chain Whole numbers, the `threads_per_chain` values
#'   to test.
#' @param chunk_size Whole numbers, the `chunk_size` values to test, or
#'   `NULL` (default) to generate them from `n_dyads`.
#' @param replicates Number of runs per (`threads_per_chain`,
#'   `chunk_size`) cell, to average over timing noise.
#' @param chunk_size_span Length-2 numeric, the low and high multipliers
#'   of `n_dyads` that bound the generated `chunk_size` values. Ignored if
#'   `chunk_size` is given.
#' @param n_chunk_size Number of `chunk_size` values to generate. Ignored
#'   if `chunk_size` is given.
#' @return A tibble with one row per run: `job_id`, `threads_per_chain`,
#'   `chunk_size`, `replicate`, and the derived `n_chunks`
#'   (`ceiling(n_dyads / chunk_size)`) and `chunks_per_thread`. Rows are
#'   ordered by `threads_per_chain`, then `chunk_size`, then `replicate`,
#'   so the rows of one `threads_per_chain` value are a contiguous block
#'   (one job array per thread count).
#' @seealso [benchmark_bilatr()], [summarize_benchmark()]
#' @examples
#' # the grid of a batch-job sweep for a 16,750-dyad panel
#' grid <- benchmark_grid(16750, threads_per_chain = c(8, 16, 32), replicates = 2)
#' grid
#' unique(grid$chunk_size)
#'
#' # explicit values
#' benchmark_grid(300, threads_per_chain = 1:2, chunk_size = c(10, 50), replicates = 1)
#' @export
benchmark_grid <- function(n_dyads,
                           threads_per_chain = c(1, 2, 4, 8, 16),
                           chunk_size = NULL,
                           replicates = 3,
                           chunk_size_span = c(1 / 200, 1 / 6),
                           n_chunk_size = 9) {
  fn <- "benchmark_grid()"
  if (!.is_count(n_dyads) || n_dyads < 2) {
    stop(fn, ": `n_dyads` must be a whole number of at least 2.", call. = FALSE)
  }
  if (!.is_count(replicates) || replicates < 1) {
    stop(fn, ": `replicates` must be a whole number of at least 1.", call. = FALSE)
  }
  .check_whole_vector(threads_per_chain, "threads_per_chain", fn)

  if (is.null(chunk_size)) {
    if (!is.numeric(chunk_size_span) || length(chunk_size_span) != 2 ||
      anyNA(chunk_size_span) || chunk_size_span[1] <= 0 ||
      chunk_size_span[1] >= chunk_size_span[2]) {
      stop(
        fn, ": `chunk_size_span` must be two increasing positive numbers.",
        call. = FALSE
      )
    }
    if (!.is_count(n_chunk_size) || n_chunk_size < 1) {
      stop(fn, ": `n_chunk_size` must be a whole number of at least 1.", call. = FALSE)
    }
    lo <- max(1, round(n_dyads * chunk_size_span[1]))
    hi <- max(lo + 1, round(n_dyads * chunk_size_span[2]))
    chunk_size <- unique(round(exp(seq(log(lo), log(hi), length.out = n_chunk_size))))
  }
  .check_whole_vector(chunk_size, "chunk_size", fn)

  tidyr::expand_grid(
    threads_per_chain = sort(unique(as.integer(threads_per_chain))),
    chunk_size = sort(unique(as.integer(chunk_size))),
    replicate = seq_len(replicates)
  ) %>%
    dplyr::mutate(
      n_chunks = as.integer(ceiling(n_dyads / .data$chunk_size)),
      chunks_per_thread = .data$n_chunks / .data$threads_per_chain,
      job_id = dplyr::row_number()
    ) %>%
    dplyr::relocate("job_id")
}

#' Check a vector of positive whole numbers
#' @param x The vector.
#' @param arg Its argument name, for the error.
#' @param fn The calling function, for the error.
#' @return `TRUE`, invisibly; errors otherwise.
#' @noRd
.check_whole_vector <- function(x, arg, fn) {
  ok <- is.numeric(x) && length(x) >= 1 && !anyNA(x) && all(is.finite(x)) &&
    all(x >= 1) && all(x == round(x))
  if (!ok) {
    stop(fn, ": `", arg, "` must be positive whole numbers.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Wrap a model's `reduce_sum` likelihood in a Stan `profile()` block
#'
#' Returns a copy of a bilatr Stan program in which the `reduce_sum`
#' likelihood statement (with the `dyad_seq` line before it) is wrapped in
#' `profile("<profile_name>") { ... }`. A fit of the returned program
#' reports the time spent in the parallel region through
#' `CmdStanMCMC$profiles()`, which separates the likelihood from the rest
#' of each gradient evaluation. [benchmark_bilatr()] uses it with
#' `profile = TRUE`; for a batch-job sweep, write the result to a `.stan`
#' file and compile that (see `vignette("tuning", package = "bilatr")`).
#'
#' The statement is found by pattern, and the same pattern matches all
#' three registered models (`"stable"`, `"stable_gamma"` and `"ou"`). It
#' errors if the statement is not found exactly once, so a changed model
#' fails loudly rather than being benchmarked unprofiled.
#'
#' @param stan_code Character scalar, the full Stan program, for example
#'   `paste(readLines(f), collapse = "\n")` of a bilatr `.stan` file.
#' @param profile_name Name of the `profile()` block.
#' @return Character scalar, the modified Stan program.
#' @seealso [benchmark_bilatr()]
#' @examples
#' stan_file <- system.file("stan", "bilatr_alphanorm.stan", package = "bilatr")
#' src <- paste(readLines(stan_file), collapse = "\n")
#' profiled <- insert_reduce_sum_profile(src)
#' # the wrapped statement
#' lines <- strsplit(profiled, "\n")[[1]]
#' i <- grep("profile(", lines, fixed = TRUE)
#' writeLines(lines[i:(i + 4)])
#' @export
insert_reduce_sum_profile <- function(stan_code, profile_name = "reduce_sum_likelihood") {
  fn <- "insert_reduce_sum_profile()"
  if (!rlang::is_string(stan_code)) {
    stop(fn, ": `stan_code` must be a single string.", call. = FALSE)
  }
  if (!rlang::is_string(profile_name) || !grepl("^[A-Za-z0-9_.-]+$", profile_name)) {
    stop(fn, ": `profile_name` must be a single name of letters, digits, `_`, `.` or `-`.", call. = FALSE)
  }
  pattern <- "array\\[D\\] int dyad_seq = linspaced_int_array\\(D, 1, D\\);\\s*\\n\\s*target \\+= reduce_sum\\([^;]*\\);"

  loc <- stringr::str_locate_all(stan_code, pattern)[[1]]
  if (nrow(loc) != 1L) {
    stop(
      fn, ": expected exactly one reduce_sum ",
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

#' Maximum tree depth of the runs of benchmark_bilatr()
#'
#' Bounds the work per iteration of a short, unadapted run (see
#' [benchmark_bilatr()], Details).
#' @noRd
.BILATR_BENCHMARK_MAX_TREEDEPTH <- 6L

#' The columns every benchmark result has
#' @noRd
.BILATR_BENCHMARK_COLUMNS <- c(
  "threads_per_chain", "chunk_size", "replicate",
  "sampling_sec", "iter_per_sec", "core_sec_per_1k_iter"
)

#' Summarize benchmark runs to one row per grid cell
#'
#' Collapses the replicates of a benchmark to one row per
#' (`threads_per_chain`, `chunk_size`) cell, within any other grouping
#' columns (for example a `scenario` column that tells data sets apart).
#'
#' @section Benchmark results:
#' A benchmark result is a data frame with one row per run (one chain).
#' [benchmark_bilatr()] returns one; in a batch-job sweep, each job writes
#' one row and the rows are bound together (see
#' `vignette("tuning", package = "bilatr")`). Every
#' benchmark function uses these column names:
#'
#' * `threads_per_chain`: the `threads_per_chain` of the run;
#' * `chunk_size`: the `chunk_size` (`reduce_sum` grainsize) of the run;
#' * `replicate`: the replicate number within the cell;
#' * `sampling_sec`: wall time of the sampling phase, in seconds (from
#'   `fit$time()`);
#' * `iter_per_sec`: sampling iterations per second,
#'   `iter_sampling / sampling_sec`;
#' * `core_sec_per_1k_iter`: core-seconds per 1,000 sampling iterations,
#'   `threads_per_chain * sampling_sec * 1000 / iter_sampling`, the cost
#'   of the throughput;
#' * `reduce_sum_wall_sec` (optional): the wall time spent in the
#'   profiled `reduce_sum` likelihood over the whole run, warmup
#'   included, when the run used a model from
#'   [insert_reduce_sum_profile()];
#' * `leapfrog_per_sec` (optional): leapfrog steps (gradient
#'   evaluations) per second of sampling, a throughput that does not
#'   depend on how long each iteration's trajectory was;
#' * any grouping columns, such as `scenario`.
#'
#' Other columns (for example `warmup_sec`, `n_chunks`) are allowed and
#' ignored.
#'
#' @param results A benchmark result (see Benchmark results), such as the
#'   output of [benchmark_bilatr()].
#' @param group_vars Columns identifying a cell. Those not in `results`
#'   are dropped, so the default works with or without a `scenario`
#'   column.
#' @return A tibble with one row per cell: the `group_vars`, `n_rep`, the
#'   median and MAD of `sampling_sec` and `iter_per_sec`, the median of
#'   `core_sec_per_1k_iter`, and the medians of `reduce_sum_wall_sec` and
#'   `leapfrog_per_sec` (`NA` when `results` has no such column).
#' @seealso [benchmark_pareto()], [plot_benchmark()]
#' @examples
#' # a benchmark result as a batch-job sweep would collect it
#' runs <- benchmark_grid(1000, threads_per_chain = c(2, 4), chunk_size = c(25, 100),
#'   replicates = 2)
#' runs$sampling_sec <- c(40, 41, 38, 39, 24, 25, 26, 26)
#' runs$iter_per_sec <- 100 / runs$sampling_sec
#' runs$core_sec_per_1k_iter <- runs$threads_per_chain * runs$sampling_sec * 1000 / 100
#' summarize_benchmark(runs)
#' @export
summarize_benchmark <- function(results,
                                group_vars = c("scenario", "threads_per_chain", "chunk_size")) {
  .check_benchmark_columns(results, "summarize_benchmark()")
  group_vars <- intersect(group_vars, names(results))
  has_profile <- "reduce_sum_wall_sec" %in% names(results)
  has_leapfrog <- "leapfrog_per_sec" %in% names(results)

  results %>%
    dplyr::ungroup() %>%
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
      leapfrog_per_sec_median = if (has_leapfrog) {
        stats::median(.data$leapfrog_per_sec, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

#' Check that a data frame has the benchmark result columns
#' @param results The data frame.
#' @param fn The calling function, for the error.
#' @return `TRUE`, invisibly; errors otherwise.
#' @noRd
.check_benchmark_columns <- function(results, fn) {
  if (!is.data.frame(results)) {
    stop(fn, ": `results` must be a data frame of benchmark runs.", call. = FALSE)
  }
  absent <- setdiff(.BILATR_BENCHMARK_COLUMNS, names(results))
  if (length(absent) > 0) {
    stop(
      fn, ": `results` lacks the benchmark column(s) ",
      paste0("`", absent, "`", collapse = ", "),
      "; see ?summarize_benchmark, section \"Benchmark results\".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Flag the Pareto-optimal cells of a benchmark
#'
#' Marks each row as on the Pareto front of (higher `throughput`, lower
#' `cost`), or dominated by another row: a cell is on the front when no
#' other cell is at least as fast and at least as cheap, and strictly
#' better in one of the two. With `by`, the front is computed within each
#' group (for example per `scenario`).
#'
#' The front lists the settings worth choosing between; which one to use
#' depends on whether wall time or core time is scarcer.
#'
#' @param df A benchmark summary from [summarize_benchmark()], or a
#'   benchmark result itself (see `?summarize_benchmark`, "Benchmark
#'   results"), such as the output of [benchmark_bilatr()].
#' @param throughput Name of the throughput column (higher is better).
#'   `NULL` (default) uses `"iter_per_sec_median"` if `df` has it, else
#'   `"iter_per_sec"`.
#' @param cost Name of the cost column (lower is better). `NULL` (default)
#'   uses `"core_sec_per_1k_iter_median"` if `df` has it, else
#'   `"core_sec_per_1k_iter"`.
#' @param by Optional character vector of grouping columns.
#' @return `df` with an added logical `pareto` column (`NA` for rows whose
#'   throughput or cost is missing).
#' @seealso [summarize_benchmark()], [plot_benchmark()]
#' @examples
#' runs <- benchmark_grid(1000, threads_per_chain = c(2, 4), chunk_size = c(25, 100),
#'   replicates = 1)
#' runs$sampling_sec <- c(40, 38, 25, 26)
#' runs$iter_per_sec <- 100 / runs$sampling_sec
#' runs$core_sec_per_1k_iter <- runs$threads_per_chain * runs$sampling_sec * 1000 / 100
#' benchmark_pareto(summarize_benchmark(runs))
#' @export
benchmark_pareto <- function(df, throughput = NULL, cost = NULL, by = NULL) {
  fn <- "benchmark_pareto()"
  if (!is.data.frame(df)) {
    stop(fn, ": `df` must be a data frame.", call. = FALSE)
  }
  pick <- function(arg, summary_col, raw_col) {
    if (!is.null(arg)) {
      return(arg)
    }
    if (summary_col %in% names(df)) summary_col else raw_col
  }
  throughput <- pick(throughput, "iter_per_sec_median", "iter_per_sec")
  cost <- pick(cost, "core_sec_per_1k_iter_median", "core_sec_per_1k_iter")
  absent <- setdiff(c(throughput, cost, by), names(df))
  if (length(absent) > 0) {
    stop(fn, ": `df` has no column(s) ", paste0("`", absent, "`", collapse = ", "), ".", call. = FALSE)
  }

  mark <- function(sub) {
    tp <- sub[[throughput]]
    co <- sub[[cost]]
    ok <- is.finite(tp) & is.finite(co)
    dominated <- rep(NA, nrow(sub))
    dominated[ok] <- vapply(which(ok), function(i) {
      any(tp[ok] >= tp[i] & co[ok] <= co[i] & (tp[ok] > tp[i] | co[ok] < co[i]))
    }, logical(1))
    sub$pareto <- !dominated
    sub
  }

  if (is.null(by)) {
    return(mark(df))
  }
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
    dplyr::group_modify(~ mark(.x)) %>%
    dplyr::ungroup()
}

#' Heatmap of a benchmark over `threads_per_chain` x `chunk_size`
#'
#' One tile per (`threads_per_chain`, `chunk_size`) cell, filled by
#' `metric`. Cells on the Pareto front (a `pareto` column, from
#' [benchmark_pareto()]) are outlined.
#'
#' @param summary_df A benchmark summary from [summarize_benchmark()],
#'   optionally passed through [benchmark_pareto()]. A
#'   `bilatr_benchmark` result from [benchmark_bilatr()] is summarized
#'   (and its Pareto front marked) first.
#' @param metric Column to map to the fill. Default
#'   `"iter_per_sec_median"`.
#' @param facet Column to facet by, or `NULL`. The default facets by
#'   `scenario` when there is such a column.
#' @return A `ggplot` object.
#' @seealso [summarize_benchmark()], [benchmark_pareto()]
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   runs <- benchmark_grid(1000, threads_per_chain = c(2, 4), chunk_size = c(25, 100),
#'     replicates = 1)
#'   runs$sampling_sec <- c(40, 38, 25, 26)
#'   runs$iter_per_sec <- 100 / runs$sampling_sec
#'   runs$core_sec_per_1k_iter <- runs$threads_per_chain * runs$sampling_sec * 1000 / 100
#'   plot_benchmark(benchmark_pareto(summarize_benchmark(runs)))
#' }
#' @export
plot_benchmark <- function(summary_df,
                           metric = "iter_per_sec_median",
                           facet = "scenario") {
  rlang::check_installed("ggplot2", "for plot_benchmark()")
  if (inherits(summary_df, "bilatr_benchmark")) {
    summary_df <- benchmark_pareto(summarize_benchmark(summary_df))
  }
  absent <- setdiff(c("threads_per_chain", "chunk_size", metric), names(summary_df))
  if (length(absent) > 0) {
    stop(
      "plot_benchmark(): `summary_df` has no column(s) ",
      paste0("`", absent, "`", collapse = ", "), ".",
      call. = FALSE
    )
  }

  p <- ggplot2::ggplot(summary_df, ggplot2::aes(
    x = factor(.data$chunk_size),
    y = factor(.data$threads_per_chain),
    fill = .data[[metric]]
  )) +
    ggplot2::geom_tile(colour = "grey92") +
    ggplot2::scale_fill_viridis_c(name = metric) +
    ggplot2::labs(
      x = "chunk_size (dyads per reduce_sum task)",
      y = "threads_per_chain",
      title = "reduce_sum tuning sweep"
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

#' Number of cores this process may use
#'
#' `parallelly::availableCores()` when parallelly is installed (it
#' respects the cores a job scheduler or container allocates), else
#' `parallel::detectCores()`.
#' @return An integer, or `NA` if it cannot be determined.
#' @noRd
.bilatr_available_cores <- function() {
  n <- if (rlang::is_installed("parallelly")) {
    parallelly::availableCores()
  } else {
    parallel::detectCores()
  }
  as.integer(unname(n))
}

#' Benchmark `threads_per_chain` and `chunk_size` on this machine
#'
#' Runs a small sweep of `threads_per_chain` x `chunk_size` on your own
#' data and machine, so you can choose the settings of [fit_panel()] by
#' measurement. Each run is one chain of the model, with the given
#' `threads_per_chain` and `chunk_size`; the result goes straight into
#' [summarize_benchmark()], [benchmark_pareto()] and [plot_benchmark()].
#' See `vignette("tuning", package = "bilatr")` for how to read it.
#'
#' Before compiling anything, `benchmark_bilatr()` checks that no
#' `threads_per_chain` value exceeds the cores available
#' (`parallelly::availableCores()` if parallelly is installed, else
#' `parallel::detectCores()`). It then compiles the model (the same cache
#' as [compile_bilatr_model()]), times a short pilot run, and reports
#' (with a message) the number of runs and an estimate of the total time,
#' before starting the sweep.
#'
#' Every run uses the same `seed` within a replicate (replicate `r` uses
#' `seed + r - 1`), so the runs of a replicate sample the same posterior
#' from the same initial values and differ only in how the work is
#' split. Every run also uses `max_treedepth = 6` (at most 63 leapfrog
#' steps per iteration): with a short warmup the step size is not yet
#' adapted, and at the default depth of 10 most iterations would take
#' 1,023 steps, making the runs slow and their work per iteration
#' uneven. Even so, a run's trajectory can drift from the others' (the
#' likelihood is summed in a different order for each `chunk_size`), and
#' then it does a different number of leapfrog steps per iteration; a
#' message names any run more than 10% off the median, whose
#' `iter_per_sec` is then not comparable. Replicates guard against this,
#' since the summaries take medians. Runs are short: their timings say
#' which settings are faster, not how long a full fit takes.
#'
#' With `profile = TRUE`, the model is a copy whose `reduce_sum`
#' likelihood is wrapped in a Stan `profile()` block (see
#' [insert_reduce_sum_profile()]). It is compiled under its own name in
#' the model cache (`bilatr_alphanorm_profile.stan` for `"stable"`, and
#' so on), next to the unprofiled model. Profiling adds a little
#' overhead to every gradient evaluation.
#'
#' @param data A dyad-period data frame (see `?bilatr_dyad_period_data`)
#'   or a Stan data list from [prepare_stan_data()], as for [fit_panel()].
#' @param threads_per_chain Whole numbers, the `threads_per_chain` values
#'   to test.
#' @param chunk_size Whole numbers, the `chunk_size` values to test, or
#'   `NULL` (default) for 5 values spaced geometrically over
#'   [benchmark_grid()]'s default span (`D / 200` to `D / 6`).
#' @param replicates Number of runs per cell.
#' @param iter_warmup,iter_sampling Warmup and sampling iterations of each
#'   run.
#' @param stan_model Which model to benchmark: `"stable"` (default),
#'   `"ou"` or `"stable_gamma"`.
#' @param profile Whether to also record the time spent in the
#'   `reduce_sum` likelihood (`reduce_sum_wall_sec`); see Details.
#' @param ... Passed to [prepare_stan_data()] when `data` is a data frame
#'   (for example `reference_category`). Not allowed when `data` is a
#'   prepared Stan data list.
#' @param seed Base seed of the runs; see Details.
#' @return A tibble of class `bilatr_benchmark`, one row per run, with the
#'   columns described in `?summarize_benchmark` ("Benchmark results"):
#'   `job_id`, `threads_per_chain`, `chunk_size`, `replicate`, `n_chunks`,
#'   `chunks_per_thread`, `warmup_sec`, `sampling_sec`, `total_sec`,
#'   `iter_per_sec`, `core_sec_per_1k_iter`, `n_leapfrog` (leapfrog steps
#'   in the sampling phase), `leapfrog_per_sec`, and, with
#'   `profile = TRUE`, `reduce_sum_wall_sec`. `attr(x, "settings")` records the model, the
#'   data dimensions, the iterations, the seed and the cores available.
#'   It has `print()` and [autoplot()] methods.
#' @seealso [summarize_benchmark()], [benchmark_pareto()],
#'   [plot_benchmark()], [benchmark_grid()] for a batch-job sweep.
#' @examplesIf bilatr_cmdstan_available()
#' # 30 simulated dyads over 10 years, aggregated by PentaClass
#' set.seed(1)
#' countries <- c("USA", "CHN", "RUS", "IND", "BRA", "ZAF")
#' days <- seq(as.Date("2011-01-01"), as.Date("2020-12-31"), by = "day")
#' pairs <- expand.grid(a1 = countries, a2 = countries, stringsAsFactors = FALSE)
#' pairs <- pairs[pairs$a1 != pairs$a2, ]
#' idx <- sample(nrow(pairs), 3000, replace = TRUE)
#' events <- recode_cameo(tibble::tibble(
#'   Actor1CountryCode = pairs$a1[idx],
#'   Actor2CountryCode = pairs$a2[idx],
#'   SQLDATE = as.integer(format(sample(days, 3000, replace = TRUE), "%Y%m%d")),
#'   EventCode = sample(cameo_lookup$CAMEOEVENTCODE, 3000, replace = TRUE)
#' ))
#' agg <- aggregate_events(events, grouping_var = "PentaClass",
#'   start = "2011-01-01", end = "2020-12-31")
#' \donttest{
#' bench <- benchmark_bilatr(
#'   agg,
#'   threads_per_chain = 1:2, chunk_size = 10,
#'   iter_warmup = 30, iter_sampling = 30,
#'   reference_category = "1"
#' )
#' bench
#' benchmark_pareto(summarize_benchmark(bench))
#' }
#' @export
benchmark_bilatr <- function(data,
                             threads_per_chain = c(1, 2, 4),
                             chunk_size = NULL,
                             replicates = 1,
                             iter_warmup = 100,
                             iter_sampling = 100,
                             stan_model = "stable",
                             profile = TRUE,
                             ...,
                             seed = 1) {
  fn <- "benchmark_bilatr()"

  # -- arguments, cores: everything that can fail before compiling --
  .check_whole_vector(threads_per_chain, "threads_per_chain", fn)
  if (!is.null(chunk_size)) .check_whole_vector(chunk_size, "chunk_size", fn)
  for (arg in c("replicates", "iter_warmup", "iter_sampling")) {
    value <- get(arg)
    if (!.is_count(value) || value < 1) {
      stop(fn, ": `", arg, "` must be a whole number of at least 1.", call. = FALSE)
    }
  }
  if (!rlang::is_bool(profile)) {
    stop(fn, ": `profile` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!.is_count(seed)) {
    stop(fn, ": `seed` must be a whole number.", call. = FALSE)
  }
  stan_model <- .canonical_stan_model(stan_model)

  n_cores <- .bilatr_available_cores()
  if (!is.na(n_cores) && max(threads_per_chain) > n_cores) {
    stop(
      fn, ": `threads_per_chain` includes ", max(threads_per_chain), ", but only ",
      n_cores, " core(s) are available",
      if (rlang::is_installed("parallelly")) " (parallelly::availableCores())" else " (parallel::detectCores())",
      ". Use values up to ", n_cores, ".",
      call. = FALSE
    )
  }

  # not list(...): that would evaluate tidyselect expressions such as
  # `action_cols = all_of(x)` outside prepare_stan_data()
  dot_names <- names(rlang::enquos(...))
  if (.is_stan_data_list(data)) {
    if (length(dot_names) > 0) {
      stop(
        fn, ": `data` is an already prepared stan_data list, so the argument(s) ",
        paste0("`", dot_names, "`", collapse = ", "),
        " cannot be passed to prepare_stan_data().",
        call. = FALSE
      )
    }
    stan_data <- data
  } else if (is.data.frame(data)) {
    stan_data <- prepare_stan_data(data, ...)
  } else {
    stop(
      fn, ": `data` must be a dyad-period data frame or a stan_data list from ",
      "prepare_stan_data().",
      call. = FALSE
    )
  }
  if (stan_data$D < 2) {
    stop(fn, ": benchmarking needs a panel (D >= 2); got D = ", stan_data$D, ".", call. = FALSE)
  }
  if (.bilatr_model_has_gamma(stan_model) && !all(c("n_countries", "ctry_a", "ctry_b", "w_send") %in% names(stan_data))) {
    stop(fn, ": stan_model = \"stable_gamma\" needs actor columns in `data`.", call. = FALSE)
  }

  grid <- if (is.null(chunk_size)) {
    benchmark_grid(stan_data$D, threads_per_chain, replicates = replicates, n_chunk_size = 5)
  } else {
    benchmark_grid(stan_data$D, threads_per_chain, chunk_size = chunk_size, replicates = replicates)
  }

  # -- compile, pilot, estimate --
  .require_cmdstan(fn)
  .bilatr_warn_if_experimental(stan_model)
  mod <- .compile_stan_model(stan_model, opt_level = 3, profile = profile)

  run_one <- function(threads, chunk, rep, warmup, sampling) {
    sd <- stan_data
    sd$C <- as.integer(chunk)
    run_seed <- seed + rep - 1L
    out_dir <- tempfile("bilatr_benchmark_")
    dir.create(out_dir)
    on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
    fit <- mod$sample(
      data = sd,
      chains = 1,
      parallel_chains = 1,
      threads_per_chain = threads,
      iter_warmup = warmup,
      iter_sampling = sampling,
      seed = run_seed,
      init = .bilatr_draw_inits(sd, stan_model, 1, run_seed),
      max_treedepth = .BILATR_BENCHMARK_MAX_TREEDEPTH,
      refresh = 0,
      show_messages = FALSE,
      show_exceptions = FALSE,
      diagnostics = NULL,
      output_dir = out_dir
    )
    times <- fit$time()$chains
    leapfrog <- posterior::extract_variable_matrix(fit$sampler_diagnostics(), "n_leapfrog__")
    row <- tibble::tibble(
      warmup_sec = sum(times$warmup),
      sampling_sec = sum(times$sampling),
      total_sec = sum(times$total),
      n_leapfrog = as.integer(sum(leapfrog))
    )
    if (profile) {
      row$reduce_sum_wall_sec <- .bilatr_profile_wall_sec(fit, "reduce_sum_likelihood")
    }
    row
  }

  pilot_iter <- 20L
  pilot <- run_one(1L, max(grid$chunk_size), 1L, pilot_iter, pilot_iter)
  sec_per_iter <- pilot$total_sec / (2 * pilot_iter)
  estimate_sec <- nrow(grid) * sec_per_iter * (iter_warmup + iter_sampling)
  message(
    fn, ": ", nrow(grid), " run(s) of one chain (", iter_warmup, " + ", iter_sampling,
    " iterations) on D = ", stan_data$D, " dyads. At the speed of a 1-thread pilot run, ",
    "that is at most about ", .format_duration(estimate_sec),
    "; more threads should make it faster."
  )

  # -- the sweep --
  rows <- purrr::pmap(
    list(grid$threads_per_chain, grid$chunk_size, grid$replicate),
    function(threads, chunk, rep) run_one(threads, chunk, rep, iter_warmup, iter_sampling)
  )
  out <- dplyr::bind_cols(grid, dplyr::bind_rows(rows)) %>%
    dplyr::mutate(
      iter_per_sec = iter_sampling / .data$sampling_sec,
      core_sec_per_1k_iter = .data$threads_per_chain * .data$sampling_sec * 1000 / iter_sampling,
      .after = "total_sec"
    ) %>%
    dplyr::mutate(leapfrog_per_sec = .data$n_leapfrog / .data$sampling_sec, .after = "n_leapfrog")

  # Same posterior, same seeds: the runs should do the same work per
  # iteration. A run whose trajectory diverged (different chunking sums in
  # a different order) can take far fewer or more leapfrog steps, and then
  # its iter_per_sec measures that, not the setting.
  typical <- stats::median(out$n_leapfrog)
  odd <- out$job_id[abs(out$n_leapfrog / typical - 1) > 0.1]
  if (length(odd) > 0) {
    message(
      fn, ": run(s) ", paste(odd, collapse = ", "), " took more than 10% more or fewer ",
      "leapfrog steps than the median run (", typical, "), so their `iter_per_sec` is not ",
      "comparable. Use replicates (the summaries take medians), or compare `leapfrog_per_sec`."
    )
  }

  attr(out, "settings") <- list(
    stan_model = stan_model, D = stan_data$D, T = stan_data$T, A = stan_data$A,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    max_treedepth = .BILATR_BENCHMARK_MAX_TREEDEPTH, seed = seed,
    profile = profile, available_cores = n_cores
  )
  class(out) <- c("bilatr_benchmark", class(out))
  out
}

#' Wall time of a Stan `profile()` block
#'
#' `CmdStanMCMC$profiles()` reports each profile's `total_time` per
#' thread that ran it; the wall time is their sum over the threads seen,
#' divided by the number of threads.
#'
#' @param fit A `CmdStanMCMC` fit of a profiled model.
#' @param name The profile's name.
#' @return Seconds, or `NA` if the profile is absent.
#' @noRd
.bilatr_profile_wall_sec <- function(fit, name) {
  pr <- tryCatch(fit$profiles(), error = function(e) NULL)
  if (is.null(pr)) {
    return(NA_real_)
  }
  pr <- if (is.data.frame(pr)) pr else dplyr::bind_rows(pr)
  rs <- pr[pr$name == name, , drop = FALSE]
  if (nrow(rs) == 0) {
    return(NA_real_)
  }
  sum(rs$total_time, na.rm = TRUE) / max(length(unique(rs$thread_id)), 1L)
}

#' Format a duration in seconds for a message
#' @param sec Seconds.
#' @return A string such as `"40 seconds"` or `"12 minutes"`.
#' @noRd
.format_duration <- function(sec) {
  if (!is.finite(sec)) {
    return("an unknown time")
  }
  if (sec < 1.5) {
    return("1 second")
  }
  if (sec < 90) {
    return(paste(round(sec), "seconds"))
  }
  if (sec < 90 * 60) {
    return(paste(round(sec / 60), "minutes"))
  }
  paste(round(sec / 3600, 1), "hours")
}

#' Print a `bilatr_benchmark` object
#'
#' @param x A `bilatr_benchmark` object, as returned by
#'   [benchmark_bilatr()].
#' @param ... Passed on to the tibble print method.
#' @return `x`, invisibly.
#' @examples
#' # a hand-made result with the columns benchmark_bilatr() returns
#' x <- benchmark_grid(300, threads_per_chain = 1:2, chunk_size = c(10, 50), replicates = 1)
#' x$sampling_sec <- c(12, 11, 7, 6.5)
#' x$iter_per_sec <- 100 / x$sampling_sec
#' x$core_sec_per_1k_iter <- x$threads_per_chain * x$sampling_sec * 1000 / 100
#' attr(x, "settings") <- list(stan_model = "stable", D = 300, T = 10, A = 5,
#'   iter_warmup = 100, iter_sampling = 100, seed = 1, profile = FALSE, available_cores = 8)
#' class(x) <- c("bilatr_benchmark", class(x))
#' x
#' @export
print.bilatr_benchmark <- function(x, ...) {
  s <- attr(x, "settings")
  if (!is.null(s)) {
    cat(sprintf(
      "<bilatr_benchmark> %s | D = %s, T = %s, A = %s | %d run(s) of 1 chain, %s + %s iterations | %s core(s) available\n\n",
      s$stan_model, s$D, s$T, s$A, nrow(x), s$iter_warmup, s$iter_sampling,
      format(s$available_cores)
    ))
  }
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

#' Plot a `bilatr_benchmark` object
#'
#' Summarizes the runs with [summarize_benchmark()], marks the Pareto
#' front with [benchmark_pareto()], and draws the heatmap of
#' [plot_benchmark()].
#'
#' @param object A `bilatr_benchmark` object, as returned by
#'   [benchmark_bilatr()].
#' @param metric Column of the summary to map to the fill; see
#'   [plot_benchmark()].
#' @param ... Ignored; present for S3 consistency.
#' @return A `ggplot` object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- benchmark_grid(300, threads_per_chain = 1:2, chunk_size = c(10, 50), replicates = 1)
#'   x$sampling_sec <- c(12, 11, 7, 6.5)
#'   x$iter_per_sec <- 100 / x$sampling_sec
#'   x$core_sec_per_1k_iter <- x$threads_per_chain * x$sampling_sec * 1000 / 100
#'   class(x) <- c("bilatr_benchmark", class(x))
#'   ggplot2::autoplot(x)
#' }
#' @exportS3Method ggplot2::autoplot
autoplot.bilatr_benchmark <- function(object, metric = "iter_per_sec_median", ...) {
  rlang::check_installed("ggplot2", "for autoplot.bilatr_benchmark()")
  plot_benchmark(object, metric = metric)
}
