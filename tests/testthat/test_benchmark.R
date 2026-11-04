fake_runs <- function(grid, iter_sampling = 100) {
  grid$sampling_sec <- 10 + grid$replicate * 0.1 + grid$chunk_size / 1000
  grid$iter_per_sec <- iter_sampling / grid$sampling_sec
  grid$core_sec_per_1k_iter <- grid$threads_per_chain * grid$sampling_sec * 1000 / iter_sampling
  grid
}

# -- benchmark_grid() ---------------------------------------------------------

test_that("benchmark_grid produces threads x chunk_size x replicate rows in order", {
  grid <- benchmark_grid(
    n_dyads = 3300,
    threads_per_chain = c(4, 16, 48),
    replicates = 2,
    n_chunk_size = 5
  )

  expect_s3_class(grid, "tbl_df")
  expect_named(grid, c(
    "job_id", "threads_per_chain", "chunk_size", "replicate", "n_chunks", "chunks_per_thread"
  ))
  expect_equal(nrow(grid), 3 * 5 * 2)
  expect_setequal(unique(grid$threads_per_chain), c(4L, 16L, 48L))
  expect_length(unique(grid$chunk_size), 5L)
  expect_equal(grid$job_id, seq_len(nrow(grid)))

  # ordered by threads, then chunk_size, then replicate, so one thread
  # count is a contiguous block (one job array per thread count)
  first_block <- dplyr::filter(grid, threads_per_chain == 4L)
  expect_equal(nrow(first_block), 5 * 2)
  expect_equal(first_block$chunk_size, rep(sort(unique(grid$chunk_size)), each = 2))
  expect_equal(first_block$job_id, 1:10)

  expect_equal(grid$n_chunks, as.integer(ceiling(3300 / grid$chunk_size)))
  expect_equal(grid$chunks_per_thread, grid$n_chunks / grid$threads_per_chain)
})

test_that("benchmark_grid's defaults are 1, 2, 4, 8, 16 threads, 9 chunk sizes, 3 replicates", {
  grid <- benchmark_grid(1000)
  expect_identical(sort(unique(grid$threads_per_chain)), c(1L, 2L, 4L, 8L, 16L))
  expect_length(unique(grid$chunk_size), 9L)
  expect_identical(sort(unique(grid$replicate)), 1:3)
  expect_equal(range(grid$chunk_size), c(5L, 167L)) # 1000 / 200, 1000 / 6
})

test_that("benchmark_grid scales the generated chunk_size range with n_dyads", {
  small <- benchmark_grid(3300, n_chunk_size = 6)
  large <- benchmark_grid(16750, n_chunk_size = 6)

  expect_lt(max(small$chunk_size), max(large$chunk_size))
  # 600, on the Pareto front of the production sweep, is inside the range
  expect_true(min(large$chunk_size) < 600 && max(large$chunk_size) > 600)
  expect_true(all(small$chunk_size >= 1))
  # a small panel still gets at least two distinct values
  expect_gte(length(unique(benchmark_grid(20, n_chunk_size = 5)$chunk_size)), 2)
})

test_that("benchmark_grid honours explicit chunk_size and validates inputs", {
  grid <- benchmark_grid(1000, chunk_size = c(300, 50, 50, 100))
  expect_equal(sort(unique(grid$chunk_size)), c(50L, 100L, 300L))

  expect_error(benchmark_grid(1), "n_dyads")
  expect_error(benchmark_grid(1000, replicates = 0), "replicates")
  expect_error(benchmark_grid(1000, threads_per_chain = 0), "threads_per_chain")
  expect_error(benchmark_grid(1000, threads_per_chain = 1.5), "threads_per_chain")
  expect_error(benchmark_grid(1000, chunk_size = c(10, NA)), "chunk_size")
  expect_error(benchmark_grid(1000, chunk_size_span = c(0.2, 0.1)), "chunk_size_span")
  expect_error(benchmark_grid(1000, n_chunk_size = 0), "n_chunk_size")
})

# -- insert_reduce_sum_profile() ------------------------------------------------

for (model in c("stable", "stable_gamma", "ou")) {
  test_that(paste0("insert_reduce_sum_profile wraps the reduce_sum call exactly once (", model, ")"), {
    src <- paste(readLines(.resolve_stan_model(model)), collapse = "\n")
    out <- insert_reduce_sum_profile(src)

    expect_match(out, 'profile("reduce_sum_likelihood") {', fixed = TRUE)
    # the reduce_sum call and its dyad_seq setup are still present, exactly once
    expect_equal(lengths(regmatches(out, gregexpr("target += reduce_sum(", out, fixed = TRUE))), 1L)
    expect_equal(lengths(regmatches(out, gregexpr("linspaced_int_array(D, 1, D)", out, fixed = TRUE))), 1L)
    # the only added non-whitespace tokens are the profile wrapper + braces
    added <- nchar(gsub("\\s", "", out)) - nchar(gsub("\\s", "", src))
    expect_equal(added, nchar('profile("reduce_sum_likelihood"){}'))
    # the reduce_sum statement sits inside the profile block
    inside <- regmatches(out, regexpr('profile\\("reduce_sum_likelihood"\\) \\{[^}]*\\}', out))
    expect_match(inside, "target += reduce_sum(", fixed = TRUE)
  })
}

test_that("insert_reduce_sum_profile takes a profile name and validates its inputs", {
  src <- paste(readLines(.resolve_stan_model("stable")), collapse = "\n")
  expect_match(insert_reduce_sum_profile(src, "lik"), 'profile("lik") {', fixed = TRUE)
  expect_error(insert_reduce_sum_profile(src, 'a"b'), "profile_name")
  expect_error(insert_reduce_sum_profile(c(src, src)), "single string")
  expect_error(
    insert_reduce_sum_profile("model { target += normal_lpdf(y | 0, 1); }"),
    "found 0"
  )
})

test_that("the profiled model is staged under its own name, and only rewritten on change", {
  dir <- withr::local_tempdir("bilatr models with spaces ")
  plain <- .bilatr_stage_stan_file("stable", dir)
  profiled <- .bilatr_stage_stan_file("stable", dir, profile = TRUE)

  expect_equal(basename(plain), "bilatr_alphanorm.stan")
  expect_equal(basename(profiled), "bilatr_alphanorm_profile.stan")
  expect_identical(
    paste(readLines(profiled), collapse = "\n"),
    insert_reduce_sum_profile(paste(readLines(plain), collapse = "\n"))
  )
  mtime <- file.mtime(profiled)
  Sys.sleep(1.1)
  .bilatr_stage_stan_file("stable", dir, profile = TRUE)
  expect_equal(file.mtime(profiled), mtime)

  expect_equal(
    basename(.bilatr_stage_stan_file("stable_gamma", dir, profile = TRUE)),
    "bilatr_alphanorm_gamma_profile.stan"
  )
})

test_that("bilatr_clear_model_cache() also deletes the profiled copies", {
  dir <- withr::local_tempdir("bilatr_models_")
  withr::local_options(bilatr.model_dir = dir)
  .bilatr_stage_stan_file("stable", dir, profile = TRUE)
  file.create(file.path(dir, c("bilatr_alphanorm_profile", "keep_me.txt")))
  bilatr_clear_model_cache()
  expect_identical(list.files(dir), "keep_me.txt")
})

# -- summarize_benchmark(), benchmark_pareto(), plot_benchmark() ---------------

test_that("summarize_benchmark collapses replicates to one row per cell", {
  raw <- tidyr::expand_grid(
    scenario = "prd_quad",
    threads_per_chain = c(8L, 16L),
    chunk_size = c(60L, 120L),
    replicate = 1:3
  )
  raw$sampling_sec <- 10 + raw$replicate * 0.1
  raw$iter_per_sec <- 75 / raw$sampling_sec
  raw$core_sec_per_1k_iter <- raw$threads_per_chain * raw$sampling_sec * 1000 / 75
  raw$reduce_sum_wall_sec <- raw$sampling_sec * 0.8

  summ <- summarize_benchmark(raw)
  expect_named(summ, c(
    "scenario", "threads_per_chain", "chunk_size", "n_rep",
    "sampling_sec_median", "sampling_sec_mad", "iter_per_sec_median", "iter_per_sec_mad",
    "core_sec_per_1k_iter_median", "reduce_sum_wall_sec_median", "leapfrog_per_sec_median"
  ))
  expect_equal(nrow(summ), 4L)
  expect_true(all(summ$n_rep == 3))
  expect_equal(
    dplyr::filter(summ, threads_per_chain == 8L, chunk_size == 60L)$sampling_sec_median,
    10.2
  )
  expect_true(all(is.finite(summ$reduce_sum_wall_sec_median)))
})

test_that("summarize_benchmark works without scenario or profile columns", {
  raw <- fake_runs(benchmark_grid(500, threads_per_chain = 1:2, chunk_size = c(5, 50), replicates = 2))
  summ <- summarize_benchmark(raw)
  expect_equal(nrow(summ), 4L)
  expect_false("scenario" %in% names(summ))
  expect_true(all(is.na(summ$reduce_sum_wall_sec_median)))
  expect_true(all(is.na(summ$leapfrog_per_sec_median)))
})

test_that("summarize_benchmark names the missing benchmark columns", {
  raw <- fake_runs(benchmark_grid(500, threads_per_chain = 1, chunk_size = 5, replicates = 1))
  old <- dplyr::rename(raw, cores = threads_per_chain, grainsize = chunk_size)
  expect_error(summarize_benchmark(old), "`threads_per_chain`, `chunk_size`")
  expect_error(summarize_benchmark(list()), "data frame")
})

test_that("benchmark_pareto marks the non-dominated cells", {
  df <- tibble::tibble(
    scenario = "a",
    cell = letters[1:4],
    tp = c(10, 12, 12, 8), # throughput, higher better
    co = c(5, 6, 9, 4) # cost, lower better
  )
  out <- benchmark_pareto(df, throughput = "tp", cost = "co")
  # b (12 @ 6) dominates c (12 @ 9); a (10 @ 5) and d (8 @ 4) survive on cost
  expect_equal(out$pareto, c(TRUE, TRUE, FALSE, TRUE))
})

test_that("benchmark_pareto computes the front within groups", {
  df <- tibble::tibble(
    scenario = rep(c("a", "b"), each = 2),
    tp = c(10, 5, 1, 2),
    co = c(1, 2, 1, 2)
  )
  out <- benchmark_pareto(df, "tp", "co", by = "scenario")
  # group a: (10@1) dominates (5@2). group b: (1@1) and (2@2) are mutually
  # non-dominated, so both are on the front.
  expect_equal(out$pareto, c(TRUE, FALSE, TRUE, TRUE))
})

test_that("benchmark_pareto picks summary or raw columns by default, and NAs missing values", {
  raw <- fake_runs(benchmark_grid(500, threads_per_chain = 1:2, chunk_size = c(5, 50), replicates = 1))
  summ <- summarize_benchmark(raw)
  expect_identical(
    benchmark_pareto(summ)$pareto,
    benchmark_pareto(summ, "iter_per_sec_median", "core_sec_per_1k_iter_median")$pareto
  )
  expect_identical(
    benchmark_pareto(raw)$pareto,
    benchmark_pareto(raw, "iter_per_sec", "core_sec_per_1k_iter")$pareto
  )
  raw$iter_per_sec[1] <- NA
  expect_true(is.na(benchmark_pareto(raw)$pareto[1]))
  expect_error(benchmark_pareto(raw, throughput = "nope"), "`nope`")
})

test_that("plot_benchmark returns a ggplot, from a summary or a bilatr_benchmark", {
  skip_if_not_installed("ggplot2")
  summ <- tibble::tibble(
    scenario = "prd_quad",
    threads_per_chain = c(8L, 16L, 8L, 16L),
    chunk_size = c(60L, 60L, 120L, 120L),
    iter_per_sec_median = c(5, 7, 6, 6.5),
    pareto = c(FALSE, TRUE, FALSE, FALSE)
  )
  expect_s3_class(plot_benchmark(summ), "ggplot")
  expect_error(plot_benchmark(summ, metric = "nope"), "`nope`")

  raw <- fake_runs(benchmark_grid(500, threads_per_chain = 1:2, chunk_size = c(5, 50), replicates = 2))
  class(raw) <- c("bilatr_benchmark", class(raw))
  expect_s3_class(plot_benchmark(raw), "ggplot")
  expect_s3_class(ggplot2::autoplot(raw), "ggplot")
})

test_that("print.bilatr_benchmark shows the settings line and the runs", {
  raw <- fake_runs(benchmark_grid(300, threads_per_chain = 1:2, chunk_size = 10, replicates = 1))
  attr(raw, "settings") <- list(
    stan_model = "stable", D = 300, T = 10, A = 5, iter_warmup = 100,
    iter_sampling = 100, seed = 1, profile = FALSE, available_cores = 8L
  )
  class(raw) <- c("bilatr_benchmark", class(raw))
  out <- capture.output(print(raw))
  expect_match(out[1], "<bilatr_benchmark> stable | D = 300, T = 10, A = 5 | 2 run(s)", fixed = TRUE)
  expect_true(any(grepl("threads_per_chain", out)))
})

# -- benchmark_bilatr() ------------------------------------------------------

test_that("benchmark_bilatr() checks threads against the available cores before compiling", {
  sd <- readRDS(file.path(dirname(bilatr_example_fit()[1]), "example_stan_data.rds"))
  local_mocked_bindings(
    .bilatr_available_cores = function() 2L,
    .compile_stan_model = function(...) stop("compiled"),
    .require_cmdstan = function(...) stop("required cmdstan")
  )
  expect_error(benchmark_bilatr(sd, threads_per_chain = c(1, 4)), "includes 4, but only 2 core")
  # within the limit, it gets as far as compiling
  expect_error(benchmark_bilatr(sd, threads_per_chain = 1:2), "required cmdstan")
})

test_that("benchmark_bilatr() validates its arguments before compiling", {
  sd <- readRDS(file.path(dirname(bilatr_example_fit()[1]), "example_stan_data.rds"))
  local_mocked_bindings(.require_cmdstan = function(...) stop("required cmdstan"))
  expect_error(benchmark_bilatr(sd, threads_per_chain = 0), "threads_per_chain")
  expect_error(benchmark_bilatr(sd, threads_per_chain = 1, chunk_size = -1), "chunk_size")
  expect_error(benchmark_bilatr(sd, threads_per_chain = 1, replicates = 0), "replicates")
  expect_error(benchmark_bilatr(sd, threads_per_chain = 1, iter_sampling = 1.5), "iter_sampling")
  expect_error(benchmark_bilatr(sd, threads_per_chain = 1, profile = NA), "profile")
  expect_error(benchmark_bilatr(sd, threads_per_chain = 1, stan_model = "nope"), "nope")
  expect_error(
    benchmark_bilatr(sd, threads_per_chain = 1, reference_category = "1"),
    "already prepared stan_data list.*`reference_category`"
  )
  expect_error(benchmark_bilatr(list(a = 1), threads_per_chain = 1), "dyad-period data frame")
})

test_that("benchmark_bilatr() passes tidyselect `...` to prepare_stan_data() unevaluated", {
  dp <- tibble::tibble(
    dyad = rep(c("a-b", "a-c", "b-c"), each = 2), period = rep(1:2, 3),
    talks = c(3, 0, 2, 1, 4, 2), threats = c(0, 1, 0, 2, 0, 1)
  )
  cols <- c("talks", "threats")
  local_mocked_bindings(.require_cmdstan = function(...) stop("required cmdstan"))
  expect_no_warning(expect_error(
    benchmark_bilatr(dp, threads_per_chain = 1, time = "period",
      action_cols = tidyselect::all_of(cols), reference_category = "talks"),
    "required cmdstan"
  ))
})

test_that("benchmark_bilatr() runs a 2 x 2 grid with real CmdStan", {
  skip_on_cran()
  skip_if_no_cmdstan()
  skip_if(isTRUE(.bilatr_available_cores() < 2), "needs 2 cores")

  agg <- readRDS(file.path(dirname(bilatr_example_fit()[1]), "example_stan_data.rds"))
  expect_message(
    bench <- benchmark_bilatr(
      agg,
      threads_per_chain = 1:2, chunk_size = c(1, 3),
      iter_warmup = 15, iter_sampling = 15
    ),
    "4 run\\(s\\) of one chain \\(15 \\+ 15 iterations\\) on D = 6 dyads"
  )

  expect_s3_class(bench, "bilatr_benchmark")
  expect_equal(nrow(bench), 4L)
  expect_true(all(.BILATR_BENCHMARK_COLUMNS %in% names(bench)))
  expect_true("reduce_sum_wall_sec" %in% names(bench))
  expect_true(all(bench$sampling_sec > 0))
  expect_true(all(bench$reduce_sum_wall_sec > 0))
  expect_equal(bench$iter_per_sec, 15 / bench$sampling_sec)
  expect_equal(bench$core_sec_per_1k_iter, bench$threads_per_chain * bench$sampling_sec * 1000 / 15)
  expect_true(all(bench$n_leapfrog >= 15))
  expect_equal(bench$leapfrog_per_sec, bench$n_leapfrog / bench$sampling_sec)
  expect_true(all(is.finite(summarize_benchmark(bench)$leapfrog_per_sec_median)))
  expect_identical(attr(bench, "settings")$stan_model, "stable")

  summ <- summarize_benchmark(bench)
  expect_equal(nrow(summ), 4L)
  expect_true("pareto" %in% names(benchmark_pareto(bench)))
  expect_true(any(benchmark_pareto(summ)$pareto))

  # profile = FALSE: no profile column, the plain model
  plain <- benchmark_bilatr(agg, threads_per_chain = 1, chunk_size = 2,
    iter_warmup = 10, iter_sampling = 10, profile = FALSE)
  expect_false("reduce_sum_wall_sec" %in% names(plain))
})
