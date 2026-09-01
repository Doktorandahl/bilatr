test_that("build_reduce_sum_grid produces cores x grainsize x replicate rows in order", {
  grid <- build_reduce_sum_grid(
    n_dyads = 3300,
    thread_levels = c(4, 16, 48),
    replicates = 2L,
    n_grainsize = 5L
  )

  expect_s3_class(grid, "tbl_df")
  expect_equal(nrow(grid), 3 * 5 * 2)
  expect_setequal(unique(grid$cores), c(4L, 16L, 48L))
  expect_length(unique(grid$grainsize), 5L)
  expect_equal(grid$job_id, seq_len(nrow(grid)))

  # ordered by cores, then grainsize, then replicate -> a SLURM array
  # filtered to one core level maps 1..(n_grain * replicates) in order
  first_block <- dplyr::filter(grid, cores == 4L)
  expect_equal(nrow(first_block), 5 * 2)
  expect_equal(first_block$grainsize, rep(sort(unique(grid$grainsize)), each = 2))

  # derived columns
  expect_equal(grid$n_chunks, as.integer(ceiling(3300 / grid$grainsize)))
  expect_equal(grid$chunks_per_core, grid$n_chunks / grid$cores)
})

test_that("build_reduce_sum_grid scales the generated grainsize range with n_dyads", {
  small <- build_reduce_sum_grid(3300, n_grainsize = 6L)
  large <- build_reduce_sum_grid(16750, n_grainsize = 6L)

  expect_lt(max(small$grainsize), max(large$grainsize))
  # ~600 (the prior full-set optimum) sits inside the large-panel sweep
  expect_true(min(large$grainsize) < 600 && max(large$grainsize) > 600)
  expect_true(all(small$grainsize >= 1))
})

test_that("build_reduce_sum_grid honours explicit grainsize_levels and validates inputs", {
  grid <- build_reduce_sum_grid(1000, grainsize_levels = c(300, 50, 50, 100))
  expect_equal(sort(unique(grid$grainsize)), c(50L, 100L, 300L))

  expect_error(build_reduce_sum_grid(1), "n_dyads")
  expect_error(build_reduce_sum_grid(1000, replicates = 0), "replicates")
})

test_that("insert_reduce_sum_profile wraps the reduce_sum call exactly once", {
  src <- paste(
    readLines(system.file("stan", "bilatr_dirmult_irt.stan", package = "bilatr")),
    collapse = "\n"
  )
  out <- insert_reduce_sum_profile(src)

  expect_match(out, 'profile("reduce_sum_likelihood") {', fixed = TRUE)
  # the reduce_sum call and its dyad_seq setup are still present, exactly once
  expect_equal(lengths(regmatches(out, gregexpr("target += reduce_sum(", out, fixed = TRUE))), 1L)
  expect_equal(lengths(regmatches(out, gregexpr("linspaced_int_array(D, 1, D)", out, fixed = TRUE))), 1L)
  # the only added non-whitespace tokens are the profile wrapper + braces
  added <- nchar(gsub("\\s", "", out)) - nchar(gsub("\\s", "", src))
  expect_equal(added, nchar('profile("reduce_sum_likelihood"){}'))
})

test_that("insert_reduce_sum_profile errors if the target statement is missing", {
  expect_error(
    insert_reduce_sum_profile("model { target += normal_lpdf(y | 0, 1); }"),
    "found 0"
  )
})

test_that("summarise_benchmark collapses replicates to one row per cell", {
  raw <- tidyr::expand_grid(
    scenario = "prd_quad",
    cores = c(8L, 16L),
    grainsize = c(60L, 120L),
    replicate = 1:3
  )
  raw$sampling_sec <- 10 + raw$replicate * 0.1
  raw$iter_per_sec <- 75 / raw$sampling_sec
  raw$core_sec_per_1k_iter <- raw$cores * raw$sampling_sec * 1000 / 75
  raw$reduce_sum_wall_sec <- raw$sampling_sec * 0.8

  summ <- summarise_benchmark(raw)
  expect_equal(nrow(summ), 4L)
  expect_true(all(summ$n_rep == 3))
  expect_equal(
    dplyr::filter(summ, cores == 8L, grainsize == 60L)$sampling_sec_median,
    10.2
  )
  expect_true(all(is.finite(summ$reduce_sum_wall_sec_median)))
})

test_that("summarise_benchmark tolerates missing profile column", {
  raw <- tibble::tibble(
    scenario = "x", cores = 8L, grainsize = 60L, replicate = 1:2,
    sampling_sec = c(10, 12),
    iter_per_sec = 75 / c(10, 12),
    core_sec_per_1k_iter = 8 * c(10, 12) * 1000 / 75
  )
  summ <- summarise_benchmark(raw)
  expect_true(is.na(summ$reduce_sum_wall_sec_median))
})

test_that("pareto_front marks the non-dominated cells", {
  df <- tibble::tibble(
    scenario = "a",
    cell = letters[1:4],
    tp = c(10, 12, 12, 8),   # throughput, higher better
    co = c(5, 6, 9, 4)       # cost, lower better
  )
  out <- pareto_front(df, throughput = "tp", cost = "co")
  # b (12 @ 6) dominates c (12 @ 9); a (10 @ 5) and d (8 @ 4) survive on cost
  expect_equal(out$pareto, c(TRUE, TRUE, FALSE, TRUE))
})

test_that("pareto_front computes the frontier within groups", {
  df <- tibble::tibble(
    scenario = rep(c("a", "b"), each = 2),
    tp = c(10, 5, 1, 2),
    co = c(1, 2, 1, 2)
  )
  out <- pareto_front(df, "tp", "co", by = "scenario")
  # group a: (10@1) dominates (5@2). group b: (1@1) and (2@2) are mutually
  # non-dominated, so both are on the frontier.
  expect_equal(out$pareto, c(TRUE, FALSE, TRUE, TRUE))
})

test_that("plot_benchmark returns a ggplot", {
  skip_if_not_installed("ggplot2")
  summ <- tibble::tibble(
    scenario = "prd_quad",
    cores = c(8L, 16L, 8L, 16L),
    grainsize = c(60L, 60L, 120L, 120L),
    iter_per_sec_median = c(5, 7, 6, 6.5),
    pareto = c(FALSE, TRUE, FALSE, FALSE)
  )
  p <- plot_benchmark(summ)
  expect_s3_class(p, "ggplot")
})
