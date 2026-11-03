# Tests for prepare_stan_data() (0.12.0) on non-event data. Equivalence
# with assemble_stan_data() on event data is in
# test_prepare_stan_data_equivalence.R. No CmdStan needed.

# Arbitrary dyad labels, no actor columns, a week_id-style time index with
# gaps (3861-3866; 3862 and 3865 never observed), and K = 3 action columns
# in non-sorted order. Counts are doubles, as user data often are.
synthetic_dp <- function() {
  tibble::tibble(
    dyad = c("north", "north", "north", "South", "South", "east"),
    week_id = c(3861L, 3863L, 3866L, 3861L, 3864L, 3866L),
    zeta = c(2, 0, 1, 0, 3, 1),
    alpha = c(0, 1, 0, 4, 0, 0),
    mid = c(1, 0, 0, 0, 1, 2)
  )
}

prep <- function(data, ...) {
  suppressMessages(prepare_stan_data(data, time = "week_id", action_cols = c(zeta, alpha, mid), ...))
}

test_that("the dyad x time grid is filled, with is_obs = 0 in the gaps", {
  sd <- prep(synthetic_dp())
  expect_identical(sd$D, 3L)
  expect_identical(sd$T, 6)
  expect_identical(sd$A, 3L)
  # dyads in C-locale order: "South" < "east" < "north"
  ids <- attr(sd, "dyad_ids")
  expect_identical(unique(ids$dyad), c("South", "east", "north"))
  expect_identical(ids$week_id, rep(3861:3866, 3))
  expect_identical(sd$is_obs, rbind(
    c(1L, 0L, 0L, 1L, 0L, 0L),
    c(0L, 0L, 0L, 0L, 0L, 1L),
    c(1L, 0L, 1L, 0L, 0L, 1L)
  ))
  expect_type(sd$Y, "integer")
  expect_identical(dim(sd$Y), c(3L, 6L, 3L))
  # north, week 3861 (t = 1): zeta 2, alpha 0, mid 1
  expect_identical(sd$Y[3, 1, ], c(2L, 0L, 1L))
  # South, week 3864 (t = 4): zeta 3, alpha 0, mid 1
  expect_identical(sd$Y[1, 4, ], c(3L, 0L, 1L))
  expect_identical(sum(sd$Y), as.integer(sum(synthetic_dp()[c("zeta", "alpha", "mid")])))
})

test_that("action columns: reference first, then the user's column order (not sorted)", {
  expect_identical(attr(prep(synthetic_dp()), "event_classes"), c("zeta", "alpha", "mid"))
  sd <- prep(synthetic_dp(), reference_category = "mid")
  expect_identical(attr(sd, "event_classes"), c("mid", "zeta", "alpha"))
  expect_identical(sd$Y[3, 1, ], c(1L, 2L, 0L))
  expect_error(prep(synthetic_dp(), reference_category = "nope"), "anchors on the first action column")
})

test_that("without actor columns, no country fields; dyad2/actor1/actor2 are NA", {
  sd <- prep(synthetic_dp())
  expect_identical(names(sd), c(
    "D", "T", "A", "C", "is_obs", "Y", "rho_prior_a", "rho_prior_b", "compute_log_lik",
    "prior_only", "compute_theta_filtered", "n_filter_dyads", "filter_dyads"
  ))
  expect_null(attr(sd, "country_codes"))
  ids <- attr(sd, "dyad_ids")
  expect_named(ids, c("dyad_id", "time_index", "dyad", "dyad2", "actor1", "actor2", "week_id"))
  expect_true(all(is.na(ids$dyad2)))
  expect_true(all(is.na(ids$actor1)))
  expect_identical(ids$dyad_id, rep(1:3, each = 6))
  expect_identical(ids$time_index, rep(1:6, 3))
})

test_that("attributes are always present", {
  sd <- prep(synthetic_dp())
  expect_identical(attr(sd, "time_col"), "week_id")
  expect_true(is.na(attr(sd, "time_origin")))
  expect_true(attr(sd, "directed"))
  expect_false(attr(sd, "w_send_imputed"))
  expect_null(attr(sd, "grouping_var"))
})

with_actors <- function() {
  dp <- synthetic_dp()
  dp$actor1 <- c("N1", "N1", "N1", "S1", "S1", "E1")
  dp$actor2 <- c("N2", "N2", "N2", "N1", "N1", "N2")
  dp
}

test_that("with actor columns, the country fields are present", {
  sd <- prep(with_actors())
  expect_true(all(c("n_countries", "ctry_a", "ctry_b", "w_send") %in% names(sd)))
  expect_identical(attr(sd, "country_codes"), c("E1", "N1", "N2", "S1"))
  # dyad order South, east, north -> (S1, N1), (E1, N2), (N1, N2)
  expect_identical(sd$ctry_a, c(4L, 1L, 2L))
  expect_identical(sd$ctry_b, c(2L, 3L, 3L))
  expect_identical(sd$w_send, c(1, 1, 1))
  ids <- attr(sd, "dyad_ids")
  expect_identical(unique(ids$dyad2), c("N1_S1", "E1_N2", "N1_N2"))
})

test_that("undirected without w_send: w_send is 0.5 and w_send_imputed is TRUE, with no warning", {
  expect_no_warning(sd <- prep(with_actors(), directed = FALSE))
  expect_identical(sd$w_send, c(0.5, 0.5, 0.5))
  expect_true(attr(sd, "w_send_imputed"))
  expect_false(attr(sd, "directed"))

  dp <- with_actors()
  dp$w_send <- c(0.2, 0.2, 0.2, 0.9, 0.9, 0)
  sd <- prep(dp, directed = FALSE)
  expect_identical(sd$w_send, c(0.9, 0, 0.2))
  expect_false(attr(sd, "w_send_imputed"))
})

test_that("directed with a w_send != 1 errors", {
  dp <- with_actors()
  dp$w_send <- c(1, 1, 1, 0.5, 0.5, 1)
  expect_error(prep(dp), "`directed = TRUE`, but `w_send` is not 1")
  expect_error(prep(dp, directed = "yes"), "`directed` must be TRUE or FALSE")
})

test_that("the dyad column is built from the actors when absent", {
  dp <- dplyr::select(with_actors(), -"dyad")
  sd <- prep(dp)
  expect_identical(unique(attr(sd, "dyad_ids")$dyad), c("E1_N2", "N1_N2", "S1_N1"))
})

test_that("time = NULL: from the attribute, from a single known name, or an error", {
  dp <- synthetic_dp()
  expect_message(
    sd <- prepare_stan_data(dp, action_cols = c(zeta, alpha, mid)),
    "using `week_id` as the time column"
  )
  expect_identical(attr(sd, "time_col"), "week_id")

  renamed <- dplyr::rename(dp, period = "week_id")
  attr(renamed, "time_col") <- "period"
  expect_identical(
    attr(prep_attr <- suppressMessages(prepare_stan_data(renamed, action_cols = c(zeta, alpha, mid))), "time_col"),
    "period"
  )
  expect_identical(prep_attr$Y, sd$Y)

  ambiguous <- dp
  ambiguous$year <- 2020L
  expect_error(
    prepare_stan_data(ambiguous, action_cols = c(zeta, alpha, mid)),
    "cannot tell which column.*found year, week_id.*set `time`"
  )
  expect_error(
    prepare_stan_data(dplyr::rename(dp, period = "week_id"), action_cols = c(zeta, alpha, mid)),
    "cannot tell which column is the time period; set `time`"
  )
})

test_that("all-zero action columns are dropped with a message, before min_n_events", {
  dp <- synthetic_dp()
  dp$empty <- 0
  expect_message(
    sd <- prepare_stan_data(dp, time = "week_id", action_cols = c(zeta, empty, alpha, mid)),
    "dropping action column\\(s\\) with zero events in every row: empty"
  )
  expect_identical(attr(sd, "event_classes"), c("zeta", "alpha", "mid"))

  # a column whose only events are in a dyad min_n_events drops is kept
  sd <- prep(synthetic_dp(), min_n_events = 4)
  expect_identical(unique(attr(sd, "dyad_ids")$dyad), c("South", "north"))
  expect_identical(attr(sd, "event_classes"), c("zeta", "alpha", "mid"))

  dp$alpha <- 0
  dp$mid <- 0
  expect_error(
    suppressMessages(prepare_stan_data(dp, time = "week_id", action_cols = c(zeta, empty, alpha, mid))),
    "fewer than two action columns have any events"
  )
})

test_that("min_n_events drops dyads, and D == 0 errors", {
  sd <- prep(synthetic_dp(), min_n_events = 5)
  expect_identical(sd$D, 2L)
  expect_identical(dim(sd$Y), c(2L, 6L, 3L))
  expect_error(prep(synthetic_dp(), min_n_events = 100), "No dyads have at least `min_n_events`")
})

test_that("filter_dyads resolves to dyad_ids and errors on unknown dyads", {
  sd <- prep(synthetic_dp(), compute_theta_filtered = 1, filter_dyads = c("north", "South"))
  expect_identical(sd$filter_dyads, c(3L, 1L))
  expect_identical(sd$n_filter_dyads, 2L)
  sd <- prep(synthetic_dp(), compute_theta_filtered = 1)
  expect_identical(sd$filter_dyads, 1:3)
  expect_error(
    prep(synthetic_dp(), compute_theta_filtered = 1, filter_dyads = "west"),
    "`filter_dyads` names not found.*west"
  )
})

test_that("time_labels are joined into dyad_ids, and the attribute may be dropped", {
  labels <- tibble::tibble(week_id = 3861:3866, week_start = as.Date("1946-01-01") + (3861:3866) * 7L)
  sd <- prep(synthetic_dp(), time_labels = labels)
  ids <- attr(sd, "dyad_ids")
  expect_named(ids, c("dyad_id", "time_index", "dyad", "dyad2", "actor1", "actor2", "week_id", "week_start"))
  expect_identical(ids$week_start[1], as.Date("2019-12-31"))
  expect_error(prep(synthetic_dp(), time_labels = tibble::tibble(year = 1)), "must be a data frame with a `week_id` column")

  # an aggregate_events() table whose attributes were dropped still works,
  # given `time` (and gives the same Y)
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "RUS"),
    Actor2CountryCode = c("RUS", "CHN", "USA"),
    SQLDATE = c(20180105L, 20190610L, 20190302L),
    PentaClass = c(0, 1, 2)
  )
  agg <- aggregate_events(events, "PentaClass")
  plain <- agg
  attributes(plain) <- attributes(agg)[c("names", "row.names", "class")]
  expect_null(attr(plain, "time_col"))
  with_attr <- suppressMessages(prepare_stan_data(agg, action_cols = c("0", "1", "2")))
  without <- suppressMessages(prepare_stan_data(plain, time = "year", action_cols = c("0", "1", "2")))
  expect_identical(without$Y, with_attr$Y)
  expect_identical(attr(with_attr, "grouping_var"), "PentaClass")
  expect_null(attr(without, "grouping_var"))
})

test_that("dyad_id order is C-locale for mixed-case and numeric labels, matching dplyr::group_by()", {
  labels <- c("usa", "USA", "Abc", "abc", "2", "10", "B", "_x", "a")
  dp <- tibble::tibble(
    dyad = rep(labels, each = 2),
    week_id = rep(1:2, length(labels)),
    k1 = 1, k2 = rep(0:1, length(labels))
  )
  sd <- suppressMessages(prepare_stan_data(dp, time = "week_id", action_cols = c(k1, k2)))
  dyads <- unique(attr(sd, "dyad_ids")$dyad)
  expect_identical(dyads, c("10", "2", "Abc", "B", "USA", "_x", "a", "abc", "usa"))
  expect_identical(dyads, sort(labels, method = "radix"))

  # dplyr::group_by() (what make_dyad_ids() uses to number dyads) gives
  # the same order under the default dplyr.legacy_locale = FALSE
  withr::local_options(dplyr.legacy_locale = NULL)
  grouped <- dplyr::group_by(tibble::tibble(dyad = labels), .data$dyad) %>%
    dplyr::summarise(id = dplyr::cur_group_id())
  expect_identical(grouped$dyad, dyads)
  expect_identical(grouped$id, seq_along(dyads))
})

test_that("non-event duplicates and bad tables fail validation through prepare_stan_data()", {
  dp <- synthetic_dp()
  dp$week_id[2] <- 3861L
  expect_error(prep(dp), "duplicate dyad x `week_id` rows.*north")
})

# --- 0.12.0b ------------------------------------------------------------

agg_fixture <- function() {
  events <- tibble::tibble(
    Actor1CountryCode = c("USA", "USA", "RUS", "CHN", "CHN", "USA"),
    Actor2CountryCode = c("RUS", "RUS", "USA", "USA", "USA", "CHN"),
    SQLDATE = c(20180105L, 20180610L, 20190302L, 20190401L, 20200101L, 20200505L),
    PentaClass = c(0, 1, 0, 1, 3, 3)
  )
  suppressMessages(aggregate_events(events, grouping_var = "PentaClass"))
}

test_that("aggregate_events() records its class columns as the action_cols attribute (D16)", {
  agg <- agg_fixture()
  expect_identical(attr(agg, "action_cols"), c("0", "1", "3"))
})

test_that("the action_cols attribute is used silently when action_cols = NULL (D16)", {
  agg <- agg_fixture()
  expect_no_warning(sd <- suppressMessages(prepare_stan_data(agg)))
  explicit <- suppressMessages(prepare_stan_data(agg, action_cols = c("0", "1", "3")))
  expect_identical(sd, explicit)
})

test_that("a verb that drops the action_cols attribute brings the default warning back (D16)", {
  agg <- agg_fixture()
  # filter()/select()/mutate()/arrange() keep custom attributes (dplyr
  # >= 1.0's dplyr_reconstruct()), so the attribute survives them ...
  filtered <- dplyr::filter(agg, year >= 2018)
  expect_identical(attr(filtered, "action_cols"), c("0", "1", "3"))
  expect_no_warning(suppressMessages(prepare_stan_data(filtered)))
  # ... but summarise() drops it.
  resummarised <- dplyr::summarise(
    agg,
    dplyr::across(c("0", "1", "3"), sum),
    .by = c("dyad", "actor1", "actor2", "year", "w_send")
  )
  expect_null(attr(resummarised, "action_cols"))
  expect_warning(
    suppressMessages(prepare_stan_data(resummarised, time = "year")),
    "`action_cols` not set"
  )
})

test_that("an action_cols attribute naming a missing column falls back to the default, with the warning (D16)", {
  agg <- dplyr::select(agg_fixture(), -"3")
  expect_identical(attr(agg, "action_cols"), c("0", "1", "3"))
  expect_warning(
    sd <- suppressMessages(prepare_stan_data(agg)),
    "`action_cols` not set; using every column except the dyad/actor/time/w_send columns as action counts: 0, 1"
  )
  expect_identical(attr(sd, "event_classes"), c("0", "1"))
})

test_that("time_range pins T and fills the grid over it (D17)", {
  sd <- prep(synthetic_dp(), time_range = c(3859, 3868))
  expect_identical(sd$T, 10)
  ids <- attr(sd, "dyad_ids")
  expect_identical(ids$week_id, rep(3859:3868, 3))
  expect_identical(sd$is_obs[, 1:2], matrix(0L, 3, 2))
  expect_identical(sd$is_obs[, 9:10], matrix(0L, 3, 2))
  # the observed part equals the default fill
  expect_identical(sd$Y[, 3:8, ], prep(synthetic_dp())$Y)
})

test_that("time_range = NULL is the data's own range", {
  expect_identical(prep(synthetic_dp()), prep(synthetic_dp(), time_range = c(3861L, 3866L)))
})

test_that("time_range drops rows outside it with a message giving the count (D17)", {
  expect_message(
    sd <- prepare_stan_data(synthetic_dp(), time = "week_id", action_cols = c(zeta, alpha, mid), time_range = c(3861, 3864)),
    "dropping 2 row\\(s\\) outside `time_range`"
  )
  expect_identical(sd$T, 4)
  expect_identical(unique(attr(sd, "dyad_ids")$week_id), 3861:3864)
  # "east" has no events in 3861:3864 and drops out under min_n_events = 1
  expect_identical(unique(attr(sd, "dyad_ids")$dyad), c("South", "north"))
})

test_that("time_range is validated (D17)", {
  expect_error(prep(synthetic_dp(), time_range = 3861), "must be two whole numbers")
  expect_error(prep(synthetic_dp(), time_range = c(3861, 3862, 3863)), "must be two whole numbers")
  expect_error(prep(synthetic_dp(), time_range = c(3861.5, 3866)), "must be two whole numbers")
  expect_error(prep(synthetic_dp(), time_range = c("3861", "3866")), "must be two whole numbers")
  expect_error(prep(synthetic_dp(), time_range = c(NA, 3866)), "must be two whole numbers")
  expect_error(prep(synthetic_dp(), time_range = c(3866, 3861)), "min \\(3866\\) is greater than its max \\(3861\\)")
  expect_error(prep(synthetic_dp(), time_range = c(1, 10)), "no rows fall within `time_range` \\(1 to 10\\); `week_id` spans 3861 to 3866")
})

test_that("event_class_labels() returns the raw column names silently for a user table (D12)", {
  sd <- prep(synthetic_dp())
  expect_null(attr(sd, "grouping_var"))
  expect_silent(labels <- event_class_labels(sd))
  expect_identical(labels, c(zeta = "zeta", alpha = "alpha", mid = "mid"))
})
