# Regression test of the Stan input: for event data,
# prepare_stan_data(aggregate_events(ev, ...), ...) must give exactly the
# Stan input that the frozen reference holds (every list field, and the
# event_classes/country_codes/grouping_var/dyad_ids attributes). The
# reference is the input the published results were fit to, so every
# comparison is expect_identical(). No CmdStan needed.
#
# The reference for every case (helper_equivalence_cases.R) is frozen in
# fixtures/assemble_stan_data_0.11.0_reference.rds, generated with an
# earlier, independent implementation of the data preparation. It is not
# regenerated: the comparison below allows only for storage-type
# differences, which do not change any value.

reference <- readRDS(test_path("fixtures", "assemble_stan_data_0.11.0_reference.rds"))
cases <- equivalence_cases()

test_that("the reference fixture covers exactly the case grid", {
  expect_length(cases, 74)
  expect_identical(names(reference), names(cases))
})

# Builds the case with aggregate_events() and prepare_stan_data() and
# compares with the stored reference. Messages are not compared (the
# function-name prefixes differ by design).
expect_matches_reference <- function(name) {
  old <- reference[[name]]
  # The action_cols attribute from aggregate_events() is used silently.
  expect_no_warning(new <- new_path_stan_data(cases[[name]]))

  expect_identical(names(new), names(old))
  # T, C and the 0/1 flags are integers now; the reference stored them as
  # doubles. Same values, so compare them as integers.
  as_int <- c("T", "C", "compute_log_lik", "prior_only", "compute_theta_filtered")
  for (field in names(old)) {
    ref <- if (field %in% as_int) as.integer(old[[field]]) else old[[field]]
    expect_identical(new[[field]], ref, label = paste0(name, ": new$", field))
  }
  for (a in c("event_classes", "country_codes", "grouping_var")) {
    expect_identical(attr(new, a), attr(old, a), label = paste0(name, ": attr(new, \"", a, "\")"))
  }
  old_ids <- dplyr::rename(attr(old, "dyad_ids"), actor1 = "actor_a", actor2 = "actor_b")
  new_ids <- attr(new, "dyad_ids")
  if (cases[[name]]$resolution == "monthly") {
    new_ids$month_id <- NULL
  }
  # dyad_ids always has `active` now (1.0.0f); without actor periods it is
  # all TRUE. The frozen reference predates it.
  expect_true(is.logical(new_ids$active) && all(new_ids$active), label = paste0(name, ": active all TRUE"))
  new_ids$active <- NULL
  expect_identical(new_ids, old_ids, label = paste0(name, ": dyad_ids"))
  invisible(list(old = old, new = new))
}

cases_in <- function(prefix) names(cases)[startsWith(names(cases), paste0(prefix, "/"))]

test_that("the full option grid on make_fake_events() is identical to the frozen reference", {
  expect_length(cases_in("grid"), 48)
  for (name in cases_in("grid")) {
    expect_matches_reference(name)
  }
})

test_that("min_n_events = 5 really drops dyads (the grid's dropping case is not vacuous)", {
  res1 <- expect_matches_reference("min_n/1")
  res5 <- expect_matches_reference("min_n/5")
  expect_lt(res5$new$D, res1$new$D)
  expect_false(any(grepl("FRA", attr(res5$new, "dyad_ids")$dyad)))
})

test_that("ModifiedRootCode and QuadClass, other chunk_size/prior settings, are identical", {
  expect_length(cases_in("schemes"), 4)
  for (name in cases_in("schemes")) {
    expect_matches_reference(name)
  }
})

test_that("a window that drops rows and drops a class is identical", {
  events <- window_events()
  expect_length(cases_in("window"), 4)
  for (name in cases_in("window")) {
    case <- cases[[name]]
    msgs <- capture_messages(aggregate_events(
      events, "action", case$resolution, case$directed,
      start = "2012-01-01", end = "2014-12-31"
    ))
    expect_match(msgs, "dropping [0-9]+ row\\(s\\) outside", all = FALSE)
    expect_match(msgs, "zero in-window events and get no column: 9", all = FALSE)
    res <- expect_matches_reference(name)
    expect_false("9" %in% attr(res$new, "event_classes"))
  }
})

test_that("the COW-code fixture (test_event_data.R) is identical, with default and custom column names", {
  expect_length(cases_in("cow"), 8)
  for (name in cases_in("cow")) {
    expect_matches_reference(name)
  }
})

test_that("the gamma fixtures (test_gamma_offset.R) are identical", {
  expect_length(cases_in("gamma"), 4)
  for (directed in c("directed", "undirected")) {
    res <- expect_matches_reference(paste0("gamma/plain/", directed))
  }
  # undirected: a non-trivial sender share made it through unchanged
  expect_true(any(res$new$w_send != 1))
  for (directed in c("directed", "undirected")) {
    res <- expect_matches_reference(paste0("gamma/iso/", directed))
    expect_false("ISO" %in% attr(res$new, "country_codes"))
  }
})

test_that("mixed-case and numeric actor codes give the same dyad_id order on both paths", {
  expect_length(cases_in("mixed"), 4)
  for (name in cases_in("mixed")) {
    res <- expect_matches_reference(name)
    dyads <- unique(attr(res$new, "dyad_ids")$dyad)
    expect_identical(dyads, sort(dyads, method = "radix"))
  }
})
