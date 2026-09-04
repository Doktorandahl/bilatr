# Keeps the shared partial_log_lik Stan function byte-identical across
# every registered model's .stan file, without depending on Stan's own
# #include mechanism at compile/sample time.
#
# cmdstanr (0.9.0, in use here, and every CRAN release as of writing) has a
# bug where an include_paths entry containing a space breaks at
# $sample()-time, not compile-time: the model still compiles (that goes
# through a shell-invoked Makefile), but $sample()'s internal
# self$variables() check calls stanc directly via processx with no shell,
# so a shell-quoted include path is passed through literally and stanc
# can't find the file. See https://github.com/stan-dev/cmdstanr/issues/820
# (fixed upstream in PR #1226, merged 2026-07-25, not yet in any release).
# This project's own devtools::load_all() working tree lives under a
# space-containing path, so this is not a CI-only concern: it would break
# real fit_*_dev()/fit_panel()/fit_dyad_ts() calls made from source.
#
# The shared likelihood's single source of truth is therefore
# inst/stan/include/partial_log_lik.stanfunctions, a plain text file that
# is never read by stanc. data-raw/sync_stan_functions.R uses the helpers
# below to splice it verbatim into a marker-delimited block in each model
# .stan file. Every registered model ships as an ordinary, fully
# self-contained .stan program with no runtime #include.
# tests/testthat/test_stan_includes.R re-derives each file's block from
# the canonical source and fails if they've drifted (i.e. the sync script
# wasn't re-run after an edit to the canonical source).

.partial_log_lik_begin_marker <- paste(
  "  // <<< BEGIN GENERATED partial_log_lik",
  "(source: inst/stan/include/partial_log_lik.stanfunctions) >>>"
)
.partial_log_lik_end_marker <- "  // <<< END GENERATED partial_log_lik >>>"

#' Render the marker-delimited partial_log_lik block for a model .stan file
#'
#' @param fragment_lines Character vector, the lines of
#'   `inst/stan/include/partial_log_lik.stanfunctions`.
#' @return Character vector of lines, indented for nesting one level inside
#'   a `functions { ... }` block, wrapped in begin/end markers.
#' @keywords internal
.render_partial_log_lik_block <- function(fragment_lines) {
  indented <- ifelse(
    nzchar(fragment_lines),
    paste0("  ", fragment_lines),
    fragment_lines
  )
  c(
    .partial_log_lik_begin_marker,
    "  // Do not hand-edit between these markers -- edit the source file",
    "  // above and rerun `Rscript data-raw/sync_stan_functions.R`",
    "  // (checked by tests/testthat/test_stan_includes.R).",
    indented,
    .partial_log_lik_end_marker
  )
}

#' Splice a freshly-rendered partial_log_lik block into a model .stan file
#'
#' @param stan_lines Character vector, the full current contents of a
#'   model `.stan` file (as read by [readLines()]), which must contain
#'   exactly one well-formed generated block (from a previous run of this
#'   function, or hand-seeded once at refactor time).
#' @param fragment_lines Character vector, the lines of
#'   `inst/stan/include/partial_log_lik.stanfunctions`.
#' @return Character vector, `stan_lines` with the block between (and
#'   including) the begin/end markers replaced by a fresh render.
#' @keywords internal
.splice_partial_log_lik_block <- function(stan_lines, fragment_lines) {
  begin_idx <- which(stan_lines == .partial_log_lik_begin_marker)
  end_idx <- which(stan_lines == .partial_log_lik_end_marker)
  if (length(begin_idx) != 1 || length(end_idx) != 1 || end_idx <= begin_idx) {
    stop(
      "Expected exactly one well-formed GENERATED partial_log_lik block ",
      "(BEGIN/END markers must each appear exactly once, in order).",
      call. = FALSE
    )
  }
  before <- if (begin_idx > 1) stan_lines[seq_len(begin_idx - 1)] else character(0)
  after <- if (end_idx < length(stan_lines)) {
    stan_lines[(end_idx + 1):length(stan_lines)]
  } else {
    character(0)
  }
  c(before, .render_partial_log_lik_block(fragment_lines), after)
}
