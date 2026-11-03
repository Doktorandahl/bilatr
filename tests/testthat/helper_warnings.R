# Evaluate `expr`, muffling only the two warnings that very short test
# fits produce by design:
#
# * posterior's "The ESS has been capped to avoid unstable estimates"
#   (a handful of draws per chain);
# * cmdstanr's "E-BFMI not computed because it is undefined for posterior
#   chains of length less than 3" (one or two sampling iterations).
#
# Every other warning passes through, so a test still sees (and the
# reporter still counts) anything unexpected.
muffle_short_fit_warnings <- function(expr) {
  patterns <- c(
    "The ESS has been capped to avoid unstable estimates",
    "E-BFMI not computed"
  )
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (any(vapply(patterns, grepl, logical(1), x = msg, fixed = TRUE))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
