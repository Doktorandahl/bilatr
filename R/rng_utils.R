#' Run code under a seeded RNG stream without disturbing the caller's own
#' random state
#'
#' Several exported functions need a reproducible-but-local random draw
#' (a stratified dyad subsample, a fallback theta-range subsample, the
#' `alpha` prior simulation) without imposing that reproducibility on
#' whatever else the caller's session is doing with randomness. Calling
#' `set.seed()` directly, as these all used to, resets the *global*
#' `.Random.seed`, so e.g. calling [check_compositional_residuals()] in
#' the middle of an unrelated simulation silently resets that
#' simulation's stream too (0.10.2, B2). This saves `.Random.seed`
#' (or notes that it was absent, e.g. before any random-number function
#' has run this session), seeds, evaluates `code`, and restores the saved
#' state -- or removes `.Random.seed` again if there wasn't one -- on
#' exit, including when `code` errors.
#'
#' @param seed Seed to set for the duration of `code`.
#' @param code An expression to evaluate under that seed.
#' @return The value of `code`.
#' @keywords internal
.with_seed <- function(seed, code) {
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  code
}
