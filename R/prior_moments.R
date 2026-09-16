.alpha_prior_moments_cache <- new.env(parent = emptyenv())

#' Simulate draws from `alpha`'s implied prior
#'
#' Mirrors `inst/stan/bilatr_alphanorm.stan`'s `transformed parameters`
#' construction of `alpha` exactly, so the two cannot silently drift apart:
#' `alpha_raw ~ std_normal()` on a `sum_to_zero_vector[A]` is equivalent in
#' distribution to an iid-normal draw with its row mean subtracted (both are
#' the unique zero-mean isotropic Gaussian confined to the sum-zero
#' hyperplane); `alpha` is then `alpha_raw` scaled to RMS 1
#' (`sqrt(A / dot_self(alpha_raw))`, see `orientation_sign()` and the
#' `alpha` assignment in that file's `transformed parameters` block), and
#' folded by the sign of the first element if `fold = TRUE` (matching
#' `orientation_sign(alpha_raw)`). If that Stan construction ever changes,
#' this function must change with it.
#'
#' @param A Number of action classes.
#' @param n_sim Number of Monte Carlo draws.
#' @param fold If `TRUE` (default), apply the same sign fold `stable`/`ou`
#'   apply to `alpha` (see `orientation_sign()`); if `FALSE`, return the
#'   unfolded (exchangeable) prior.
#' @param seed Optional seed for reproducibility.
#' @return An `n_sim x A` numeric matrix of simulated `alpha` draws.
#' @keywords internal
.simulate_alpha_prior <- function(A, n_sim = 2e6, fold = TRUE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  raw <- matrix(stats::rnorm(n_sim * A), nrow = n_sim, ncol = A)
  raw <- raw - rowMeans(raw)
  rms <- sqrt(rowSums(raw^2) / A)
  alpha <- raw / rms
  if (fold) {
    s <- ifelse(alpha[, 1] >= 0, 1, -1)
    alpha <- alpha * s
  }
  alpha
}

#' Prior mean/sd of `alpha`, per action index
#'
#' `alpha` is not exchangeable across action classes even though its prior
#' is: `orientation_sign()`'s fold breaks the symmetry for the reference
#' class specifically (`action_index == 1`), since it is defined by the
#' sign of that class's own raw coordinate. `E[alpha_k^2] = 1` exactly for
#' every `k` and any `A` (the fold does not change any coordinate's square),
#' but folding a sum-to-zero vector by its first element's sign gives that
#' element a materially different mean/sd than the rest -- e.g. at `A = 10`,
#' `alpha[1]` has prior mean 0.820 and sd 0.572 while the other classes have
#' mean -0.091 (`= -E[alpha[1]] / (A - 1)`, an identity that holds in every
#' single draw, not just on average, since `alpha` sums to exactly 0 both
#' before and after folding) and sd 0.996. Treating the prior as `N(0, 1)`
#' for `k = 1` gets both the contraction and z-score materially wrong; use
#' this function's output instead wherever a per-class prior sd/mean is
#' needed (see [diagnose_category_merges()]).
#'
#' Computed by Monte Carlo (see [.simulate_alpha_prior()]) rather than a
#' closed form, so it stays correct if the underlying Stan construction
#' ever changes; results are memoised per `(A, n_sim, fold, seed)`.
#'
#' @inheritParams .simulate_alpha_prior
#' @return A tibble with columns `action_index`, `prior_mean`, `prior_sd`.
#' @export
#' @examples
#' alpha_prior_moments(4)
alpha_prior_moments <- function(A, n_sim = 2e6, fold = TRUE, seed = NULL) {
  cache_key <- paste(A, n_sim, fold, seed)
  if (exists(cache_key, envir = .alpha_prior_moments_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .alpha_prior_moments_cache, inherits = FALSE))
  }

  alpha <- .simulate_alpha_prior(A, n_sim = n_sim, fold = fold, seed = seed)
  out <- tibble::tibble(
    action_index = seq_len(A),
    prior_mean = colMeans(alpha),
    prior_sd = apply(alpha, 2, stats::sd)
  )

  assign(cache_key, out, envir = .alpha_prior_moments_cache)
  out
}
