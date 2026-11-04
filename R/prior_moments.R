.alpha_prior_moments_cache <- new.env(parent = emptyenv())

#' Number of batches [alpha_prior_moments()] splits its simulation into
#'
#' Not a user-facing argument: `alpha_prior_moments()`'s cache key is
#' `(A, n_sim, fold, seed)`, so a batch count that *could* vary
#' independently of that key would let a later call with the same key
#' but a different batch count silently return a stale cached result
#' computed under the other one -- the moments are unchanged in
#' expectation but not bit-for-bit (different batching consumes the same
#' seeded RNG stream in a different shape), so that mismatch would be
#' invisible. Fixed at 10, a tenth of the memory of an unbatched
#' simulation.
#' @noRd
.BILATR_ALPHA_PRIOR_N_BATCHES <- 10L

#' Simulate draws from `alpha`'s implied prior
#'
#' Mirrors the Stan programs' `transformed parameters`
#' construction of `alpha` exactly, so the two cannot silently drift apart:
#' `alpha_raw ~ std_normal()` on a `sum_to_zero_vector[A]` is equivalent in
#' distribution to an iid-normal draw with its row mean subtracted (both are
#' the unique zero-mean isotropic Gaussian confined to the sum-zero
#' hyperplane); `alpha` is then `alpha_raw` scaled to RMS 1
#' (`sqrt(A / dot_self(alpha_raw))`, see `orientation_sign()` and the
#' `alpha` assignment in the `transformed parameters` block), and
#' folded by the sign of the first element if `fold = TRUE` (matching
#' `orientation_sign(alpha_raw)`). If that Stan construction ever changes,
#' this function must change with it.
#'
#' @param A Number of action classes.
#' @param n_sim Number of Monte Carlo draws.
#' @param fold If `TRUE` (default), apply the same sign fold the models
#'   apply to `alpha` (see `orientation_sign()`); if `FALSE`, return the
#'   unfolded (exchangeable) prior.
#' @param seed Optional seed for reproducibility. When supplied, uses
#'   [.with_seed()], so this never disturbs the caller's own
#'   global RNG stream; when `NULL`, draws are taken from whatever stream
#'   is already running (used internally by [alpha_prior_moments()]'s own
#'   batched simulation, which seeds once around the whole batch loop
#'   rather than once per batch).
#' @return An `n_sim x A` numeric matrix of simulated `alpha` draws.
#' @noRd
.simulate_alpha_prior <- function(A, n_sim = 2e6, fold = TRUE, seed = NULL) {
  sim <- function() {
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
  if (!is.null(seed)) .with_seed(seed, sim()) else sim()
}

#' Split `n_sim` Monte Carlo draws into batch sizes summing to `n_sim`
#'
#' @param n_sim Total draws.
#' @param n_batches Target number of batches.
#' @return Integer vector of batch sizes (length `<= n_batches`; any
#'   remainder from an uneven split is folded into the last batch).
#' @noRd
.alpha_prior_batch_sizes <- function(n_sim, n_batches) {
  base <- n_sim %/% n_batches
  remainder <- n_sim %% n_batches
  sizes <- rep(base, n_batches)
  if (remainder > 0) sizes[n_batches] <- sizes[n_batches] + remainder
  sizes[sizes > 0]
}

#' Sum and sum-of-squares of `.simulate_alpha_prior()`'s columns, batched
#'
#' Accumulates `colSums()`/`colSums(.^2)` across batches instead of
#' materializing the full `n_sim x A` matrix at once (at `A = 18` and
#' the default `n_sim = 2e6`, the full matrix plus its row-mean/RMS/fold
#' copies is over 1 GB of transient memory). Peak
#' memory is one batch's matrix. The batch loop draws from whatever RNG
#' stream is already running (`seed = NULL` on each
#' [.simulate_alpha_prior()] call) -- the caller seeds once, via
#' [.with_seed()], around the whole loop.
#'
#' @inheritParams .simulate_alpha_prior
#' @param n_batches Number of batches to split `n_sim` into.
#' @return A list with `n` (total draws), `sum` and `sumsq` (length-`A`
#'   numeric vectors).
#' @noRd
.batched_alpha_prior_sums <- function(A, n_sim, fold, n_batches) {
  batch_sizes <- .alpha_prior_batch_sizes(n_sim, n_batches)
  sum_vec <- numeric(A)
  sumsq_vec <- numeric(A)
  for (b in batch_sizes) {
    batch <- .simulate_alpha_prior(A, n_sim = b, fold = fold, seed = NULL)
    sum_vec <- sum_vec + colSums(batch)
    sumsq_vec <- sumsq_vec + colSums(batch^2)
  }
  list(n = sum(batch_sizes), sum = sum_vec, sumsq = sumsq_vec)
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
#' Computed by Monte Carlo, with the same construction of `alpha` as the
#' Stan programs (a sum-to-zero standard normal vector, scaled to
#' root-mean-square 1 and, with `fold = TRUE`, multiplied by the sign of
#' its first element), simulated in batches to bound memory.
#'
#' Results are reproducible: with a seed (the default is `1L`), the
#' simulation runs under that seed without changing the caller's random
#' number stream, and the result is cached for the session per `(A, n_sim,
#' fold, seed)`. `seed = NULL` draws from the current random number stream
#' and is never cached, so repeated calls give different Monte Carlo
#' noise.
#'
#' @param A Number of action classes.
#' @param n_sim Number of Monte Carlo draws. The default, `2e6`, gives
#'   moments accurate to about three decimal places.
#' @param fold If `TRUE` (default), apply the models' orientation fold
#'   (the sign of `alpha[1]`, the reference class, is made non-negative);
#'   if `FALSE`, return the unfolded, exchangeable prior.
#' @param seed Seed for reproducibility, or `NULL` for an unseeded,
#'   uncached draw. Defaults to `1L`.
#' @return A tibble with columns `action_index`, `prior_mean`, `prior_sd`.
#' @export
#' @examples
#' alpha_prior_moments(4)
alpha_prior_moments <- function(A, n_sim = 2e6, fold = TRUE, seed = 1L) {
  cache_key <- if (!is.null(seed)) paste(A, n_sim, fold, seed) else NULL
  if (!is.null(cache_key) && exists(cache_key, envir = .alpha_prior_moments_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .alpha_prior_moments_cache, inherits = FALSE))
  }

  compute <- function() .batched_alpha_prior_sums(A, n_sim, fold, .BILATR_ALPHA_PRIOR_N_BATCHES)
  sums <- if (!is.null(seed)) .with_seed(seed, compute()) else compute()

  mean_vec <- sums$sum / sums$n
  var_vec <- (sums$sumsq - sums$n * mean_vec^2) / (sums$n - 1)
  out <- tibble::tibble(
    action_index = seq_len(A),
    prior_mean = mean_vec,
    prior_sd = sqrt(var_vec)
  )

  if (!is.null(cache_key)) assign(cache_key, out, envir = .alpha_prior_moments_cache)
  out
}
