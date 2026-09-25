# Compositional residual check for a missing dyad-level offset beta_d
# (0.6.0). See dev/claude_code_prompt_0.6.0_residuals_and_icc.md for the
# full derivation this file implements.
#
# The model is eta_dt = alpha * theta_dt - mu_intercept. A dyad-specific
# compositional offset beta_d, if it existed, would have to be orthogonal
# to both alpha (else confounded with theta) and to 1 (softmax is
# shift-invariant), so it lives in the (A-2)-dimensional subspace
# perpendicular to both. clr(p_dt) = eta_dt - mean(eta_dt) exactly for
# softmax p, and since alpha/mu_intercept both sum to zero this is just
# alpha*theta_dt - mu_intercept -- residuals in clr space are therefore on
# the linear-predictor scale beta would live on. Decomposing an observed
# (or replicated) clr residual into its along-alpha and orthogonal-to-alpha
# components isolates exactly the part theta cannot already absorb.

#' Row-wise softmax
#'
#' @param eta An `n x A` matrix; softmax is applied to each row.
#' @return An `n x A` matrix of row-stochastic probabilities.
#' @keywords internal
.softmax_rows <- function(eta) {
  m <- apply(eta, 1, max)
  z <- exp(eta - m)
  z / rowSums(z)
}

#' Row-wise centred-log-ratio transform
#'
#' `clr(p) = log(p) - mean(log(p))`. Accepts either a plain vector (a
#' single fixed composition) or an `n x A` matrix (one composition per
#' row, e.g. per posterior draw), matching how [.var_pi_alpha()] and
#' friends in `R/category_merges.R` accept either shape for `shares`.
#'
#' @param x A positive numeric vector or `n x A` matrix.
#' @return `clr(x)`, same shape as `x`.
#' @keywords internal
.clr <- function(x) {
  # A component that underflowed to exactly 0 (possible for pbar_d, a
  # weighted mean of softmax() rows, when every observed period's eta_dt
  # for that category is extreme enough that exp() underflows) would
  # otherwise give log(0) = -Inf, poisoning the row mean and hence every
  # component of the row via NaN. Flooring at the smallest positive
  # double is a no-op for any component that isn't already zero.
  x <- pmax(x, .Machine$double.xmin)
  if (is.null(dim(x))) {
    log(x) - mean(log(x))
  } else {
    log(x) - rowMeans(log(x))
  }
}

#' Recycle a length-`A` vector into every row of an `n x A` matrix
#'
#' @keywords internal
.broadcast_rows <- function(vec, n) {
  matrix(vec, nrow = n, ncol = length(vec), byrow = TRUE)
}

#' Strip a `posterior::draws_matrix`'s class down to a plain base matrix
#'
#' `posterior::draws_matrix` objects define their own arithmetic (`Ops`)
#' methods, which do not follow plain R's row-recycling rules -- the
#' `alpha_mat * theta_dt` (an `n_draws x A` matrix times a length-
#' `n_draws` vector, relying on ordinary column-wise recycling) used
#' throughout [check_compositional_residuals()]'s per-dyad-period loop
#' needs a genuine base matrix, not a `draws_matrix`, or it errors with
#' "non-conformable arrays". Every draws matrix is stripped via this
#' immediately after subsetting to rows/columns of interest.
#'
#' @keywords internal
.as_plain_matrix <- function(m) {
  matrix(as.numeric(m), nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))
}

#' Vectorized multinomial sampling with a FIXED size and per-row
#' probabilities
#'
#' `n` (the multinomial size) is a single scalar shared across every row
#' of `p_mat` -- true here because the Dirichlet-multinomial replicate is
#' always simulated at the OBSERVED `n_dt` (see the module header), which
#' does not vary by posterior draw, only the concentration/probability
#' does. This lets the whole `nrow(p_mat)`-row multinomial draw be
#' vectorized as `A - 1` sequential [stats::rbinom()] calls (stick-
#' breaking) instead of `nrow(p_mat)` separate [stats::rmultinom()] calls
#' -- the difference that keeps [check_compositional_residuals()]'s
#' per-dyad-period simulation loop fast at the scale 1d describes (up to
#' ~1500 sampled dyads x ~35 periods).
#'
#' @param n A single non-negative integer, the multinomial size (assumed
#'   identical across rows).
#' @param p_mat An `n_draws x A` matrix of row-stochastic probabilities.
#' @return An `n_draws x A` integer matrix of multinomial draws, each row
#'   summing to `n`.
#' @keywords internal
.rmultinom_rows <- function(n, p_mat) {
  n_draws <- nrow(p_mat)
  A <- ncol(p_mat)
  y <- matrix(0L, n_draws, A)
  remaining_n <- rep(n, n_draws)
  remaining_p <- rep(1, n_draws)
  if (A > 1) {
    for (k in seq_len(A - 1L)) {
      # remaining_p can drift to a tiny negative value by the last few
      # categories from float accumulation even though it is
      # mathematically >= p_mat[, k] always; guard against a 0/~0 (or
      # negative-denominator) division producing NaN/Inf, which would
      # otherwise propagate through rbinom() as NA.
      pk <- ifelse(remaining_p > 1e-12, p_mat[, k] / remaining_p, 0)
      pk <- pmin(pmax(pk, 0), 1)
      yk <- stats::rbinom(n_draws, size = remaining_n, prob = pk)
      y[, k] <- yk
      remaining_n <- remaining_n - yk
      remaining_p <- pmax(remaining_p - p_mat[, k], 0)
    }
  }
  y[, A] <- remaining_n
  y
}

#' Vectorized Dirichlet sampling, one draw per row of a concentration
#' matrix
#'
#' @param conc_mat An `n_draws x A` matrix of positive Dirichlet
#'   concentrations, one row per draw.
#' @return An `n_draws x A` matrix of row-stochastic Dirichlet draws.
#' @keywords internal
.rdirichlet_rows <- function(conc_mat) {
  a <- as.vector(conc_mat)
  if (any(a <= 0)) {
    bad <- which(conc_mat <= 0, arr.ind = TRUE)[1, ]
    stop(sprintf(
      "`.rdirichlet_rows()` got a non-positive concentration (%.3g) at draw row %d, category column %d. conc = phi_d * p_dt is positive by construction, so this points to a bug upstream, not a data issue.",
      conc_mat[bad["row"], bad["col"]], bad["row"], bad["col"]
    ), call. = FALSE)
  }
  # Marsaglia-Tsang boost in log space: shape a+1 >= 1 never underflows, and
  # log(U)/a carries the small-shape behaviour that a direct rgamma(a) loses
  # to underflow for a << 1 (which is routine here: conc = phi_d * p_dt, and
  # phi_d is barely identified for sparse dyads). Exact, not approximate --
  # if X ~ Gamma(a + 1) and U ~ Uniform(0, 1) then X * U^(1/a) ~ Gamma(a).
  log_g <- log(stats::rgamma(length(a), shape = a + 1, rate = 1)) +
    log(stats::runif(length(a))) / a
  log_g <- matrix(log_g, nrow = nrow(conc_mat))
  log_g <- log_g - apply(log_g, 1, max) # stabilise; max component -> exp(0) = 1
  g <- exp(log_g)
  g / rowSums(g)
}

#' Simulate `yrep ~ DirMult(n, conc)`, one draw per row of `conc_mat`
#'
#' `yrep_k = Multinomial(n, q)` with `q ~ Dirichlet(conc)` -- the standard
#' Dirichlet-multinomial construction, matching
#' `dirichlet_multinomial_lpmf(Y[d,t] | conc)` in
#' `inst/stan/bilatr_alphanorm.stan` with `conc = phi[d] * p`.
#'
#' @param n A single non-negative integer (the observed `n_dt`; see
#'   [.rmultinom_rows()] for why this can be a scalar shared across rows).
#' @param conc_mat An `n_draws x A` matrix of concentrations
#'   (`phi_d * p_dt`, one row per draw).
#' @return An `n_draws x A` integer matrix.
#' @keywords internal
.dirichlet_multinomial_rows <- function(n, conc_mat) {
  .rmultinom_rows(n, .rdirichlet_rows(conc_mat))
}

#' Drop posterior draws with non-finite compositional-residual accumulators
#'
#' A defensive backstop for [check_compositional_residuals()]: both known
#' sources of non-finite `clr()` residuals are fixed at the source
#' (`.rdirichlet_rows()`'s Dirichlet-gamma underflow, `.clr()`'s
#' zero-component log; see `NEWS.md` 0.6.1), but a multi-hour production
#' run should not die on its last step if some future, unanticipated
#' source slips a non-finite value through. Rather than `na.rm = TRUE`
#' (which would hide the problem silently -- how the original crash
#' reached production undetected), whole draws with any non-finite entry
#' are dropped from every pooled statistic, with a single warning naming
#' how many, and a hard `stop()` if too few draws remain to trust the
#' result.
#'
#' @param T_obs_acc,T_rep_acc `n_draws x n_eps` pooled-statistic
#'   accumulators.
#' @param cat_contrib_obs,cat_contrib_rep `A x n_draws` per-category
#'   accumulators.
#' @param sum_theta_bar_sq,sum_theta_sq Length-`n_draws` accumulators.
#' @param min_draws Minimum number of surviving draws before this
#'   `stop()`s instead of warning. Default `20`.
#' @return A list of the same objects, each subset to the surviving
#'   draws, plus `n_draws_dropped`.
#' @keywords internal
.drop_nonfinite_draws <- function(T_obs_acc, T_rep_acc, cat_contrib_obs, cat_contrib_rep,
                                   sum_theta_bar_sq, sum_theta_sq, min_draws = 20L) {
  n_draws_used <- nrow(T_obs_acc)
  bad <- apply(T_obs_acc, 1, function(r) any(!is.finite(r))) |
    apply(T_rep_acc, 1, function(r) any(!is.finite(r))) |
    apply(cat_contrib_obs, 2, function(cc) any(!is.finite(cc))) |
    apply(cat_contrib_rep, 2, function(cc) any(!is.finite(cc)))
  n_dropped <- sum(bad)

  if (n_dropped == 0L) {
    return(list(
      T_obs_acc = T_obs_acc, T_rep_acc = T_rep_acc,
      cat_contrib_obs = cat_contrib_obs, cat_contrib_rep = cat_contrib_rep,
      sum_theta_bar_sq = sum_theta_bar_sq, sum_theta_sq = sum_theta_sq,
      n_draws_dropped = 0L
    ))
  }

  n_remaining <- n_draws_used - n_dropped
  if (n_remaining < min_draws) {
    stop(sprintf(
      paste(
        "%d of %d posterior draws produced non-finite compositional residuals,",
        "leaving only %d usable draws (fewer than %d); aborting rather than",
        "reporting an unreliable statistic. Check dyads$phi_min -- a phi_d",
        "posterior near 1e-6 is the known cause (see NEWS.md 0.6.1)."
      ),
      n_dropped, n_draws_used, n_remaining, min_draws
    ), call. = FALSE)
  }
  warning(sprintf(
    paste(
      "%d of %d posterior draws produced non-finite compositional residuals",
      "and were dropped from the pooled statistics (see global$n_draws_dropped)."
    ),
    n_dropped, n_draws_used
  ), call. = FALSE)

  good <- !bad
  list(
    T_obs_acc = T_obs_acc[good, , drop = FALSE],
    T_rep_acc = T_rep_acc[good, , drop = FALSE],
    cat_contrib_obs = cat_contrib_obs[, good, drop = FALSE],
    cat_contrib_rep = cat_contrib_rep[, good, drop = FALSE],
    sum_theta_bar_sq = sum_theta_bar_sq[good],
    sum_theta_sq = sum_theta_sq[good],
    n_draws_dropped = n_dropped
  )
}

#' Stratified sample of dyad indices by a volume statistic
#'
#' Splits `1:length(n_d)` into `n_strata` roughly equal-sized strata by
#' [dplyr::ntile()] on `n_d` (tie-robust, unlike a quantile-breaks
#' approach, which can collapse strata when `n_d` has many repeated
#' values -- common for event counts), then samples `n_target` indices
#' total, spread as evenly as possible across strata, so the volume
#' gradient stays represented in the subsample (see
#' [check_compositional_residuals()]'s `n_dyads`/`n_strata` args).
#'
#' @param n_d Numeric vector, one volume statistic per dyad (the pool
#'   being sampled from, already restricted to eligible dyads).
#' @param n_target Target sample size (capped at `length(n_d)`).
#' @param n_strata Number of strata.
#' @param seed Seed for reproducibility.
#' @return An integer vector of sampled positions into `n_d`, sorted
#'   ascending, length `min(n_target, length(n_d))`.
#' @keywords internal
.stratified_dyad_sample <- function(n_d, n_target, n_strata, seed) {
  D <- length(n_d)
  n_target <- min(n_target, D)
  n_strata <- max(1L, min(n_strata, D))
  strat <- dplyr::ntile(n_d, n_strata)

  set.seed(seed)
  per_stratum <- floor(n_target / n_strata)
  remainder <- n_target - per_stratum * n_strata
  sampled <- integer(0)
  for (s in seq_len(n_strata)) {
    idx <- which(strat == s)
    k <- min(per_stratum + as.integer(s <= remainder), length(idx))
    if (k > 0) sampled <- c(sampled, sample(idx, k))
  }
  if (length(sampled) < n_target) {
    pool <- setdiff(seq_len(D), sampled)
    extra_n <- min(n_target - length(sampled), length(pool))
    if (extra_n > 0) sampled <- c(sampled, sample(pool, extra_n))
  }
  sort(unique(sampled))
}

#' Posterior predictive check for a missing dyad-level compositional
#' offset
#'
#' The model has no dyad-specific compositional offset: `eta_dt = alpha *
#' theta_dt - mu_intercept`. An identified offset `beta_d` would have to
#' be orthogonal to both `alpha` (else confounded with `theta`) and to
#' `1` (softmax is shift-invariant), so it lives in an `A - 2`-dimensional
#' subspace. This function checks whether the observed dyad-level
#' compositional residuals in that subspace exceed what the
#' Dirichlet-multinomial likelihood itself generates.
#'
#' **The statistic** (per dyad `d`, posterior draw `s`, pooling over
#' *observed* periods only):
#' ```
#' pbar_d^s = sum_t n_dt * p_dt^s / n_d
#' r_obs^s  = clr(y_d + eps)      - clr(pbar_d^s)     # y_d = sum_t Y[d, t, ]
#' r_rep^s  = clr(yrep_d^s + eps) - clr(pbar_d^s)      # yrep_dt^s ~ DirMult(phi_d^s * p_dt^s) at the OBSERVED n_dt
#' along^s  = <r, alpha_c> / <alpha_c, alpha_c>        # theta-absorbable misfit
#' perp^s   = r - along^s * alpha_c                    # beta-absorbable misfit
#' ```
#' `alpha_c = alpha - mean(alpha)` is computed defensively every draw
#' even though it is currently a no-op (`alpha` is a `sum_to_zero_vector`
#' in every registered model) -- keep it, so a future model whose `alpha`
#' does not sum to zero stays correct.
#'
#' **The pooled statistic is the one with power; the per-dyad check is
#' calibrated but nearly powerless** (see the prompt's own simulation:
#' at `D = 400`, a per-dyad `beta` of RMS 0.40 flags only a fifth of
#' dyads, while the pooled statistic `T = sum_d ||perp_d||^2` detects
#' decisively from about `beta` RMS 0.15). Accordingly:
#' - the **pooled PPP** (`mean(T_rep >= T_obs)` across draws) is the
#'   headline test, in `global$pooled_ppp`;
#' - `dyads$ppp_dyad` is a *secondary* diagnostic for *which* dyads,
#'   never the test itself;
#' - at production scale the pooled p-value saturates (rejects for
#'   almost any non-zero `beta`), so `global$implied_beta_rms` -- a
#'   **rough moment estimator**, `sqrt(max(0, (T_obs - mean(T_rep)) / D) /
#'   A)`, good to a few tens of percent and biased slightly high, not a
#'   calibrated estimate -- is the number that actually carries decision
#'   content.
#'
#' **The comparator that makes `implied_beta_rms` interpretable** is
#' `global$theta_between_rms`: because `alpha`/`mu_intercept` both sum to
#' zero, `clr(p_dt) = alpha * theta_dt - mu_intercept` exactly, and
#' `implied_beta_rms` is a per-*component* RMS (the excess is divided by
#' `A`), so its like-for-like comparator is also per-component. The
#' per-component RMS of the theta-driven part of `clr` is `|theta_dt| *
#' sqrt(mean(alpha_k^2)) = |theta_dt|`, since `mean(alpha_k^2) = 1`
#' exactly under the RMS-1 constraint -- so the comparator reduces to the
#' RMS of `theta` itself, no outer product needed. Two versions are
#' reported: `theta_between_rms = sqrt(mean(theta_bar_d^2))` (dyad means;
#' the **headline** comparator, since `beta_d` is dyad-constant, so the
#' like-for-like quantity is the dyad-constant offset along `alpha` the
#' model already grants) and `theta_total_rms = sqrt(mean(theta_dt^2))`
#' (also includes within-dyad temporal drift, which `theta` already
#' handles and `beta` could not -- context only, typically ~10% larger).
#' Both are computed **per draw, then averaged** (not from posterior-mean
#' `theta`): per-draw RMS carries posterior noise, which slightly
#' inflates the comparator and shrinks `global$beta_signal_ratio =
#' implied_beta_rms / theta_between_rms` -- the conservative direction, a
#' deliberate bias against proposing `D x (A - 2)` new parameters. A
#' ratio near 0.05 is a curiosity; approaching 0.3-0.5 is a real case for
#' `beta_d`. (The equivalence this reduction relies on --
#' `sd(outer(theta_dt, alpha)) == sqrt(mean(theta_dt^2))` -- is asserted
#' as a test of `alpha`'s RMS-1 identification, not recomputed here; see
#' `tests/testthat/test_residual_check.R`.)
#'
#' **Two further interpretation points:**
#' - The `along` component is *conservative by construction*: `theta_dt`
#'   is fitted per dyad-period from the same data, so the model has
#'   already absorbed the along-`alpha` direction and `along ~ 0` is
#'   expected, not evidence of good fit. The exception runs the other
#'   way: for sparse dyads `theta` is shrunk toward the prior, so `along`
#'   can be non-zero there without indicating misfit -- `dyads$along_mean`
#'   is reported alongside `dyads$n_d` so this shrinkage pattern is
#'   visible rather than confusing.
#' - `mu_intercept` pins the cross-dyad mean, so the share-weighted mean
#'   of `perp` across dyads is ~0 by construction; the signal is the
#'   *dispersion*, which is why the pooled statistic sums squared norms
#'   rather than averaging a residual vector. There is deliberately no
#'   "mean perp residual" test -- it would test something the model has
#'   already fitted. `categories` instead decomposes the pooled statistic
#'   itself by category (`sum_d perp_dk^2`, observed vs. replicate), to
#'   show which categories drive the dispersion.
#'
#' **Scale and subsampling.** At production scale (`D` ~ 16750, `T` ~ 35)
#' a full sweep needs ~586k Dirichlet-multinomial simulations per draw --
#' not feasible, and not necessary for a PPC. Dyads are subsampled
#' (stratified by volume, `n_dyads`/`n_strata`, reproducible from `seed`;
#' the volume gradient is scientifically interesting in its own right,
#' since whether `beta` is needed may depend on it) and draws are
#' subsampled (`n_draws`; a PPC does not need the full posterior). Only
#' `theta[d,t]` for the sampled dyads' *observed* periods is ever read
#' (via [.get_draws()], the same reader [extract_theta()] uses), and the
#' pooled statistic is reported per sampled dyad (`D_sample`, in
#' `settings`) so results are comparable across subsample sizes.
#'
#' @param fit A `CmdStanMCMC` fit object, or a character vector of raw
#'   CmdStan CSV file paths (one per chain) -- the same forms
#'   [extract_alpha()]/[.get_draws()] accept.
#' @param stan_data The Stan data list used to produce `fit`, as returned
#'   by [assemble_stan_data()] (must still carry its `dyad_ids`
#'   attribute). Supplies `Y`, `is_obs`, `D`, `T`, `A`.
#' @param stan_model Name registered in `.bilatr_stan_models`; see
#'   [.canonical_stan_model()].
#' @param n_dyads Target size of the stratified dyad subsample (capped at
#'   the number of eligible -- i.e. `n_d > 0` -- dyads). Default `1500`.
#' @param n_strata Number of volume strata for [.stratified_dyad_sample()].
#'   Default `8`.
#' @param n_draws Target number of posterior draws to subsample (capped
#'   at the number available). Default `200`.
#' @param eps Pseudo-count added before `clr()` (needed since `clr`
#'   requires positive entries; does not bias the check -- replicates go
#'   through the identical pipeline, so `eps` is absorbed into the
#'   reference distribution -- but does affect power). Default `0.5`.
#' @param eps_sensitivity Additional `eps` values reported alongside the
#'   default in `global$eps_sensitivity`, so the choice is visible rather
#'   than assumed. Default `c(0.1, 1.0)`.
#' @param seed Seed for the dyad/draw subsampling (reproducible; the
#'   sampled dyad ids are returned in `settings$sampled_dyad_ids`).
#' @param probs Posterior interval bounds reported alongside means/
#'   medians (`dyads$along_lower`/`along_upper`,
#'   `categories$contribution_*_lower`/`_upper`). Default `c(0.05, 0.95)`.
#' @param event_classes Optional character vector of event-class labels,
#'   in `stan_data`'s action-dimension order. Defaults to `stan_data`'s
#'   `"event_classes"` attribute.
#' @param class_label_fn Optional function mapping an integer
#'   `action_index` vector to pretty labels, matching
#'   [diagnose_category_merges()]'s argument of the same name.
#' @return A `bilatr_residual_check` object: `dyads` (one row per sampled
#'   dyad: `dyad_id`, `dyad`, `dyad2`, `n_d`, `n_obs_t`, `along_mean`,
#'   `along_lower`/`along_upper`, `perp_norm2_mean`, `ppp_dyad`,
#'   `phi_min` -- the smallest `phi_d` posterior draw seen for that dyad,
#'   a modelling signal worth seeing in its own right when a dyad's `phi`
#'   is barely identified),
#'   `categories` (one row per action category: `action_index`,
#'   `event_class`, `class_label`, `contribution_obs`/`contribution_rep`
#'   and their intervals, `contribution_ratio`), `global` (`pooled_ppp`,
#'   `T_obs_mean`, `T_rep_mean`, `implied_beta_rms`, `theta_between_rms`,
#'   `theta_total_rms`, `beta_signal_ratio`, `eps_sensitivity`,
#'   `D_sample`, `n_draws_used`, `n_draws_dropped` -- posterior draws
#'   excluded from the pooled statistics because they produced a
#'   non-finite compositional residual, see [.drop_nonfinite_draws()] --
#'   and `phi_min`, the minimum of `dyads$phi_min`), and `settings` (the
#'   resolved arguments, including `sampled_dyad_ids` and
#'   `n_draws_dropped`).
#' @export
check_compositional_residuals <- function(
  fit,
  stan_data,
  stan_model = .BILATR_DEFAULT_MODEL,
  n_dyads = 1500,
  n_strata = 8,
  n_draws = 200,
  eps = 0.5,
  eps_sensitivity = c(0.1, 1.0),
  seed = 1,
  probs = c(0.05, 0.95),
  event_classes = attr(stan_data, "event_classes"),
  class_label_fn = NULL
) {
  stan_model <- .canonical_stan_model(stan_model)

  dyad_ids <- attr(stan_data, "dyad_ids")
  if (is.null(dyad_ids)) {
    stop(
      "`stan_data` must be the output of assemble_stan_data() ",
      "(missing the 'dyad_ids' attribute).",
      call. = FALSE
    )
  }

  # Model-aware, via the registry (R/model_registry.R's
  # .bilatr_model_has_gamma()) rather than a user-facing flag: `eta` for
  # a `stable_gamma` fit is alpha*theta - mu_intercept - g_d, and this
  # function is the one that motivated the offset in the first place, so
  # it must subtract g_d or it is simply testing the wrong model (see
  # inst/stan/bilatr_alphanorm_gamma.stan's header).
  has_gamma <- .bilatr_model_has_gamma(stan_model)
  if (has_gamma && is.null(stan_data$ctry_a)) {
    stop(
      "`stan_data` is missing `ctry_a`/`ctry_b`/`w_send`, needed to ",
      "compute stable_gamma's country offset g_d. This means `stan_data` ",
      "was assembled by a pre-0.7.0 assemble_stan_data() -- re-assemble ",
      "it (>= 0.7.0) before checking a stable_gamma fit; see NEWS.md.",
      call. = FALSE
    )
  }

  D <- stan_data$D
  Tn <- stan_data$T
  A <- stan_data$A
  Y <- stan_data$Y
  is_obs <- stan_data$is_obs

  # --- 1. dyad subsample (stratified by volume; free -- Y/is_obs are
  # already in memory as data, not draws) ---
  n_d_all <- vapply(seq_len(D), function(d) {
    Y_d <- Y[d, , ]
    sum(Y_d[is_obs[d, ] == 1, , drop = FALSE])
  }, numeric(1))
  eligible <- which(n_d_all > 0 & rowSums(is_obs) > 0)
  if (length(eligible) == 0) {
    stop("No dyads have any observed periods with events; nothing to check.", call. = FALSE)
  }
  sampled_rel <- .stratified_dyad_sample(n_d_all[eligible], n_target = n_dyads, n_strata = n_strata, seed = seed)
  sampled_dyad_ids <- sort(eligible[sampled_rel])

  # --- 2. read alpha/mu_intercept (Tier 1, cheap), determine the draw
  # subsample from their total draw count ---
  alpha_vars <- paste0("alpha[", seq_len(A), "]")
  mu_vars <- paste0("mu_intercept[", seq_len(A), "]")
  am_draws <- .get_draws(fit, c(alpha_vars, mu_vars))
  am_mat <- posterior::as_draws_matrix(am_draws)
  total_draws <- nrow(am_mat)
  n_draws_used <- min(n_draws, total_draws)
  set.seed(seed + 1L)
  draw_idx <- sort(sample(total_draws, n_draws_used))

  alpha_mat <- .as_plain_matrix(am_mat[draw_idx, alpha_vars, drop = FALSE])
  mu_mat <- .as_plain_matrix(am_mat[draw_idx, mu_vars, drop = FALSE])

  # --- 2b. read gamma too, for stable_gamma fits (Tier 1, cheap: A x
  # n_countries). gamma is orientation-FREE by construction (see
  # inst/stan/bilatr_alphanorm_gamma.stan's header). Row i of `gamma_mat`
  # corresponds to the SAME posterior draw as row i of
  # `alpha_mat`/`mu_mat` (draw_idx, already determined above, is reused
  # rather than re-derived -- same reasoning the module's existing
  # comment gives for phi/theta below).
  gamma_col <- NULL
  if (has_gamma) {
    gamma_draws <- .get_draws(fit, "gamma")
    gamma_mat_full <- posterior::as_draws_matrix(gamma_draws)
    gamma_mat <- .as_plain_matrix(gamma_mat_full[draw_idx, , drop = FALSE])
    gamma_col <- function(c) gamma_mat[, paste0("gamma[", seq_len(A), ",", c, "]"), drop = FALSE]
  }

  # --- 3. read phi for the sampled dyads only ---
  phi_vars <- paste0("phi[", sampled_dyad_ids, "]")
  phi_draws <- .get_draws(fit, phi_vars)
  phi_mat_full <- posterior::as_draws_matrix(phi_draws)[, phi_vars, drop = FALSE]
  phi_mat <- .as_plain_matrix(phi_mat_full[draw_idx, , drop = FALSE])

  # --- 4. read theta[d, t] for the sampled dyads' OBSERVED periods only
  # (1d: never all of Tier 3) ---
  obs_list <- lapply(sampled_dyad_ids, function(d) which(is_obs[d, ] == 1))
  n_obs_t <- vapply(obs_list, length, integer(1))
  keep <- n_obs_t > 0
  if (!all(keep)) {
    sampled_dyad_ids <- sampled_dyad_ids[keep]
    obs_list <- obs_list[keep]
    n_obs_t <- n_obs_t[keep]
    phi_mat <- phi_mat[, keep, drop = FALSE]
  }
  D_sample <- length(sampled_dyad_ids)

  theta_vars <- unlist(purrr::map2(sampled_dyad_ids, obs_list, function(d, ts) paste0("theta[", d, ",", ts, "]")))
  theta_draws <- .get_draws(fit, theta_vars)
  # Row i of every matrix read above corresponds to the same posterior
  # draw: posterior::as_draws_matrix() flattens chains/iterations in a
  # fixed, deterministic order given the same underlying fit, whether
  # read from an in-memory CmdStanMCMC or from raw CSVs via the same
  # `fit` value -- so draw_idx (computed once, from the alpha/mu read)
  # is valid to reuse for phi and theta without re-deriving it.
  theta_mat_full <- posterior::as_draws_matrix(theta_draws)[, theta_vars, drop = FALSE]
  theta_mat <- .as_plain_matrix(theta_mat_full[draw_idx, , drop = FALSE])

  # --- 5. per-dyad loop: pbar_d, y_d, yrep_d, then eps-cheap
  # r/along/perp at eps and eps_sensitivity ---
  alpha_c <- alpha_mat - rowMeans(alpha_mat) # defensive; a no-op here since alpha sums to zero exactly by construction (sum_to_zero_vector in every registered model) -- keep for a future model that does not
  alpha_c_ss <- rowSums(alpha_c^2)

  eps_all <- c(eps, setdiff(eps_sensitivity, eps))
  n_eps <- length(eps_all)

  T_obs_acc <- matrix(0, n_draws_used, n_eps)
  T_rep_acc <- matrix(0, n_draws_used, n_eps)
  sum_theta_bar_sq <- numeric(n_draws_used)
  sum_theta_sq <- numeric(n_draws_used)
  n_dt_total <- 0L
  cat_contrib_obs <- matrix(0, A, n_draws_used)
  cat_contrib_rep <- matrix(0, A, n_draws_used)
  dyad_rows <- vector("list", D_sample)

  col_offset <- 0L
  for (i in seq_len(D_sample)) {
    d <- sampled_dyad_ids[i]
    ts <- obs_list[[i]]
    k_t <- length(ts)
    cols <- (col_offset + 1L):(col_offset + k_t)
    col_offset <- col_offset + k_t
    theta_d <- theta_mat[, cols, drop = FALSE] # n_draws_used x k_t

    Y_d <- Y[d, , ]
    Y_d_obs <- Y_d[ts, , drop = FALSE] # k_t x A
    n_dt_vec <- rowSums(Y_d_obs)
    y_d <- colSums(Y_d_obs)
    n_d <- sum(n_dt_vec)

    theta_bar_d <- rowMeans(theta_d)
    sum_theta_bar_sq <- sum_theta_bar_sq + theta_bar_d^2
    sum_theta_sq <- sum_theta_sq + rowSums(theta_d^2)
    n_dt_total <- n_dt_total + k_t

    pbar_num <- matrix(0, n_draws_used, A)
    pbar_den <- 0
    yrep_d <- matrix(0, n_draws_used, A)
    phi_vec_d <- phi_mat[, i]

    # g_d, once per dyad (reused across every observed period below), not
    # once per dyad-period -- matching stable_gamma's own Stan likelihood
    # (see inst/stan/include/partial_log_lik.stanfunctions'
    # partial_log_lik_offset(), which builds g the same way, once per
    # dyad, outside its own t loop).
    g_d <- if (has_gamma) {
      ca <- stan_data$ctry_a[d]
      cb <- stan_data$ctry_b[d]
      wsend <- stan_data$w_send[d]
      if (wsend == 1) gamma_col(ca) else wsend * gamma_col(ca) + (1 - wsend) * gamma_col(cb)
    } else {
      0
    }

    for (j in seq_len(k_t)) {
      theta_dt <- theta_d[, j]
      eta <- alpha_mat * theta_dt - mu_mat - g_d
      p_dt <- .softmax_rows(eta)
      n_dt <- n_dt_vec[j]
      pbar_num <- pbar_num + n_dt * p_dt
      pbar_den <- pbar_den + n_dt
      conc <- phi_vec_d * p_dt
      yrep_d <- yrep_d + .dirichlet_multinomial_rows(n_dt, conc)
    }
    pbar_d <- pbar_num / pbar_den
    clr_pbar <- .clr(pbar_d)

    along_default <- NULL
    ppp_dyad <- NA_real_
    perp_norm2_mean <- NA_real_

    for (e_idx in seq_len(n_eps)) {
      e <- eps_all[e_idx]
      r_obs <- .broadcast_rows(.clr(y_d + e), n_draws_used) - clr_pbar
      r_rep <- .clr(yrep_d + e) - clr_pbar
      along_obs <- rowSums(r_obs * alpha_c) / alpha_c_ss
      along_rep <- rowSums(r_rep * alpha_c) / alpha_c_ss
      perp_obs <- r_obs - along_obs * alpha_c
      perp_rep <- r_rep - along_rep * alpha_c
      perp_obs_norm2 <- rowSums(perp_obs^2)
      perp_rep_norm2 <- rowSums(perp_rep^2)

      T_obs_acc[, e_idx] <- T_obs_acc[, e_idx] + perp_obs_norm2
      T_rep_acc[, e_idx] <- T_rep_acc[, e_idx] + perp_rep_norm2

      if (e_idx == 1L) {
        cat_contrib_obs <- cat_contrib_obs + t(perp_obs^2)
        cat_contrib_rep <- cat_contrib_rep + t(perp_rep^2)
        along_default <- along_obs
        ppp_dyad <- mean(perp_rep_norm2 >= perp_obs_norm2)
        perp_norm2_mean <- mean(perp_obs_norm2)
      }
    }

    along_q <- stats::quantile(along_default, probs = probs)
    dyad_rows[[i]] <- tibble::tibble(
      dyad_id = d, n_d = n_d, n_obs_t = k_t,
      along_mean = mean(along_default),
      along_lower = unname(along_q[1]), along_upper = unname(along_q[length(along_q)]),
      perp_norm2_mean = perp_norm2_mean,
      ppp_dyad = ppp_dyad,
      phi_min = min(phi_vec_d)
    )
  }

  dyads <- dplyr::bind_rows(dyad_rows) %>%
    dplyr::left_join(dplyr::distinct(dyad_ids, dyad_id, dyad, dyad2), by = "dyad_id")

  # --- 6. drop any draws with non-finite accumulators, then summarize
  # (a backstop, not the expected path -- see .drop_nonfinite_draws()) ---
  dropped <- .drop_nonfinite_draws(
    T_obs_acc, T_rep_acc, cat_contrib_obs, cat_contrib_rep, sum_theta_bar_sq, sum_theta_sq
  )
  T_obs_acc <- dropped$T_obs_acc
  T_rep_acc <- dropped$T_rep_acc
  cat_contrib_obs <- dropped$cat_contrib_obs
  cat_contrib_rep <- dropped$cat_contrib_rep
  sum_theta_bar_sq <- dropped$sum_theta_bar_sq
  sum_theta_sq <- dropped$sum_theta_sq
  n_draws_dropped <- dropped$n_draws_dropped

  theta_between_rms <- mean(sqrt(sum_theta_bar_sq / D_sample))
  theta_total_rms <- mean(sqrt(sum_theta_sq / n_dt_total))

  T_obs_mean <- colMeans(T_obs_acc)
  T_rep_mean <- colMeans(T_rep_acc)
  pooled_ppp <- colMeans(T_rep_acc >= T_obs_acc)
  implied_beta_rms <- sqrt(pmax(0, (T_obs_mean - T_rep_mean) / D_sample) / A)
  eps_sensitivity_tbl <- tibble::tibble(eps = eps_all, pooled_ppp = pooled_ppp, implied_beta_rms = implied_beta_rms)

  if (is.null(event_classes)) event_classes <- as.character(seq_len(A))
  labels <- if (!is.null(class_label_fn)) as.character(class_label_fn(seq_len(A))) else event_classes

  cat_obs_mean <- rowMeans(cat_contrib_obs)
  cat_rep_mean <- rowMeans(cat_contrib_rep)
  cat_obs_q <- t(apply(cat_contrib_obs, 1, stats::quantile, probs = probs))
  cat_rep_q <- t(apply(cat_contrib_rep, 1, stats::quantile, probs = probs))

  categories <- tibble::tibble(
    action_index = seq_len(A),
    event_class = event_classes,
    class_label = labels,
    contribution_obs = cat_obs_mean,
    contribution_obs_lower = cat_obs_q[, 1], contribution_obs_upper = cat_obs_q[, ncol(cat_obs_q)],
    contribution_rep = cat_rep_mean,
    contribution_rep_lower = cat_rep_q[, 1], contribution_rep_upper = cat_rep_q[, ncol(cat_rep_q)],
    contribution_ratio = contribution_obs / contribution_rep
  ) %>%
    dplyr::arrange(dplyr::desc(contribution_obs))

  global <- list(
    pooled_ppp = pooled_ppp[1],
    T_obs_mean = T_obs_mean[1],
    T_rep_mean = T_rep_mean[1],
    implied_beta_rms = implied_beta_rms[1],
    theta_between_rms = theta_between_rms,
    theta_total_rms = theta_total_rms,
    beta_signal_ratio = implied_beta_rms[1] / theta_between_rms,
    eps_sensitivity = eps_sensitivity_tbl,
    D_sample = D_sample,
    n_draws_used = n_draws_used,
    n_draws_dropped = n_draws_dropped,
    phi_min = min(dyads$phi_min)
  )

  settings <- list(
    seed = seed, n_dyads_requested = n_dyads, n_strata = n_strata,
    n_draws_requested = n_draws, n_draws_used = n_draws_used,
    n_draws_dropped = n_draws_dropped,
    eps = eps, eps_sensitivity = eps_sensitivity, probs = probs,
    stan_model = stan_model, D_sample = D_sample,
    sampled_dyad_ids = sampled_dyad_ids
  )

  structure(
    list(dyads = dyads, categories = categories, global = global, settings = settings),
    class = "bilatr_residual_check"
  )
}

#' Print a `bilatr_residual_check` object
#'
#' @param x A `bilatr_residual_check` object, as returned by
#'   [check_compositional_residuals()].
#' @param n_categories Maximum number of category rows to print.
#' @param ... Ignored; present for S3 consistency.
#' @return `x`, invisibly.
#' @export
print.bilatr_residual_check <- function(x, n_categories = 10, ...) {
  g <- x$global
  extreme <- g$pooled_ppp < 0.05 || g$pooled_ppp > 0.95
  cat("<bilatr_residual_check>\n\n")
  cat(sprintf(
    "dyads sampled: %d (seed = %d) | posterior draws used: %d\n\n",
    x$settings$D_sample, x$settings$seed, x$settings$n_draws_used
  ))
  cat(sprintf("pooled PPP: %.3f%s\n", g$pooled_ppp, if (extreme) " (extreme)" else ""))
  cat(sprintf(
    "the model already gives each dyad a constant offset along alpha of typical size\ntheta_between_rms = %.3f; the residuals imply each dyad also wants an offset\northogonal to alpha of typical size implied_beta_rms = %.3f.\n",
    g$theta_between_rms, g$implied_beta_rms
  ))
  cat(sprintf(
    "beta_signal_ratio = implied_beta_rms / theta_between_rms = %.3f\n(near 0.05 is a curiosity; 0.3-0.5 is a real case for beta_d; theta_total_rms = %.3f for context)\n\n",
    g$beta_signal_ratio, g$theta_total_rms
  ))
  cat(sprintf(
    "minimum phi_d draw seen across sampled dyads: %.3g%s\n",
    g$phi_min, if (g$phi_min < 1e-3) " (very small; that dyad's phi posterior is barely identified)" else ""
  ))
  if (g$n_draws_dropped > 0) {
    cat(sprintf(
      "%d of %d posterior draws were dropped from the pooled statistics (non-finite compositional residuals)\n",
      g$n_draws_dropped, x$settings$n_draws_used
    ))
  }
  cat("\n== eps sensitivity ==\n")
  print(g$eps_sensitivity)
  cat("\n== categories (top contributors to the pooled statistic, obs vs. replicate) ==\n")
  print(utils::head(dplyr::select(x$categories, class_label, contribution_obs, contribution_rep, contribution_ratio), n_categories), n = Inf)
  invisible(x)
}

.plot_residual_perp_vs_n <- function(x) {
  d <- dplyr::mutate(x$dyads, flagged = .data$ppp_dyad < 0.05)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$n_d, y = sqrt(.data$perp_norm2_mean), colour = .data$flagged)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(name = "dyad PPP < 0.05", values = c(`FALSE` = "grey60", `TRUE` = "#440154")) +
    ggplot2::labs(
      x = "n_d (total observed events)", y = "perp_obs norm (posterior mean)",
      title = "Per-dyad compositional residual magnitude vs. dyad volume"
    ) +
    ggplot2::theme_minimal()
}

.plot_residual_ppp_hist <- function(x) {
  n_bins <- 20
  expected <- nrow(x$dyads) / n_bins
  ggplot2::ggplot(x$dyads, ggplot2::aes(x = .data$ppp_dyad)) +
    ggplot2::geom_histogram(bins = n_bins, boundary = 0, fill = "#31688E", colour = "white") +
    ggplot2::geom_hline(yintercept = expected, linetype = "dashed", colour = "grey30") +
    ggplot2::labs(
      x = "per-dyad posterior predictive p-value", y = "count",
      title = "Per-dyad PPP distribution against Uniform(0,1)"
    ) +
    ggplot2::theme_minimal()
}

.plot_residual_category_contrib <- function(x) {
  d <- x$categories %>%
    dplyr::mutate(class_label = factor(.data$class_label, levels = rev(.data$class_label))) %>%
    dplyr::select("class_label", Observed = "contribution_obs", Replicate = "contribution_rep") %>%
    tidyr::pivot_longer(c("Observed", "Replicate"), names_to = "source", values_to = "contribution")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$class_label, y = .data$contribution, fill = .data$source)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::scale_fill_viridis_d(name = NULL, begin = 0.15, end = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL, y = "mean contribution to pooled T (sum of squared perp)",
      title = "Per-category contribution to the pooled statistic"
    ) +
    ggplot2::theme_minimal()
}

#' Diagnostic plots for a `bilatr_residual_check` object
#'
#' Three diagnostics: per-dyad `perp` magnitude versus dyad volume `n_d`;
#' the per-dyad PPP distribution against Uniform(0,1); and each
#' category's contribution to the pooled statistic, observed versus
#' replicate.
#'
#' @param object A `bilatr_residual_check` object.
#' @param which `"all"` (default, a named list of the three `ggplot`
#'   objects) or one of `"perp_vs_n"`, `"ppp_hist"`, `"category"` (a
#'   single `ggplot` object).
#' @param ... Ignored; present for S3 consistency.
#' @return A `ggplot` object, or (for `which = "all"`) a named list of
#'   three.
#' @exportS3Method ggplot2::autoplot
autoplot.bilatr_residual_check <- function(object, which = c("all", "perp_vs_n", "ppp_hist", "category"), ...) {
  rlang::check_installed("ggplot2", "for autoplot.bilatr_residual_check()")
  which <- match.arg(which)
  plots <- list(
    perp_vs_n = .plot_residual_perp_vs_n(object),
    ppp_hist = .plot_residual_ppp_hist(object),
    category = .plot_residual_category_contrib(object)
  )
  if (which == "all") plots else plots[[which]]
}

#' @rdname autoplot.bilatr_residual_check
#' @param x A `bilatr_residual_check` object.
#' @export
plot.bilatr_residual_check <- function(x, which = c("all", "perp_vs_n", "ppp_hist", "category"), ...) {
  rlang::check_installed("ggplot2", "for plot.bilatr_residual_check()")
  which <- match.arg(which)
  p <- ggplot2::autoplot(x, which = which)
  if (which == "all") {
    for (nm in names(p)) print(p[[nm]])
  } else {
    print(p)
  }
  invisible(x)
}

#' Write a `bilatr_residual_check` object to disk
#'
#' The one-call convenience for a thin runscript: [check_compositional_residuals()]
#' itself never touches disk (a function that writes as a side effect is
#' untestable and awkward to compose), so this is a separate,
#' non-exported-computation wrapper writing `dyads`/`categories`/
#' `global`/`eps_sensitivity` as CSVs and the three [autoplot()]
#' diagnostics as PNGs into `dir`.
#'
#' @param x A `bilatr_residual_check` object.
#' @param dir Output directory (created if it does not exist).
#' @return `x`, invisibly.
#' @export
write_residual_check <- function(x, dir) {
  stopifnot(inherits(x, "bilatr_residual_check"))
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  readr::write_csv(x$dyads, file.path(dir, "residual_check_dyads.csv"))
  readr::write_csv(x$categories, file.path(dir, "residual_check_categories.csv"))
  readr::write_csv(x$global$eps_sensitivity, file.path(dir, "residual_check_eps_sensitivity.csv"))
  global_scalars <- tibble::as_tibble(x$global[setdiff(names(x$global), "eps_sensitivity")])
  readr::write_csv(global_scalars, file.path(dir, "residual_check_global.csv"))

  rlang::check_installed("ggplot2", "for write_residual_check()'s plots")
  plots <- ggplot2::autoplot(x, which = "all")
  ggplot2::ggsave(file.path(dir, "residual_check_perp_vs_n.png"), plots$perp_vs_n, width = 7, height = 5)
  ggplot2::ggsave(file.path(dir, "residual_check_ppp_hist.png"), plots$ppp_hist, width = 7, height = 5)
  ggplot2::ggsave(file.path(dir, "residual_check_category_contrib.png"), plots$category, width = 7, height = 5)

  invisible(x)
}
