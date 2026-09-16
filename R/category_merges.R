.softmax <- function(x) {
  z <- exp(x - max(x))
  z / sum(z)
}

#' Read a vector variable's draws as an `n_draws x length` matrix, columns
#' ordered by their bracket index (`1, 2, ..., length`) regardless of how
#' `posterior` happened to order them.
#'
#' @keywords internal
.as_ordered_matrix <- function(draws, base_name) {
  mat <- posterior::as_draws_matrix(posterior::subset_draws(draws, variable = base_name))
  idx <- as.integer(stringr::str_extract(colnames(mat), "(?<=\\[)\\d+(?=\\])"))
  mat[, order(idx), drop = FALSE]
}

#' Per-draw `Var_pi(alpha)` and its per-category decomposition
#'
#' `Var_pi(alpha) = sum_k pi_k (alpha_k - alpha_bar)^2`, the Fisher
#' information about `theta` from one unit-weight multinomial observation
#' (see [diagnose_category_merges()]). `shares` may be a fixed length-`A`
#' vector (recycled across draws) or an `n_draws x A` matrix (one share
#' vector per draw, e.g. `shares = "theta0"`).
#'
#' @return A list: `alpha_bar` (length `n_draws`), `var_pi` (length
#'   `n_draws`), `info_share` (`n_draws x A`, each row summing to 1).
#' @keywords internal
.var_pi_alpha <- function(alpha_mat, shares) {
  if (is.null(dim(shares))) {
    shares <- matrix(shares, nrow = nrow(alpha_mat), ncol = length(shares), byrow = TRUE)
  }
  alpha_bar <- rowSums(shares * alpha_mat)
  centered_sq <- (alpha_mat - alpha_bar)^2
  contrib <- shares * centered_sq
  var_pi <- rowSums(contrib)
  info_share <- contrib / var_pi
  list(alpha_bar = alpha_bar, var_pi = var_pi, info_share = info_share)
}

#' Per-draw pairwise merge loss, `loss(j, k) = pi_j * pi_k * (alpha_j -
#' alpha_k)^2 / (pi_j + pi_k)` -- see [diagnose_category_merges()] for the
#' derivation.
#'
#' @keywords internal
.pairwise_merge_loss <- function(alpha_mat, shares) {
  if (is.null(dim(shares))) {
    shares <- matrix(shares, nrow = nrow(alpha_mat), ncol = length(shares), byrow = TRUE)
  }
  A <- ncol(alpha_mat)
  pairs <- utils::combn(A, 2)
  loss <- matrix(NA_real_, nrow(alpha_mat), ncol(pairs))
  for (p in seq_len(ncol(pairs))) {
    j <- pairs[1, p]
    k <- pairs[2, p]
    pi_j <- shares[, j]
    pi_k <- shares[, k]
    loss[, p] <- pi_j * pi_k * (alpha_mat[, j] - alpha_mat[, k])^2 / (pi_j + pi_k)
  }
  list(pairs = pairs, loss = loss)
}

#' Loss from merging a general group `S` of categories (columns), relative
#' to the same fixed `alpha_bar` used for the ungrouped `Var_pi(alpha)`.
#'
#' @keywords internal
.group_merge_loss <- function(alpha_mat, shares, alpha_bar, cols) {
  if (is.null(dim(shares))) {
    shares <- matrix(shares, nrow = nrow(alpha_mat), ncol = length(shares), byrow = TRUE)
  }
  pi_s <- rowSums(shares[, cols, drop = FALSE])
  alpha_s <- rowSums(shares[, cols, drop = FALSE] * alpha_mat[, cols, drop = FALSE]) / pi_s
  within <- rowSums(shares[, cols, drop = FALSE] * (alpha_mat[, cols, drop = FALSE] - alpha_bar)^2)
  between <- pi_s * (alpha_s - alpha_bar)^2
  within - between
}

#' Greedy agglomerative merge ladder, ranked by posterior-mean
#' `pct_info_lost` at each step (see [diagnose_category_merges()] for why
#' the ranking uses the mean while every reported loss stays a full
#' per-draw distribution).
#'
#' @return A tibble, one row per step: `step`, `merged_a`, `merged_b`
#'   (labels), `loss_mean`, `loss_lower`, `loss_upper` (posterior
#'   mean/interval of that step's `pct_info_lost`), `cumulative_loss_mean`,
#'   `n_remaining`.
#' @keywords internal
.greedy_merge_ladder <- function(alpha_mat, shares, labels, probs) {
  if (is.null(dim(shares))) {
    shares <- matrix(shares, nrow = nrow(alpha_mat), ncol = length(shares), byrow = TRUE)
  }
  A <- ncol(alpha_mat)
  n_draws <- nrow(alpha_mat)
  vp <- .var_pi_alpha(alpha_mat, shares)
  alpha_bar <- vp$alpha_bar
  total_var <- vp$var_pi

  # working state: one column per surviving group, a list mapping group
  # index -> character vector of original labels it contains
  work_alpha <- alpha_mat
  work_shares <- shares
  groups <- as.list(labels)

  steps <- vector("list", A - 1)
  cumulative_mean <- 0

  for (step in seq_len(A - 1)) {
    n_groups <- length(groups)
    pairs <- utils::combn(n_groups, 2)
    loss <- matrix(NA_real_, n_draws, ncol(pairs))
    for (p in seq_len(ncol(pairs))) {
      j <- pairs[1, p]
      k <- pairs[2, p]
      pi_j <- work_shares[, j]
      pi_k <- work_shares[, k]
      loss[, p] <- pi_j * pi_k * (work_alpha[, j] - work_alpha[, k])^2 / (pi_j + pi_k)
    }
    pct_loss <- loss / total_var
    mean_pct <- colMeans(pct_loss)
    best <- which.min(mean_pct)
    j <- pairs[1, best]
    k <- pairs[2, best]

    q <- stats::quantile(pct_loss[, best], probs = probs)
    cumulative_mean <- cumulative_mean + mean_pct[best]

    steps[[step]] <- tibble::tibble(
      step = step,
      merged_a = paste(groups[[j]], collapse = "+"),
      merged_b = paste(groups[[k]], collapse = "+"),
      loss_mean = mean(pct_loss[, best]),
      loss_lower = unname(q[1]),
      loss_upper = unname(q[length(q)]),
      cumulative_loss_mean = cumulative_mean,
      n_remaining = n_groups - 1L
    )

    # merge j and k into a new working group, drop the old two columns
    pi_jk <- work_shares[, j] + work_shares[, k]
    alpha_jk <- (work_shares[, j] * work_alpha[, j] + work_shares[, k] * work_alpha[, k]) / pi_jk
    keep <- setdiff(seq_len(n_groups), c(j, k))
    work_alpha <- cbind(work_alpha[, keep, drop = FALSE], alpha_jk)
    work_shares <- cbind(work_shares[, keep, drop = FALSE], pi_jk)
    groups <- c(groups[keep], list(c(groups[[j]], groups[[k]])))
  }

  dplyr::bind_rows(steps)
}

#' Which event categories can be merged, and what it would cost
#'
#' Two categories are redundant *for the latent scale* exactly when their
#' discriminations (`alpha`) are equal: the Fisher information about
#' `theta` from one dyad-period's multinomial is `n * Var_pi(alpha)`, with
#' `Var_pi(alpha) = sum_k pi_k (alpha_k - alpha_bar)^2` and `alpha_bar =
#' sum_k pi_k alpha_k`. Merging categories `j` and `k` gives the merged
#' category the share-weighted mean discrimination, and the drop in
#' `Var_pi(alpha)` is exactly
#'
#' ```
#' loss(j, k) = pi_j * pi_k * (alpha_j - alpha_k)^2 / (pi_j + pi_k)
#' ```
#'
#' Derivation: before merging, `j`/`k` contribute `pi_j(alpha_j-alpha_bar)^2
#' + pi_k(alpha_k-alpha_bar)^2` to `Var_pi(alpha)`; merging replaces this
#' with `pi_jk(alpha_jk-alpha_bar)^2` where `pi_jk = pi_j+pi_k` and
#' `alpha_jk` is the share-weighted mean (chosen so `pi_jk*alpha_jk =
#' pi_j*alpha_j + pi_k*alpha_k`, which leaves the OVERALL `alpha_bar`
#' unchanged by the merge). Expanding the difference, every `alpha_bar`
#' term cancels and what remains is exactly `pi_j*pi_k*(alpha_j-alpha_k)^2
#' / (pi_j+pi_k)`. The general-`S` formula is the same algebra applied to
#' a group. **This is exactly Ward's linkage** on the one-dimensional
#' `alpha` values with weights `pi` (`stats::hclust(method = "ward.D")`
#' reproduces the greedy ladder below, used as a test, not the
#' implementation).
#'
#' Losses are reported as a **fraction of `Var_pi(alpha)`**, so they read
#' directly as "merging these costs x% of the information about theta."
#' The Dirichlet-multinomial overdispersion correction (`phi`) scales
#' every term equally, so it does not affect this ranking at all -- see
#' `summary$effective_info` for the one absolute (phi-corrected) figure
#' this function reports, using `phi` and the mean per-dyad-period event
#' count from `stan_data$Y`.
#'
#' **`alpha` is a discrimination, not a severity.** A category can have a
#' strongly negative `alpha` while representing objectively severe
#' conduct: rare extreme acts discriminate poorly precisely because they
#' are rare and can occur in dyads that are not systematically hostile,
#' whereas pervasive low-level behavior tracks the latent state closely.
#' Do not read the `alpha` ordering as a severity ranking.
#'
#' **Cheapness is a price, not a recommendation.** This function reports
#' what merging costs in information about `theta`; it does not say
#' whether to merge, and it never emits a recommended grouping. When the
#' cost is negligible, the decision belongs to interpretability,
#' identification hygiene, and whether a category is so sparse that its
#' Dirichlet concentration is near zero and the cell is effectively
#' unconstrained.
#'
#' **Compute per draw, then summarise.** The loss is nonlinear in `alpha`,
#' so its value at the posterior mean is not its posterior mean; every
#' quantity here is computed for each draw and then summarised, which
#' also handles `alpha`'s correlation structure for free (the sum-to-zero
#' and RMS-1 constraints make the `alpha_k` strongly negatively
#' correlated, so `sd(alpha_j - alpha_k)` is emphatically not
#' `sqrt(sd_j^2 + sd_k^2)`).
#'
#' @param fit A `CmdStanMCMC`-like fit object, or a character vector of
#'   CmdStan CSV file paths (see [.get_draws()]).
#' @param stan_data The Stan data list used to produce `fit`, as returned
#'   by [assemble_stan_data()] (used for the default empirical `shares`
#'   and the `effective_info` figure's mean event count).
#' @param stan_model Name registered in `.bilatr_stan_models`, or a
#'   recognized pre-0.4.0 alias; see [.canonical_stan_model()].
#' @param event_classes Character vector of event-class labels, length
#'   `A`, in the same order as `stan_data`'s `"event_classes"` attribute.
#'   Defaults to that attribute if present, else `"1", "2", ...`.
#' @param class_label_fn Optional function mapping an integer
#'   `action_index` vector to pretty labels (e.g. `bilatr:::
#'   eventrootcode3_name`), matching the scheme `event_classes` was built
#'   under. There is no scheme registry in this package (each CAMEO
#'   regrouping has its own standalone namer in `R/cameo_recode.R`), so
#'   the caller supplies the right one for their scheme, exactly as the
#'   `runscripts/` do. Defaults to `NULL` (raw `event_classes` used as
#'   the label too).
#' @param shares Reference category shares: `NULL` (default) uses the
#'   empirical shares from `stan_data$Y`, summed over dyads/periods and
#'   normalised (a single, FIXED reference point); `"theta0"` evaluates at
#'   `softmax(-mu_intercept)` instead, per draw; or a user-supplied
#'   numeric vector of length `A` (fixed). `pi` is a function of `theta`,
#'   so any single share vector makes this a *local* measure at one point
#'   on the scale -- a category that discriminates only at extreme
#'   `theta` will be undervalued regardless of which reference is used.
#'   The reference actually used is recorded in `$summary$shares_used`
#'   and printed.
#' @param phi Dirichlet-multinomial dispersion used only for the
#'   `effective_info` summary figure (the merge ranking is exactly
#'   invariant to it): `NULL` (default) uses, per draw, the posterior
#'   median of `phi` across dyads; or a user-supplied scalar/per-dyad
#'   vector (per-dyad vectors are averaged, since this feeds one
#'   representative figure, not a per-dyad-period computation).
#' @param probs Posterior interval bounds to report alongside the mean
#'   (first and last are used as the interval; a third, e.g. the median,
#'   is ignored for `pairwise`/`ladder` but not `categories`).
#' @return A `bilatr_category_merges` object: `categories` (one row per
#'   class: `action_index`, `event_class`, `class_label`, `post_mean`,
#'   `post_sd`, `share`, `info_share`, `contraction`, `z_score`,
#'   `precision_gain`), `pairwise` (every pair: indices, labels, `d_alpha`,
#'   `pct_info_lost` mean/interval), `ladder` (the greedy agglomerative
#'   path from `A` down to 2 -- **read this one**: it gives the whole
#'   curve so a knee can be picked rather than pricing one grouping at a
#'   time; greedy is not guaranteed globally optimal for a fixed target
#'   `A`, standard for agglomerative methods), and `summary` (scalars for
#'   `print()`).
#' @export
diagnose_category_merges <- function(
  fit,
  stan_data,
  stan_model = .BILATR_DEFAULT_MODEL,
  event_classes = NULL,
  class_label_fn = NULL,
  shares = NULL,
  phi = NULL,
  probs = c(0.05, 0.5, 0.95)
) {
  stan_model <- .canonical_stan_model(stan_model)

  need_theta0 <- identical(shares, "theta0")
  need_phi_draws <- is.null(phi)
  vars <- c("alpha", if (need_theta0) "mu_intercept", if (need_phi_draws) "phi")

  draws <- .get_draws(fit, vars)
  if (length(.bilatr_flip_variables(stan_model)) > 0) {
    draws <- bilatr_orient(draws, stan_model = stan_model, variables = vars)
  }

  alpha_mat <- .as_ordered_matrix(draws, "alpha")
  A <- ncol(alpha_mat)
  n_draws <- nrow(alpha_mat)

  if (is.null(event_classes)) {
    event_classes <- attr(stan_data, "event_classes")
  }
  if (is.null(event_classes)) {
    event_classes <- as.character(seq_len(A))
  }
  labels <- if (!is.null(class_label_fn)) {
    as.character(class_label_fn(seq_len(A)))
  } else {
    event_classes
  }

  # shares
  if (is.null(shares)) {
    totals <- apply(stan_data$Y, 3, sum)
    shares_used <- totals / sum(totals)
    shares_ref <- "empirical (stan_data$Y, summed over dyads/periods)"
  } else if (identical(shares, "theta0")) {
    mu_mat <- .as_ordered_matrix(draws, "mu_intercept")
    shares_used <- t(apply(-mu_mat, 1, .softmax))
    shares_ref <- "theta0 (softmax(-mu_intercept), per draw)"
  } else {
    if (!is.numeric(shares) || length(shares) != A) {
      stop("`shares` must be NULL, \"theta0\", or a numeric vector of length A (", A, ").", call. = FALSE)
    }
    shares_used <- shares / sum(shares)
    shares_ref <- "user-supplied"
  }

  # phi (representative scalar per draw, used only for the summary's
  # effective_info figure -- the ranking is exactly invariant to it)
  if (is.null(phi)) {
    phi_mat <- .as_ordered_matrix(draws, "phi")
    phi_vec <- apply(phi_mat, 1, stats::median)
    phi_ref <- "posterior median of phi across dyads"
  } else if (length(phi) == 1) {
    phi_vec <- rep(phi, n_draws)
    phi_ref <- "user-supplied scalar"
  } else {
    phi_vec <- rep(mean(phi), n_draws)
    phi_ref <- "user-supplied per-dyad vector (averaged)"
  }

  vp <- .var_pi_alpha(alpha_mat, shares_used)

  # categories table
  post_mean <- colMeans(alpha_mat)
  post_sd <- apply(alpha_mat, 2, stats::sd)
  prior <- alpha_prior_moments(A)
  contraction <- 1 - (post_sd^2) / (prior$prior_sd^2)
  z_score <- (post_mean - prior$prior_mean) / prior$prior_sd
  precision_gain <- prior$prior_sd / post_sd
  share_mean <- if (is.null(dim(shares_used))) shares_used else colMeans(shares_used)
  info_share_mean <- colMeans(vp$info_share)

  categories <- tibble::tibble(
    action_index = seq_len(A),
    event_class = event_classes,
    class_label = labels,
    post_mean = post_mean,
    post_sd = post_sd,
    share = share_mean,
    info_share = info_share_mean,
    contraction = contraction,
    z_score = z_score,
    precision_gain = precision_gain
  ) %>%
    dplyr::arrange(dplyr::desc(post_mean))

  # pairwise table
  pw <- .pairwise_merge_loss(alpha_mat, shares_used)
  pct_loss <- pw$loss / vp$var_pi
  q <- t(apply(pct_loss, 2, stats::quantile, probs = probs))
  pairwise <- tibble::tibble(
    index_a = pw$pairs[1, ],
    index_b = pw$pairs[2, ],
    label_a = labels[pw$pairs[1, ]],
    label_b = labels[pw$pairs[2, ]],
    d_alpha = post_mean[pw$pairs[1, ]] - post_mean[pw$pairs[2, ]],
    pct_info_lost = colMeans(pct_loss),
    pct_info_lost_lower = q[, 1],
    pct_info_lost_upper = q[, ncol(q)]
  ) %>%
    dplyr::arrange(pct_info_lost)

  ladder <- .greedy_merge_ladder(alpha_mat, shares_used, labels, probs)

  mean_n <- mean(apply(stan_data$Y, c(1, 2), sum)[stan_data$is_obs == 1])
  effective_info_draws <- mean_n * (1 + phi_vec) / (mean_n + phi_vec) * vp$var_pi

  summary_info <- list(
    A = A,
    shares_used = shares_ref,
    phi_used = phi_ref,
    mean_n_per_dyad_period = mean_n,
    var_pi_alpha_mean = mean(vp$var_pi),
    effective_info_mean = mean(effective_info_draws)
  )

  structure(
    list(
      categories = categories,
      pairwise = pairwise,
      ladder = ladder,
      summary = summary_info,
      alpha_draws = alpha_mat,
      shares = shares_used,
      event_classes = event_classes,
      class_labels = labels
    ),
    class = "bilatr_category_merges"
  )
}

#' Price an explicit grouping of event categories
#'
#' Sums [diagnose_category_merges()]'s per-draw group-merge loss over an
#' explicit list of groupings (categories not mentioned stay singleton),
#' reusing `x`'s stored per-draw `alpha`/shares rather than re-reading the
#' fit -- so a scheme already in mind can be priced directly, without
#' walking the greedy ladder.
#'
#' @param x A `bilatr_category_merges` object, as returned by
#'   [diagnose_category_merges()].
#' @param groups A list of character vectors of `event_class` labels (or
#'   integer `action_index` values), one per group to merge.
#' @return A tibble with one row: `pct_info_lost` (mean), `pct_info_lost_lower`,
#'   `pct_info_lost_upper` (a posterior interval, `probs`' outer bounds as
#'   used by `x`).
#' @export
merge_cost <- function(x, groups) {
  stopifnot(inherits(x, "bilatr_category_merges"))

  alpha_mat <- x$alpha_draws
  shares <- x$shares
  vp <- .var_pi_alpha(alpha_mat, shares)

  cols_list <- lapply(groups, function(g) {
    if (is.character(g)) {
      match(g, x$event_classes)
    } else {
      as.integer(g)
    }
  })
  if (any(vapply(cols_list, anyNA, logical(1)))) {
    stop("`groups` contains a label/index not found in `x$event_classes`.", call. = FALSE)
  }

  total_loss <- Reduce(`+`, lapply(cols_list, function(cols) {
    .group_merge_loss(alpha_mat, shares, vp$alpha_bar, cols)
  }))
  pct_loss <- total_loss / vp$var_pi
  q <- stats::quantile(pct_loss, probs = c(0.05, 0.5, 0.95))

  tibble::tibble(
    pct_info_lost = mean(pct_loss),
    pct_info_lost_lower = unname(q[1]),
    pct_info_lost_upper = unname(q[3])
  )
}

#' Print a `bilatr_category_merges` object
#'
#' @param x A `bilatr_category_merges` object, as returned by
#'   [diagnose_category_merges()].
#' @param n_pairwise Maximum number of pairwise merges to print.
#' @param n_ladder Maximum number of ladder steps to print.
#' @param ... Ignored; present for S3 consistency.
#' @return `x`, invisibly.
#' @export
print.bilatr_category_merges <- function(x, n_pairwise = 10, n_ladder = 10, ...) {
  cat("<bilatr_category_merges>\n\n")
  cat(
    "alpha is a DISCRIMINATION, not a severity -- do not read the ordering",
    "below as a severity ranking (see ?diagnose_category_merges).\n",
    "This reports what merging costs in information about theta; it does",
    "not say whether to merge.\n\n"
  )
  cat(sprintf("shares: %s\n", x$summary$shares_used))
  cat(sprintf("phi: %s\n", x$summary$phi_used))
  cat(sprintf(
    "Var_pi(alpha) (posterior mean): %.4f | effective info (DM-corrected, mean n = %.1f): %.4f\n\n",
    x$summary$var_pi_alpha_mean, x$summary$mean_n_per_dyad_period, x$summary$effective_info_mean
  ))

  cat("== Categories (sorted by alpha) ==\n")
  print(dplyr::select(x$categories, class_label, post_mean, post_sd, share, info_share, contraction, z_score), n = Inf)

  cat("\n== Cheapest pairwise merges ==\n")
  print(utils::head(dplyr::select(x$pairwise, label_a, label_b, d_alpha, pct_info_lost), n_pairwise), n = Inf)

  cat("\n== Merge ladder (greedy, cumulative) ==\n")
  print(utils::head(dplyr::select(x$ladder, step, merged_a, merged_b, loss_mean, cumulative_loss_mean, n_remaining), n_ladder), n = Inf)

  invisible(x)
}
