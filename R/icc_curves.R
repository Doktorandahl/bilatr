# Category response curves (0.6.0). p_k(theta) = softmax(alpha * theta -
# mu_intercept)_k over a range of theta is exactly Bock's (1972) nominal
# response model's category response function -- this likelihood is that
# model with a random-walk latent trait and Dirichlet-multinomial
# overdispersion layered on top. See
# dev/claude_code_prompt_0.6.0_residuals_and_icc.md for the design.

#' Resolve the theta grid range for [icc_curves()], cheapest path first
#'
#' Never reads all of `theta`: an explicit `theta_range` is used as-is;
#' otherwise a `theta_summary` tibble (e.g. [extract_theta()]'s output,
#' which is saved to disk per-spec as `<id>_theta.rds` in this project's
#' own runscripts, so this is normally a free, already-in-memory path)
#' has quantiles of its `mean` column taken; only as a last resort is a
#' draw-and-dyad subsample read directly off `fit`, with a `message()`
#' noting that fallback was used.
#'
#' @keywords internal
.resolve_icc_theta_range <- function(fit, stan_data, theta_range, theta_summary, quantile_probs,
                                      needs_orient, stan_model, n_dyads_fallback, seed) {
  if (!is.null(theta_range)) {
    return(theta_range)
  }
  if (!is.null(theta_summary)) {
    return(unname(stats::quantile(theta_summary$mean, probs = quantile_probs)))
  }
  if (is.null(stan_data)) {
    stop(
      "icc_curves() needs one of `theta_range`, `theta_summary`, or ",
      "`stan_data` (for a last-resort subsample) to determine the theta grid.",
      call. = FALSE
    )
  }

  D <- stan_data$D
  Tn <- stan_data$T
  n_sample <- min(n_dyads_fallback, D)
  set.seed(seed)
  sampled_d <- sort(sample(D, n_sample))
  theta_vars <- unlist(lapply(sampled_d, function(d) paste0("theta[", d, ",", seq_len(Tn), "]")))

  if (needs_orient) {
    draws <- .get_draws(fit, c("alpha[1]", theta_vars))
    draws <- bilatr_orient(draws, stan_model = stan_model, variables = "theta")
  } else {
    draws <- .get_draws(fit, theta_vars)
  }
  theta_mat <- posterior::as_draws_matrix(draws)[, theta_vars, drop = FALSE]
  theta_means <- colMeans(theta_mat)

  message(
    "icc_curves(): no `theta_range`/`theta_summary` supplied; deriving the ",
    "range from a ", n_sample, "-dyad posterior-mean subsample (last resort, ",
    "seed = ", seed, "). Pass `theta_range` or `theta_summary` to skip this read."
  )
  unname(stats::quantile(theta_means, probs = quantile_probs))
}

#' Category response curves for the bilatr dyadic IRT model
#'
#' `p_k(theta) = softmax(alpha * theta - mu_intercept)_k` over a range of
#' `theta`, for all action categories or a subset -- exactly the category
#' response curves of Bock's (1972) nominal response model, which is what
#' this likelihood is once a random-walk latent trait and
#' Dirichlet-multinomial overdispersion are layered on top.
#'
#' **Cheap, except for the default `theta_range`.** Computing the curves
#' themselves needs only `alpha` and `mu_intercept`, both length-`A`
#' vectors, read once in full (never all of `theta`). Resolving a
#' *default* `theta_range`, however, would naively need `theta`'s
#' observed spread -- see [.resolve_icc_theta_range()] for the
#' cheapest-first resolution: an explicit `theta_range`, then quantiles
#' of a `theta_summary` tibble's `mean` column (e.g. [extract_theta()]'s
#' output -- already in memory or on disk in most workflows), and only as
#' a last resort a small draw-and-dyad subsample read directly off `fit`
#' (with a `message()` noting the fallback). The default range is the
#' 1st-99th percentile of observed `theta` (`quantile_probs`), not the
#' full range, so one extreme dyad does not flatten the plot.
#'
#' **`type = "information"`** gives the per-category Fisher information
#' contribution at `theta`, `p_k(theta) * (alpha_k - alpha_bar(theta))^2`
#' with `alpha_bar(theta) = sum_j p_j(theta) alpha_j` -- exactly the
#' quantity [diagnose_category_merges()] computes at a single reference
#' share vector (`R/category_merges.R`'s [.var_pi_alpha()], reused here
#' rather than reimplemented), evaluated across a whole range of `theta`
#' instead of one point. This directly addresses that function's own
#' documented caveat that its measure is local in `theta`: seeing where
#' each category is informative across the scale is the natural
#' companion to the merge ladder. The total `Var_p(theta)(alpha)` --
#' the scale's test information function -- is attached as
#' `attr(x, "total_information")` (a `theta`/`median`/`lower`/`upper`
#' tibble), and shown as a reference curve by [autoplot()]/[plot()].
#'
#' **Bands are pointwise, not simultaneous**: at each `theta` grid point,
#' the curve is computed per posterior draw and summarised by
#' `stats::quantile(probs)` independently of every other grid point.
#'
#' **These are *mean* share curves.** The Dirichlet-multinomial adds
#' dispersion around them, so observed shares scatter more widely than
#' the bands suggest -- the bands are posterior uncertainty in the curve
#' itself, not a predictive interval for an observation.
#'
#' @param fit A `CmdStanMCMC` fit object, or a character vector of raw
#'   CmdStan CSV file paths -- the same forms [.get_draws()] accepts.
#' @param stan_data The Stan data list used to produce `fit`, as returned
#'   by [assemble_stan_data()]. Only used for the last-resort
#'   `theta_range` fallback (see Details); `NULL` is fine if `theta_range`
#'   or `theta_summary` is supplied.
#' @param theta_range Explicit `c(min, max)` for the theta grid. If
#'   `NULL` (default), resolved from `theta_summary` or, failing that,
#'   `stan_data` (see Details).
#' @param theta_summary A tibble with a `mean` column (e.g.
#'   [extract_theta()]'s output), used to derive `theta_range` via
#'   `quantile_probs` when `theta_range` is not supplied directly.
#' @param quantile_probs The two quantiles of observed `theta` bounding
#'   the default range. Default `c(0.01, 0.99)`.
#' @param categories Optional subset of categories to compute curves for:
#'   a character vector of `event_class` labels, or an integer vector of
#'   `action_index` values. `NULL` (default) computes every category.
#' @param event_classes Optional character vector of event-class labels,
#'   in `stan_data`'s action-dimension order. Defaults to `stan_data`'s
#'   `"event_classes"` attribute if `stan_data` is supplied.
#' @param class_label_fn Optional function mapping an integer
#'   `action_index` vector to pretty labels, matching
#'   [diagnose_category_merges()]'s argument of the same name.
#' @param n_grid Number of theta grid points. Default `201`.
#' @param probs The lower/upper posterior interval bounds (a third,
#'   e.g. the median, is always included). Default `c(0.05, 0.95)`.
#' @param type `"probability"` (default) for `p_k(theta)`, or
#'   `"information"` for the per-category Fisher information
#'   contribution (see Details).
#' @param stan_model Name registered in `.bilatr_stan_models`, or a
#'   recognized pre-0.4.0 alias; see [.canonical_stan_model()].
#' @param n_dyads_fallback Number of dyads to subsample for the
#'   last-resort `theta_range` fallback (see Details). Default `300`.
#' @param seed Seed for the last-resort fallback's dyad subsample.
#' @return A tibble (subclassed `bilatr_icc_curves`, with an
#'   [autoplot()]/[plot()] method and usable directly by
#'   [icc_crossings()]) of `theta`, `action_index`, `event_class`,
#'   `class_label`, `median`, `lower`, `upper`.
#' @export
icc_curves <- function(
  fit,
  stan_data = NULL,
  theta_range = NULL,
  theta_summary = NULL,
  quantile_probs = c(0.01, 0.99),
  categories = NULL,
  event_classes = NULL,
  class_label_fn = NULL,
  n_grid = 201,
  probs = c(0.05, 0.95),
  type = c("probability", "information"),
  stan_model = .BILATR_DEFAULT_MODEL,
  n_dyads_fallback = 300,
  seed = 1
) {
  type <- match.arg(type)
  stan_model <- .canonical_stan_model(stan_model)
  needs_orient <- length(.bilatr_flip_variables(stan_model)) > 0

  am_draws <- .get_draws(fit, c("alpha", "mu_intercept"))
  if (needs_orient) {
    am_draws <- bilatr_orient(am_draws, stan_model = stan_model, variables = c("alpha", "mu_intercept"))
  }
  alpha_mat <- .as_plain_matrix(.as_ordered_matrix(am_draws, "alpha"))
  mu_mat <- .as_plain_matrix(.as_ordered_matrix(am_draws, "mu_intercept"))
  A <- ncol(alpha_mat)
  n_draws <- nrow(alpha_mat)

  if (is.null(event_classes) && !is.null(stan_data)) {
    event_classes <- attr(stan_data, "event_classes")
  }
  event_classes_used <- event_classes %||% as.character(seq_len(A))
  labels <- if (!is.null(class_label_fn)) as.character(class_label_fn(seq_len(A))) else event_classes_used

  if (is.null(categories)) {
    cat_idx <- seq_len(A)
  } else if (is.character(categories)) {
    cat_idx <- match(categories, event_classes_used)
    if (anyNA(cat_idx)) {
      stop("`categories` contains a label not found in `event_classes`.", call. = FALSE)
    }
  } else {
    cat_idx <- as.integer(categories)
  }
  A_used <- length(cat_idx)

  theta_range <- .resolve_icc_theta_range(
    fit, stan_data, theta_range, theta_summary, quantile_probs,
    needs_orient, stan_model, n_dyads_fallback, seed
  )
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = n_grid)
  q_probs <- c(probs[1], 0.5, probs[length(probs)])

  median_mat <- matrix(NA_real_, n_grid, A_used)
  lower_mat <- matrix(NA_real_, n_grid, A_used)
  upper_mat <- matrix(NA_real_, n_grid, A_used)
  if (type == "information") {
    total_median <- numeric(n_grid)
    total_lower <- numeric(n_grid)
    total_upper <- numeric(n_grid)
  }

  for (g in seq_len(n_grid)) {
    theta_g <- theta_grid[g]
    eta <- alpha_mat * theta_g - mu_mat
    p_mat <- .softmax_rows(eta)

    if (type == "probability") {
      vals <- p_mat[, cat_idx, drop = FALSE]
    } else {
      vp <- .var_pi_alpha(alpha_mat, p_mat)
      contrib <- vp$info_share * vp$var_pi
      vals <- contrib[, cat_idx, drop = FALSE]
      tq <- stats::quantile(vp$var_pi, probs = q_probs)
      total_lower[g] <- unname(tq[1])
      total_median[g] <- unname(tq[2])
      total_upper[g] <- unname(tq[3])
    }

    q <- apply(vals, 2, stats::quantile, probs = q_probs)
    median_mat[g, ] <- q[2, ]
    lower_mat[g, ] <- q[1, ]
    upper_mat[g, ] <- q[3, ]
  }

  result <- tibble::tibble(
    theta = rep(theta_grid, times = A_used),
    action_index = rep(cat_idx, each = n_grid),
    median = as.vector(median_mat),
    lower = as.vector(lower_mat),
    upper = as.vector(upper_mat)
  ) %>%
    dplyr::mutate(
      event_class = event_classes_used[.data$action_index],
      class_label = labels[.data$action_index]
    ) %>%
    dplyr::select("theta", "action_index", "event_class", "class_label", "median", "lower", "upper")

  class(result) <- c("bilatr_icc_curves", class(result))
  attr(result, "type") <- type
  attr(result, "theta_range") <- theta_range
  attr(result, "probs") <- probs
  attr(result, "labels") <- labels
  attr(result, "event_classes") <- event_classes_used
  attr(result, "alpha_draws") <- alpha_mat
  attr(result, "mu_draws") <- mu_mat
  if (type == "information") {
    attr(result, "total_information") <- tibble::tibble(
      theta = theta_grid, median = total_median, lower = total_lower, upper = total_upper
    )
  }
  attr(result, "settings") <- list(
    n_grid = n_grid, quantile_probs = quantile_probs, stan_model = stan_model,
    categories = categories, n_draws = n_draws
  )

  result
}

#' Print a `bilatr_icc_curves` object
#'
#' @param x A `bilatr_icc_curves` object.
#' @param ... Passed on to the underlying tibble print method.
#' @return `x`, invisibly.
#' @export
print.bilatr_icc_curves <- function(x, ...) {
  cat(sprintf(
    "<bilatr_icc_curves> type = %s | theta range = [%.3f, %.3f]\n\n",
    attr(x, "type"), attr(x, "theta_range")[1], attr(x, "theta_range")[2]
  ))
  print(tibble::as_tibble(x), ...)
  invisible(x)
}

.plot_icc_curves <- function(x, show_total = TRUE) {
  type <- attr(x, "type")
  total <- if (type == "information") attr(x, "total_information") else NULL

  df <- dplyr::mutate(x, class_label = factor(.data$class_label, levels = unique(.data$class_label)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$theta, y = .data$median, colour = .data$class_label)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper, fill = .data$class_label),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_colour_viridis_d(name = NULL) +
    ggplot2::scale_fill_viridis_d(name = NULL) +
    ggplot2::theme_minimal()

  if (type == "probability") {
    p <- p + ggplot2::labs(
      x = "theta", y = "p_k(theta)",
      title = "Category response curves (Bock's nominal response model)"
    )
  } else {
    p <- p + ggplot2::labs(
      x = "theta", y = "per-category information contribution",
      title = "Category information curves"
    )
    if (show_total && !is.null(total)) {
      p <- p +
        ggplot2::geom_line(
          data = total, ggplot2::aes(x = .data$theta, y = .data$median),
          inherit.aes = FALSE, colour = "black", linetype = "dashed", linewidth = 0.9
        ) +
        ggplot2::labs(caption = "dashed black: total test information Var_p(theta)(alpha)")
    }
  }
  p
}

#' Plot `icc_curves()` output
#'
#' A line + pointwise-interval-ribbon plot, one series per category. For
#' `type = "information"`, the total test information curve
#' (`attr(x, "total_information")`) is overlaid as a dashed reference
#' line unless `show_total = FALSE`.
#'
#' @param object,x A `bilatr_icc_curves` object, as returned by
#'   [icc_curves()].
#' @param show_total For `type = "information"` objects, whether to
#'   overlay the total test information curve. Ignored for
#'   `type = "probability"`.
#' @param ... Ignored; present for S3 consistency.
#' @return A `ggplot` object.
#' @exportS3Method ggplot2::autoplot
autoplot.bilatr_icc_curves <- function(object, show_total = TRUE, ...) {
  rlang::check_installed("ggplot2", "for autoplot.bilatr_icc_curves()")
  .plot_icc_curves(object, show_total = show_total)
}

#' @rdname autoplot.bilatr_icc_curves
#' @export
plot.bilatr_icc_curves <- function(x, show_total = TRUE, ...) {
  rlang::check_installed("ggplot2", "for plot.bilatr_icc_curves()")
  print(ggplot2::autoplot(x, show_total = show_total))
  invisible(x)
}

#' Pairwise category-crossing points of an `icc_curves()` object
#'
#' The theta at which two categories' response curves cross is available
#' in closed form: `theta = (mu_j - mu_k) / (alpha_j - alpha_k)`, since
#' `p_j(theta) == p_k(theta)` iff `eta_j(theta) == eta_k(theta)` (softmax
#' is strictly increasing, so equal probabilities need equal linear
#' predictors) iff `alpha_j * theta - mu_j == alpha_k * theta - mu_k`.
#' Computed per posterior draw (reusing `x`'s stored `alpha`/
#' `mu_intercept` draws -- no re-read of `fit`) and summarised to a
#' median and posterior interval; undefined (`NA`) for draws where
#' `alpha_j == alpha_k` exactly.
#'
#' This ties back to [diagnose_category_merges()]'s merge diagnostics: a
#' pair whose crossing falls outside the observed `theta_range`
#' (`in_range = FALSE`) is a pair the scale never actually distinguishes
#' over the range dyads occupy, regardless of how different their
#' discriminations look in the abstract.
#'
#' @param x A `bilatr_icc_curves` object, as returned by [icc_curves()].
#' @return A tibble, one row per category pair: `index_a`, `index_b`,
#'   `label_a`, `label_b`, `crossing_median`, `crossing_lower`,
#'   `crossing_upper`, `in_range` (whether `crossing_median` falls inside
#'   `x`'s `theta_range`).
#' @export
icc_crossings <- function(x) {
  stopifnot(inherits(x, "bilatr_icc_curves"))
  alpha_mat <- attr(x, "alpha_draws")
  mu_mat <- attr(x, "mu_draws")
  labels <- attr(x, "labels")
  theta_range <- attr(x, "theta_range")
  probs <- attr(x, "probs")
  q_probs <- c(probs[1], 0.5, probs[length(probs)])

  A <- ncol(alpha_mat)
  pairs <- utils::combn(A, 2)
  n_pairs <- ncol(pairs)
  crossing_med <- numeric(n_pairs)
  crossing_low <- numeric(n_pairs)
  crossing_up <- numeric(n_pairs)

  for (p in seq_len(n_pairs)) {
    j <- pairs[1, p]
    k <- pairs[2, p]
    denom <- alpha_mat[, j] - alpha_mat[, k]
    valid <- abs(denom) > 1e-10
    if (!any(valid)) {
      crossing_med[p] <- NA_real_
      crossing_low[p] <- NA_real_
      crossing_up[p] <- NA_real_
      next
    }
    theta_jk <- (mu_mat[valid, j] - mu_mat[valid, k]) / denom[valid]
    q <- stats::quantile(theta_jk, probs = q_probs)
    crossing_low[p] <- unname(q[1])
    crossing_med[p] <- unname(q[2])
    crossing_up[p] <- unname(q[3])
  }

  in_range <- !is.na(crossing_med) & crossing_med >= theta_range[1] & crossing_med <= theta_range[2]

  tibble::tibble(
    index_a = pairs[1, ], index_b = pairs[2, ],
    label_a = labels[pairs[1, ]], label_b = labels[pairs[2, ]],
    crossing_median = crossing_med, crossing_lower = crossing_low, crossing_upper = crossing_up,
    in_range = in_range
  )
}
