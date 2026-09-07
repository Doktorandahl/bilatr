# Post-hoc sign relabelling for alphanorm/alphanorm_ou's alpha/theta
# reflection symmetry.
#
# alphanorm/alphanorm_ou normalize alpha via a sum_to_zero_vector with no
# fixed element, which leaves an exact reflection symmetry: negating alpha
# together with theta and their shared upstream raw parameters leaves the
# likelihood, every prior, and the sum_to_zero_vector Jacobian unchanged
# (see each model's .stan file header, "REFLECTION SYMMETRY"). Those
# models' soft sign anchor (anchor_scale data field) makes the TARGET
# correctly specified, but does not make a chain visit the correct basin:
# the two modes are separated by a likelihood barrier of thousands of
# nats, so whichever basin a chain's init happened to land in is the one
# it reports. bilatr_init_fn() (R/fit.R) biases inits toward the
# alpha[1] > 0 basin, and fit_bilatr() warns post-sampling if a fit still
# came back in the wrong one -- but neither is a hard guarantee for every
# init/data/seed combination. bilatr_orient() is the deterministic
# fallback that always works regardless of which basin the sampler found:
# given ANY posterior draws from one of these two models, it checks the
# posterior median of alpha[1] and negates every quantity the reflection
# symmetry ties to alpha's sign if that median is negative. The
# relabelling is exact -- a genuine symmetry of the posterior -- not an
# approximation.

#' Variables the alphanorm/alphanorm_ou reflection symmetry ties to
#' alpha's sign
#'
#' `character(0)` for every other registered model: their identification
#' already fixes `alpha[1]`'s sign (`alpha[1] = 1`, hard), so they have no
#' reflection symmetry to correct.
#'
#' @param stan_model Name registered in `.bilatr_stan_models`.
#' @return Character vector of variable base names (matched against
#'   `posterior` draws column names via `.bilatr_match_draws_columns()`).
#' @keywords internal
.bilatr_flip_variables <- function(stan_model) {
  switch(
    stan_model,
    alphanorm = c("alpha", "alpha_raw", "theta", "theta0", "z_theta0", "theta_raw"),
    alphanorm_ou = c("alpha", "alpha_raw", "theta", "mu_dyad", "mu_dyad_raw", "theta_raw"),
    character(0)
  )
}

#' Match posterior draws column names against variable base names
#'
#' A "base name" like `"alpha"` should match `"alpha[1]"`, `"alpha[2]"`,
#' ... but not `"alpha_raw[1]"`; `"theta"` should match `"theta[1,1]"` but
#' not `"theta0[1]"`. Anchors each base to either the end of the name or
#' an immediately following `[`.
#'
#' @param all_names Character vector of column names (e.g.
#'   `names(posterior::as_draws_df(draws))`).
#' @param bases Character vector of variable base names to match.
#' @return The subset of `all_names` matching one of `bases`.
#' @keywords internal
.bilatr_match_draws_columns <- function(all_names, bases) {
  if (length(bases) == 0) {
    return(character(0))
  }
  pattern <- paste0("^(", paste(bases, collapse = "|"), ")(\\[|$)")
  grep(pattern, all_names, value = TRUE)
}

#' Reorient posterior draws so `alpha[1]` has the canonical (positive) sign
#'
#' `mu_intercept`, `phi`, and every scale/dispersion/ratio quantity
#' (`sigma_theta0`/`sigma_mu`, `process_noise`, `sd_stat`, `rho`,
#' `mu_log_*`, `sigma_log_*`, `within_between_ratio`) are always left
#' untouched: `alpha .* theta` is invariant under the joint negation, so
#' `mu_intercept` doesn't need to flip, and every scale/dispersion/ratio
#' quantity is a positive quantity, not a location, so none of them need
#' to either.
#'
#' @param draws A `posterior` draws object (e.g. from `fit$draws()`) that
#'   includes `"alpha[1]"` -- required to determine orientation even if
#'   `"alpha[1]"` itself is not requested via `variables`.
#' @param stan_model Name registered in `.bilatr_stan_models`. Only
#'   `"alphanorm"`/`"alphanorm_ou"` have a reflection symmetry to correct
#'   (see [.bilatr_flip_variables()]); for any other name, `draws` is
#'   returned unmodified (restricted to `variables`, if supplied), since
#'   there is nothing to fix.
#' @param variables Which variables to return. Defaults to every variable
#'   [.bilatr_flip_variables()] lists for `stan_model`; pass a subset
#'   (e.g. `"theta"`) when only that one is needed downstream.
#' @return A `posterior::draws_df` restricted to `variables`, with every
#'   variable in `intersect(variables, .bilatr_flip_variables(stan_model))`
#'   negated if the posterior median of `alpha[1]` in `draws` was
#'   negative.
#' @keywords internal
bilatr_orient <- function(draws, stan_model, variables = NULL) {
  flip_vars <- .bilatr_flip_variables(stan_model)
  requested <- variables %||% flip_vars

  df <- as.data.frame(posterior::as_draws_df(draws))
  if (!("alpha[1]" %in% names(df))) {
    stop(
      "bilatr_orient() needs \"alpha[1]\" present in `draws` to determine ",
      "orientation; it was not found. Include it in the `variables` ",
      "passed to fit$draws().",
      call. = FALSE
    )
  }

  out_cols <- .bilatr_match_draws_columns(names(df), requested)
  out <- df[c(out_cols, ".chain", ".iteration", ".draw")]

  if (length(flip_vars) > 0 && stats::median(df[["alpha[1]"]]) < 0) {
    flip_cols <- .bilatr_match_draws_columns(names(df), intersect(requested, flip_vars))
    out[flip_cols] <- -out[flip_cols]
  }

  posterior::as_draws_df(out)
}
