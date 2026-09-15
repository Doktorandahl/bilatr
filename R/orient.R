# Sign-symmetric quantities of the alpha/theta reflection symmetry
# (promoted from alphanorm/alphanorm_ou in 0.4.0; see NEWS.md), and where
# bilatr_orient() -- LEGACY-ONLY since 0.4.2 -- still applies.
#
# Every registered model's alpha_raw is a free sum_to_zero_vector with no
# fixed element, which leaves an exact reflection symmetry: negating
# alpha_raw together with theta's own raw parameters leaves the
# likelihood, every prior, and the sum_to_zero_vector Jacobian unchanged
# (see each model's .stan file header). The retired stable_soft_anchor/
# ou_soft_anchor programs (what stable/ou were named 0.4.0-0.4.1) fix
# this with a soft sign anchor (anchor_scale data field) that makes the
# TARGET correctly specified but does not make a chain visit the correct
# basin -- the two modes are separated by a likelihood barrier of
# thousands of nats, so whichever basin a chain's init happened to land
# in is the one it reports in EVERY reported quantity. bilatr_orient() is
# the deterministic fallback that always works regardless of which basin
# the sampler found: given ANY posterior draws from one of these two
# LEGACY models, it checks the posterior median of alpha[1] and negates
# every quantity the reflection symmetry ties to alpha's sign if that
# median is negative. The relabelling is exact -- a genuine symmetry of
# the posterior -- not an approximation.
#
# The current stable/ou programs (0.4.2b+) instead fold the sign into
# what's REPORTED (each .stan file's header, "IDENTIFICATION:
# ORIENTATION FOLD"): alpha and theta are multiplied by
# orientation_sign(alpha_raw) in transformed parameters, so alpha[1] >= 0
# always in the reported draws regardless of which raw-space basin a
# chain occupies. This is NOT the same claim as 0.4.2's (reverted) hard
# constraint made -- nothing is excluded from the parameter space, and
# the RAW parameters (alpha_raw itself, and the theta-side raws that
# share its sign) remain genuinely sign-ambiguous; only the reported
# alpha/theta are fixed. There is nothing left for bilatr_orient() to fix
# in the REPORTED quantities for stable/ou, so .bilatr_flip_variables()
# returns character(0) for them and every extract_*() call site skips
# bilatr_orient() rather than calling it as a no-op. This file, and
# bilatr_orient() itself, are retained only to read/re-derive CmdStan
# output produced before 0.4.2 (registered as stable_soft_anchor/
# ou_soft_anchor), and are slated for removal once those runs are gone.
#
# .bilatr_sign_tied_names() below is the single source of truth for BOTH
# this file's .bilatr_flip_variables() (legacy models: everything tied to
# alpha's sign, raw and reported together, since nothing else fixes it)
# AND R/diagnose_convergence.R's .bilatr_sign_ambiguous_raw_names()
# (every model, current and legacy alike: just the RAW parameters, which
# are sign-ambiguous regardless of whether something downstream folds or
# relabels them) -- split, not duplicated, so the two lists cannot drift
# apart.

#' Raw and derived variable names the alpha/theta reflection symmetry
#' ties to alpha's sign, per registered model
#'
#' `raw` is the set of actually-sampled parameters whose sign is
#' ambiguous under the symmetry (present, and sign-ambiguous, in every
#' registered model regardless of what -- if anything -- fixes the
#' REPORTED quantities downstream). `derived` is the set of transformed
#' parameters built from them that inherit the same ambiguity UNLESS
#' something corrects it: for the legacy `_soft_anchor` models, nothing
#' does, so `derived` is exactly what still needs `bilatr_orient()`; for
#' `stable`/`ou`, the orientation fold corrects `derived` automatically
#' (see each `.stan` file's header, "IDENTIFICATION: ORIENTATION FOLD"),
#' so it is never combined with `raw` for them (see
#' [.bilatr_flip_variables()]).
#'
#' @param stan_model Name registered in `.bilatr_stan_models`, or a
#'   recognized pre-0.4.0 alias (see [.canonical_stan_model()], called
#'   here first).
#' @return A list with elements `raw` and `derived`, each a character
#'   vector of variable base names (matched against `posterior` draws
#'   column names via [.bilatr_match_draws_columns()]); both
#'   `character(0)` for any model with no such symmetry.
#' @keywords internal
.bilatr_sign_tied_names <- function(stan_model) {
  stan_model <- .canonical_stan_model(stan_model)
  switch(
    stan_model,
    stable = ,
    stable_soft_anchor = list(
      raw = c("alpha_raw", "z_theta0", "theta_raw"),
      derived = c("alpha", "theta", "theta0")
    ),
    ou = ,
    ou_soft_anchor = list(
      raw = c("alpha_raw", "mu_dyad_raw", "theta_raw"),
      derived = c("alpha", "theta", "mu_dyad")
    ),
    list(raw = character(0), derived = character(0))
  )
}

#' Variables the legacy `stable_soft_anchor`/`ou_soft_anchor` reflection
#' symmetry ties to alpha's sign
#'
#' `character(0)` for `stable`/`ou` (0.4.2b+): those programs fold
#' `alpha[1]`'s sign into the REPORTED `alpha`/`theta` (each `.stan`
#' file's header, "IDENTIFICATION: ORIENTATION FOLD"), so nothing is left
#' for post-hoc relabelling to do to the reported quantities -- see
#' [.bilatr_sign_tied_names()] for why their RAW parameters are still
#' sign-ambiguous regardless (handled separately, in
#' `R/diagnose_convergence.R`, not here). Only the retired
#' `stable_soft_anchor`/`ou_soft_anchor` programs (registered
#' `status = "legacy"`; what `stable`/`ou` were named 0.4.0-0.4.1) still
#' need `bilatr_orient()`, since their free `sum_to_zero_vector` alpha
#' with only a soft sign anchor cannot guarantee a chain lands in the
#' `alpha[1] > 0` basin, and nothing folds the sign for them.
#'
#' @param stan_model Name registered in `.bilatr_stan_models`, or a
#'   recognized pre-0.4.0 alias (see [.canonical_stan_model()], called
#'   here first -- an unrecognized name errors rather than silently
#'   returning `character(0)`, per B1).
#' @return Character vector of variable base names (matched against
#'   `posterior` draws column names via `.bilatr_match_draws_columns()`).
#' @keywords internal
.bilatr_flip_variables <- function(stan_model) {
  stan_model <- .canonical_stan_model(stan_model)
  if (stan_model %in% c("stable_soft_anchor", "ou_soft_anchor")) {
    tied <- .bilatr_sign_tied_names(stan_model)
    c(tied$raw, tied$derived)
  } else {
    character(0)
  }
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
#' **Deprecated as of 0.4.2**, in favor of doing nothing: `stable`/`ou`
#' fold `alpha[1]`'s sign into the REPORTED `alpha`/`theta` (see each
#' `.stan` file's header, "IDENTIFICATION: ORIENTATION FOLD"), leaving no
#' reflection symmetry in those reported quantities for this function to
#' correct. It is retained only for reading/re-deriving CmdStan output
#' produced by the pre-0.4.2 `stable_soft_anchor`/`ou_soft_anchor`
#' programs, and is slated for removal once those runs are gone -- this
#' is not enforced with [base::.Deprecated()], which would fire warnings
#' on that legitimate legacy use; every current call site (`extract_*()`,
#' `R/diagnose_and_extract.R`) already checks [.bilatr_flip_variables()]
#' first and skips calling this function entirely for `stable`/`ou`.
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
#' @param stan_model Name registered in `.bilatr_stan_models`, or a
#'   recognized pre-0.4.0 alias (see [.canonical_stan_model()]). Only
#'   `"stable_soft_anchor"`/`"ou_soft_anchor"` have a reflection symmetry
#'   to correct (see [.bilatr_flip_variables()]); for any other
#'   registered name (including the current `"stable"`/`"ou"`), `draws`
#'   is returned unmodified (restricted to `variables`, if supplied),
#'   since there is nothing to fix. An unrecognized, non-alias name
#'   errors instead.
#' @param variables Which variables to return. Defaults to every variable
#'   [.bilatr_flip_variables()] lists for `stan_model`; pass a subset
#'   (e.g. `"theta"`) when only that one is needed downstream.
#' @return A `posterior::draws_df` restricted to `variables`, with every
#'   variable in `intersect(variables, .bilatr_flip_variables(stan_model))`
#'   negated if the posterior median of `alpha[1]` in `draws` was
#'   negative.
#' @keywords internal
bilatr_orient <- function(draws, stan_model, variables = NULL) {
  stan_model <- .canonical_stan_model(stan_model)
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
