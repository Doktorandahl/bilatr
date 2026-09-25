# Raw, sampled parameter names that are sign-ambiguous under the
# alpha/theta reflection symmetry, per registered model.
#
# Every registered model's alpha_raw is a free sum_to_zero_vector with no
# fixed element, which leaves an exact reflection symmetry: negating
# alpha_raw together with theta's own raw parameters leaves the
# likelihood, every prior, and the sum_to_zero_vector Jacobian unchanged
# (see each model's .stan file header). `stable`/`ou` (0.4.2b+) fold the
# sign into what's REPORTED instead (each .stan file's header,
# "IDENTIFICATION: ORIENTATION FOLD"): alpha and theta are both multiplied
# by orientation_sign(alpha_raw) in transformed parameters, so alpha[1] >=
# 0 always in the reported draws regardless of which raw-space basin a
# chain occupies, and there is no post-hoc relabelling left to do (0.10.0
# retired the pre-0.4.2 stable_soft_anchor/ou_soft_anchor programs and the
# bilatr_orient() machinery that did that relabelling for them; see
# NEWS.md).
#
# What the fold does NOT fix is the RAW parameters it's built from
# (alpha_raw itself, plus whichever of z_theta0/mu_dyad_raw and theta_raw
# a given model has, and gamma_z for stable_gamma, for a different
# reason -- see .bilatr_sign_tied_names()'s docs): if two chains land in
# opposite raw-space basins, those variables' own cross-chain Rhat is
# meaningless even though everything reported is fine.
# .bilatr_sign_tied_names() below is the single source of truth
# R/diagnose_convergence.R's .bilatr_sign_ambiguous_raw_names() draws on to
# keep those names out of the tiered diagnostics tables (tier = NA)
# entirely, rather than let a meaningless Rhat surface as a false Tier 1
# alarm. Its name and output must not change: the runscripts call
# `bilatr:::.bilatr_sign_ambiguous_raw_names()` directly.

#' Raw, sampled parameter names the alpha/theta reflection symmetry ties
#' to alpha's sign, per registered model
#'
#' Present, and sign-ambiguous, in every registered model regardless of
#' whether the current orientation fold fixes the REPORTED quantities
#' built from them -- it does, for every model registered since 0.10.0,
#' but the raw parameters themselves stay ambiguous.
#'
#' @param stan_model Name registered in `.bilatr_stan_models`; see
#'   [.canonical_stan_model()].
#' @return A list with element `raw`, a character vector of variable base
#'   names (matched against `posterior` draws column names); `character(0)`
#'   for any model with no such symmetry.
#' @keywords internal
.bilatr_sign_tied_names <- function(stan_model) {
  stan_model <- .canonical_stan_model(stan_model)
  switch(
    stan_model,
    stable = list(raw = c("alpha_raw", "z_theta0", "theta_raw")),
    # stable_gamma (0.7.0) shares `stable`'s alpha_raw/z_theta0/theta_raw
    # reflection symmetry unchanged (it is built from a copy of `stable`;
    # see inst/stan/bilatr_alphanorm_gamma.stan's header) -- plus
    # `gamma_z`, which is listed here too even though it is NOT tied to
    # alpha's sign (the bilatr_alphanorm_gamma.stan projection that builds
    # `gamma` from `gamma_z` is invariant to alpha_raw -> -alpha_raw; see
    # that file's header). `gamma_z` is unidentified for a different
    # reason: its projected-out directions (the alpha direction, the
    # all-ones direction, the across-country mean) are pinned only by
    # `gamma_z ~ std_normal()`, not by the likelihood, so its cross-chain
    # Rhat is just as uninformative as a genuinely sign-ambiguous raw
    # parameter's -- listing it here reuses
    # R/diagnose_convergence.R's .bilatr_sign_ambiguous_raw_names()/
    # .classify_bilatr_tier() `tier = NA` exclusion for that different
    # reason, rather than inventing a second exclusion mechanism. `gamma`
    # itself is NOT included -- it is fully identified and must stay in
    # the tiered tables, as Tier 1; see .bilatr_tier1_names in
    # R/diagnose_convergence.R.
    stable_gamma = list(raw = c("alpha_raw", "z_theta0", "theta_raw", "gamma_z")),
    ou = list(raw = c("alpha_raw", "mu_dyad_raw", "theta_raw")),
    list(raw = character(0))
  )
}
