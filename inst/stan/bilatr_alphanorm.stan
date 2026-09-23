// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// Registered as `stable` (R/model_registry.R, status = "stable") since
// 0.4.0, promoted from the experimental `alphanorm` variant it was
// developed under -- see NEWS.md. This is now the model
// fit_dyad_ts()/fit_panel() fit by default; the file itself is
// unchanged/unrenamed across that promotion, only the registry entry
// pointing at it. NOTE: every bare "stable"/"ou" reference in the
// comparative/historical discussion BELOW this point (written before the
// promotion) means the ORIGINAL stable/ou models this one replaced,
// retired to inst/stan/legacy/bilatr_dirmult_irt_pre_0.4.0.stan and
// bilatr_ou_pre_0.4.0.stan respectively -- not this file.
//
// 0.4.6: the dyad_weight/period_weight/action_weight likelihood-weighting
// data fields are retired from this file (and from partial_log_lik() /
// dyad_period_log_lik() / assemble_stan_data()) -- never used in
// production, and action_weight in particular made the Dirichlet-
// multinomial concentration depend on theta, complicating the
// hand-differentiated forward filter built on top of this likelihood.
// At unit weights (the only weights ever used) this is an exact,
// bit-identical no-op; see NEWS.md. The legacy stable_soft_anchor
// program (inst/stan/legacy/bilatr_stable_soft_anchor.stan) still
// declares and applies all three, unchanged, since it must keep
// matching the fits that produced its output.
//
// 0.4.2: alpha_raw's sum_to_zero_vector[A] (with only a SOFT sign anchor
// on alpha[1]) was briefly replaced with a hand-built construction whose
// first element was positive BY DECLARATION (real<lower=0> alpha_raw_1).
// A production cluster run falsified that approach: a hard constraint
// removes one of the two reflection-symmetric modes from the parameter
// space, but not the likelihood barrier between them, so a chain that
// starts in the excluded orientation cannot cross to the feasible one --
// it slides to the constraint boundary and parks there instead (observed:
// alpha_raw_1 = 0.00046, 5-95% interval [0.000025, 0.00127], wandering on
// the log scale, lp__ materially worse than a healthy chain; confirmed to
// be the mirrored solution projected onto alpha[1] = 0 and re-fitted
// around that boundary, not a repairable mirror image -- correlation
// 0.994 with the healthy chain's negation, but NOT identical to it). See
// "IDENTIFICATION: ORIENTATION FOLD" below for what replaced it (0.4.2b)
// and why no hard constraint of any kind can work here. The soft-anchor
// version of this file (with its `anchor_scale` data field, the version
// before EITHER of these) is retired to
// inst/stan/legacy/bilatr_stable_soft_anchor.stan, registered as
// `stable_soft_anchor` -- fits made under this file are not
// parameter-comparable to fits made under that one (see NEWS.md).
//
// Motivation: the model this replaced (retired to
// inst/stan/legacy/bilatr_dirmult_irt_pre_0.4.0.stan) had a residual
// affine ridge in its identification. With alpha[1] = 1 and
// mu_intercept[1] = 0, the softmax only sees
// eta_a - eta_1 = (alpha_a - 1) * theta - mu_a. Writing A_a = alpha_a - 1,
// the map theta' = (theta - b)/c, A'_a = c*A_a, mu'_a = mu_a - A_a*b
// leaves the likelihood exactly invariant and preserves both constraints.
// Location and scale are therefore broken only softly, by
// mu_theta0 ~ normal(0,1) and alpha_raw ~ std_normal() -- this inflates
// posterior intervals and is a likely contributor to treedepth pathology
// in the stable model, and is part of why cross-dyad theta levels are not
// reliably comparable.
//
// This variant closes the ridge:
//   - location has no separate free parameter: mu_theta0 is removed
//     entirely (theta0 is sigma_theta0 * z_theta0). This removes the
//     ADDITIONAL soft anchor mu_theta0 ~ normal(0,1) stable puts on top
//     of z_theta0's own prior -- but z_theta0's own population mean is
//     itself still only softly (if stiffly) pinned toward 0 by
//     z_theta0 ~ std_normal(), not pinned exactly; see "Identification"
//     below for the correction to this point and why it doesn't
//     practically matter
//   - scale is pinned on the ALPHA side, not the theta side: alpha_raw is
//     a sum_to_zero_vector[A] (mean exactly 0 by construction, so no
//     mu_intercept-style location dof leaks into alpha), and alpha is
//     alpha_raw normalized to RMS 1 in transformed parameters (RMS == SD
//     here because alpha_raw sums to zero). This makes theta's unit a
//     property of the coding scheme (the discrimination profile across
//     action classes), not of whatever dyad sample happens to be fit --
//     which matters because the dyad-set restriction is still being
//     varied across runs; under stable's theta-side normalization
//     (sigma_theta0 as the scale anchor), fits on different dyad samples
//     are not on a common scale, but under alpha-side normalization they
//     are (the action-class coding scheme doesn't change between runs).
//   - mu_intercept is also a sum_to_zero_vector[A] rather than
//     fixed-first-to-0, for the same reason: no residual location dof
//     hiding in the softmax level-shift.
//
// Identification:
//   - alpha has RMS (== population SD, since sum-to-zero) exactly 1, by
//     construction in transformed parameters -- not softly shrunk by a
//     prior
//   - mu_intercept sums to exactly 0, by construction (sum_to_zero_vector)
//   - theta0's population mean is softly but STIFFLY pinned toward 0 --
//     NOT exact, and this is a correction to earlier text in this
//     header (and to this variant's registry description): z_theta0 is
//     a plain vector[D] with z_theta0 ~ std_normal(), a prior, not a
//     sum-to-zero constraint, so a common shift b in z_theta0 costs
//     D * b^2 / 2 in log density rather than being forbidden outright.
//     The pin is nonetheless very stiff in practice (D is in the
//     thousands for this project's dyad sets), so the mean is
//     effectively 0 for any purpose that matters -- but "effectively"
//     is doing real work in that sentence, and it belongs in a
//     different category than alpha's/mu_intercept's genuinely exact
//     constraints above. sum_to_zero_vector[D] z_theta0 would make it
//     exact too, if that distinction ever mattered enough to act on.
//   - no dyad-specific intercept, as in stable: cross-dyad level
//     differences are forced into theta via mu_intercept
//
// PRIOR UNITS HAVE CHANGED, NOT JUST BEEN PORTED. This has NOT yet been
// checked against prior predictive simulation:
//   - sigma_theta0 ~ normal(0, 2): under stable this was
//     sigma_theta0 ~ normal(0, 0.5), interpreted loosely (theta0's scale
//     was also softly anchored by mu_theta0 ~ normal(0,1), so
//     sigma_theta0 alone didn't fully determine theta's units). Here,
//     with alpha carrying the scale anchor instead, sigma_theta0 is
//     GENUINELY the cross-dyad SD of theta0 in the units alpha defines --
//     a materially different and more consequential prior than in
//     stable. The normal(0, 2) here is a placeholder guess, not a
//     recalibrated choice.
//   - mu_log_noise ~ normal(log(0.2), 0.5) is UNCHANGED IN TEXT from
//     stable but its MEANING has changed: process_noise here is built as
//     sigma_theta0 * exp(mu_log_noise + sigma_log_noise * ...), i.e. it
//     is now a RATIO to sigma_theta0 (innovation SD as a fraction of the
//     cross-dyad SD), not an absolute theta-unit quantity. Read as
//     "annual innovation SD is ~20% of the cross-dyad SD" -- the
//     between/within ratio that actually drives cross-dyad rank
//     inversion under detente. This reinterpretation is deliberate (see
//     task rationale), but the specific normal(log(0.2), 0.5) prior on
//     that ratio has not been separately justified/calibrated; it is
//     carried over textually from stable's absolute-units prior only
//     because 0.2 happens to be a plausible ratio too.
//
// RADIAL DEGENERACY: hard-normalizing alpha in transformed parameters
// (alpha = alpha_raw * sqrt(A / dot_self(alpha_raw))) leaves
// ||alpha_raw|| itself unidentified -- only its direction matters, since
// any rescaling of alpha_raw is undone by the normalization. This makes
// alpha_raw ~ std_normal() load-bearing (it is what identifies the radial
// component), not merely regularizing, and dot_self(alpha_raw) near zero
// is a numerical hazard (division blows up). This can produce a mild
// funnel along the radial direction. The diagnostic script for this
// variant should report divergence count AND the posterior of
// dot_self(alpha_raw) (it should concentrate near A - 1, i.e. near where
// std_normal() puts most of a sum-to-zero (A-1)-dimensional vector's
// squared norm -- not near 0).
//
// ESCAPE HATCH if the radial degeneracy misbehaves in practice: keep
// alpha[1] = 1 (stable's fixed-reference identification) and pin the
// scale instead via sigma_theta0 ~ lognormal(0, 0.2) (a tight prior
// directly on theta's cross-dyad SD, playing the scale-anchor role alpha
// plays here). Not implemented; recorded here as the documented fallback
// per task instructions, left for a future decision.
//
// alpha[1] IS the reference/neutral action class here, not an arbitrary
// index: assemble_stan_data() (via grouped_events_to_dyad_period()'s
// `reference_category` argument) already reorders that class to be
// first in the action dimension before this model ever sees the data
// (verified directly: with grouping_var = "BilatrClass2",
// reference_category = 2, `attr(stan_data, "event_classes")[1]` is
// `"2"`), so no separate index needs to be threaded through as new
// data -- `alpha[1]` already IS the anchor position.
//
// IDENTIFICATION: ORIENTATION FOLD (0.4.2b). eta = alpha .* theta -
// mu_intercept is invariant under the JOINT negation alpha -> -alpha,
// theta -> -theta (i.e. theta0 -> -theta0, z_theta0 -> -z_theta0, and
// theta_raw -> -theta_raw at every t), because alpha .* theta is a
// product of two negations while mu_intercept is untouched. With a free
// sum_to_zero_vector[A] alpha_raw, every prior on the flipped quantities
// is symmetric about 0, so the two mirror modes carry exactly equal
// posterior mass, and whichever one a chain's init happens to land in is
// the one it reports -- making cross-chain Rhat on alpha/theta
// uninterpretable for independently-initialized chains. That is the
// actual problem; two different fixes for it have been tried and
// rejected before this one (both recorded here so neither is tried
// again):
//   - a SOFT sign anchor (target += log_inv_logit(alpha[1] *
//     inv(anchor_scale)), the pre-0.4.2 approach, retired to
//     inst/stan/legacy/bilatr_stable_soft_anchor.stan / registered as
//     `stable_soft_anchor`): makes the target correctly specified but
//     cannot move a chain across the likelihood barrier (thousands of
//     nats) separating the two modes once warmup has landed it in one.
//   - a HARD constraint (real<lower=0> alpha_raw_1, the 0.4.2 approach,
//     briefly shipped and reverted): removes one mode from the parameter
//     space entirely, but not the barrier -- a chain that would have
//     landed in the excluded mode instead slides to the constraint
//     boundary and parks there (see the top-of-file note above for the
//     cluster numbers). A constraint identifies by removing ambiguity;
//     it cannot move a chain across a barrier. Constraining a different
//     element, the largest-magnitude element, or strengthening the soft
//     anchor all have the identical failure mode -- they only relocate
//     where a trapped chain parks, they do not free it. Do not try any
//     of them.
//
// The fix that actually works changes what is REPORTED, not what is
// REACHABLE. Let s = sign(alpha_raw[1]) (see orientation_sign(), in the
// functions block above) and multiply BOTH alpha and theta by it.
// Because the likelihood depends on those two only through the
// elementwise product alpha .* theta, and both flip together, the
// likelihood is completely unchanged -- so a chain in EITHER raw-space
// basin reports identical alpha and theta, with alpha[1] >= 0 always. No
// barrier needs crossing, no region is excluded, nothing is truncated.
// Two properties matter, and must survive any future edit to this file:
//   - The target stays smooth. The discontinuity is in the reported
//     transformed parameters, not in the density HMC differentiates.
//     Crossing alpha_raw[1] = 0 flips s, alpha and theta simultaneously,
//     leaving alpha .* theta -- and hence every term of the likelihood --
//     continuous and differentiable; the priors are on the raw
//     parameters (alpha_raw, z_theta0, theta_raw) and are untouched by
//     s. There is no Jacobian adjustment to make, because this is not a
//     change of variables: it is the same model, reported under a
//     canonical labelling.
//   - It is exact, not approximate. Both raw-space modes map to the same
//     reported values, so the pushforward posterior on (alpha, theta) is
//     correct whichever basin a chain occupies.
// orientation_sign() is a functions-block helper, deliberately NOT cached
// into a local `real s` here: any bare top-level declaration in
// transformed parameters is written to the output CSV, and a raw sign
// scalar would show the same meaningless cross-chain Rhat this fold
// exists to remove (see R/diagnose_convergence.R's
// .bilatr_sign_ambiguous_raw_names() / .classify_bilatr_tier() for where
// that exclusion is actually enforced, for alpha_raw/z_theta0/theta_raw
// themselves).
//
// alpha_raw ~ std_normal() is on the free parameter directly (no
// constructed-variable Jacobian to reason about, unlike the reverted
// 0.4.2 approach) -- see "RADIAL DEGENERACY" above for why it is
// load-bearing, not merely regularizing.
//
// ORIENTATION: positive alpha[1] means higher theta corresponds to
// better (less hostile) relations at the reference/neutral action class,
// matching stable/ou and the package's stated quantity (bilateral
// relationship quality) -- and this holds for every fit made under this
// file via the fold above, not just on average across many. `alpha`,
// `theta`, `mu_intercept`, `phi`, and every scale/dispersion parameter
// are always in their final orientation as REPORTED;
// bilatr_orient()/.bilatr_flip_variables() return `character(0)` for
// `stable`/`ou` accordingly (R/orient.R) -- there is nothing left for
// post-hoc relabeling to do to the reported quantities. The RAW
// parameters the fold consumes (alpha_raw, z_theta0, theta_raw) remain
// genuinely sign-ambiguous themselves: if two chains land in opposite
// raw-space basins, those three variables' own cross-chain Rhat is
// meaningless even though everything reported is fine -- this is why
// they are excluded from the tiered diagnostics tables rather than left
// to be (mis)read as a convergence failure (R/diagnose_convergence.R).
// Only fits made under the retired `stable_soft_anchor`/`ou_soft_anchor`
// programs still need bilatr_orient()'s post-hoc relabeling of the
// reported quantities themselves.
functions {
  // <<< BEGIN GENERATED partial_log_lik (source: inst/stan/include/partial_log_lik.stanfunctions) >>>
  // Do not hand-edit between these markers -- edit the source file
  // above and rerun `Rscript data-raw/sync_stan_functions.R`
  // (checked by tests/testthat/test_stan_includes.R).
  // Shared reduce_sum likelihood for all bilatr Stan model variants.
  //
  // Indifferent to how `alpha` and `mu_intercept` were constructed upstream
  // (fixed-reference vs. sum-to-zero, static vs. OU-derived theta, ...): it
  // only consumes the already-built `theta`/`alpha`/`mu_intercept`/`phi`
  // vectors.
  //
  // This is the single canonical source. It is NOT included at compile/
  // sample time via Stan's `#include` (cmdstanr breaks `#include` resolution
  // at $sample()-time whenever include_paths contains a space -- see
  // https://github.com/stan-dev/cmdstanr/issues/820 -- which bites this
  // project's own devtools::load_all() working tree). Instead,
  // `data-raw/sync_stan_functions.R` splices this file's contents verbatim
  // into a marker-delimited block in each registered model's `.stan` file;
  // see `R/stan_includes.R`. Edit only this file, then rerun the sync
  // script -- do not hand-edit the generated blocks, and do not fork this
  // file per variant.
  real partial_log_lik(array[] int slice_d,
                        int start, int end,
                        int T, int A,
                        array[,] int is_obs,
                        array[,,] int Y,
                        array[,] real theta,
                        vector mu_intercept,
                        vector phi,
                        vector alpha) {
    real lp = 0;
    for (d in start:end) {
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_intercept;
          vector[A] p = softmax(eta);
          vector[A] conc = phi[d] * p;
          lp += dirichlet_multinomial_lpmf(Y[d, t] | conc);
        }
      }
    }
    return lp;
  }

  // Country-offset variant of partial_log_lik(), added alongside it (0.7.0,
  // bilatr_alphanorm_gamma.stan / registered `stable_gamma`) -- not a fork,
  // per this file's own header: partial_log_lik() itself is untouched, and
  // every other registered program (stable, ou, and the legacy soft-anchor
  // pair) gets this spliced in too via data-raw/sync_stan_functions.R, as
  // an unused extra function, so their .stan text changes (one-time
  // recompile) but their data/parameter blocks -- and hence fits already on
  // disk -- are untouched.
  //
  // g_d, the per-dyad country offset, is built ONCE PER DYAD, outside the
  // `t` loop -- not stored anywhere (see bilatr_alphanorm_gamma.stan's
  // header, "Part 0b": a per-dyad-period or per-dyad matrix in transformed
  // parameters would be written to the output CSV at prohibitive cost).
  // w_send[d] == 1 covers directed data (g_d = gamma[, ctry_a[d]] exactly,
  // ctry_b[d] unused) and the directed-as-degenerate-undirected case alike;
  // any other value mixes ctry_a's and ctry_b's columns by event share (see
  // bilatr_alphanorm_gamma.stan's header for the directed/undirected
  // design this implements).
  real partial_log_lik_offset(array[] int slice_d,
                               int start, int end,
                               int T, int A,
                               array[,] int is_obs,
                               array[,,] int Y,
                               array[,] real theta,
                               vector mu_intercept,
                               vector phi,
                               vector alpha,
                               matrix gamma,
                               array[] int ctry_a,
                               array[] int ctry_b,
                               vector w_send) {
    real lp = 0;
    for (d in start:end) {
      // once per dyad, not per dyad-period
      vector[A] g = w_send[d] == 1.0
                    ? gamma[, ctry_a[d]]
                    : w_send[d] * gamma[, ctry_a[d]] + (1 - w_send[d]) * gamma[, ctry_b[d]];
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_intercept - g;
          vector[A] p = softmax(eta);
          vector[A] conc = phi[d] * p;
          lp += dirichlet_multinomial_lpmf(Y[d, t] | conc);
        }
      }
    }
    return lp;
  }
  // <<< END GENERATED partial_log_lik >>>

  // Per-dyad-period log-likelihood, not reduced/summed -- used only by the
  // compute_log_lik generated quantities block below. Must stay
  // numerically identical to the per-cell term inside partial_log_lik
  // above (same eta/softmax/conc/dirichlet_multinomial_lpmf). Duplicated
  // rather than shared via the GENERATED mechanism because it is
  // GQ-only, needed by the three experimental variants only (not
  // stable), and small enough that hand-verified identity across the
  // three files is lower-risk than extending the sync tooling for it.
  real dyad_period_log_lik(int obs_dt, array[] int y_dt, real theta_dt,
                            int A, vector mu_intercept, real phi_d,
                            vector alpha) {
    if (obs_dt == 0) {
      return 0;
    }
    vector[A] eta = alpha .* rep_vector(theta_dt, A) - mu_intercept;
    vector[A] p = softmax(eta);
    vector[A] conc = phi_d * p;
    return dirichlet_multinomial_lpmf(y_dt | conc);
  }

  // +1/-1 orientation of a draw, from the reference class's sign. Applied
  // to BOTH alpha and theta below; see header, "IDENTIFICATION:
  // ORIENTATION FOLD". Hand-duplicated (not spliced via the GENERATED
  // mechanism) in bilatr_alphanorm_ou.stan -- identical there, same
  // reasoning as dyad_period_log_lik() above (tiny, not worth extending
  // the sync tooling for).
  real orientation_sign(vector alpha_raw) {
    return alpha_raw[1] >= 0 ? 1.0 : -1.0;
  }
}
data {
  int<lower=1> T;                            // number of time points
  int<lower=1> D;                            // number of dyads
  int<lower=1> A;                            // number of action types
  int<lower=1> C;                            // reduce_sum grainsize
  array[D, T] int<lower=0, upper=1> is_obs;  // observed indicator
  array[D, T, A] int<lower=0> Y;             // event counts
  int<lower=0, upper=1> compute_log_lik;     // 1 = also compute per-dyad-period
                                              // log_lik in generated quantities
                                              // (D x T x draws; default 0/off)
  int<lower=0, upper=1> prior_only;          // 1 = skip the likelihood entirely
                                              // (fit the prior only; see
                                              // alpha_prior_moments())
  int<lower=0, upper=1> compute_theta_filtered; // 1 = also compute
                                              // theta_filtered/theta_filtered_sd
                                              // in generated quantities, for the
                                              // dyads in filter_dyads
  int<lower=0> n_filter_dyads;               // length of filter_dyads; 0 if
                                              // compute_theta_filtered is 0
  array[n_filter_dyads] int<lower=1, upper=D> filter_dyads; // which dyads to filter
}
parameters {
  // Latent states per dyad
  array[D, T] real theta_raw;
  sum_to_zero_vector[A] mu_intercept;   // softmax level-shift; sums to 0 exactly
  sum_to_zero_vector[A] alpha_raw;      // pre-normalization discrimination; sums to 0 exactly
  real<lower=0> sigma_theta0;
  vector[D] z_theta0;                   // no mu_theta0: theta0 has mean exactly 0

  // Hierarchical process/dispersion parameters
  vector[D] log_process_noise_raw;  // non-centered: process_noise built in TP
  real mu_log_noise;                // NOTE: now a ratio to sigma_theta0, see header
  real<lower=0> sigma_log_noise;
  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;
}
transformed parameters {
  // Normalize alpha to RMS (== population SD, sum-to-zero) 1. This is the
  // scale anchor for the whole model; see header for the radial
  // degeneracy this creates. Sign is NOT fixed here -- see
  // orientation_sign() below and header, "IDENTIFICATION: ORIENTATION
  // FOLD".
  vector[A] alpha = orientation_sign(alpha_raw)
                    * alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  // theta0 folds the same sign as alpha, so alpha .* theta (all the
  // likelihood ever sees) is invariant to which raw-space basin the
  // sampler is in -- see header.
  vector[D] theta0 = orientation_sign(alpha_raw) * sigma_theta0 * z_theta0;

  // Per-dyad process noise, now a RATIO to sigma_theta0 (non-centered
  // lognormal hierarchy on the ratio, not on an absolute theta-unit
  // quantity -- see header).
  vector<lower=0>[D] process_noise =
    sigma_theta0 * exp(mu_log_noise + sigma_log_noise * log_process_noise_raw);

  array[D, T] real theta;
  for (d in 1:D) {
    theta[d, 1] = theta0[d] + orientation_sign(alpha_raw) * process_noise[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = theta[d, t - 1] + orientation_sign(alpha_raw) * process_noise[d] * theta_raw[d, t];
    }
  }
}
model {
  // hyperpriors for pooling (lognormal hierarchy on the process_noise /
  // sigma_theta0 ratio)
  mu_log_noise ~ normal(log(0.2), 0.5);
  sigma_log_noise ~ normal(0, 0.5);
  log_process_noise_raw ~ std_normal();

  mu_log_phi ~ normal(0, 1);
  sigma_log_phi ~ normal(0, 1);
  phi ~ lognormal(mu_log_phi, sigma_log_phi);

  // theta random-walk innovations
  for (d in 1:D) {
    for (t in 1:T) {
      theta_raw[d, t] ~ std_normal();
    }
  }

  sigma_theta0 ~ normal(0, 2);
  z_theta0 ~ std_normal();

  mu_intercept ~ std_normal();

  // load-bearing, not merely regularizing: identifies the radial
  // direction of alpha_raw (see header, "RADIAL DEGENERACY"). The
  // orientation fold (transformed parameters, above) leaves this prior
  // untouched -- it is symmetric about 0 and orientation_sign() is
  // applied only to the REPORTED alpha/theta0/theta, not to alpha_raw
  // itself.
  alpha_raw ~ std_normal();

  // likelihood, chunked via reduce_sum -- skipped entirely if prior_only,
  // gating this exactly the way compute_log_lik gates its own generated
  // quantity below
  if (!prior_only) {
    array[D] int dyad_seq = linspaced_int_array(D, 1, D);
    target += reduce_sum(partial_log_lik, dyad_seq, C,
                          T, A, is_obs, Y, theta, mu_intercept, phi, alpha);
  }
}
generated quantities {
  array[compute_log_lik ? D : 0, compute_log_lik ? T : 0] real log_lik;

  if (compute_log_lik) {
    for (d in 1:D) {
      for (t in 1:T) {
        log_lik[d, t] = dyad_period_log_lik(
          is_obs[d, t], Y[d, t], theta[d, t], A, mu_intercept, phi[d], alpha
        );
      }
    }
  }

  // Forward-filtered theta: conditional on THIS draw's hyperparameters,
  // a Fisher-scoring (West-Harrison linear-Bayes) filter over the
  // observations -- an approximation to p(theta_t | y_1:t), not the exact
  // marginal (the hyperparameters here were themselves fit on all T
  // periods). No autodiff in generated quantities, so the score is
  // hand-derived (see data-raw or dev notes for the derivation): with
  // conc_0 = phi[d] constant in theta (true unconditionally since the
  // weights were retired in 0.4.6), d/dtheta log P(y|conc) reduces to
  // sum_k phi[d]*p_k*(alpha_k - a_bar) * (digamma(y_k + phi[d]*p_k) -
  // digamma(phi[d]*p_k)); the Fisher information uses the standard DM
  // overdispersion correction n*(1+phi)/(n+phi) applied to
  // Var_pi(alpha) = dot_product(p, square(alpha - a_bar)) (see
  // diagnose_category_merges()'s use of the same quantity). No
  // orientation_sign() call anywhere here: alpha, mu_intercept,
  // sigma_theta0, and process_noise are already-oriented quantities, and
  // the initial state (0) is orientation-free, so the recursion's output
  // is automatically on the same oriented scale as theta.
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered;
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered_sd;

  if (compute_theta_filtered) {
    for (i in 1:n_filter_dyads) {
      int d = filter_dyads[i];
      real m = 0;
      real p_var = square(sigma_theta0);
      for (t in 1:T) {
        real m_pred = m;
        real p_pred = p_var + square(process_noise[d]);
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(m_pred, A) - mu_intercept;
          vector[A] p = softmax(eta);
          real a_bar = dot_product(p, alpha);
          int n = sum(Y[d, t]);
          real info = n * (1 + phi[d]) / (n + phi[d]) * dot_product(p, square(alpha - a_bar));
          real g = 0;
          for (k in 1:A) {
            real conc_k = phi[d] * p[k];
            g += conc_k * (alpha[k] - a_bar) * (digamma(Y[d, t, k] + conc_k) - digamma(conc_k));
          }
          p_var = 1 / (1 / p_pred + info);
          m = m_pred + p_var * g;
        } else {
          m = m_pred;
          p_var = p_pred;
        }
        theta_filtered[i, t] = m;
        theta_filtered_sd[i, t] = sqrt(p_var);
      }
    }
  }
}
