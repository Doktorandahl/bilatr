// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// EXPERIMENTAL variant `alphanorm` (registered in R/model_registry.R,
// status = "experimental"). Not the default model; fit only via
// fit_dyad_ts_dev()/fit_panel_dev(stan_model = "alphanorm").
//
// Motivation: the stable model (bilatr_dirmult_irt.stan) has a residual
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
//   - location is pinned HARD: mu_theta0 is removed entirely (theta0 is
//     sigma_theta0 * z_theta0, mean exactly 0 across the population by
//     construction of z_theta0 ~ std_normal(), not merely softly shrunk
//     toward 0 by a prior)
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
//   - theta0 has population mean exactly 0, by construction (no
//     mu_theta0)
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
// REFLECTION SYMMETRY -- a SECOND, DISCRETE degeneracy, distinct from the
// radial one above and exact regardless of alpha_raw's realized value:
// eta = alpha .* theta - mu_intercept is invariant under the JOINT
// negation alpha -> -alpha, theta -> -theta (i.e. theta0 -> -theta0,
// z_theta0 -> -z_theta0, and theta_raw -> -theta_raw at every t),
// because alpha .* theta is a product of two negations while
// mu_intercept is untouched. Every prior on the flipped quantities
// (alpha_raw ~ std_normal(), z_theta0 ~ std_normal(), theta_raw ~
// std_normal()) is symmetric about 0, so the two mirror modes have
// exactly equal posterior mass. With mass split 50/50 across two modes,
// different chains can land in either: Rhat on alpha/theta becomes
// uninterpretable, and pooled posterior means get pulled toward 0. This
// is also the explanation for the sign difference previously observed
// between this model's output and stable's/ou's: those pin
// alpha[1] = 1, which selects a mode directly; this model removed that
// pin (see "This variant closes the ridge" above) without replacing
// what it was doing for sign identification, and evidently landed in
// the other mode.
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
// Fixed with a SOFT sign anchor on alpha[1], not a hard one:
//   target += log_inv_logit(alpha[1] * inv(anchor_scale))
// (anchor_scale is a NEW data field, default 0.1). This is ~0 when
// alpha[1] is comfortably positive and ~ -|alpha[1]| / anchor_scale when
// negative -- it penalizes SIGN only, not magnitude. At alpha[1] ~ 0.786
// (this model's current fit, in the wrong-sign mode) that's ~7.9 nats of
// penalty, ample to make the negative-alpha[1] mode posterior-negligible
// without pulling alpha[1]'s estimated magnitude toward any particular
// value. Two things deliberately NOT done here, both considered and
// rejected:
//   - a hard `alpha[1] <lower=0>` constraint: a boundary the sampler
//     must approach whenever the true posterior mass sits near 0, which
//     this one plausibly does
//   - an informative location prior like alpha[1] ~ normal(0.8, 0.3):
//     confounds sign-breaking with a magnitude belief, and needs
//     re-tuning every time the action-class coding scheme changes
// Also deliberately not done: a second anchor on a hostile class. The
// symmetry is a discrete two-element group; one constraint removes it
// entirely. A second constraint would cut a region of parameter space
// unrelated to this symmetry and would distort alpha, the same failure
// mode as the two hard alpha anchors tried previously (see stable's
// history / this model's own "This variant closes the ridge" section).
//
// ORIENTATION: positive alpha[1] means higher theta corresponds to
// better (less hostile) relations at the reference/neutral action class
// -- matching stable/ou and the package's stated quantity (bilateral
// relationship quality). Runs from before this anchor was added may be
// sign-flipped relative to runs after it. To compare them:
//   FLIP sign:  alpha, theta, theta0, z_theta0, theta_raw
//   UNCHANGED:  mu_intercept, phi, sigma_theta0, process_noise,
//               mu_log_noise, sigma_log_noise, mu_log_phi, sigma_log_phi
// mu_intercept does not flip: alpha .* theta is invariant under the
// joint negation, so only the theta-side and alpha-side terms change;
// every scale/dispersion parameter is a positive quantity, not a
// location, so none of them flip either.
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
                        vector alpha,
                        vector dyad_weight,
                        vector period_weight,
                        vector action_weight) {
    real lp = 0;
    for (d in start:end) {
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_intercept;
          vector[A] p = softmax(eta);
          // action_weight rescales concentration per action type
          vector[A] conc = phi[d] * (action_weight .* p);
          // dyad_weight/period_weight scale each dyad-period's contribution
          lp += dyad_weight[d] * period_weight[t] *
            dirichlet_multinomial_lpmf(Y[d, t] | conc);
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
                            vector alpha, real dyad_weight_d,
                            real period_weight_t, vector action_weight) {
    if (obs_dt == 0) {
      return 0;
    }
    vector[A] eta = alpha .* rep_vector(theta_dt, A) - mu_intercept;
    vector[A] p = softmax(eta);
    vector[A] conc = phi_d * (action_weight .* p);
    return dyad_weight_d * period_weight_t * dirichlet_multinomial_lpmf(y_dt | conc);
  }
}
data {
  int<lower=1> T;                            // number of time points
  int<lower=1> D;                            // number of dyads
  int<lower=1> A;                            // number of action types
  int<lower=1> C;                            // reduce_sum grainsize
  array[D, T] int<lower=0, upper=1> is_obs;  // observed indicator
  array[D, T, A] int<lower=0> Y;             // event counts
  vector<lower=0>[D] dyad_weight;            // per-dyad reweighting, default 1s
  vector<lower=0>[T] period_weight;          // per-period reweighting, default 1s
  vector<lower=0>[A] action_weight;          // per-action-type reweighting, default 1s
  int<lower=0, upper=1> compute_log_lik;     // 1 = also compute per-dyad-period
                                              // log_lik in generated quantities
                                              // (D x T x draws; default 0/off)
  real<lower=0> anchor_scale;                // soft sign-anchor scale on alpha[1];
                                              // default 0.1 (see header,
                                              // "REFLECTION SYMMETRY")
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
  // degeneracy this creates.
  vector[A] alpha = alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  vector[D] theta0 = sigma_theta0 * z_theta0;

  // Per-dyad process noise, now a RATIO to sigma_theta0 (non-centered
  // lognormal hierarchy on the ratio, not on an absolute theta-unit
  // quantity -- see header).
  vector<lower=0>[D] process_noise =
    sigma_theta0 * exp(mu_log_noise + sigma_log_noise * log_process_noise_raw);

  array[D, T] real theta;
  for (d in 1:D) {
    theta[d, 1] = theta0[d] + process_noise[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = theta[d, t - 1] + process_noise[d] * theta_raw[d, t];
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
  // direction of alpha_raw (see header, "RADIAL DEGENERACY")
  alpha_raw ~ std_normal();

  // soft sign anchor: breaks the alpha/theta reflection symmetry by
  // penalizing alpha[1] < 0, not by constraining or centering it (see
  // header, "REFLECTION SYMMETRY")
  target += log_inv_logit(alpha[1] * inv(anchor_scale));

  // likelihood, chunked via reduce_sum
  array[D] int dyad_seq = linspaced_int_array(D, 1, D);
  target += reduce_sum(partial_log_lik, dyad_seq, C,
                        T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                        dyad_weight, period_weight, action_weight);
}
generated quantities {
  array[compute_log_lik ? D : 0, compute_log_lik ? T : 0] real log_lik;

  if (compute_log_lik) {
    for (d in 1:D) {
      for (t in 1:T) {
        log_lik[d, t] = dyad_period_log_lik(
          is_obs[d, t], Y[d, t], theta[d, t], A, mu_intercept, phi[d], alpha,
          dyad_weight[d], period_weight[t], action_weight
        );
      }
    }
  }
}
