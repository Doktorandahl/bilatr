// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// Registered as `ou` (R/model_registry.R, status = "experimental") since
// 0.4.0, promoted (still experimental) from the `alphanorm_ou` variant it
// was developed under, once `alphanorm` itself was promoted to `stable`
// -- see NEWS.md. Fit via fit_dyad_ts_dev()/fit_panel_dev(stan_model =
// "ou"); the file itself is unchanged/unrenamed across that promotion.
// The ORIGINAL `ou` model this replaces (an OU/AR(1) variant of the
// pre-0.4.0 stable model, sharing ITS alpha[1] = 1 identification, not
// this file's) was retired to
// inst/stan/legacy/bilatr_ou_pre_0.4.0.stan -- every bare "ou" reference
// below this point that isn't clearly about this file means that
// retired model.
//
// 0.4.6: the dyad_weight/period_weight/action_weight likelihood-weighting
// data fields are retired from this file, same change and rationale as
// bilatr_alphanorm.stan's header -- see there for the full explanation.
// The legacy ou_soft_anchor program
// (inst/stan/legacy/bilatr_ou_soft_anchor.stan) still declares and
// applies all three, unchanged.
//
// 0.4.2: alpha_raw's sum_to_zero_vector[A] (with only a SOFT sign anchor
// on alpha[1]) was briefly replaced with the hand-built,
// first-element-positive construction described (and since reverted) in
// bilatr_alphanorm.stan's header -- identical change, same rationale,
// same failure, applied here. 0.4.2b's orientation fold (see
// "IDENTIFICATION: ORIENTATION FOLD" below) replaces it; see
// bilatr_alphanorm.stan's header for the full derivation, the cluster
// numbers that falsified the hard constraint, and why no constraint of
// any kind can fix this. The soft-anchor version of this file is retired
// to inst/stan/legacy/bilatr_ou_soft_anchor.stan, registered as
// `ou_soft_anchor` -- fits made under this file are not
// parameter-comparable to fits made under that one (see NEWS.md).
//
// Combines `alphanorm`'s identification (bilatr_alphanorm.stan, now
// registered as `stable`) with the retired `ou`'s OU/AR(1) dynamics --
// see those two files for the full rationale behind each half; this
// header covers only how the two combine and what's genuinely new in
// the combination.
//
// From `alphanorm` (closing the affine ridge):
//   - alpha_raw and mu_intercept are both sum_to_zero_vector[A]; alpha is
//     alpha_raw normalized to RMS 1 in transformed parameters (same
//     radial-degeneracy caveat as `alphanorm`: alpha_raw ~ std_normal()
//     is load-bearing, not merely regularizing; see that file's header
//     for the diagnostic/escape-hatch discussion, which applies
//     unchanged here)
//   - location has no separate free parameter: there is no mu_theta_bar
//     here (unlike `ou`, which keeps stable's mu_theta0-equivalent
//     location anchor). mu_dyad = sigma_mu * mu_dyad_raw; mu_dyad_raw is
//     a plain vector[D] with mu_dyad_raw ~ std_normal(), so mu_dyad's
//     population mean is softly but STIFFLY pinned toward 0 (a prior, not
//     a constraint -- D is in the thousands for this project's dyad
//     sets, so a common shift b costs D * b^2 / 2 in log density, which
//     makes it effectively 0 in practice but not exactly 0 by
//     construction the way alpha's/mu_intercept's sum_to_zero_vector
//     constraints are; see bilatr_alphanorm.stan's header for the same
//     correction applied there to theta0/z_theta0).
//     sum_to_zero_vector[D] mu_dyad_raw would make it exact, if that
//     distinction ever mattered enough to act on.
//
// From `ou` (restoring force):
//   - theta follows the same OU/AR(1) process with dyad-specific
//     equilibria (mu_dyad) and global persistence rho, stationary
//     initialization at t = 1, and the stationary-sd-first
//     parameterization (process_noise derived from sd_stat and rho, not
//     the reverse) -- all for the reasons given in bilatr_ou.stan's
//     header.
//
// What's NEW in the combination (not just a splice of the two):
//   - sd_stat is made RELATIVE TO sigma_mu, not given its own free
//     hierarchy: sd_stat[d] = sigma_mu * exp(mu_log_sd_stat +
//     sigma_log_sd_stat * log_sd_stat_raw[d]). Because alpha (not theta)
//     carries the scale anchor here (as in `alphanorm`), tying sd_stat to
//     sigma_mu multiplicatively means exp(mu_log_sd_stat) becomes
//     DIRECTLY the within/between-dyad sd ratio -- the single number the
//     whole cross-dyad-comparability problem reduces to -- rather than a
//     ratio that has to be reconstructed post hoc (as in `ou`, where it's
//     sigma_mu / mean(sd_stat) computed in generated quantities).
//
// PRIOR UNITS -- carried over but NOT yet prior-predictive calibrated:
//   - sigma_mu ~ normal(0, 2): same placeholder as `alphanorm`'s
//     sigma_theta0 ~ normal(0, 2) (alpha now carries the scale anchor, so
//     this is genuinely the cross-dyad equilibrium SD in alpha's units --
//     a materially different prior than stable's sigma_theta0, not a
//     recalibrated one).
//   - mu_log_sd_stat ~ normal(log(1), 0.5): centers the within/between
//     ratio at 1 (innovation-scale stationary sd equal to the cross-dyad
//     equilibrium sd), a neutral-seeming but unvalidated starting point.
//   - rho_prior_a/rho_prior_b (default Beta(8, 2)) and
//     sigma_log_sd_stat ~ normal(0, 0.5): same as `ou`, same caveat.
//
// alpha[1] IS the reference/neutral action class here (same mechanism as
// `alphanorm`: assemble_stan_data()'s `reference_category` argument
// already reorders it to be first, so no new index needs to be passed
// as data).
//
// IDENTIFICATION: ORIENTATION FOLD -- inherited from `alphanorm` (see
// bilatr_alphanorm.stan's header for the full derivation, the rejected
// soft-anchor and hard-constraint alternatives, and why neither can work;
// this section only restates the mechanism in this model's own variable
// names, since mu_dyad/mu_dyad_raw take the role theta0/z_theta0 play
// there). Under a free sum_to_zero_vector[A] alpha_raw, eta = alpha .*
// theta - mu_intercept is invariant under the JOINT negation alpha ->
// -alpha, theta -> -theta (i.e. mu_dyad -> -mu_dyad, mu_dyad_raw ->
// -mu_dyad_raw, and theta_raw -> -theta_raw at every t; the OU recursion
// theta[t] = mu_dyad + rho*(theta[t-1] - mu_dyad) + process_noise*
// theta_raw[t] is self-consistent under this joint negation, since rho
// and process_noise are untouched positive/unsigned quantities), with the
// same soft-anchor/hard-constraint failure modes described in
// `alphanorm`'s header.
//
// The fold: s = orientation_sign(alpha_raw) (functions block above)
// multiplies alpha AND theta (via mu_dyad and each dyad's own innovation
// term -- see transformed parameters), so alpha .* theta is unchanged
// regardless of which raw-space basin a chain occupies, and alpha[1] >=
// 0 always in the REPORTED quantities. This is exact and Jacobian-free
// for the same reason given in `alphanorm`'s header: it changes what is
// reported, not the target HMC differentiates, and both raw-space modes
// map to the same reported values. Induction through the OU recursion
// (checked separately from the plain random-walk case, since it is not
// identical -- see the comment at the theta[d, t] build in transformed
// parameters): mu_dyad already carries s, so theta[d, t-1] - mu_dyad[d]
// carries s once theta[d, t-1] does, and multiplying by the unsigned rho
// preserves it, so every theta[d, t] carries the same s as theta[d, 1].
//
// ORIENTATION: positive alpha[1] means higher theta corresponds to
// better (less hostile) relations at the reference/neutral action class
// -- matching stable/ou, holding for every fit made under this file via
// the fold. `alpha`, `theta`, `mu_dyad`, `mu_intercept`, `phi`, and every
// scale/dispersion/ratio quantity (`sigma_mu`, `sd_stat`, `process_noise`,
// `rho`, `within_between_ratio`) are always in their final orientation as
// REPORTED; bilatr_orient()/.bilatr_flip_variables() return `character(0)`
// for `stable`/`ou` accordingly (R/orient.R). The RAW parameters the fold
// consumes (`alpha_raw`, `mu_dyad_raw`, `theta_raw`) remain genuinely
// sign-ambiguous themselves and are excluded from the tiered diagnostics
// tables for exactly that reason (R/diagnose_convergence.R) -- see
// `alphanorm`'s header for the full statement of this. Only fits made
// under the retired `stable_soft_anchor`/`ou_soft_anchor` programs still
// need bilatr_orient()'s post-hoc relabeling of the reported quantities
// themselves.
//
// GENERATED QUANTITIES: within_between_ratio = exp(mu_log_sd_stat)
// directly (no mean() over dyads needed, unlike `ou`), computed
// unconditionally (cheap). log_lik gated behind compute_log_lik (D x T x
// draws, expensive -- default off).
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
  //
  // 0.7.1: g_d is folded into the intercept ONCE PER DYAD (mu_eff =
  // mu_intercept + g_d), not subtracted a second time inside the `t`
  // loop's per-cell `eta` -- g_d is dyad-constant, so re-subtracting it at
  // every observed cell was A extra autodiff nodes per cell for nothing
  // (~3-8% more nodes in eta's part of the graph at production scale, all
  // avoidable). The per-cell body below (`eta = alpha .* theta - mu_eff;
  // ...`) is now the same shape, and the same autodiff cost, as
  // partial_log_lik()'s -- the only added work is A new vars per dyad, not
  // per dyad-period. This makes stable_gamma NOT bit-identical to 0.6.x/
  // pre-0.7.1 output at n_countries > 1 ((x - mu) - g and x - (mu + g)
  // differ in floating-point association) -- but bit-identical at
  // n_countries = 1, where gamma is identically zero and mu_intercept + 0
  // is exact in IEEE 754, which is why the exact-nesting test stays valid
  // evidence for the shared-gamma design after this change.
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
      // fold the country offset into the intercept ONCE PER DYAD: the
      // per-cell expression below is then character-for-character the
      // same shape, and the same autodiff cost, as partial_log_lik()'s.
      vector[A] mu_eff = mu_intercept + (w_send[d] == 1.0
                         ? gamma[, ctry_a[d]]
                         : w_send[d] * gamma[, ctry_a[d]] + (1 - w_send[d]) * gamma[, ctry_b[d]]);
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          vector[A] eta = alpha .* rep_vector(theta[d, t], A) - mu_eff;
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
  // to alpha/mu_dyad/theta below; see header, "IDENTIFICATION:
  // ORIENTATION FOLD". Hand-duplicated (not spliced via the GENERATED
  // mechanism) from bilatr_alphanorm.stan -- identical there, same
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
  real<lower=0> rho_prior_a;                 // Beta(rho_prior_a, rho_prior_b) on
  real<lower=0> rho_prior_b;                 // rho; default 8, 2 (weighted toward
                                              // strong persistence)
}
parameters {
  array[D, T] real theta_raw;
  sum_to_zero_vector[A] mu_intercept;   // softmax level-shift; sums to 0 exactly
  sum_to_zero_vector[A] alpha_raw;      // pre-normalization discrimination; sums to 0 exactly

  real<lower=0> sigma_mu;             // takes the role of stable's sigma_theta0;
                                       // no mu_theta_bar -- location pinned hard
  vector[D] mu_dyad_raw;
  real<lower=0, upper=1> rho;         // global OU/AR(1) persistence

  real mu_log_sd_stat;                 // exp(mu_log_sd_stat) IS the within/between
  real<lower=0> sigma_log_sd_stat;     // sd ratio directly -- see header
  vector[D] log_sd_stat_raw;

  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;
}
transformed parameters {
  // Sign is NOT fixed here -- see orientation_sign() (functions block
  // above) and header, "IDENTIFICATION: ORIENTATION FOLD".
  vector[A] alpha = orientation_sign(alpha_raw)
                    * alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  // mu_dyad folds the same sign as alpha, so alpha .* theta (all the
  // likelihood ever sees) is invariant to which raw-space basin the
  // sampler is in -- see header.
  vector[D] mu_dyad = orientation_sign(alpha_raw) * sigma_mu * mu_dyad_raw;

  // sd_stat relative to sigma_mu: exp(mu_log_sd_stat) is directly the
  // within/between-dyad sd ratio (see header). Unsigned -- a positive
  // scale, not a location; untouched by the fold.
  vector<lower=0>[D] sd_stat =
    sigma_mu * exp(mu_log_sd_stat + sigma_log_sd_stat * log_sd_stat_raw);
  vector<lower=0>[D] process_noise = sd_stat * sqrt(1 - square(rho));

  array[D, T] real theta;
  for (d in 1:D) {
    // stationary initialization at t = 1 (sd_stat, not process_noise).
    // Fold applies to this dyad's own innovation term only -- mu_dyad is
    // already signed above, and rho is an unsigned persistence, so the
    // recursion below preserves the same sign through every t (verified
    // by induction: theta[d, t-1] - mu_dyad[d] carries orientation_sign()
    // once mu_dyad does, and multiplying by the unsigned rho preserves
    // it).
    theta[d, 1] = mu_dyad[d] + orientation_sign(alpha_raw) * sd_stat[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = mu_dyad[d] + rho * (theta[d, t - 1] - mu_dyad[d])
                    + orientation_sign(alpha_raw) * process_noise[d] * theta_raw[d, t];
    }
  }
}
model {
  rho ~ beta(rho_prior_a, rho_prior_b);
  sigma_mu ~ normal(0, 2);
  mu_dyad_raw ~ std_normal();

  mu_log_sd_stat ~ normal(log(1), 0.5);
  sigma_log_sd_stat ~ normal(0, 0.5);
  log_sd_stat_raw ~ std_normal();

  mu_log_phi ~ normal(0, 1);
  sigma_log_phi ~ normal(0, 1);
  phi ~ lognormal(mu_log_phi, sigma_log_phi);

  // theta OU/AR(1) innovations
  for (d in 1:D) {
    for (t in 1:T) {
      theta_raw[d, t] ~ std_normal();
    }
  }

  mu_intercept ~ std_normal();

  // load-bearing, not merely regularizing: identifies the radial
  // direction of alpha_raw (see bilatr_alphanorm.stan header, "RADIAL
  // DEGENERACY"). The orientation fold (transformed parameters, above)
  // leaves this prior untouched -- it is symmetric about 0 and
  // orientation_sign() is applied only to the REPORTED alpha/mu_dyad/
  // theta, not to alpha_raw itself.
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
  // the single number the cross-dyad-comparability problem reduces to
  real within_between_ratio = exp(mu_log_sd_stat);

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

  // Forward-filtered theta -- see bilatr_alphanorm.stan's generated
  // quantities block for the full derivation (identical score/Fisher-
  // information math; only the state-transition step differs here).
  // The prediction step for a mean-reverting OU/AR(1) transition is
  // m_pred = mu_dyad + rho*(m - mu_dyad), p_pred = rho^2*p_var +
  // process_noise^2 -- NOT p_var + process_noise^2 (that unscaled form is
  // only correct for stable, where rho is implicitly 1). Sanity check:
  // at the stationary variance p_var = sd_stat[d]^2 (this filter's own
  // initial value, and process_noise[d]^2 = sd_stat[d]^2*(1-rho^2) from
  // transformed parameters), rho^2*p_var + process_noise[d]^2 =
  // sd_stat[d]^2 exactly -- so an unobserved dyad's filtered variance
  // stays at the stationary level forever, as it should, rather than
  // growing without bound. No orientation_sign() call here either: mu_dyad,
  // sd_stat, process_noise, and rho are already-oriented/unsigned
  // quantities.
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered;
  array[compute_theta_filtered ? n_filter_dyads : 0,
        compute_theta_filtered ? T : 0] real theta_filtered_sd;

  if (compute_theta_filtered) {
    for (i in 1:n_filter_dyads) {
      int d = filter_dyads[i];
      real m = mu_dyad[d];
      real p_var = square(sd_stat[d]);
      for (t in 1:T) {
        real m_pred = mu_dyad[d] + rho * (m - mu_dyad[d]);
        real p_pred = square(rho) * p_var + square(process_noise[d]);
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
