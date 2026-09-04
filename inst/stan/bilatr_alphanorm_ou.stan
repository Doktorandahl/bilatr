// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// EXPERIMENTAL variant `alphanorm_ou` (registered in R/model_registry.R,
// status = "experimental"). Not the default model; fit only via
// fit_dyad_ts_dev()/fit_panel_dev(stan_model = "alphanorm_ou"). Combines
// `alphanorm`'s identification with `ou`'s dynamics -- see those two
// files for the full rationale behind each half; this header covers only
// how the two combine and what's genuinely new in the combination.
//
// From `alphanorm` (closing the affine ridge):
//   - alpha_raw and mu_intercept are both sum_to_zero_vector[A]; alpha is
//     alpha_raw normalized to RMS 1 in transformed parameters (same
//     radial-degeneracy caveat as `alphanorm`: alpha_raw ~ std_normal()
//     is load-bearing, not merely regularizing; see that file's header
//     for the diagnostic/escape-hatch discussion, which applies
//     unchanged here)
//   - location is pinned HARD: there is no mu_theta_bar here (unlike
//     `ou`, which keeps stable's mu_theta0-equivalent location anchor).
//     mu_dyad = sigma_mu * mu_dyad_raw, mean exactly 0 across dyads by
//     construction.
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
  vector[A] alpha = alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  vector[D] mu_dyad = sigma_mu * mu_dyad_raw;

  // sd_stat relative to sigma_mu: exp(mu_log_sd_stat) is directly the
  // within/between-dyad sd ratio (see header).
  vector<lower=0>[D] sd_stat =
    sigma_mu * exp(mu_log_sd_stat + sigma_log_sd_stat * log_sd_stat_raw);
  vector<lower=0>[D] process_noise = sd_stat * sqrt(1 - square(rho));

  array[D, T] real theta;
  for (d in 1:D) {
    // stationary initialization at t = 1 (sd_stat, not process_noise)
    theta[d, 1] = mu_dyad[d] + sd_stat[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = mu_dyad[d] + rho * (theta[d, t - 1] - mu_dyad[d])
                    + process_noise[d] * theta_raw[d, t];
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
  // direction of alpha_raw (see bilatr_alphanorm.stan header,
  // "RADIAL DEGENERACY")
  alpha_raw ~ std_normal();

  // likelihood, chunked via reduce_sum
  array[D] int dyad_seq = linspaced_int_array(D, 1, D);
  target += reduce_sum(partial_log_lik, dyad_seq, C,
                        T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                        dyad_weight, period_weight, action_weight);
}
generated quantities {
  // the single number the cross-dyad-comparability problem reduces to
  real within_between_ratio = exp(mu_log_sd_stat);

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
