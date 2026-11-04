// bilatr model `ou` (experimental): hierarchical dynamic IRT model for
// latent bilateral relations with mean-reverting (Ornstein-Uhlenbeck /
// AR(1)) latent states. Fit with fit_panel(stan_model = "ou").
//
// The measurement model (likelihood, alpha, mu_intercept, phi) and its
// identification are those of the `stable` model; see
// bilatr_alphanorm.stan's header. This header covers the latent process.
//
// DATA
//   As in `stable`, plus rho_prior_a, rho_prior_b: the Beta prior on rho
//   (default Beta(8, 2), weighted toward strong persistence).
//
// LIKELIHOOD
//   As in `stable`: for is_obs[d, t] = 1,
//   Y[d, t] ~ DirichletMultinomial(phi[d] * softmax(alpha * theta[d, t] - mu_intercept)).
//
// LATENT PROCESS (non-centered; the sign factor of the orientation fold,
// below, is left out here)
//   mu_dyad[d]       = sigma_mu * mu_dyad_raw[d]      dyad equilibrium
//   sd_stat[d]       = sigma_mu * exp(mu_log_sd_stat
//                        + sigma_log_sd_stat * log_sd_stat_raw[d])
//   process_noise[d] = sd_stat[d] * sqrt(1 - rho^2)
//   theta[d, 1]      = mu_dyad[d] + sd_stat[d] * theta_raw[d, 1]   (stationary start)
//   theta[d, t]      = mu_dyad[d] + rho * (theta[d, t - 1] - mu_dyad[d])
//                        + process_noise[d] * theta_raw[d, t]
//   theta reverts toward its dyad's equilibrium mu_dyad[d] with a global
//   persistence rho. The stationary SD sd_stat is the primary parameter
//   and process_noise is derived from it, so the process is stationary
//   from t = 1. sd_stat is a multiple of sigma_mu, so
//   exp(mu_log_sd_stat) is directly the ratio of within-dyad to
//   between-dyad SD.
//
// PRIORS
//   theta_raw, mu_dyad_raw, log_sd_stat_raw ~ std_normal()
//   rho               ~ beta(rho_prior_a, rho_prior_b)
//   sigma_mu          ~ half-normal(0, 2)
//   mu_log_sd_stat    ~ normal(log(1), 0.5)   (ratio centered at 1)
//   sigma_log_sd_stat ~ half-normal(0, 0.5)
//   phi, mu_intercept, alpha_raw: as in `stable`.
//   These priors have not been calibrated by prior-predictive simulation.
//
// IDENTIFICATION
//   Scale, radial direction and softmax location: as in `stable` (alpha
//     normalized to RMS 1; alpha_raw ~ std_normal() load-bearing;
//     alpha_raw and mu_intercept sum to zero).
//   Location of theta: mu_dyad has no free mean parameter. Its cross-dyad
//     mean is pinned at 0 by mu_dyad_raw ~ std_normal() -- a stiff prior
//     for large D (a common shift b costs D * b^2 / 2), not an exact
//     constraint.
//   Sign: the ORIENTATION FOLD, as in `stable`, with mu_dyad and
//     mu_dyad_raw in the roles of theta0 and z_theta0. The likelihood is
//     invariant under negating alpha and theta together (mu_dyad,
//     mu_dyad_raw and theta_raw at every t); the OU recursion is
//     consistent with that negation because rho and process_noise are
//     unsigned. With s = orientation_sign(alpha_raw), the reported alpha
//     and mu_dyad, and each dyad's innovation terms, are multiplied by s.
//     By induction every theta[d, t] then carries the same s:
//     theta[d, t - 1] - mu_dyad[d] carries s once theta[d, t - 1] does,
//     and multiplying by the unsigned rho preserves it.
//
// ORIENTATION: the reported alpha[1] is never negative, so higher theta
// makes the reference action class more likely relative to the other
// action types (more cooperative relations, with a cooperative reference
// class) -- matching stable, and holding for every fit made under this
// file via the fold. `alpha`, `theta`, `mu_dyad`, `mu_intercept`, `phi`,
// and every scale/dispersion/ratio quantity (`sigma_mu`, `sd_stat`,
// `process_noise`, `rho`, `within_between_ratio`) are always in their
// final orientation as REPORTED; no post-hoc relabeling is needed. The
// RAW parameters the fold consumes (`alpha_raw`, `mu_dyad_raw`,
// `theta_raw`) remain sign-ambiguous themselves and are left out of
// bilatr's tiered convergence diagnostics for that reason (see
// ?diagnose_convergence).
//
// GENERATED QUANTITIES
//   within_between_ratio = exp(mu_log_sd_stat), always computed (cheap).
//   log_lik (when compute_log_lik = 1) and theta_filtered /
//   theta_filtered_sd (when compute_theta_filtered = 1), as in `stable`,
//   with the OU prediction step in the forward filter.
functions {
  // <<< BEGIN GENERATED partial_log_lik (source: inst/stan/include/partial_log_lik.stanfunctions) >>>
  // Do not hand-edit between these markers -- edit the source file
  // above and rerun `Rscript data-raw/sync_stan_functions.R`
  // (checked by tests/testthat/test_stan_includes.R).
  // Shared reduce_sum likelihood for all bilatr Stan model variants.
  //
  // Indifferent to how `theta`, `alpha` and `mu_intercept` were constructed
  // upstream (random walk or OU process, with or without country offsets):
  // it only consumes the already-built `theta`/`alpha`/`mu_intercept`/`phi`
  // values.
  //
  // This is the single canonical source. It is not included via Stan's
  // `#include`, which cmdstanr cannot resolve at sampling time when the
  // include path contains a space
  // (https://github.com/stan-dev/cmdstanr/issues/820). Instead, its contents
  // are copied verbatim into a marker-delimited block in each model's
  // `.stan` file, and the package tests check that the copies match.
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

  // Country-offset variant of partial_log_lik(), used by the `stable_gamma`
  // model (and present, unused, in the other models' copies of this file).
  //
  // The per-dyad country offset g_d is built once per dyad, outside the `t`
  // loop, and folded into the intercept (mu_eff = mu_intercept + g_d), so
  // the per-cell expression has the same shape and autodiff cost as
  // partial_log_lik()'s. g_d is not stored: a per-dyad matrix in
  // transformed parameters would be written to the output for every draw.
  // w_send[d] == 1 covers directed data (g_d = gamma[, ctry_a[d]] exactly,
  // ctry_b[d] unused); any other value mixes ctry_a's and ctry_b's columns
  // by event share (see the `stable_gamma` model's header).
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

  // Per-dyad-period log-likelihood, not summed -- used only by the
  // compute_log_lik generated quantities below. Must stay numerically
  // identical to the per-cell term inside partial_log_lik() above (same
  // eta/softmax/conc/dirichlet_multinomial_lpmf), and to the copy in
  // bilatr_alphanorm.stan.
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
  // to alpha/mu_dyad/theta below; see header, "Sign: the ORIENTATION
  // FOLD". The same function is written out in each model file.
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
                                       // no free location parameter (see header)
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
  // The sign comes from orientation_sign() (functions block above); see
  // header, "Sign: the ORIENTATION FOLD".
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
  // direction of alpha_raw (see bilatr_alphanorm.stan's header, "Radial
  // direction"). The orientation fold (transformed parameters, above)
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
  // the ratio of within-dyad to between-dyad SD (see header)
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
  // quantities block for the derivation (identical score/Fisher-
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
