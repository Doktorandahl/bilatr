// bilatr model `stable`: hierarchical dynamic IRT model for latent
// bilateral relations, with a Gaussian random walk for each dyad's latent
// state. The model fit_panel() and fit_dyad_ts() fit by default.
//
// PURPOSE
//   Measures the latent state theta[d, t] of each dyad d (a directed or
//   undirected pair of actors) in each period t from the counts of its
//   events in A action classes. Each class a has a discrimination
//   alpha[a] and an intercept mu_intercept[a]; higher theta shifts a
//   dyad's event mix toward the classes with larger alpha.
//
// DATA
//   T, D, A        number of periods, dyads and action classes
//   C              reduce_sum grainsize: dyads per task (chunk_size in R)
//   is_obs[d, t]   1 if dyad d has at least one event in period t, else 0
//   Y[d, t, a]     event counts
//   compute_log_lik, prior_only, compute_theta_filtered   0/1 switches
//   n_filter_dyads, filter_dyads   the dyads to forward-filter
//   Action class 1 is the reference class: prepare_stan_data() puts its
//   `reference_category` first.
//
// LIKELIHOOD (dyad-periods with is_obs[d, t] = 1 only)
//   eta[d, t] = alpha * theta[d, t] - mu_intercept           (length A)
//   Y[d, t]   ~ DirichletMultinomial(phi[d] * softmax(eta[d, t]))
//   phi[d] is the dyad's concentration (smaller phi, more overdispersion).
//   The sum over dyads is parallelized with reduce_sum (partial_log_lik()).
//   Periods without events add nothing; theta there is informed by the
//   latent process alone. prior_only = 1 drops the likelihood.
//
// LATENT PROCESS (non-centered; the sign factor of the orientation fold,
// below, is left out here)
//   theta0[d]        = sigma_theta0 * z_theta0[d]
//   theta[d, 1]      = theta0[d] + process_noise[d] * theta_raw[d, 1]
//   theta[d, t]      = theta[d, t - 1] + process_noise[d] * theta_raw[d, t]
//   process_noise[d] = sigma_theta0 * exp(mu_log_noise
//                        + sigma_log_noise * log_process_noise_raw[d])
//   process_noise is a ratio to sigma_theta0: exp(mu_log_noise) is the
//   typical per-period innovation SD as a fraction of the cross-dyad SD
//   of the initial states.
//
// PRIORS
//   theta_raw, z_theta0, log_process_noise_raw ~ std_normal()
//   sigma_theta0    ~ half-normal(0, 2)
//   mu_log_noise    ~ normal(log(0.2), 0.5)
//   sigma_log_noise ~ half-normal(0, 0.5)
//   phi             ~ lognormal(mu_log_phi, sigma_log_phi),
//                     mu_log_phi ~ normal(0, 1), sigma_log_phi ~ half-normal(0, 1)
//   mu_intercept    ~ std_normal()
//   alpha_raw       ~ std_normal()   (load-bearing; see "Radial direction")
//   These priors have not been calibrated by prior-predictive simulation.
//
// IDENTIFICATION
//   Scale. alpha = alpha_raw * sqrt(A / dot_self(alpha_raw)), so alpha has
//     root-mean-square 1 (equal to its SD, since it sums to zero). The unit
//     of theta is set by the discrimination profile of the action classes,
//     not by the dyads in the sample, so fits of the same class scheme on
//     different dyad sets share a scale.
//   Radial direction. The normalization leaves the length of alpha_raw
//     unidentified by the likelihood; alpha_raw ~ std_normal() identifies
//     it (dot_self(alpha_raw) concentrates near A - 1). A posterior of
//     dot_self(alpha_raw) near 0 signals a numerical problem.
//   Softmax location. alpha_raw and mu_intercept are sum_to_zero_vector[A],
//     so neither carries a free level shift.
//   Location of theta. theta0 has no free mean parameter: its cross-dyad
//     mean is pinned at 0 by z_theta0 ~ std_normal(). A common shift b
//     costs D * b^2 / 2 in log density, a stiff pin for large D, but a
//     prior rather than an exact constraint.
//   Sign: the ORIENTATION FOLD. eta depends on alpha and theta only through
//     alpha .* theta, which is unchanged when both are negated (theta0,
//     z_theta0 and theta_raw at every t negated together), and every prior
//     on the negated quantities is symmetric about 0. The posterior
//     therefore has two mirror-image modes of equal mass, separated by a
//     large likelihood barrier. A sign constraint cannot resolve this: a
//     chain initialized on the excluded side cannot cross the barrier and
//     settles on the constraint boundary instead. The fold changes what is
//     reported, not what the sampler can reach. With
//     s = orientation_sign(alpha_raw), the sign of alpha_raw[1], the
//     reported alpha, theta0 and theta are all multiplied by s:
//     - the likelihood is unchanged, since alpha .* theta is;
//     - the target stays smooth: the priors are on the raw parameters,
//       which s does not touch, and crossing alpha_raw[1] = 0 flips alpha
//       and theta together;
//     - no Jacobian adjustment is needed: it is the same model, reported
//       under a canonical sign;
//     - it is exact: both raw-space modes map to the same reported values.
//     s is computed by a function rather than stored, because a stored
//     transformed parameter would be written to the output, where its
//     cross-chain R-hat would be meaningless.
//
// ORIENTATION: the reported alpha[1] is never negative, so higher theta
// makes the reference action class more likely relative to the other
// action types; with a cooperative reference class, higher theta means
// more cooperative relations (see ?prepare_stan_data, reference_category).
// This holds for every fit made under this file via the fold above, not
// just on average across many. `alpha`, `theta`, `mu_intercept`, `phi`,
// and every scale/dispersion parameter are always in their final
// orientation as REPORTED; no post-hoc relabeling is needed. The RAW
// parameters the fold consumes (alpha_raw, z_theta0, theta_raw) remain
// sign-ambiguous themselves: if two chains land in opposite raw-space
// basins, those three variables' own cross-chain R-hat is meaningless even
// though everything reported is fine, which is why bilatr's tiered
// convergence diagnostics leave them out (see ?diagnose_convergence).
//
// GENERATED QUANTITIES
//   log_lik[d, t]   per-dyad-period log-likelihood (0 where is_obs is 0),
//                   when compute_log_lik = 1 (D x T values per draw).
//   theta_filtered[i, t], theta_filtered_sd[i, t]
//                   for dyad filter_dyads[i], when compute_theta_filtered
//                   = 1: the mean and SD of theta[d, t] given events up to
//                   and including period t, from a linear-Bayes forward
//                   filter run conditional on the draw's hyperparameters
//                   (see the comment in the block).
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
  // eta/softmax/conc/dirichlet_multinomial_lpmf). Written out here rather
  // than generated from the shared source because it is small and used
  // only in generated quantities.
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
  // to BOTH alpha and theta below; see header, "Sign: the ORIENTATION
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
  real mu_log_noise;                // log ratio to sigma_theta0, see header
  real<lower=0> sigma_log_noise;
  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;
}
transformed parameters {
  // Normalize alpha to RMS (== population SD, sum-to-zero) 1. This is the
  // scale anchor for the whole model; see header, "Scale" and "Radial
  // direction". The sign comes from orientation_sign(); see header,
  // "Sign: the ORIENTATION FOLD".
  vector[A] alpha = orientation_sign(alpha_raw)
                    * alpha_raw * sqrt((1.0 * A) / dot_self(alpha_raw));

  // theta0 folds the same sign as alpha, so alpha .* theta (all the
  // likelihood ever sees) is invariant to which raw-space basin the
  // sampler is in -- see header.
  vector[D] theta0 = orientation_sign(alpha_raw) * sigma_theta0 * z_theta0;

  // Per-dyad process noise, a RATIO to sigma_theta0 (non-centered
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
  // direction of alpha_raw (see header, "Radial direction"). The
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
  // derived by hand: with the total concentration phi[d] constant in
  // theta, d/dtheta log P(y | conc) reduces to
  // sum_k phi[d]*p_k*(alpha_k - a_bar) * (digamma(y_k + phi[d]*p_k) -
  // digamma(phi[d]*p_k)); the Fisher information uses the standard DM
  // overdispersion correction n*(1+phi)/(n+phi) applied to
  // Var_pi(alpha) = dot_product(p, square(alpha - a_bar)) (the same
  // quantity diagnose_category_merges() uses). No
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
