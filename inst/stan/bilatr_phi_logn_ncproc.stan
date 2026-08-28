// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// Variant: identical to bilatr_phi_logn.stan EXCEPT that the per-dyad
// process_noise hierarchy is also sampled in non-centered form. Instead of
// process_noise[d] ~ lognormal(mu_log_noise, sigma_log_noise) directly, we
// sample log_process_noise_raw[d] ~ std_normal() and set
// process_noise[d] = exp(mu_log_noise + sigma_log_noise * log_process_noise_raw[d])
// in transformed parameters -- the same reparameterization already applied to
// the log_phi0 hierarchy in this file, now extended to process_noise to break
// its funnel (mu_log_noise / sigma_log_noise tracked lp__ Rhat/ESS degradation
// in the centered form). The phi-as-log(n) logic is otherwise unchanged.
//
// phi (Dirichlet-multinomial concentration) is a per-dyad-period quantity,
// modeled as log(phi[d,t]) = log_phi0[d] + beta_logn * centered_log_n[d,t],
// where log_phi0[d] retains the original per-dyad hierarchical structure
// (mu_log_phi, sigma_log_phi) and beta_logn is a single global slope capturing
// how overdispersion scales with observed count magnitude. Setting beta_logn's
// posterior near 0 recovers the original per-dyad-constant-phi model.
//
// dyad_weight / period_weight / action_weight retain their previous role and
// are unaffected by this change; setting all three to 1s plus beta_logn's
// prior alone does NOT recover the exact old model bit-for-bit (phi is now
// time-varying via log_phi0 + beta_logn*0 at the mean), but is equivalent in
// spirit and nests it as beta_logn -> 0.
//
// Identification:
//   - alpha[1]  = 1  (scale reference, neutral action)
//   - mu_intercept[1] = 0  (softmax level-shift reference)
//   - alpha[A]  = -alpha_hostile, alpha_hostile > 0 (known hostile action
//     anchors the negative end of the discrimination scale)
//   - theta0 is fully hierarchical (mu_theta0, sigma_theta0), not supplied
//     as data
//   - no dyad-specific intercept: cross-dyad level differences are forced
//     into theta via the global mu_intercept, keeping theta comparable
//     across dyads
functions {
  real partial_log_lik(array[] int slice_d,
                        int start, int end,
                        int T, int A,
                        array[,] int is_obs,
                        array[,,] int Y,
                        array[,] real theta,
                        vector mu_intercept,
                        array[,] real phi,          // CHANGED: was vector[D]
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
          vector[A] conc = phi[d, t] * (action_weight .* p);   // CHANGED: phi[d,t]
          // dyad_weight/period_weight scale each dyad-period's contribution
          lp += dyad_weight[d] * period_weight[t] *
            dirichlet_multinomial_lpmf(Y[d, t] | conc);
        }
      }
    }
    return lp;
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
}
transformed data {
  // Precompute log dyad-period totals once; fixed function of the data,
  // so no need to touch this inside the sampler loop.
  array[D, T] real log_n;
  real mean_log_n;
  {
    real sum_log_n = 0;
    int n_obs_cells = 0;
    for (d in 1:D) {
      for (t in 1:T) {
        if (is_obs[d, t] == 1) {
          int n_dt = sum(Y[d, t]);
          // guard log(0) in case an observed cell has all-zero counts
          log_n[d, t] = log(n_dt > 0 ? n_dt : 1);
          sum_log_n += log_n[d, t];
          n_obs_cells += 1;
        } else {
          log_n[d, t] = 0; // unused; likelihood is gated by is_obs
        }
      }
    }
    mean_log_n = sum_log_n / n_obs_cells;
  }
}
parameters {
  // Latent states per dyad
  array[D, T] real theta_raw;
  vector[A - 1] mu_intercept_raw;   // global action baseline (index 2:A)
  real mu_theta0;
  real<lower=0> sigma_theta0;
  vector[D] z_theta0;

  // Hierarchical process/dispersion parameters
  vector[D] log_process_noise_raw;  // non-centered: process_noise built in TP
  real mu_log_noise;
  real<lower=0> sigma_log_noise;

  // Dispersion: per-dyad baseline (non-centered) + global log(n) slope
  vector[D] log_phi0_raw;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;
  real beta_logn;                   // NEW: slope of log(phi) on centered log(n)

  // Global discrimination parameters
  vector[A - 2] alpha_raw;
  real<lower=0> alpha_hostile;      // known positive, used as negative anchor
}
transformed parameters {
  vector[A] alpha;
  alpha[1] = 1;
  alpha[A] = -alpha_hostile;
  for (i in 2:(A - 1)) {
    alpha[i] = alpha_raw[i - 1];
  }

  vector[A] mu_intercept;
  mu_intercept[1] = 0;
  for (i in 2:A) {
    mu_intercept[i] = mu_intercept_raw[i - 1];
  }

  vector[D] theta0;
  for (d in 1:D) {
    theta0[d] = mu_theta0 + sigma_theta0 * z_theta0[d];
  }

  // Per-dyad process noise (non-centered lognormal hierarchy)
  vector<lower=0>[D] process_noise =
    exp(mu_log_noise + sigma_log_noise * log_process_noise_raw);

  array[D, T] real theta;
  for (d in 1:D) {
    theta[d, 1] = theta0[d] + process_noise[d] * theta_raw[d, 1];
    for (t in 2:T) {
      theta[d, t] = theta[d, t - 1] + process_noise[d] * theta_raw[d, t];
    }
  }

  // Per-dyad baseline dispersion (non-centered), then per-dyad-period phi
  // as a function of centered log(n_dt).
  vector[D] log_phi0;
  for (d in 1:D) {
    log_phi0[d] = mu_log_phi + sigma_log_phi * log_phi0_raw[d];
  }

  array[D, T] real<lower=0> phi;
  for (d in 1:D) {
    for (t in 1:T) {
      phi[d, t] = exp(log_phi0[d] + beta_logn * (log_n[d, t] - mean_log_n));
    }
  }
}
model {
  // hyperpriors for pooling (lognormal hierarchy)
  mu_log_noise ~ normal(log(0.2), 0.5);
  sigma_log_noise ~ normal(0, 0.5);
  log_process_noise_raw ~ std_normal();

  mu_log_phi ~ normal(0, 1);
  sigma_log_phi ~ normal(0, 1);
  log_phi0_raw ~ std_normal();
  beta_logn ~ normal(0, 0.5);   // shrinks toward the old constant-phi model

  // theta random-walk innovations
  for (d in 1:D) {
    for (t in 1:T) {
      theta_raw[d, t] ~ std_normal();
    }
  }

  mu_theta0 ~ normal(0, 1);
  sigma_theta0 ~ normal(0, 0.5);
  z_theta0 ~ std_normal();

  mu_intercept_raw ~ std_normal();

  alpha_raw ~ std_normal();
  alpha_hostile ~ normal(0, 1);

  // likelihood, chunked via reduce_sum
  array[D] int dyad_seq = linspaced_int_array(D, 1, D);
  target += reduce_sum(partial_log_lik, dyad_seq, C,
                        T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                        dyad_weight, period_weight, action_weight);
}
