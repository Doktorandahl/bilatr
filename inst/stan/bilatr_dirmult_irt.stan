// bilatr: hierarchical dynamic IRT model for dyadic conflict trajectories.
//
// The stable model. Single collapsed program covering all three
// previously-separate variants (base / event-weighted /
// dyad+period+event-weighted) via the dyad_weight, period_weight,
// action_weight data vectors. Setting all three to vectors of 1s recovers
// the unweighted base model exactly; setting only dyad_weight and
// period_weight to 1s recovers the event-weighted-only model exactly. See
// package vignettes for details.
//
// The per-dyad process_noise hierarchy is sampled in non-centered form:
// log_process_noise_raw[d] ~ std_normal(), with
// process_noise[d] = exp(mu_log_noise + sigma_log_noise * log_process_noise_raw[d])
// built in transformed parameters. This broke the process_noise funnel
// (mu_log_noise / sigma_log_noise previously tracked lp__ Rhat/ESS
// degradation in the centered form). The pre-0.3.0 centered version is kept
// under inst/stan/legacy/ for reference.
//
// Identification:
//   - alpha[1]  = 1  (scale/sign reference, neutral action); alpha[2:A] are
//     all freely estimated (no known-hostile anchor)
//   - mu_intercept[1] = 0  (softmax level-shift reference)
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
  vector<lower=0>[D] phi;
  real mu_log_phi;
  real<lower=0> sigma_log_phi;

  // Global discrimination parameters (alpha[1] fixed to 1; the rest free)
  vector[A - 1] alpha_raw;
}
transformed parameters {
  vector[A] alpha;
  alpha[1] = 1;
  for (i in 2:A) {
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
}
model {
  // hyperpriors for pooling (lognormal hierarchy)
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

  mu_theta0 ~ normal(0, 1);
  sigma_theta0 ~ normal(0, 0.5);
  z_theta0 ~ std_normal();

  mu_intercept_raw ~ std_normal();

  alpha_raw ~ std_normal();

  // likelihood, chunked via reduce_sum
  array[D] int dyad_seq = linspaced_int_array(D, 1, D);
  target += reduce_sum(partial_log_lik, dyad_seq, C,
                        T, A, is_obs, Y, theta, mu_intercept, phi, alpha,
                        dyad_weight, period_weight, action_weight);
}
